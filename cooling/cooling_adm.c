#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"
#include "./adm_cooling_functions.h"

/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM,
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by Volker Springel. The code has been modified heavily by
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; essentially everything has been re-written at this point */


// Note that I have commented out all the FIRE flags i.e.STELLAREVOLUTION, RT_HIIHEATING, UVHEATING, RT_USE_GRAVTREE.... This should be changed in the future when ADM feedback is implemented.

#ifdef COOLING
#ifdef ADM

/* these are variables of the cooling tables. they are static but this shouldnt be a problem for shared-memory structure because
    they are only defined once in a global operation, then locked for particle-by-particle operations */
/* requires the cooling table TREECOOL, which is included in the GIZMO source in the cooling directory */
#define NCOOLTAB_ADM  2000 /* defines size of cooling table */

#if !defined(CHIMES)
static double Tmin_adm = -1.0, Tmax_adm = 9.0, deltaT_adm; /* minimum/maximum temp, in log10(T/K) and temperature gridding: will be appropriately set in make_cooling_tables subroutine below */
static double *BetaH0_adm, *BetaHep_adm, *Betaff_adm, *AlphaHpRate_adm, *AlphaHp_adm, *AlphaHep_adm, *Alphad_adm, *AlphaHepp_adm, *GammaeH0_adm, *GammaeHe0_adm, *GammaeHep_adm; // UV background parameters
#ifdef COOL_METAL_LINES_BY_SPECIES
/* if this is enabled, the cooling table files should be in a folder named 'spcool_tables' in the run directory.
 cooling tables can be downloaded at: http://www.tapir.caltech.edu/~phopkins/public/spcool_tables.tgz or on the Bitbucket site (downloads section) */
static float *SpCoolTable0_adm, *SpCoolTable1_adm;
#endif
/* these are constants of the UV background at a given redshift: they are interpolated from TREECOOL but then not modified particle-by-particle */
static double J_UV_adm = 0, gJH0_adm = 0, gJHep_adm = 0, gJHe0_adm = 0, epsH0_adm = 0, epsHep_adm = 0, epsHe0_adm = 0;
#endif

#if defined(CHIMES)
int ChimesEqmMode, ChimesUVBMode, ChimesInitIonState, N_chimes_full_output_freq, Chimes_incl_full_output = 1;
double chimes_rad_field_norm_factor, shielding_length_factor, cr_rate;
char ChimesDataPath[256], ChimesEqAbundanceTable[196], ChimesPhotoIonTable[196];
struct gasVariables *ChimesGasVars;
struct globalVariables ChimesGlobalVars;
#ifdef CHIMES_METAL_DEPLETION
struct Chimes_depletion_data_structure *ChimesDepletionData;
#endif
#endif




/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle_adm(int i)
{
    printf("ADM Alert! cooling_adm.c We have a problem!\n");
    double unew, dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);

    if((dtime>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction_adm(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif
        double uold = DMAX(All.MinEgySpec_adm, SphP[i].InternalEnergy);
//////////////////////
// STELLAREVOLUTION //
// ///////////////////

//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
//        double uion=HIIRegion_Temp/(0.59*(5./3.-1.)*U_TO_TEMP_UNITS); if(SphP[i].DelayTimeHII>0) {if(uold<uion) {uold=uion;}} /* u_old should be >= ionized temp if used here [unless using newer model] */
//#endif

#ifndef COOLING_OPERATOR_SPLIT
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        double grav_acc; int k;
        for(k = 0; k < 3; k++)
        {
            grav_acc = All.cf_a2inv * P[i].GravAccel[k];
#ifdef PMGRID
            grav_acc += All.cf_a2inv * P[i].GravPM[k];
#endif
            SphP[i].DtInternalEnergy -= SphP[i].GravWorkTerm[k] * All.cf_atime * grav_acc;
        }
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
        SphP[i].DtInternalEnergy = DMAX(SphP[i].DtInternalEnergy , -0.99*SphP[i].InternalEnergy/dtime ); // equivalent to saying this wouldn't lower internal energy to below 1% in one timestep
        SphP[i].DtInternalEnergy = DMIN(SphP[i].DtInternalEnergy ,  1.e4*SphP[i].InternalEnergy/dtime ); // equivalent to saying we cant massively enhance internal energy in a single timestep from the hydro work terms: should be big, since just numerical [shocks are real!]
        /* and convert to cgs before use in the cooling sub-routine */
        SphP[i].DtInternalEnergy *= (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (All.ADM_ProtonMass/HYDROGEN_MASSFRAC_ADM);
#endif


#ifndef RT_COOLING_PHOTOHEATING_OLDFORMAT
        /* Call the actual COOLING subroutine! */
#ifdef CHIMES
        double dummy_ne = 0.0;
        unew = DoCooling_adm(uold, SphP[i].Density * All.cf_a3inv, dtime, dummy_ne, i);
#else
        unew = DoCooling_adm(uold, SphP[i].Density * All.cf_a3inv, dtime, SphP[i].Ne, i);
#endif
#else
        unew = uold + dtime * (rt_DoHeating(i, dtime) + rt_DoCooling(i, dtime));
#endif

////////////////////////
//// STELLAREVOLUTION //
//// ///////////////////

//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING) /* for older model, set internal energy to minimum level if marked as ionized by stars */
//        if(SphP[i].DelayTimeHII > 0)
//        {
//            if(unew<uion) {unew=uion; if(SphP[i].DtInternalEnergy<0) SphP[i].DtInternalEnergy=0;}
//#ifndef CHIMES
//            SphP[i].Ne = 1.0 + 2.0*yhelium(i); /* fully ionized. note that this gives Ne as free electron fraction per H */
//#endif
//        }
//#endif


#if defined(BH_THERMALFEEDBACK)
        if(SphP[i].Injected_BH_Energy) {unew += SphP[i].Injected_BH_Energy / P[i].Mass; SphP[i].Injected_BH_Energy = 0;}
#endif


#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_ALT_DISABLE_LOSSES)
        CR_cooling_and_losses(i, SphP[i].Ne, SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, dtime*UNIT_TIME_IN_CGS );
#endif


#ifdef RT_INFRARED /* assume (for now) that all radiated/absorbed energy comes from the IR bin [not really correct, this should just be the dust term] */
        double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS;	/* hydrogen number dens in cgs units */
        double ratefact = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * nHcgs * nHcgs / (SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS); /* need to account for RSOL factors in emission/absorption rates */
        double de_u = -SphP[i].LambdaDust * ratefact * (dtime*UNIT_TIME_IN_CGS) / (UNIT_SPECEGY_IN_CGS) * P[i].Mass; /* energy gained by gas needs to be subtracted from radiation */
        if(de_u<=-0.99*SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]) {de_u=-0.99*SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]; unew=DMAX(0.01*SphP[i].InternalEnergy , SphP[i].InternalEnergy-de_u/P[i].Mass);}
        SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] += de_u; /* energy gained by gas is lost here */
        SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED] = SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]; /* updated drifted */
#if defined(RT_EVOLVE_INTENSITIES)
        int k_tmp; for(k_tmp=0;k_tmp<N_RT_INTENSITY_BINS;k_tmp++) {SphP[i].Rad_Intensity[RT_FREQ_BIN_INFRARED][k_tmp] += de_u/RT_INTENSITY_BINS_DOMEGA; SphP[i].Rad_Intensity_Pred[RT_FREQ_BIN_INFRARED][k_tmp] += de_u/RT_INTENSITY_BINS_DOMEGA;}
#endif
        int kv; // add leading-order relativistic corrections here, accounting for gas motion in the addition/subtraction to the flux:
#if defined(RT_EVOLVE_FLUX)
        for(kv=0;kv<3;kv++) {double fluxfac = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE)*SphP[i].VelPred[kv]/All.cf_atime * de_u;
            SphP[i].Rad_Flux[RT_FREQ_BIN_INFRARED][kv] += fluxfac; SphP[i].Rad_Flux_Pred[RT_FREQ_BIN_INFRARED][kv] += fluxfac;}
#endif
        double momfac = 1. - de_u / (P[i].Mass * C_LIGHT_CODE*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission [note peculiar units here, its b/c of how we fold in the existing value of v and tilde[u] in our derivation - one rsol factor in denominator needed]
        for(kv=0;kv<3;kv++) {P[i].Vel[kv] *= momfac; SphP[i].VelPred[kv] *= momfac;}
#endif


        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated.
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        SphP[i].InternalEnergy = unew;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
#ifndef COOLING_OPERATOR_SPLIT
        SphP[i].DtInternalEnergy = 0;
#endif

#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction_adm(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif
//#if defined(GALSF_FB_FIRE_RT_HIIHEATING) /* count off time which has passed since ionization 'clock' */
//        if(SphP[i].DelayTimeHII > 0) {SphP[i].DelayTimeHII -= dtime;}
//        if(SphP[i].DelayTimeHII < 0) {SphP[i].DelayTimeHII = 0;}
//#endif

    } // closes if((dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}




/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling_adm(double u_old, double rho, double dt, double ne_guess, int target)
{
    double u, du; u=0; du=0;
    printf("ADM Alert! cooling_adm, doCooling_adm\n");

#ifdef COOL_GRACKLE
#ifndef COOLING_OPERATOR_SPLIT
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below to split the step */
    du = dt * SphP[target].DtInternalEnergy / ( (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (All.ADM_ProtonMass/HYDROGEN_MASSFRAC_ADM));
    u_old += 0.5*du;
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
    /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
     term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
     the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
     to which we then add the heating term from hydro and compute the solution over a full timestep */
    double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(r-1)/log(r);} r=DMAX(0,DMIN(r,1));
    du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
    return DMAX(u,All.MinEgySpec);
#endif

#ifdef CHIMES
    chimes_update_gas_vars(target);

    /* Call CHIMES to evolve the chemistry and temperature over
     * the hydro timestep. */
    chimes_network(&(ChimesGasVars[target]), &ChimesGlobalVars);

    // Compute updated internal energy
    u = (double) ChimesGasVars[target].temperature * BOLTZMANN / ((GAMMA(target)-1) * PROTONMASS * calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars));
    u /= UNIT_SPECEGY_IN_CGS;  // code units

#ifdef CHIMES_TURB_DIFF_IONS
    chimes_update_turbulent_abundances(target, 1);
#endif

    return DMAX(u, All.MinEgySpec);

#else // CHIMES

    int iter=0, iter_upper=0, iter_lower=0, iter_condition = 0; double LambdaNet, ratefact, u_upper, u_lower;
#ifdef RT_INFRARED
    double LambdaDust;
#endif
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC_ADM * rho / All.ADM_ProtonMass;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;

    u = u_old; u_lower = u; u_upper = u; /* initialize values */
    LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);

    /* bracketing */
    if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
        u_upper *= sqrt(1.1); u_lower /= sqrt(1.1);
        while((iter_upper<MAXITER)&&(u_upper - u_old - ratefact * CoolingRateFromU_adm(u_upper, rho, ne_guess, target) * dt < 0))
        {
            u_upper *= 1.1; u_lower *= 1.1; iter_upper++;
        }

    }

    if(u - u_old - ratefact * LambdaNet * dt > 0) /* cooling */
    {
        u_lower /= sqrt(1.1); u_upper *= sqrt(1.1);
        while((iter_lower<MAXITER)&&(u_lower - u_old - ratefact * CoolingRateFromU_adm(u_lower, rho, ne_guess, target) * dt > 0))
        {
            u_upper /= 1.1; u_lower /= 1.1; iter_lower++;
        }
    }

    /* core iteration to convergence */
    do
    {
        u = 0.5 * (u_lower + u_upper);
#ifdef RT_INFRARED
        LambdaDust = SphP[target].LambdaDust;
#endif
        LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);
        if(u - u_old - ratefact * LambdaNet * dt > 0) {u_upper = u;} else {u_lower = u;}
        du = u_upper - u_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("u=%g u_old=%g u_upper=%g u_lower=%g ne_guess=%g dt=%g iter=%d \n", u,u_old,u_upper,u_lower,ne_guess,dt,iter);}

        iter_condition = ((fabs(du/u) > 3.0e-2)||((fabs(du/u) > 3.0e-4)&&(iter < 10)));
#ifdef RT_INFRARED
        iter_condition = iter_condition || (((fabs(LambdaDust - SphP[target].LambdaDust) > 1e-2*fabs(LambdaDust)) || (fabs(u - u_old - ratefact * LambdaNet * dt) > 0.01*fabs(u-u_old)))  && (iter < MAXITER-11));
#endif
        iter_condition = iter_condition &&  (iter < MAXITER); // make sure we don't iterate more than MAXITER times

    }
    while(iter_condition); /* iteration condition */
    /* crash condition */
    if(iter >= MAXITER) {printf("failed to converge in DoCooling_adm(): u_in=%g rho_in=%g dt=%g ne_in=%g target=%d \n",u_old,rho,dt,ne_guess,target); endrun(10);}
    double specific_energy_codeunits_toreturn = u / UNIT_SPECEGY_IN_CGS;    /* in internal units */

#ifdef RT_CHEM_PHOTOION
    /* set variables used by RT routines; this must be set only -outside- of iteration, since this is the key chemistry update */
    double u_in=specific_energy_codeunits_toreturn, rho_in=SphP[target].Density*All.cf_a3inv, mu=1, temp, ne=1, nHI=SphP[target].HI, nHII=SphP[target].HII, nHeI=1, nHeII=0, nHeIII=0;
    temp = ThermalProperties_adm(u_in, rho_in, target, &mu, &ne, &nHI, &nHII, &nHeI, &nHeII, &nHeIII);
    SphP[target].HI = nHI; SphP[target].HII = nHII;
#ifdef RT_CHEM_PHOTOION_HE
    SphP[target].HeI = nHeI; SphP[target].HeII = nHeII; SphP[target].HeIII = nHeIII;
#endif
#endif

    /* safe return */
    return specific_energy_codeunits_toreturn;
#endif // CHIMES
}



#ifndef CHIMES
/* returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime_adm(double u_old, double rho, double ne_guess, int target)
{
printf("ADM Alert! cooling_adm, getcoolingtime\n");
#if defined(COOL_GRACKLE) && !defined(GALSF_EFFECTIVE_EQS)
    double LambdaNet = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(LambdaNet >= 0) LambdaNet = 0.0;
    return LambdaNet / UNIT_TIME_IN_CGS;
#else
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC_ADM * rho / All.ADM_ProtonMass;	/* hydrogen number dens in cgs units */
    double LambdaNet = CoolingRateFromU_adm(u_old, rho, ne_guess, target);
    if(LambdaNet >= 0) {return 0;} /* net heating due to UV background */
    return u_old / (-(nHcgs * nHcgs / rho) * LambdaNet) / UNIT_TIME_IN_CGS;
#endif
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoInstabilityCooling_adm(double m_old, double u, double rho, double dt, double fac, double ne_guess, int target)
{
    printf("ADM Alert! cooling_adm, instabilitycooling\n");
    if(fac <= 0) {return 0.01*m_old;} /* the hot phase is actually colder than the cold reservoir! */
    double m, dm, m_lower, m_upper, ratefact, LambdaNet;
    int iter = 0;

    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    fac /= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC_ADM * rho / All.ADM_ProtonMass;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho * fac;
    m = m_old; m_lower = m; m_upper = m;
    LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);

    /* bracketing */
    if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
        m_upper *= sqrt(1.1); m_lower /= sqrt(1.1);
        while(m_upper - m_old - m_upper * m_upper / m_old * ratefact * CoolingRateFromU_adm(u, rho * m_upper / m_old, ne_guess, target) * dt < 0)
        {
            m_upper *= 1.1; m_lower *= 1.1;
        }
    }
    if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
        m_lower /= sqrt(1.1); m_upper *= sqrt(1.1);
        while(m_lower - m_old - m_lower * m_lower / m_old * ratefact * CoolingRateFromU_adm(u, rho * m_lower / m_old, ne_guess, target) * dt > 0)
        {
            m_upper /= 1.1; m_lower /= 1.1;
        }
    }

    do
    {
        m = 0.5 * (m_lower + m_upper);
        LambdaNet = CoolingRateFromU_adm(u, rho * m / m_old, ne_guess, target);
        if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0) {m_upper = m;} else {m_lower = m;}
        dm = m_upper - m_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("->m= %g\n", m);}
    }
    while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);
    if(iter >= MAXITER) {printf("failed to converge in DoInstabilityCooling_adm(): m_in=%g u_in=%g rho=%g dt=%g fac=%g ne_in=%g target=%d \n",m_old,u,rho,dt,fac,ne_guess,target); endrun(11);}
    return m;
}

#endif // !(CHIMES)






#ifdef CHIMES
/* This function converts thermal energy to temperature, using the mean molecular weight computed from the non-equilibrium CHIMES abundances. */
double chimes_convert_u_to_temp(double u, double rho, int target)
{
  return u * (GAMMA(target)-1) * PROTONMASS * ((double) calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars)) / BOLTZMANN;
}

#else  // CHIMES

/* this function determines the electron fraction, and hence the mean molecular weight. With it arrives at a self-consistent temperature.
 * Ionization abundances and the rates for the emission are also computed */
double convert_u_to_temp_adm(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess)
{
    printf("ADM Alert! cooling_adm, convert u to temp\n");
    int iter = 0;
    double temp, temp_old, temp_old_old = 0, temp_new, prefac_fun_old, prefac_fun, fac, err_old, err_new, T_bracket_errneg = 0, T_bracket_errpos = 0, T_bracket_min = 0, T_bracket_max = 1.e20, bracket_sign = 0; // double max = 0;
    double u_input = u, rho_input = rho, temp_guess;
    double T_0 = u * All.ADM_ProtonMass / BOLTZMANN; // this is the dimensional temperature, which since u is fixed is -frozen- in this calculation: we can work dimensionlessly below
    temp_guess = (GAMMA(target)-1) * T_0; // begin assuming mu ~ 1
    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp_guess, rho, nH0_guess, ne_guess, 0., target); // get mu with that temp
    prefac_fun = (GAMMA(target)-1) * (*mu_guess); // dimensionless pre-factor determining the temperature
    err_new = prefac_fun - temp_guess / T_0; // define initial error from this iteration
    if(err_new < 0) {T_bracket_errneg = temp_guess;} else {T_bracket_errpos = temp_guess;}
    temp = prefac_fun * T_0; // re-calculate temo with the new mu

    do
    {
        //qfun_old = *ne_guess; // guess for ne
        //qfun_old = *mu_guess; // guess for mu
        prefac_fun_old = prefac_fun;
        err_old = err_new; // error from previous timestep
        find_abundances_and_rates_adm(log10(temp), rho, target, -1, 0, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess); // all the thermo variables for this T
        prefac_fun = (GAMMA(target)-1) * (*mu_guess); // new value of the dimensionless pre-factor we need to solve
        temp_old = temp; // guess for T we just used
        temp_new = prefac_fun * T_0; // updated temp using the new values from the iteration of find_abundances_and_rates_adm above
        err_new = (temp_new - temp_old) / T_0; // new error
        if(T_bracket_errpos == 0) {if(err_new > 0) {T_bracket_errpos = temp_old;} else {T_bracket_errneg = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg == 0) {if(err_new < 0) {T_bracket_errneg = temp_old;} else {T_bracket_errpos = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg > 0 && T_bracket_errpos > 0)
        {
            if(bracket_sign == 0) {if(T_bracket_errpos > T_bracket_errneg) {bracket_sign=1;} else {bracket_sign=-1;}}
            if(err_new > 0) {
                if(bracket_sign > 0) {T_bracket_errpos = DMIN(T_bracket_errpos, temp_old); /* Tpos>Tneg */} else {T_bracket_errpos = DMAX(T_bracket_errpos, temp_old); /* Tpos<Tneg */}
            } else {
                if(bracket_sign > 0) {T_bracket_errneg = DMAX(T_bracket_errneg, temp_old); /* Tpos>Tneg */} else {T_bracket_errneg = DMIN(T_bracket_errneg, temp_old); /* Tpos<Tneg */}
            } /* update bracket values if we can */
            if(bracket_sign > 0) {T_bracket_max=T_bracket_errpos; T_bracket_min=T_bracket_errneg;} else {T_bracket_max=T_bracket_errneg; T_bracket_min=T_bracket_errpos;}
        }

        //max = DMAX(max, temp_new * (*mu_guess) * HYDROGEN_MASSFRAC * fabs((*ne_guess - qfun_old) / (temp_new - temp_old + 1.0))); // old iteration: hardwired assumption that ne is only varying quanity in mu, and that Tmin_adm ~ 1e4 or so
        //max = DMAX(max , temp_new / (*mu_guess) * fabs(*mu_guess - qfun_old) / (fabs(temp_new - temp_old) + 1.e-4*(All.MinGasTemp+0.1))); // newer - more flexible mu, and dimensionless T dependence
        //temp = temp_old + (temp_new - temp_old) / (1 + max);

        if(fabs(prefac_fun-prefac_fun_old) < 1.e-4 || fabs(temp_new-temp_old)/(temp_new+temp_old) < 1.e-4) {break;} // break pre-emptively if we'll trigger a nan below
        fac = (prefac_fun-prefac_fun_old)*T_0 / (temp_old-temp_old_old); // numerical derivative factor: want to use this to limit for convergence

        if(fac > 1) {fac = 1;} // don't allow us to move in the opposite direction from the new evaluation (should 'guess' in the direction of T_new-T_old) -- this tells us Newton-Raphson/Secant-type method fails here, so we simply follow the iteration to t_new
        if(fac > 0.9) {fac=0.9;} // don't allow us to 'jump' by a factor >10 times the temperature difference (arbitrary choice, slows convergence a bit but helps limit bad overshoot)
        if(fac < -9999.) {fac=-9999.;} // don't allow smaller step than 1e-4 times the temperature difference (since that's below our error tolerance anyways)

        temp = temp_old + (temp_new - temp_old) * 1./(1. - fac); // standard Newton-Raphson-type (technically Secant method) iteration
        if(temp < 0.5*temp_old) {temp = 0.5*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established
        if(temp > 3.0*temp_old) {temp = 3.0*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established

        temp = temp_old + (temp_new - temp_old) * 1./(1. - fac); // standard Newton-Raphson iteration

        if(T_bracket_errneg > 0 && T_bracket_errpos > 0) // if have bracketing and this wants to go outside brackets, revert to bisection
        {
            if(temp >= T_bracket_max || temp <= T_bracket_min) {temp = sqrt(T_bracket_min*T_bracket_max);} // bisect (in log-space)
        }
#ifndef RT_INFRARED
        if(fabs(temp-temp_old_old)/(temp+temp_old_old) < 1.e-3) {double wt=get_random_number(12*iter+340*ThisTask+5435*target); temp=(wt*temp_old + (1.-wt)*temp_new);}
#endif
        temp_old_old = temp_old;
        iter++;
        if(iter > (MAXITER - 10)) {printf("-> temp_next/new/old/oldold=%g/%g/%g/%g ne=%g mu=%g rho=%g iter=%d target=%d err_new/prev=%g/%g gamma_minus_1_mu_new/prev=%g/%g Brackets: Error_bracket_positive=%g Error_bracket_negative=%g T_bracket_Min/Max=%g/%g fac_for_SecantDT=%g \n", temp,temp_new,temp_old,temp_old_old,*ne_guess, (*mu_guess) ,rho,iter,target,err_new,err_old,prefac_fun,prefac_fun_old,T_bracket_errpos,T_bracket_errneg,T_bracket_min,T_bracket_max,fac); fflush(stdout);}
    }
    while(
#ifdef RT_INFRARED
        (fabs(temp - temp_old) > 1e-3 * temp) && iter < MAXITER);
#else
          ((fabs(temp - temp_old) > 0.25 * temp) ||
           ((fabs(temp - temp_old) > 0.1 * temp) && (temp > 20.)) ||
           ((fabs(temp - temp_old) > 0.05 * temp) && (temp > 200.)) ||
           ((fabs(temp - temp_old) > 0.01 * temp) && (temp > 200.) && (iter<100)) ||
           ((fabs(temp - temp_old) > 1.0e-3 * temp) && (temp > 200.) && (iter<10))) && iter < MAXITER);
#endif
    if(iter >= MAXITER) {printf("failed to converge in convert_u_to_temp_adm(): u_input= %g rho_input=%g n_elec_input=%g target=%d\n", u_input, rho_input, *ne_guess, target); endrun(12);}

    if(temp<=0) temp=pow(10.0,Tmin_adm);
    if(log10(temp)<Tmin_adm) temp=pow(10.0,Tmin_adm);
    return temp;
}
#endif // CHIMES




#ifndef CHIMES
/* this function computes the actual ionization states, relative abundances, and returns the ionization/recombination rates if needed */
double find_abundances_and_rates_adm(double logT, double rho, int target, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess,
                                 double *mu_guess)
{
    printf("ADM Alert! cooling_adm, abundances and rates\n");
    int j, niter;
    double Tlow, Thi, flow, fhi, t, gJH0ne, gJHe0ne, gJHepne, logT_input, rho_input, ne_input, neold, nenew;
    double bH0, bHep, bff, aHpRate, aHp, aHep, aHepp, ad, geH0, geHe0, geHep, EPSILON_SMALL=1.e-40;
    double n_elec, nH0, nHe0, nHp, nHep, nHepp; /* ionization states */
    logT_input = logT; rho_input = rho; ne_input = *ne_guess; /* save inputs (in case of failed convergence below) */
    if(!isfinite(logT)) {logT=Tmin_adm;}    /* nan trap (just in case) */
    if(!isfinite(rho)) {logT=Tmin_adm;}

    if(logT <= Tmin_adm)		/* everything neutral */
    {
        nH0 = 1.0; nHe0 = yhelium(target); nHp = 0; nHep = 0; nHepp = 0; n_elec = 1.e-22;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 0, target);
        return 0;
    }
    if(logT >= Tmax_adm)		/* everything is ionized */
    {
        nH0 = 0; nHe0 = 0; nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 1.e3, target);
        return 0;
    }

    /* initialize quantities needed for iteration below */
    t = (logT - Tmin_adm) / deltaT_adm;
    j = (int) t;
    if(j<0){j=0;}
    if(j>NCOOLTAB_ADM){
        PRINT_WARNING("j>NCOOLTAB_ADM : j=%d t %g Tlow %g Thi %g logT %g Tmin_adm %g deltaT_adm %g \n",j,t,Tmin_adm+deltaT_adm*j,Tmin_adm+deltaT_adm*(j+1),logT,Tmin_adm,deltaT_adm);fflush(stdout);
        j=NCOOLTAB_ADM;
    }
    Tlow = Tmin_adm + deltaT_adm * j;
    Thi = Tlow + deltaT_adm;
    fhi = t - j;
    flow = 1 - fhi;
    if(*ne_guess == 0) /* no guess provided, try to start from something sensible */
    {
        *ne_guess = 1.0;
        if(logT < 3.8) {*ne_guess = 0.1;}
        if(logT < 2) {*ne_guess = 1.e-10;}
    }
    /* CAFG: this is the density that we should use for UV background threshold */
    double local_gammamultiplier = return_local_gammamultiplier_adm(target); // account for local UVB terms in some expressions below
    double nHcgs = HYDROGEN_MASSFRAC_ADM * rho / All.ADM_ProtonMass;	/* hydrogen number dens in cgs units */
    if(shieldfac < 0) {shieldfac = return_uvb_shieldfac_adm(target, local_gammamultiplier*gJH0_adm/1.0e-12, nHcgs, logT);} // if < 0, that's a key to tell us this needs to be recalculated
    n_elec = *ne_guess; if(!isfinite(n_elec)) {n_elec=1;}
    neold = n_elec; niter = 0;
    double dt = 0, fac_noneq_cgs = 0, necgs = n_elec * nHcgs; /* more initialized quantities */
    if(target >= 0) {dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);} // dtime [code units]
    fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * necgs; // factor needed below to asses whether timestep is larger/smaller than recombination time

#if defined(RT_CHEM_PHOTOION)
    double c_light_ne=0;
    if(target >= 0)
    {
        nH0 = SphP[target].HI; // need to initialize a value for the iteration below
#ifdef RT_CHEM_PHOTOION_HE
        nHe0 = SphP[target].HeI; nHep = SphP[target].HeII; // need to intialize a value for the iteration below
#endif
    }
#endif

    /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
    do
    {
        niter++;

        aHp = flow * AlphaHp_adm[j] + fhi * AlphaHp_adm[j + 1];
        aHpRate = flow*AlphaHpRate_adm[j] + fhi*AlphaHpRate_adm[j+1];
	aHep = flow * AlphaHep_adm[j] + fhi * AlphaHep_adm[j + 1];
        aHepp = flow * AlphaHepp_adm[j] + fhi * AlphaHepp_adm[j + 1];
        ad = flow * Alphad_adm[j] + fhi * Alphad_adm[j + 1];
        geH0 = flow * GammaeH0_adm[j] + fhi * GammaeH0_adm[j + 1];
        geH0 = DMAX(geH0, EPSILON_SMALL);
        geHe0 = flow * GammaeHe0_adm[j] + fhi * GammaeHe0_adm[j + 1];
        geHe0 = DMAX(geHe0, EPSILON_SMALL);
        geHep = flow * GammaeHep_adm[j] + fhi * GammaeHep_adm[j + 1];
        geHep = DMAX(geHep, EPSILON_SMALL);
        fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * necgs; // factor needed below to asses whether timestep is larger/smaller than recombination time
        if(necgs <= 1.e-25 || J_UV_adm == 0)
        {
            gJH0ne = gJHe0ne = gJHepne = 0;
        }
        else
        {
            /* account for self-shielding in calculating UV background effects */
            gJH0ne = gJH0_adm * local_gammamultiplier / necgs * shieldfac; // check units, should be = c_light * n_photons_vol * rt_ion_sigma_HI[0] / necgs;
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = gJHe0_adm * local_gammamultiplier / necgs * shieldfac;
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = gJHep_adm * local_gammamultiplier / necgs * shieldfac;
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
        }
#if defined(RT_DISABLE_UV_BACKGROUND)
        gJH0ne = gJHe0ne = gJHepne = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if(target >= 0)
        {
            int k;
            c_light_ne = C_LIGHT / ((MIN_REAL_NUMBER + necgs) * UNIT_LENGTH_IN_CGS); // want physical cgs units for quantities below
            double gJH0ne_0=gJH0_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHe0ne_0=gJHe0_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHepne_0=gJHep_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs); // need a baseline, so we don't over-shoot below
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
#if defined(RT_DISABLE_UV_BACKGROUND)
            gJH0ne_0=gJHe0ne_0=gJHepne_0=MAX_REAL_NUMBER;
#endif
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if((k==RT_FREQ_BIN_H0)||(k==RT_FREQ_BIN_He0)||(k==RT_FREQ_BIN_He1)||(k==RT_FREQ_BIN_He2))
                {
                    double c_ne_time_n_photons_vol = c_light_ne * rt_return_photon_number_density(target,k); // gives photon flux
                    double cross_section_ion, dummy, thold=1.0e20;
#ifdef GALSF
                    if(All.ComovingIntegrationOn) {thold=1.0e10;}
#endif
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        dummy = rt_ion_sigma_HI[k] * c_ne_time_n_photons_vol;// egy per photon x cross section x photon flux (w attenuation factors already included in flux/energy update:) * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        if(dummy > thold*gJH0ne_0) {dummy = thold*gJH0ne_0;}
                        gJH0ne += dummy;
                    }
#ifdef RT_CHEM_PHOTOION_HE
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        dummy = rt_ion_sigma_HeI[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        if(dummy > thold*gJHe0ne_0) {dummy = thold*gJHe0ne_0;}
                        gJHe0ne += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        dummy = rt_ion_sigma_HeII[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        if(dummy > thold*gJHepne_0) {dummy = thold*gJHepne_0;}
                        gJHepne += dummy;
                    }
#endif
                }
            }
        }
#endif


        nH0 = aHp / (MIN_REAL_NUMBER + aHp + geH0 + gJH0ne);	/* eqn (33) */
#ifdef RT_CHEM_PHOTOION
        if(target >= 0) {nH0 = (SphP[target].HI + fac_noneq_cgs * aHp) / (1 + fac_noneq_cgs * (aHp + geH0 + gJH0ne));} // slightly more general formulation that gives linear update but interpolates to equilibrium solution when dt >> dt_recombination
#endif
        nHp = 1.0 - nH0;		/* eqn (34) */

        if( ((gJHe0ne + geHe0) <= MIN_REAL_NUMBER) || (aHepp <= MIN_REAL_NUMBER) ) 	/* no ionization at all */
        {
            nHep = 0.0;
            nHepp = 0.0;
            nHe0 = yhelium(target);
        }
        else
        {
            nHep = yhelium(target) / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
            nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
            nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
        }
#if defined(RT_CHEM_PHOTOION) && defined(RT_CHEM_PHOTOION_HE)
        if(target >= 0)
        {
            double yHe = yhelium(target); // will use helium fraction below
            nHep = SphP[target].HeII + yHe * fac_noneq_cgs * (geHe0 + gJHe0ne) - SphP[target].HeIII * (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp));
            nHep /= 1.0 + fac_noneq_cgs*(geHe0 + gJHe0ne + aHep + ad + geHep + gJHepne) + (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp)) * fac_noneq_cgs*(geHep + gJHepne);
            if(nHep < 0) {nHep=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHep > yHe) {nHep=yHe;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHepp = (SphP[target].HeIII + SphP[target].HeII * fac_noneq_cgs*(geHep + gJHepne)) / (1. + fac_noneq_cgs*aHepp);
            if(nHepp < 0) {nHepp=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHepp > yHe-nHep) {nHepp=yHe-nHep;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHe0 = yHe - (nHep + nHepp); // remainder is neutral
        }
#endif
        if(!isfinite(n_elec)) {printf("target=%d niter=%d logT=%g n_elec/old=%g/%g nHp/nHep/nHepp=%g/%g/%g nHcgs=%g yHe=%g dt=%g shieldfac/local_gammamult=%g/%g aHp/aHep/aHepp=%g/%g/%g geH0/geHe0/geHep=%g/%g/%g gJH0ne/gJHe0ne/gJHepne=%g/%g/%g \n",target,niter,logT,n_elec,neold,nHp,nHep,nHepp,nHcgs,yhelium(target),dt,shieldfac,local_gammamultiplier,aHp,aHep,aHepp,geH0,geHe0,geHep,gJH0ne,gJHe0ne,gJHepne);}

        neold = n_elec;
        n_elec = nHp + nHep + 2 * nHepp;	/* eqn (38) */

/* ////////////////////////////////////// 
 * COMMENTED OUT COOL_LOW_TEMPERATURES //
 *///////////////////////////////////////

/* #ifdef COOL_LOW_TEMPERATURES
        n_elec += return_electron_fraction_from_heavy_ions_adm(target, pow(10.,logT), rho, n_elec);
#endif */

        necgs = n_elec * nHcgs;

        if(J_UV_adm == 0) break;

        nenew = 0.5 * (n_elec + neold);
        n_elec = nenew;
        if(!isfinite(n_elec)) {n_elec=1;}
        necgs = n_elec * nHcgs;

        double dneTHhold = DMAX(n_elec*0.01 , 1.0e-4);
        if(fabs(n_elec - neold) < dneTHhold) break;

        if(niter > (MAXITER - 10)) {printf("n_elec= %g/%g/%g yh=%g nHcgs=%g niter=%d\n", n_elec,neold,nenew, yhelium(target), nHcgs, niter);}
    }
    while(niter < MAXITER);

    if(niter >= MAXITER) {printf("failed to converge in find_abundances_and_rates_adm(): logT_input=%g  rho_input=%g  ne_input=%g target=%d shieldfac=%g cooling_return=%d", logT_input, rho_input, ne_input, target, shieldfac, return_cooling_mode); endrun(13);}

    bH0 = flow * BetaH0_adm[j] + fhi * BetaH0_adm[j + 1];
    bHep = flow * BetaHep_adm[j] + fhi * BetaHep_adm[j + 1];
    bff = flow * Betaff_adm[j] + fhi * Betaff_adm[j + 1];
    *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec; /* write to send back */
    *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, sqrt(shieldfac)*(gJH0_adm/2.29e-10), target);
    if(target >= 0) /* if this is a cell, update some of its thermodynamic stored quantities */
    {
        SphP[target].Ne = n_elec;
#if defined(OUTPUT_MOLECULAR_FRACTION)
        SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, pow(10.,logT), nH0, n_elec, sqrt(shieldfac)*(gJH0_adm/2.29e-10));
#endif
    }

    /* now check if we want to return the ionization/recombination heating/cooling rates calculated with all the above quantities */
    if(return_cooling_mode==1)
    {
        /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        double LambdaExcH0 = bH0 * n_elec * nH0;
        double LambdaExcHep = bHep * n_elec * nHep;
        double LambdaExc =  LambdaExcH0 + LambdaExcHep;	/* collisional excitation */

        //double LambdaIonH0 = 2.18e-11 * geH0 * n_elec * nH0; // GIZMO Rate
        double LambdaIonH0 = 0.5 * All.ADM_ElectronMass * pow(C_LIGHT,2.0) * pow(All.ADM_FineStructure,2.0) * geH0 * n_elec * nH0; // ADM Rate
	double LambdaIonHe0 = 3.94e-11 * geHe0 * n_elec * nHe0;
        double LambdaIonHep = 8.72e-11 * geHep * n_elec * nHep;
        double LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* collisional ionization */

        double T_lin = pow(10.0, logT);
        //double LambdaRecHp = 1.036e-16 * T_lin * n_elec * (aHp * nHp); // GIZMO rate
        double LambdaRecHp = aHpRate * n_elec * nHp; // ADM Rate
	double LambdaRecHep = 1.036e-16 * T_lin * n_elec * (aHep * nHep);
        double LambdaRecHepp = 1.036e-16 * T_lin * n_elec * (aHepp * nHepp);
        double LambdaRecHepd = 6.526e-11 * ad * n_elec * nHep;
        double LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd; /* recombination */

        double LambdaFF = bff * (nHp + nHep + 4 * nHepp) * n_elec; /* free-free (Bremsstrahlung) */	

        double Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF; /* sum all of the above */
        return Lambda; /* send it back */
    }
    return 0;
} // end of find_abundances_and_rates_adm() //



/*  this function first computes the self-consistent temperature and abundance ratios, and then it calculates (heating rate-cooling rate)/n_h^2 in cgs units */
double CoolingRateFromU_adm(double u, double rho, double ne_guess, int target)
{
    double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu; nH0_guess = DMAX(0,DMIN(1,1.-ne_guess/1.2));
    double temp = convert_u_to_temp_adm(u, rho, target, &ne_guess, &nH0_guess, &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess, &mu);
    return CoolingRate_adm(log10(temp), rho, ne_guess, target);
}


#endif // !(CHIMES)



extern FILE *fd;



#ifndef CHIMES
////////////////////////////////////////////
// DISABLED COMPTON HEATING. CHANGE THIS! //
// /////////////////////////////////////////


/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units
 */
double CoolingRate_adm(double logT, double rho, double n_elec_guess, int target)
{
    printf("ADM Alert! cooling_adm, coolingrate\n");
    double n_elec=n_elec_guess, nH0, nHe0, nHp, nHep, nHepp, mu; /* ionization states [computed below] */
    double Lambda, Heat, LambdaFF, LambdaCompton, LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
    double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd, T, shieldfac, LambdaMol, LambdaMetal;
    double nHcgs = HYDROGEN_MASSFRAC_ADM * rho / All.ADM_ProtonMass;	/* hydrogen number dens in cgs units */
    LambdaMol=0; LambdaMetal=0; LambdaCompton=0;
    if(logT <= Tmin_adm) {logT = Tmin_adm + 0.5 * deltaT_adm;}	/* floor at Tmin_adm */
    if(!isfinite(rho)) {return 0;}
    T = pow(10.0, logT);

    /* some blocks below to define useful variables before calculation of cooling rates: */

#ifdef COOL_METAL_LINES_BY_SPECIES
    double *Z;
    if(target>=0)
    {
        Z = P[target].Metallicity;
    } else { /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
        int k; double Zsol[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) {Zsol[k]=All.SolarAbundances[k];}
        Z = Zsol;
    }
#endif
    double local_gammamultiplier = return_local_gammamultiplier_adm(target);
    shieldfac = return_uvb_shieldfac_adm(target, local_gammamultiplier*gJH0_adm/1.0e-12, nHcgs, logT);

////////////////////////////////////////
// COMMENTED OUT COOL_LOW_TEMPERATURE //
// /////////////////////////////////////

/*
#if defined(COOL_LOW_TEMPERATURES)
    double Tdust = 30., LambdaDust = 0.;*/ /* set variables needed for dust heating/cooling. if dust cooling not calculated, default to 0 */
/*#if (defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate_adm(target, shieldfac);
#endif
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(SINGLE_STAR_SINK_DYNAMICS) && !defined(SINGLE_STAR_FB_RT_HEATING)
    Tdust = DMIN(DMAX(10., 2.73/All.cf_atime),300.); // runs looking at colder clouds, use a colder default dust temp [floored at CMB temperature] //
#endif
#endif
*/

#if defined(RT_CHEM_PHOTOION) || defined(RT_PHOTOELECTRIC)
    double cx_to_kappa = HYDROGEN_MASSFRAC / PROTONMASS * UNIT_MASS_IN_CGS; // pre-factor for converting cross sections into opacities
#endif
    if(logT < Tmax_adm)
    {
        /* get ionization states for H and He with associated ionization, collision, recombination, and free-free heating/cooling */
        Lambda = find_abundances_and_rates_adm(logT, rho, target, shieldfac, 1, &n_elec, &nH0, &nHp, &nHe0, &nHep, &nHepp, &mu);

        LambdaCompton = evaluate_Compton_heating_cooling_rate_adm(target,T,nHcgs,n_elec,shieldfac); /* note this can have either sign: heating or cooling */
	if(LambdaCompton > 0) {Lambda += LambdaCompton;}

#ifdef COOL_METAL_LINES_BY_SPECIES
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
//////////////////////
//// STELLAREVOLUTION //
//// ///////////////////
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
//        if(J_UV_adm != 0)
//#else
        if((J_UV_adm != 0)&&(logT > 4.00))
//#endif
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies_adm(nHcgs, logT, Z); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= n_elec;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted; ni not used b/c it should vary species-to-species */
//////////////////////
//// STELLAREVOLUTION //
//// ///////////////////
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
//            if(logT<2) {LambdaMetal *= exp(-DMIN((2.-logT)*(2.-logT)/0.1,40.));}
//            if(LambdaMetal > 0) {Lambda += LambdaMetal;}
//#else
            Lambda += LambdaMetal;
//#endif
#if defined(OUTPUT_COOLRATE_DETAIL)
            if(target >= 0) {SphP[target].MetalCoolingRate = LambdaMetal;}
#endif
        }
#endif

/////////////////////////////////
// COMMENTED OUT COOL_LOW_TEMP //
// //////////////////////////////

/*
#ifdef COOL_LOW_TEMPERATURES
        if(logT <= 5.3)
        {
*/            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
/*            LambdaMol = 2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));
            LambdaMol *= (1-shieldfac) / (1. + nHcgs/700.); // above the critical density, cooling rate suppressed by ~1/n; use critical density of CO[J(1-0)] as a proxy for this
            double Z_sol=1, truncation_factor=1;*/ /* if don't have actual metallicities, we'll assume solar */
/*            if(logT>4.5) {double dx=(logT-4.5)/0.20; truncation_factor *= exp(-DMIN(dx*dx,40.));}*/ /* continuous cutoff here just to avoid introducing artificial features in temperature-density */
/*#ifdef COOL_METAL_LINES_BY_SPECIES
            Z_sol = Z[0] / All.SolarAbundances[0];*/  /* use actual metallicity for this */
/*#endif
            LambdaMol *= (1+Z_sol)*(0.001 + 0.1*nHcgs/(1.+nHcgs) + 0.09*nHcgs/(1.+0.1*nHcgs) + Z_sol*Z_sol/(1.0+nHcgs)); // gives very crude estimate of metal-dependent terms //
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
            double column = evaluate_NH_from_GradRho(P[target].GradRho,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs
            double Z_C = DMAX(1.e-6, Z[2]/All.SolarAbundances[2]), sqrt_T=sqrt(T), ncrit_CO=1.9e4*sqrt_T, Sigma_crit_CO=3.0e-5*T/Z_C, T3=T/1.e3, EXPmax=90.; // carbon abundance (relative to solar), critical density and column
            double f_Cplus_CCO=1./(1.+nHcgs/3.e3), photoelec=0; // very crude estimate used to transition between C+ cooling curve and C/CO [nearly-identical] cooling curves above C+ critical density, where C+ rate rapidly declines
#ifdef GALSF_FB_FIRE_RT_UVHEATING
            photoelec = SphP[target].Rad_Flux_UV; if(gJH0_adm>0 && shieldfac>0) {photoelec+=sqrt(shieldfac) * (gJH0_adm/2.29e-10);} // uvb contribution //
#endif
#ifdef RT_PHOTOELECTRIC
            photoelec = SphP[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; photoelec=DMAX(DMIN(photoelec,1.e4),0); // convert to Habing field //
#endif
#ifdef RT_ISRF_BACKGROUND
            photoelec += RT_ISRF_BACKGROUND * 1.7 * exp(-DMAX(P[target].Metallicity[0]/All.SolarAbundances[0],1e-4) * column * 500.); // RT_ISRF_BACKGROUND rescales the overal ISRF, factor of 1.7 gives Draine 1978 field in Habing units, extinction factor assumes the same FUV band-integrated dust opacity as RT module
#endif
#if !(defined(GALSF_FB_RT_UVHEATING) || defined(RT_PHOTOELECTRIC) || defined(RT_ISRF_BACKGROUND))
            photoelec = 1; // if no explicit modeling of FUV, just assume Habing
#endif
            f_Cplus_CCO = (nHcgs/(340.*DMAX(0.1,photoelec))); f_Cplus_CCO=1./(1.+f_Cplus_CCO*f_Cplus_CCO/sqrt_T); // fco/(1-fco) ~ 0.0022 * ((n/50 cm^-3)/G0)^2 * (100K/T)^(1/2) from Tielens
            double Lambda_Cplus = Z_C * (4.7e-28 * (pow(T,0.15) + 1.04e4*n_elec/sqrt_T) * exp(-DMIN(91.211/T,EXPmax)) + 2.08e-29*exp(-DMIN(23.6/T,EXPmax))); // fit from Barinovs et al., ApJ, 620, 537, 2005, and Wilson & Bell MNRAS 337 1027 2002; assuming factor of 0.5 depletion factor in ISM; rate per C+ relative to solar; + plus [CI]-609 µm line cooling from Hocuk⋆ et al. 2016MNRAS.456.2586H
            double Lambda_CCO = Z_C * T*sqrt_T * 2.73e-31 / (1. + (nHcgs/ncrit_CO)*(1.+1.*DMAX(column,0.017)/Sigma_crit_CO)); // fit from Hollenbach & McKee 1979 for CO (+CH/OH/HCN/OH/HCl/H20/etc., but those don't matter), with slight re-calibration of normalization (factor ~1.4 or so) to better fit the results from the full Glover+Clark network. As Glover+Clark show, if you shift gas out of CO into C+ and O, you have almost no effect on the integrated cooling rate, so this is a surprisingly good approximation without knowing anything about the detailed chemical/molecular state of the gas. uncertainties in e.g. ambient radiation are -much- larger. also note this rate is really carbon-dominated as the limiting abundance, so should probably use that.
            double Lambda_Metals = f_Cplus_CCO * Lambda_Cplus + (1.-f_Cplus_CCO) * Lambda_CCO; // interpolate between both regimes //
          */  /* in the above Lambda_Metals expression, the column density expression attempts to account for the optically-thick correction in a slab. this is largely redundant (not exactly, b/c this is specific for CO-type molecules) with our optically-thick cooling module already included below, so we will not double-count it here [coefficient set to zero]. But it's included so you can easily turn it back on, if desired, instead of using the module below. */
/*            double Lambda_H2_thick = (6.7e-19*exp(-DMIN(5.86/T3,EXPmax)) + 1.6e-18*exp(-DMIN(11.7/T3,EXPmax)) + 3.e-24*exp(-DMIN(0.51/T3,EXPmax)) + 9.5e-22*pow(T3,3.76)*exp(-DMIN(0.0022/(T3*T3*T3),EXPmax))/(1.+0.12*pow(T3,2.1))) / nHcgs; // super-critical H2-H cooling rate [per H2 molecule]
            double Lambda_HD_thin = ((1.555e-25 + 1.272e-26*pow(T,0.77))*exp(-DMIN(128./T,EXPmax)) + (2.406e-25 + 1.232e-26*pow(T,0.92))*exp(-DMIN(255./T,EXPmax))) * exp(-DMIN(T3*T3/25.,EXPmax)); // optically-thin HD cooling rate [assuming all D locked into HD at temperatures where this is relevant], per molecule
            double f_molec = 0.5 * Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(gJH0_adm/2.29e-10)); // [0.5*f_molec for H2/HD cooling b/c cooling rates above are per molecule, not per nucleon]

            double q = logT - 3., Y_Hefrac=DMAX(0.,DMIN(1.,Z[1])), X_Hfrac=DMAX(0.,DMIN(1.,1.-Y_Hefrac-Z[0])); // variable used below
            double Lambda_H2_thin = DMAX(nH0-2.*f_molec,0) * X_Hfrac * pow(10., DMAX(-103. + 97.59*logT - 48.05*logT*logT + 10.8*logT*logT*logT - 0.9032*logT*logT*logT*logT , -50.)); // sub-critical H2 cooling rate from H2-H collisions [per H2 molecule]; this from Galli & Palla 1998
            Lambda_H2_thin += Y_Hefrac * pow(10., DMAX(-23.6892 + 2.18924*q -0.815204*q*q + 0.290363*q*q*q -0.165962*q*q*q*q + 0.191914*q*q*q*q*q, -50.)); // H2-He; often more efficient than H2-H at very low temperatures (<100 K); this and other H2-x terms below from Glover & Abel 2008
            Lambda_H2_thin += f_molec * X_Hfrac * pow(10., DMAX(-23.9621 + 2.09434*q -0.771514*q*q + 0.436934*q*q*q -0.149132*q*q*q*q -0.0336383*q*q*q*q*q, -50.)); // H2-H2; can be more efficient than H2-H when H2 fraction is order-unity
            Lambda_H2_thin += nHp * X_Hfrac * pow(10., DMAX(-21.7167 + 1.38658*q -0.379153*q*q + 0.114537*q*q*q -0.232142*q*q*q*q + 0.0585389*q*q*q*q*q, -50.)); // H2-H+; very efficient if somehow appreciable H+ fraction remains
            double logLambdaH2_e = -34.2862 -48.5372*q -77.1212*q*q -51.3525*q*q*q -15.1692*q*q*q*q -0.981203*q*q*q*q*q; // H2-e [generally sub-dominant to H2-H+; can dominate if free e- largely from other sources (e.g. Mg, etc.), but in those conditions essentially impossible for H2 cooling to dominate
            if(logT>2.30103) {logLambdaH2_e = -22.1903 + 1.5729*q -0.213351*q*q + 0.961498*q*q*q -0.910232*q*q*q*q + 0.137497*q*q*q*q*q;}
            Lambda_H2_thin += n_elec * X_Hfrac * pow(10., DMAX(logLambdaH2_e, -50.));

            double f_HD = DMIN(0.00126*f_molec , 4.0e-5*nH0); // ratio of HD molecules to H2 molecules: in low limit, HD easier to form so saturates at about 0.13% of H2 molecules, following Galli & Palla 1998, but obviously cannot exceed the cosmic ratio of D/H=4e-5
            double nH_over_ncrit = Lambda_H2_thin / Lambda_H2_thick , Lambda_HD = f_HD * Lambda_HD_thin / (1. + (f_HD/(f_molec+MIN_REAL_NUMBER))*nH_over_ncrit), Lambda_H2 = f_molec * Lambda_H2_thin / (1. + nH_over_ncrit); // correct cooling rates for densities above critical
            double Lambda_Metals_Neutral = nH0 * Lambda_Metals; // finally note our metal terms here are all for atomic or molecular, not ionized (handled in tables above)
            if(!isfinite(Lambda_Metals_Neutral) || Lambda_Metals_Neutral < 0) {Lambda_Metals_Neutral=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_H2) || Lambda_H2 < 0) {Lambda_H2=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_HD) || Lambda_HD < 0) {Lambda_HD=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            LambdaMol = Lambda_Metals_Neutral + Lambda_H2 + Lambda_HD; // combine to get total cooling rate
#endif
            LambdaMol *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            if(!isfinite(LambdaMol)) {LambdaMol=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            Lambda += LambdaMol;

*/            /* now add the dust cooling/heating terms */
/*            LambdaDust = 1.116e-32 * (Tdust-T) * sqrt(T)*(1.-0.8*exp(-75./T)) * Z_sol;  // Meijerink & Spaans 2005; Hollenbach & McKee 1979,1989 //
#ifdef RT_INFRARED
            if(target >= 0) {LambdaDust = get_rt_ir_lambdadust_effective(T, rho, &nH0, &n_elec, target);} // call our specialized subroutine, because radiation and gas energy fields are co-evolving and tightly-coupled here //
#endif
            if(T>3.e5) {double dx=(T-3.e5)/2.e5; LambdaDust *= exp(-DMIN(dx*dx,40.));} *//* needs to truncate at high temperatures b/c of dust destruction */
/*            LambdaDust *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
#ifdef RT_INFRARED
            SphP[target].LambdaDust = LambdaDust;
#endif
            if(!isfinite(LambdaDust)) {LambdaDust=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(LambdaDust<0) {Lambda -= LambdaDust;} *//* add the -positive- Lambda-dust associated with cooling */
/*        }
#endif
*/


        Heat = 0;  /* Now, collect heating terms */

//////////////////////
//// STELLAREVOLUTION //
//// ///////////////////
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
        // here we account for the fact that the local spectrum is softer than the UVB which includes AGN and is hardened by absorption within galaxies. we do this by simply lowering the effective heating rate [mean photon energy absorbed per ionization], which captures the leading-order effect //
//        if(J_UV_adm != 0) {Heat += (1. + (local_gammamultiplier-1.) * 0.33) * (nH0 * epsH0_adm + nHe0 * epsHe0_adm + nHep * epsHep_adm) / nHcgs * shieldfac;} // shieldfac allows for self-shielding from background
//#else
        if(J_UV_adm != 0) {Heat += local_gammamultiplier * (nH0 * epsH0_adm + nHe0 * epsHe0_adm + nHep * epsHep_adm) / nHcgs * shieldfac;} // shieldfac allows for self-shielding from background
//#endif
#if defined(RT_DISABLE_UV_BACKGROUND)
        Heat = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if((target >= 0) && (nHcgs > MIN_REAL_NUMBER))
        {
            int k; double c_light_nH = C_LIGHT / (nHcgs * UNIT_LENGTH_IN_CGS) * UNIT_ENERGY_IN_CGS; // want physical cgs units for quantities below
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if((k==RT_FREQ_BIN_H0)||(k==RT_FREQ_BIN_He0)||(k==RT_FREQ_BIN_He1)||(k==RT_FREQ_BIN_He2))
                {
                    double c_nH_time_n_photons_vol = c_light_nH * rt_return_photon_number_density(target,k); // gives photon flux
                    double cross_section_ion, kappa_ion, dummy;
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HI[k] * cross_section_ion * c_nH_time_n_photons_vol;// (egy per photon x cross section x photon flux) :: attenuation factors [already in flux/energy update]: * slab_averaging_function(kappa_ion * Sigma_particle); // egy per photon x cross section x photon flux (w attenuation factors) // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeI[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion * Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeII[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion*Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                }
            }
        }
#endif


#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_ALT_DISABLE_LOSSES)
        Heat += CR_gas_heating(target, n_elec, nHcgs);
#else
///////////////////////////////////////
// COMMENTED OUT COOL_LOW_TEMPERATURE// 
// ////////////////////////////////////

/*#ifdef COOL_LOW_TEMPERATURES
*/        /* if COSMIC_RAYS is not enabled, but low-temperature cooling is on, we account for the CRs as a heating source using
         a more approximate expression (assuming the mean background of the Milky Way clouds) */
/*        if(logT <= 5.2)
        {
            double prefac_CR=1.; if(All.ComovingIntegrationOn) {
                double rhofac = rho / (1000.*COSMIC_BARYON_DENSITY_CGS);
                if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
            double cr_zeta=1.e-16, e_per_cr_ioniz=8.8e-12; // either high background (zeta=1e-16), with softer spectrum (5.5eV per ionization), following e.g. van Dishoeck & Black (1986); or equivalently lower rate with higher ~20eV per ionization per Goldsmith & Langer (1978); this is formally degenerate here. however both produce ~3-10x higher rates than more modern estimates (basically in both cases, assuming a CR energy density of ~2-5 eV/cm^3, instead of more modern ~0.5-2 eV/cm^3
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
            cr_zeta=1.e-17; e_per_cr_ioniz=3.0e-11; // follow e.g. Glover+Jappsen 2007, Le Teuff et al. 2000, gives about a 3x lower CR heating rate compared to older numbers above
#endif
            Heat += prefac_CR * cr_zeta * (1. + 1.68*n_elec*HYDROGEN_MASSFRAC) / (1.e-2 + nHcgs) * e_per_cr_ioniz; // final result
        }
#endif */
#endif

/*#if defined(COOL_LOW_TEMPERATURES)
        if(LambdaDust>0) {Heat += LambdaDust;}*/ /* Dust collisional heating (Tdust > Tgas) */
/*#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2)
        if(LambdaMetal<0) {Heat -= LambdaMetal;} // potential net heating from low-temperature gas-phase metal line absorption //
#endif
#endif
*/
        if(LambdaCompton<0) {Heat -= LambdaCompton;} /* Compton heating rather than cooling */

#if defined(GALSF_FB_FIRE_RT_UVHEATING) || defined(RT_PHOTOELECTRIC) || defined(RT_ISRF_BACKGROUND)
        /* Photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995); now with 'update' from Wolfire 2005 for PAH [fudge factor 0.5 below] */
        if((target >= 0) && (T < 1.0e6))
        {
            double photoelec = 0;
#ifdef GALSF_FB_FIRE_RT_UVHEATING
            photoelec += SphP[target].Rad_Flux_UV;
	    //printf("ADM UV Flux: %.4e, %.4e\n",SphP[target].Rad_Flux_UV,SphP[target].Rad_Flux_EUV);
#ifdef COOL_UVB_SELFSHIELD_RAHMATI
            if(gJH0_adm>0 && shieldfac>0) {photoelec += sqrt(shieldfac) * (gJH0_adm / 2.29e-10);} // uvb contribution //
#endif
#endif
#ifdef RT_PHOTOELECTRIC
            photoelec += SphP[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
#ifdef RT_ISRF_BACKGROUND // add a constant assumed FUV background, for isolated ISM simulations that don't get FUV from local sources self-consistently
            double column = evaluate_NH_from_GradRho(P[target].GradRho,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs
            photoelec += RT_ISRF_BACKGROUND * 1.7 * exp(-DMAX(P[target].Metallicity[0]/All.SolarAbundances[0],1e-4) * column * 500.); // RT_ISRF_BACKGROUND rescales the overal ISRF, factor of 1.7 gives Draine 1978 field in Habing units, extinction factor assumes the same FUV band-integrated dust opacity as RT module
#endif
            if(photoelec > 0) {if(photoelec > 1.e4) {photoelec = 1.e4;}}

            if(photoelec > 0)
            {
                double LambdaPElec = 1.3e-24 * photoelec / nHcgs * (P[target].Metallicity[0]/All.SolarAbundances[0]);
                double x_photoelec = photoelec * sqrt(T) / (0.5 * (1.0e-12+n_elec) * nHcgs);
                LambdaPElec *= 0.049/(1+pow(x_photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+x_photoelec/5000.);
                Heat += LambdaPElec;
            }
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>Tmax_adm K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.  Assumes no heating. */
      Heat = LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;
      nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp; /* very hot: H and He both fully ionized */

      //LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * n_elec; // free-free GIZMO Rate
      LambdaFF = 3.7e-27 * (pow(All.ADM_FineStructure/0.01, 3.0) / pow(All.ADM_ElectronMass/ELECTRONMASS, 1.5)) * sqrt(T) * (nHp + 4*nHepp)*n_elec; // ADM Rate

      LambdaCompton = evaluate_Compton_heating_cooling_rate_adm(target,T,nHcgs,n_elec,shieldfac); // Compton

      Lambda = LambdaFF + DMAX(LambdaCompton,0);
    }

    double Q = Heat - Lambda;
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0) {SphP[target].CoolingRate = Lambda; SphP[target].HeatingRate = Heat;}
#endif
/////////////////////////////////
// COMMENTED OUT COOL_LOW_TEMP //
// //////////////////////////////

//#if defined(COOL_LOW_TEMPERATURES) && !defined(COOL_LOWTEMP_THIN_ONLY)
    /* if we are in the optically thick limit, we need to modify the cooling/heating rates according to the appropriate limits;
        this flag does so by using a simple approximation. we consider the element as if it were a slab, with a column density
        calculated from the simulation properties and the Sobolev approximation. we then assume it develops an equilibrium internal
        temperature structure on a radiative diffusion timescale much faster than the dynamical time, and so the surface radiation
        from a photosphere can be simply related to the local density by the optical depth to infinity. the equations here follow
        Rafikov, 2007 (ApJ, 662, 642):
            denergy/dt/dArea = sigma*T^4 / fc(tau)
            fc(tau) = tau^eta + 1/tau (taking chi, phi~1; the second term describes the optically thin limit, which is calculated above
                more accurately anyways - that was just Kirchoff's Law; so we only need to worry about the first term)
            eta = 4*(gamma-1) / [gamma*(1+alpha+beta*(gamma-1)/gamma)], where gamma=real polytropic index, and alpha/beta follow
                an opacity law kappa=kappa_0 * P^alpha * T^beta. for almost all the regimes of interest, however, eta~1, which is also
                what is obtained for a convectively stable slab. so we will use this.
            now, this gives sigma*T^4/tau * Area_eff / nHcgs as the 'effective' cooling rate in our units of Heat or Lambda above.
                the nHcgs just puts it in the same volumetric terms. The Area_eff must be defined as ~m_particle/surface_density
                to have the same meaning for a slab as assumed in Rafikov (and to integrate correctly over all particles in the slab,
                if/when the slab is resolved). We estimate this in our usual fashion with the Sobolev-type column density
            tau = kappa * surface_density; we estimate kappa ~ 5 cm^2/g * (0.001+Z/Z_solar), as the frequency-integrated kappa for warm
                dust radiation (~150K), weighted by the dust-to-gas ratio (with a floor for molecular absorption). we could make this
                temperature-dependent, though, fairly easily - for this particular problem it won't make much difference
        This rate then acts as an upper limit to the net heating/cooling calculated above (restricts absolute value)
     */
/*    if( (nHcgs > 0.1) && (target >= 0) ) */ /* don't bother at very low densities, since youre not optically thick, and protect from target=-1 with GALSF_EFFECTIVE_EQS */
/*    {
        double surface_density = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target);
        surface_density *= 0.2 * UNIT_SURFDEN_IN_CGS; // converts to cgs; 0.2 is a tuning factor so that the Masunaga & Inutsuka 2000 solution is reproduced
        double effective_area = 2.3 * PROTONMASS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
        if(T < 1500.)
        {
            if(T < 150.) {kappa_eff=0.0027*T*sqrt(T);} else {kappa_eff=5.;}
            kappa_eff *= P[target].Metallicity[0]/All.SolarAbundances[0];
            if(kappa_eff < 0.1) {kappa_eff=0.1;}
        } else {
*/            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
//            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * n_elec; /* Thompson scattering (non-relativistic) */
//            double k_molecular = 0.1 * P[target].Metallicity[0]; /* molecular line opacities */
//            double k_Hminus = 1.1e-25 * sqrt(P[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
//            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
//            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
//            double k_conductive = 2.6e-7 * n_elec * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
//            kappa_eff = 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
/*        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif
*/

#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0){SphP[target].NetHeatingRateQ = Q;}
#endif
#ifdef OUTPUT_MOLECULAR_FRACTION
    if(target>0) {SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(gJH0_adm/2.29e-10));}
#endif

#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) {Q += SphP[target].DtInternalEnergy / nHcgs;}
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target >= 0) {SphP[target].HydroHeatingRate = SphP[target].DtInternalEnergy / nHcgs;}
#endif
#endif

  return Q;
} // ends CoolingRate_adm





void InitCoolMemory_adm(void)
{
    //printf("ADM Alert! cooling_adm, initcoolmemory\n");
    BetaH0_adm = (double *) mymalloc("BetaH0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    BetaHep_adm = (double *) mymalloc("BetaHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHp_adm = (double *) mymalloc("AlphaHp_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHpRate_adm = (double *) mymalloc("AlphaHpRate_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHep_adm = (double *) mymalloc("AlphaHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    Alphad_adm = (double *) mymalloc("Alphad_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHepp_adm = (double *) mymalloc("AlphaHepp_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeH0_adm = (double *) mymalloc("GammaeH0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeHe0_adm = (double *) mymalloc("GammaeHe0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeHep_adm = (double *) mymalloc("GammaeHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    Betaff_adm = (double *) mymalloc("Betaff_adm", (NCOOLTAB_ADM + 1) * sizeof(double));

#ifdef COOL_METAL_LINES_BY_SPECIES
    long i_nH=41; long i_T=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES;
    SpCoolTable0_adm = (float *) mymalloc("SpCoolTable0_adm",(kspecies*i_nH*i_T)*sizeof(float));
    if(All.ComovingIntegrationOn) {SpCoolTable1_adm = (float *) mymalloc("SpCoolTable1_adm",(kspecies*i_nH*i_T)*sizeof(float));}
#endif
}



void MakeCoolingTable_adm(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    //printf("ADM Alert! cooling_adm, make cooling table\n");
    int i; double T,Tfact;
    if(All.MinGasTemp > 0.0) {Tmin_adm = log10(All.MinGasTemp);} else {Tmin_adm=-1.0;} // set minimum temperature in this table to some very low value if zero, where none of the cooling approximations above make sense
    deltaT_adm = (Tmax_adm - Tmin_adm) / NCOOLTAB_ADM;
    double y2;    
    double gtoGeV = 1.0 / (1.79e-24);
    double TtoGeV = 1.0 / (1.16e13);
    double InGeVcm = 1.9733e-14;
    double GeVtoergs = 0.001602;    

    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB_ADM; i++)
    {
        BetaH0_adm[i] = BetaHep_adm[i] = Betaff_adm[i] = AlphaHp_adm[i] = AlphaHpRate_adm[i] = AlphaHep_adm[i] = AlphaHepp_adm[i] = Alphad_adm[i] = GammaeH0_adm[i] = GammaeHe0_adm[i] = GammaeHep_adm[i] = 0;
        T = pow(10.0, Tmin_adm + deltaT_adm * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
	y2 = All.ADM_ElectronMass * pow(All.ADM_FineStructure, 2.0) * pow(C_LIGHT,2.0) / (2 * BOLTZMANN * T);
	BetaH0_adm[i] = 7.4e-18 * pow(All.ADM_FineStructure/0.01,2.0) * sqrt(ELECTRONMASS / All.ADM_ElectronMass) *  sqrt(1.0e5/T) * g_integral(y2);
	//printf("BetaH0_i = %.10e\n", BetaH0_adm[i]);
	//if(118348. / T < 70.) {BetaH0_adm[i] = 7.5e-19 * exp(-118348 / T) * Tfact;}
        
	if(473638. / T < 70.) {BetaHep_adm[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;} // UNCOMMENT LATER

        //Betaff_adm[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3)); //UNCOMMENT LATER
        Betaff_adm[i] = 3.7e-27 * (pow(All.ADM_FineStructure/0.01, 3.0) / pow(All.ADM_ElectronMass/ELECTRONMASS, 1.5)) * sqrt(T); 

        //AlphaHp_adm[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */ //UNCOMMENT LATER 
	AlphaHp_adm[i] = (pow(All.ADM_FineStructure, 5.0) / sqrt(pow(T,3.0) * All.ADM_ElectronMass)) * sqrt(pow(2.0, 11.0) * M_PI / 27.0) * (C_LIGHT * pow(InGeVcm, 2.0) / sqrt(gtoGeV * pow(TtoGeV,3.0)))*recomb_rate_integral(y2);
	AlphaHpRate_adm[i] = (pow(All.ADM_FineStructure,5.0)/sqrt(T*All.ADM_ElectronMass))*sqrt(pow(2.0,11.0)*M_PI/27.0)*(C_LIGHT*pow(InGeVcm,2.0)*GeVtoergs/sqrt(gtoGeV*TtoGeV))*recomb_cool_integral(y2);
	AlphaHep_adm[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */ //UNCOMMENT LATER
        AlphaHepp_adm[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*AlphaHp_adm[1,T/Z^2] */ //UNCOMMENT LATER

        
	if(470000.0 / T < 70) {Alphad_adm[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));} // GIZMO rate
        //if(157809.1 / T < 70) {GammaeH0_adm[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;} // GIZMO rate
        GammaeH0_adm[i] = (2.2e-7)*pow(ELECTRONMASS / All.ADM_ElectronMass,1.5)*sqrt(1.0e5 / T) * f_integral(y2);
	if(285335.4 / T < 70) {GammaeHe0_adm[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;} // GIZMO rate
        if(631515.0 / T < 70) {GammaeHep_adm[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;} // GIZMO rate
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables_adm(void)
{
    //printf("ADM Alert! cooling_adm, multispeciestable_adm\n");
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables_adm(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables_adm(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables_adm(0);
    }
}

void ReadMultiSpeciesTables_adm(int iT)
{
    //printf("ADM Alert! cooling_adm, readmultispeciestables\n");
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES; long i,j,k,r;
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;

    fname=GetMultiSpeciesFilename_adm(iT,0);
    if(ThisTask == 0) printf(" ..opening ADM Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read ADM species cooling table in file `%s'\n", fname); endrun(456);}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&SpCoolTable0_adm[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached ADM Cooling EOF! \n");
                }
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename_adm(iT,&fname,1);
     if(!(fdcool = fopen(fname, "r"))) {
     printf(" Cannot read species (He) cooling table in file `%s'\n", fname); endrun(456);}
     for(i=0;i<2;i++)
     for(j=0;j<i_nH;j++)
     for(k=0;k<i_Temp;k++)
     for(l=0;l<i_He;l++)
     fread(&SpCoolTable0_He[i][j][k][l],sizeof(float),1,fdcool);
     fclose(fdcool);
     */
    if (All.ComovingIntegrationOn && i<48) {
        fname=GetMultiSpeciesFilename_adm(iT+1,0);
        if(ThisTask == 0) printf(" ..opening (z+) ADM Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read ADM species 1 cooling table in file `%s'\n", fname); endrun(456);}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&SpCoolTable1_adm[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached ADM Cooling EOF! \n");
                    }
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename_adm(iT+1,&fname,1);
         if(!(fdcool = fopen(fname, "r"))) {
         printf(" Cannot read species 1 (He) cooling table in file `%s'\n", fname); endrun(456);}
         for(i=0;i<2;i++)
         for(j=0;j<i_nH;j++)
         for(k=0;k<i_Temp;k++)
         for(l=0;l<i_He;l++)
         fread(&SpCoolTable1_He[i][j][k][l],sizeof(float),1,fdcool);
         fclose(fdcool);
         */
    }
}

char *GetMultiSpeciesFilename_adm(int i, int hk)
{
    //printf("ADM Alert! cooling_adm, multispeciesfilename\n");
    static char fname[100];
    if(i<0) i=0; if(i>48) i=48;
    if(hk==0) {
        sprintf(fname,"./spcool_tables/spcool_%d",i);
    } else {
        sprintf(fname,"./spcool_tables/spcool_He_%d",i);
    }
    return fname;
}

#endif


////////////////////////////////////////////////
// CHANGE JAMPL_ADM IN FUTURE IF NECESSARY!!! //
// /////////////////////////////////////////////

/* table input (from file TREECOOL) for ionizing parameters */
#define JAMPL_ADM	0.0		/* amplitude factor relative to input table */
#define TABLESIZE_ADM 250		/* Max # of lines in TREECOOL */
static float inlogz_adm[TABLESIZE_ADM];
static double gH0_adm[TABLESIZE_ADM], gHe_adm[TABLESIZE_ADM], gHep_adm[TABLESIZE_ADM]; // upgrade from float to double, should read fine
static double eH0_adm[TABLESIZE_ADM], eHe_adm[TABLESIZE_ADM], eHep_adm[TABLESIZE_ADM]; // upgrade from float to double, should read fine
static int nheattab_adm;		/* length of table */


void ReadIonizeParams_adm(char *fname)
{
    //printf("ADM Alert! cooling_adm, read ionize params\n");
    int i; FILE *fdcool;
    if(!(fdcool = fopen(fname, "r"))) {printf(" Cannot read ionization table in file `%s'. Make sure the correct ADM TREECOOL file is placed in the code run-time directory, and that any leading comments (e.g. lines preceded by ##) are deleted from the file.\n", fname); endrun(456);}
    for(i=0; i<TABLESIZE_ADM; i++) {gH0_adm[i]=0;}
    for(i=0; i<TABLESIZE_ADM; i++) {if(fscanf(fdcool, "%g %lg %lg %lg %lg %lg %lg", &inlogz_adm[i], &gH0_adm[i], &gHe_adm[i], &gHep_adm[i], &eH0_adm[i], &eHe_adm[i], &eHep_adm[i]) == EOF) {break;}}
    fclose(fdcool);
    for(i=0, nheattab_adm=0; i<TABLESIZE_ADM; i++) {if(gH0_adm[i] != 0.0) {nheattab_adm++;} else {break;}} /*  nheattab_adm is the number of entries in the table */
    if(ThisTask == 0) printf(" ..read ADM ionization table [TREECOOL] with %d non-zero UVB entries in file `%s'. Make sure to cite the authors from which the UV background was compiled! (See user guide for the correct references).\n", nheattab_adm, fname);
}


void IonizeParams_adm(void)
{
    IonizeParamsTable_adm();
}



void IonizeParamsTable_adm(void)
{
    //printf("ADM Alert! cooling_adm, ionise params table\n");
    int i, ilow;
    double logz, dzlow, dzhi;
    double redshift;

    if(All.ComovingIntegrationOn)
        {redshift = 1 / All.Time - 1;}
    else
    {
        /* in non-cosmological mode, still use, but adopt z=0 background */
        redshift = 0;
        /*
         gJHe0_adm = gJHep_adm = gJH0_adm = epsHe0_adm = epsHep_adm = epsH0_adm = J_UV_adm = 0;
         return;
         */
    }

    logz = log10(redshift + 1.0);
    ilow = 0;
    for(i=0; i<nheattab_adm; i++) {if(inlogz_adm[i] < logz) {ilow = i;} else {break;}}
    dzlow = logz - inlogz_adm[ilow];
    dzhi = inlogz_adm[ilow + 1] - logz;

    if(logz > inlogz_adm[nheattab_adm - 1] || gH0_adm[ilow] == 0 || gH0_adm[ilow + 1] == 0 || nheattab_adm == 0)
    {
        gJHe0_adm = gJHep_adm = gJH0_adm = 0; epsHe0_adm = epsHep_adm = epsH0_adm = 0; J_UV_adm = 0;
        return;
    }
    else {J_UV_adm = 1.e-21;}		/* irrelevant as long as it's not 0 */

    gJH0_adm = JAMPL_ADM * pow(10., (dzhi * log10(gH0_adm[ilow]) + dzlow * log10(gH0_adm[ilow + 1])) / (dzlow + dzhi));
    gJHe0_adm = JAMPL_ADM * pow(10., (dzhi * log10(gHe_adm[ilow]) + dzlow * log10(gHe_adm[ilow + 1])) / (dzlow + dzhi));
    gJHep_adm = JAMPL_ADM * pow(10., (dzhi * log10(gHep_adm[ilow]) + dzlow * log10(gHep_adm[ilow + 1])) / (dzlow + dzhi));
    epsH0_adm = JAMPL_ADM * pow(10., (dzhi * log10(eH0_adm[ilow]) + dzlow * log10(eH0_adm[ilow + 1])) / (dzlow + dzhi));
    epsHe0_adm = JAMPL_ADM * pow(10., (dzhi * log10(eHe_adm[ilow]) + dzlow * log10(eHe_adm[ilow + 1])) / (dzlow + dzhi));
    epsHep_adm = JAMPL_ADM * pow(10., (dzhi * log10(eHep_adm[ilow]) + dzlow * log10(eHep_adm[ilow + 1])) / (dzlow + dzhi));

    return;
}


void SetZeroIonization_adm(void)
{
    //printf("ADM Alert! cooling_adm, set zero ionisation\n");
    gJHe0_adm = gJHep_adm = gJH0_adm = 0; epsHe0_adm = epsHep_adm = epsH0_adm = 0; J_UV_adm = 0;
}


void IonizeParamsFunction_adm(void)
{
    //printf("ADM Alert! cooling_adm, ionise params function\n");
    int i, nint;
    double a0, planck, ev, e0_H, e0_He, e0_Hep;
    double gint, eint, t, tinv, fac, eps;
    double at, beta, s;
    double pi;

#define UVALPHA_ADM         1.0
    double Jold = -1.0;
    double redshift;

    J_UV_adm = 0.;
    gJHe0_adm = gJHep_adm = gJH0_adm = 0.;
    epsHe0_adm = epsHep_adm = epsH0_adm = 0.;


    if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
        redshift = 1 / All.Time - 1;

        if(redshift >= 6) {J_UV_adm = 0.;}
        else
        {
            if(redshift >= 3) {J_UV_adm = 4e-22 / (1 + redshift);}
            else
            {
                if(redshift >= 2) {J_UV_adm = 1e-22;}
                else {J_UV_adm = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);}
            }
        }
        if(J_UV_adm == Jold) {return;}
        Jold = J_UV_adm;
        if(J_UV_adm == 0) {return;}


        a0 = 6.30e-18;
        planck = 6.6262e-27;
        ev = 1.6022e-12;
        e0_H = 13.6058 * ev;
        e0_He = 24.59 * ev;
        e0_Hep = 54.4232 * ev;

        gint = 0.0;
        eint = 0.0;
        nint = 5000;
        at = 1. / ((double) nint);

        for(i = 1; i <= nint; i++)
        {
            t = (double) i;
            t = (t - 0.5) * at;
            tinv = 1. / t;
            eps = sqrt(tinv - 1.);
            fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA_ADM + 3.);
            gint += fac * at;
            eint += fac * (tinv - 1.) * at;
        }

        gJH0_adm = a0 * gint / planck;
        epsH0_adm = a0 * eint * (e0_H / planck);
        gJHep_adm = gJH0_adm * pow(e0_H / e0_Hep, UVALPHA_ADM) / 4.0;
        epsHep_adm = epsH0_adm * pow((e0_H / e0_Hep), UVALPHA_ADM - 1.) / 4.0;

        at = 7.83e-18;
        beta = 1.66;
        s = 2.05;

        gJHe0_adm = (at / planck) * pow((e0_H / e0_He), UVALPHA_ADM) *
        (beta / (UVALPHA_ADM + s) + (1. - beta) / (UVALPHA_ADM + s + 1));
        epsHe0_adm = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA_ADM) *
        (beta / (UVALPHA_ADM + s - 1) + (1 - 2 * beta) / (UVALPHA_ADM + s) - (1 - beta) / (UVALPHA_ADM + s + 1));

        pi = M_PI;
        gJH0_adm *= 4. * pi * J_UV_adm;
        gJHep_adm *= 4. * pi * J_UV_adm;
        gJHe0_adm *= 4. * pi * J_UV_adm;
        epsH0_adm *= 4. * pi * J_UV_adm;
        epsHep_adm *= 4. * pi * J_UV_adm;
        epsHe0_adm *= 4. * pi * J_UV_adm;
    }
}
#endif // !(CHIMES)




void InitCool_adm(void)
{
    //printf("ADM Alert! cooling_adm, init cool adm\n");
    if(ThisTask == 0) printf("Initializing ADM cooling ...\n");

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

#ifdef COOL_GRACKLE
    InitGrackle();
#endif

#ifdef CHIMES
    sprintf(ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", ChimesDataPath);
    sprintf(ChimesGlobalVars.EqAbundanceTablePath, "%s/EqAbundancesTables/%s", ChimesDataPath, ChimesEqAbundanceTable);
    sprintf(ChimesGlobalVars.PhotoIonTablePath[0], "%s/%s", ChimesDataPath, ChimesPhotoIonTable);

    // By default, use 1 spectrum, unless
    // stellar fluxes are enabled.
    ChimesGlobalVars.N_spectra = 1;

#ifdef CHIMES_STELLAR_FLUXES
    ChimesGlobalVars.N_spectra += CHIMES_LOCAL_UV_NBINS;

    int spectrum_idx;
    for (spectrum_idx = 1; spectrum_idx < ChimesGlobalVars.N_spectra; spectrum_idx++)
      sprintf(ChimesGlobalVars.PhotoIonTablePath[spectrum_idx], "%s/starburstCrossSections/cross_sections_SB%d.hdf5", ChimesDataPath, spectrum_idx);
#endif

    if (ChimesUVBMode > 0)
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = 0;
	if (ChimesUVBMode == 1)
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
	else
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 1;
      }
    else
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = -1;
	ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
      }

    // Hybrid cooling has not yet been
    // implemented in Gizmo. Switch
    // it off for now.
    ChimesGlobalVars.hybrid_cooling_mode = 0;

    // Set the chimes_exit() function
    // to use the Gizmo-specific
    // chimes_gizmo_exit().
    chimes_exit = &chimes_gizmo_exit;

    // Initialise the CHIMES module.
    init_chimes(&ChimesGlobalVars);

#ifdef CHIMES_METAL_DEPLETION
    chimes_init_depletion_data();
#endif

#else // CHIMES
    InitCoolMemory_adm();
    MakeCoolingTable_adm();
    ReadIonizeParams_adm("TREECOOL");
    IonizeParams_adm();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables_adm();
#endif
#endif // CHIMES
}



#ifndef CHIMES
#ifdef COOL_METAL_LINES_BY_SPECIES
double GetCoolingRateWSpecies_adm(double nHcgs, double logT, double *Z)
{
    printf("ADM Alert! cooling_adm, GetCoolingRateWSpecies_adm\n");
    double ne_over_nh_tbl=1, Lambda=0;
    int k, N_species_active = (int)NUM_LIVE_SPECIES_FOR_COOLTABLES;

    /* pre-calculate the indices for density and temperature, then we just need to call the tables by species */
    int ixmax=40, iymax=175;
    int ix0, iy0, ix1, iy1;
    double dx, dy, dz, mdz;
    long i_T=iymax+1, inHT=i_T*(ixmax+1);
    if(All.ComovingIntegrationOn && All.SpeciesTableInUse<48) {dz=log10(1/All.Time)*48; dz=dz-(int)dz; mdz=1-dz;} else {dz=0; mdz=1;}

    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    if(dx<0) {dx=0;} else {if(dx>ixmax) {dx=ixmax;}}
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) {ix1=ixmax;}
    dx=dx-ix0;
    if(dy<0) {dy=0;} else {if(dy>iymax) {dy=iymax;}}
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) {iy1=iymax;}
    dy=dy-iy0;
    long index_x0y0=iy0+ix0*i_T, index_x0y1=iy1+ix0*i_T, index_x1y0=iy0+ix1*i_T, index_x1y1=iy1+ix1*i_T;

    ne_over_nh_tbl = GetLambdaSpecies_adm(0,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz);
    if(ne_over_nh_tbl > 0)
    {
        double zfac = 0.0127 / All.SolarAbundances[0];
        for (k=1; k<N_species_active; k++)
        {
            long k_index = k * inHT;
            Lambda += GetLambdaSpecies_adm(k_index,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz) * Z[k+1]/(All.SolarAbundances[k+1]*zfac);
        }
        Lambda /= ne_over_nh_tbl;
    }
    return Lambda;
}


double GetLambdaSpecies_adm(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz)
{
    printf("ADM Alert! cooling_adm, getlambdaspecies_adm\n");
    long x0y0 = index_x0y0 + k_index;
    long x0y1 = index_x0y1 + k_index;
    long x1y0 = index_x1y0 + k_index;
    long x1y1 = index_x1y1 + k_index;
    double i1, i2, j1, j2, w1, w2, u1;
    i1 = SpCoolTable0_adm[x0y0];
    i2 = SpCoolTable0_adm[x0y1];
    j1 = SpCoolTable0_adm[x1y0];
    j2 = SpCoolTable0_adm[x1y1];
    if(dz > 0)
    {
        i1 = mdz * i1 + dz * SpCoolTable1_adm[x0y0];
        i2 = mdz * i2 + dz * SpCoolTable1_adm[x0y1];
        j1 = mdz * j1 + dz * SpCoolTable1_adm[x1y0];
        j2 = mdz * j2 + dz * SpCoolTable1_adm[x1y1];
    }
    w1 = i1*(1-dy) + i2*dy;
    w2 = j1*(1-dy) + j2*dy;
    u1 = w1*(1-dx) + w2*dx;
    return u1;
}

#endif // COOL_METAL_LINES_BY_SPECIES
#endif // !(CHIMES)


#ifdef GALSF_FB_FIRE_RT_UVHEATING
void selfshield_local_incident_uv_flux_adm(void)
{   /* include local self-shielding with the following */
printf("ADM Alert! cooling_adm, selfshield_local_uv_flux\n");
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type==0)
        {
	    //printf("ADM UV Flux selfshield: %.4e, %.4e\n",SphP[i].Rad_Flux_UV,SphP[i].Rad_Flux_EUV);
            if((SphP[i].Rad_Flux_UV>0) && (PPP[i].Hsml>0) && (SphP[i].Density>0) && (P[i].Mass>0) && (All.Time>0))
            {
                SphP[i].Rad_Flux_UV *= UNIT_FLUX_IN_CGS * 1276.19; SphP[i].Rad_Flux_EUV *= UNIT_FLUX_IN_CGS * 1276.19; // convert to Habing units (normalize strength to local MW field in this [narrow] band, so not the 'full' Habing flux)
                double surfdensity = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i); // in CGS
                double tau_nuv = rt_kappa(i,RT_FREQ_BIN_FIRE_UV) * surfdensity; // optical depth: this part is attenuated by dust //
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
//                tau_nuv *= (1.0e-3 + P[i].Metallicity[0]/All.SolarAbundances[0]); // if using older FIRE defaults, this was manually added instead of rolled into rt_kappa -- annoying but here for completeness //
//#endif
                double tau_euv = 3.7e6 * surfdensity * UNIT_SURFDEN_IN_CGS; // optical depth: 912 angstrom kappa_euv: opacity from neutral gas //
                SphP[i].Rad_Flux_UV *= exp(-DMIN(tau_nuv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                //SphP[i].Rad_Flux_UV *= 0.01 + 0.99/(1.0 + 0.8*tau_nuv + 0.85*tau_nuv*tau_nuv); // attenuate (for clumpy medium with hard-minimum 1% scattering)
                //SphP[i].Rad_Flux_EUV *= exp(-DMIN(tau_euv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                SphP[i].Rad_Flux_EUV *= 0.01 + 0.99/(1.0 + 0.8*tau_euv + 0.85*tau_euv*tau_euv); // attenuate (for clumpy medium with hard-minimum 1% scattering) //
            } else {SphP[i].Rad_Flux_UV = SphP[i].Rad_Flux_EUV = 0;}
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING) && !defined(CHIMES_HII_REGIONS)
//            if(SphP[i].DelayTimeHII > 0)
//            {   /* assign typical strong HII region flux + enough flux to maintain cell fully-ionized, regardless (x'safety-factor') */
//                double n1000 = SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / 1000.; // density in 1000 cm^-3
//                double flux_compactHII = DMAX(0.85*pow(n1000,1./3.) , 1) * 2.6e5*n1000; // set to typical value in HII region or minimum needed to maintain f_neutral < 1e-5-ish, whichever is larger
//                SphP[i].Rad_Flux_UV += flux_compactHII; SphP[i].Rad_Flux_EUV += flux_compactHII;
//           }
//#endif
        }
    }
}
#endif




/* subroutine to update the molecular fraction using our implicit solver for a simple --single-species-- network (just H2) 
 * This function does nothing unless COOL_MOLEFRAC_NONEQM is enabled.
 * */
void update_explicit_molecular_fraction_adm(int i, double dtime_cgs)
{
    printf("ADM Alert! cooling_adm, update_explicit_molecular_fraction\n");
    if(dtime_cgs <= 0) {return;}
#ifdef COOL_MOLECFRAC_NONEQM
    // first define a number of environmental variables that are fixed over this update step
    double fH2_initial = SphP[i].MolecularMassFraction_perNeutralH; // initial molecular fraction per H atom, entering this subroutine, needed for update below
    double xn_e=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties_adm(u0, rho, i, &mu_meanwt, &xn_e, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties [will assume fixed value of fH2 at previous update value]
    double T=1, Z_Zsol=1, urad_G0=1, xH0=1, x_e=0, nH_cgs=rho*UNIT_DENSITY_IN_NHCGS; // initialize definitions of some variables used below to prevent compiler warnings
    Z_Zsol=1; urad_G0=1; // initialize metal and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
//#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
//    if(SphP[i].DelayTimeHII > 0) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} // force gas flagged as in HII regions to have zero molecular fraction
//#endif
    if(temperature > 3.e5) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(nh0, 0.),1.); // get neutral fraction [given by call to this program]
    if(xH0 <= MIN_REAL_NUMBER) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} // no neutral gas, no molecules!
    x_e = DMIN(DMAX(xn_e, 0.),2.); // get free electron ratio [number per H nucleon]
    double gamma_12=return_local_gammamultiplier_adm(i)*gJH0_adm/1.0e-12, shieldfac=return_uvb_shieldfac_adm(i,gamma_12,nH_cgs,log10(T)), urad_from_uvb_in_G0=sqrt(shieldfac)*(gJH0_adm/2.29e-10); // estimate UVB contribution if we have partial shielding, to full photo-dissociation rates //
#ifdef METALS
    Z_Zsol = P[i].Metallicity[0]/All.SolarAbundances[0]; // metallicity in solar units [scale to total Z, since this mixes dust and C opacity], and enforce a low-Z floor to prevent totally unphysical behaviors at super-low Z [where there is still finite opacity in reality; e.g. Kramer's type and other opacities enforce floor around ~1e-3]
#endif
    /* get incident radiation field from whatever module we are using to track it */

///////////////////////////////
//  Commented out UV_HEATING //
//  ///////////////////////////

//#ifdef GALSF_FB_FIRE_RT_UVHEATING
//    urad_G0 = DMAX(SphP[i].Rad_Flux_UV, 1.e-10); // note this is ALREADY self-shielded by dust, so we need to be careful about 2x-counting the self-shielding approximation below; hence limit this to a rather sizeable value  //
//#endif
#if defined(RT_PHOTOELECTRIC) || defined(RT_LYMAN_WERNER)
    int whichbin = RT_FREQ_BIN_LYMAN_WERNER;
#if !defined(RT_LYMAN_WERNER)
    whichbin = RT_FREQ_BIN_PHOTOELECTRIC; // use photo-electric bin as proxy (very close) if don't evolve LW explicitly
#endif
    urad_G0 = SphP[i].Rad_E_gamma[whichbin] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
    urad_G0 += urad_from_uvb_in_G0; // include whatever is contributed from the meta-galactic background, fed into this routine
    urad_G0 = DMIN(DMAX( urad_G0 , 1.e-10 ) , 1.e10 ); // limit values, because otherwise exponential self-shielding approximation easily artificially gives 0 incident field
    // define a number of variables needed in the shielding module
    double dx_cell = Get_Particle_Size(i) * All.cf_atime; // cell size
    double surface_density_H2_0 = 5.e14 * PROTONMASS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
    w0 = 0.035; // actual calibration from Drain, Gnedin, Richings, others: 0.2 is more appropriate as a re-calibration for sims doing local eqm without ability to resolve shielding at higher columns
    //double surface_density_local = xH0 * SphP[i].Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. note G0 is -already- attenuated in the pre-processing by dust.
    double surface_density_local = xH0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth to infinity with our Sobolev-type approximation. Note G0 is already attenuated by dust, but we need to include H2 self-shielding, for which it is appropriate to integrate to infinity.
    double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
    double dv2=0; int j,k; for(j=0;j<3;j++) {for(k=0;k<3;k++) {double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
        if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
        dv2 += vt*vt;}} // calculate magnitude of the velocity shear across cell from || grad -otimes- v ||^(1/2)
    double dv_turb=sqrt(dv2)*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
    double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, fH2_min, Q_max, Q_min, Q_initial, Q_0, Q_1, fH2_0, fH2_1, fH2_new; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]
    double b_time_Mach = 0.5 * dv_turb / (v_thermal_rms/sqrt(3.)); // cs_thermal for molecular [=rms v_thermal / sqrt(3)], dv_turb to full inside dx, assume "b" prefactor for compressive-to-solenoidal ratio corresponding to the 'natural mix' = 0.5. could further multiply by 1.58 if really needed to by extended dvturb to 2h = H, and vthermal from molecular to atomic for the generating field, but not as well-justified
    double clumping_factor = 1. + b_time_Mach*b_time_Mach; // this is the exact clumping factor for a standard lognormal PDF with S=ln[1+b^2 Mach^2] //
    double clumping_factor_3 = clumping_factor*clumping_factor*clumping_factor; // clumping factor N for <rho^n>/<rho>^n = clumping factor^(N*(N-1)/2) //

    /* evolve dot[nH2]/nH0 = d_dt[fH2[neutral]] = (1/nH0) * (a_H2*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] ) */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, n_e=x_e*nH_cgs, EXPmax=40.; int iter=0; // define some variables for below, including neutral H number density, free electron number, etc.
    double x_p = DMIN(DMAX(nhp , x_e/10.), 2.); // get free H+ fraction [cap because irrelevant to below in very low regime //
    double a_Z  = (9.e-19 * T / (1. + 0.04*sqrt_T + 0.002*T + 8.e-6*T*T)) * (0.5*Z_Zsol) * nH_cgs * clumping_factor; // dust formation
    double a_GP = (1.833e-18 * pow(T,0.88)) / (1. + x_p*1846.*(1.+T/20000.)/sqrt(T)) * n_e * clumping_factor; // gas-phase formation [Glover & Abel 2008, using fitting functions slightly more convenient and assuming H-->H2 much more rapid than other reactions, from Krumholz & McKee 2010; denominator factor accounts for p+H- -> H + H, instead of H2]
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * clumping_factor_3; // 3-body collisional formation
    double b_H2HI = (7.073e-19 * pow(T,2.012) * exp(-DMIN(5.179e4/T,EXPmax)) / pow(1. + 2.130e-5*T , 3.512)) * (nH0/2.) * clumping_factor; // collisional H2-H dissociation
    b_H2HI += 4.49e-9 * pow(T,0.11) * exp(-DMIN(101858./T,EXPmax)) * (n_e) * clumping_factor; // collisional H2-e- dissociation [note assuming ground-state optically thin dissociation here as thats where this is most relevant, see Glover+Abel 2008)
    double b_H2H2 = (5.996e-30 * pow(T,4.1881) * exp(-DMIN(5.466e4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881)) * (nH0/2.) * (1./2.) * clumping_factor; // collisional H2-H2 dissociation
    double G_LW = 3.3e-11 * urad_G0 * (1./2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0
    double xi_cr_H2 = (7.525e-16) * (1./2.), prefac_CR=1.;; // CR dissociation (+ionization)
#if defined(COSMIC_RAYS) // scale ionization+dissociation rates with local CR energy density
    prefac_CR=0; {int kcr; for(kcr=0;kcr<N_CR_PARTICLE_BINS;kcr++) {prefac_CR += SphP[i].CosmicRayEnergyPred[kcr];}} // add up CR energy
    prefac_CR *= (SphP[i].Density * All.cf_a3inv / P[i].Mass) * UNIT_PRESSURE_IN_EV; // convert to CR volume energy density in eV/cm^3
    prefac_CR /= 3.0; // scale ionization rate relative to the CR energy density [in eV/cm3] assumed to give rise to this level of ionization [from Indriolo], for a universal spectrum
#else
    if(All.ComovingIntegrationOn) {double rhofac = (rho*UNIT_DENSITY_IN_CGS)/(1000.*COSMIC_BARYON_DENSITY_CGS);
        if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
#endif
    xi_cr_H2 *= prefac_CR;

    // want to solve the implicit equation: f_f = f_0 + g[f_f]*dt, where g[f_f] = df_dt evaluated at f=f_f, so root-find: dt*g[f_f] + f_0-f_f = 0
    // can write this as a quadtratic: x_a*f^2 - x_b_0*f - xb_LW*f + x_c = 0, where xb_LW is a non-linear function of f accounting for the H2 self-shielding terms
    double G_LW_dt_unshielded = G_LW * dtime_cgs; // LW term without shielding, multiplied by timestep for dimensions needed below
    double x_a = (b_3B + b_H2HI - b_H2H2) * dtime_cgs; // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b_0 = (a_GP + a_Z + 2.*b_3B + b_H2HI + xi_cr_H2) * dtime_cgs + 1.; // terms linear in f [note sign, pulling the -sign out here] -- positive-definite
    double x_c = (a_GP + a_Z + b_3B) * dtime_cgs + fH2_initial; // terms independent of f -- positive-definite
    double y_a = x_a / (x_c + MIN_REAL_NUMBER), x_b, y_b, z_a; // convenient to convert to dimensionless variable needed for checking definite-ness
    // use the previous-timestep value of fH2 to guess the shielding term and then compute the resulting fH2
    fH2_tmp=fH2_initial; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
    Q_initial = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the previous value of fH2
    z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant

    /* now comes the tricky bit -- need to account for non-linear part of the solution for the molecular line self-shielding [depends on fH2, not just the dust external shielding already accounted for */
    if((fH2 > 1.e-10) && (fH2 < 1) && (G_LW_dt_unshielded > 0.1*x_b_0)) // fH2 is non-trivial, and the radiation term is significant, so to get an accurate update we need to invoke a non-linear solver here
    {
        // set updated guess values
        fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        fH2_1 = fH2; Q_1 = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        // set lower values for bracketing
        x_b=x_b_0+G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // recalculate all terms that depend on the shielding
        fH2_min = DMAX(0,DMIN(1,fH2)); // this serves as a lower-limit for fH2
        fH2_tmp=fH2_min; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        Q_min = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        // set upper values for bracketing
        fH2_tmp=1.; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for fH2
        fH2_tmp=fH2_max; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        Q_max = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        if(fH2_1 < fH2_min) {fH2 = fH2_min;} // hitting lower bound already, set to that value and exit
        else if(fH2_1 > fH2_max) {fH2 = fH2_max;} // hitting upper bound already, set to that value and exit
        else if(Q_min*Q_max >= 0) // bracketing indicates that in this timestep, we will move fully to one or the other limit -- so do that, and don't need to iterate!
            {if(fabs(Q_min) < fabs(Q_max)) {if(fabs(Q_min) < fabs(Q_1)) {fH2 = fH2_min;}} else {if(fabs(Q_max) < fabs(Q_1)) {fH2 = fH2_max;}}} // decide if Qmin/max corresponds more closely to desired zero, so move to fH2 matching that value
        else if(fH2_max > 1.01*fH2_min) // worth attempting to iterate
        {
            Q_0 = Q_initial; fH2_0 = fH2_initial; // first iteration step is already done in all the above
            fH2 = exp( (log(fH2_0)*Q_1 - log(fH2_1)*Q_0) / (Q_1-Q_0) ); // do a Newton-Raphson step in log[f_H2] space now that we have good initial brackets
            if(fH2 > fH2_max) // ok we overshot the upper limit, test if bracketing works between fH2 and fH2_max, otherwise just use fH2_max
                {if(Q_1*Q_max < 0) {fH2 = exp( (log(fH2_max)*Q_1 - log(fH2_1)*Q_max) / (Q_1-Q_max) );} else {fH2=fH2_max;}}
            else if(fH2 < fH2_min)
                {if(Q_1*Q_min < 0) {fH2 = exp( (log(fH2_min)*Q_1 - log(fH2_1)*Q_min) / (Q_1-Q_min) );} else {fH2=fH2_min;}}
            else while(1)
            {
                    fH2_tmp=fH2_1; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
                    Q_1 = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2
                    //if(Q_1*Q_0 >= 0) {break;} // no longer bracketing, end while loop
                    fH2_new = exp( (log(fH2_0)*Q_1 - log(fH2_1)*Q_0) / (Q_1-Q_0) ); fH2_0=fH2_1; Q_0=Q_1; fH2_1=fH2_new; // update guess and previous values //
                    iter++; // count iterations
                    if(fabs(fH2_1-fH2_0) < 0.01*(0.5*(fH2_1+fH2_0))) {break;} // converged well enough for our purposes!
                    if((y_ss > 0.95) || (y_ss*G_LW_dt_unshielded < 0.1*x_b)) {break;} // negligible shielding, or converged to point where external LW is not dominant dissociator so no further iteration needed
                    if((fH2 > 0.95*fH2_max) || (fH2 > 0.99) || (fH2 < 1.e-10) || (fH2 < 1.05*fH2_min) || (iter > 10)) {break;} // approached physical limits or bounds of validity, or end of iteration cycle
            } // end of convergence iteration to find solution for fmol
        } // opened plausible iteration clause
    } // opened self-shielding clause [attempting to bracket]
    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    SphP[i].MolecularMassFraction_perNeutralH = fH2; // record variable -- this will be used for the next update, meanwhile the total fraction will be used in various routines through the code
    SphP[i].MolecularMassFraction = xH0 * SphP[i].MolecularMassFraction_perNeutralH; // record variable -- this is largely what is needed below
#endif
}





/* simple subroutine to estimate the dust temperatures in our runs without detailed tracking of these individually [more detailed chemistry models do this] */
double get_equilibrium_dust_temperature_estimate_adm(int i, double shielding_factor_for_exgalbg)
{   /* simple three-component model [can do fancier] with cmb, dust, high-energy photons */
printf("ADM Alert! cooling_adm, get_equilibrium_dust_temperature_estimate_adm\n");
#if defined(RT_INFRARED)
    if(i >= 0) {return SphP[i].Dust_Temperature;} // this is pre-computed -- simply return it
#endif
    double e_CMB=0.262*All.cf_a3inv/All.cf_atime, T_cmb=2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    double e_IR=0.31, Tdust_ext=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_HiEgy=0.66, T_hiegy=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
#ifdef RT_ISRF_BACKGROUND
    e_IR *= RT_ISRF_BACKGROUND; e_HiEgy *= RT_ISRF_BACKGROUND; // need to re-scale the assumed ISRF components
#endif
    if(i >= 0)
    {
//////////////////////////////////////////////////////////////////////////
// COMMENTED OUT ALL RT_USE_GRAVTREE LINES HERE. EDIT FOR FUTURE USE!!! //
// ///////////////////////////////////////////////////////////////////////
//#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
/*        e_HiEgy=0; e_IR = 0; int k; double E_tot_to_evol_eVcgs = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++) {e_HiEgy+=SphP[i].Rad_E_gamma_Pred[k];}
#if defined(GALSF_FB_FIRE_RT_LONGRANGE)
        e_IR += SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_FIRE_IR]; // note IR
#endif
#if defined(RT_INFRARED)
        e_IR += SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED]; Tdust_ext = SphP[i].Radiation_Temperature; // note IR [irrelevant b/c of call above, but we'll keep this as a demo]
#endif
        e_HiEgy -= e_IR; // don't double-count the IR component flagged above //
        e_IR *= E_tot_to_evol_eVcgs; e_HiEgy *= E_tot_to_evol_eVcgs;
#endif */
    }
    e_HiEgy += shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    double Tdust_eqm = 10.; // arbitrary initial value //
    if(Tdust_ext*e_IR < 1.e-10 * (T_cmb*e_CMB + T_hiegy*e_HiEgy)) { // IR term is totally negligible [or zero exactly], use simpler expression assuming constant temperature for it to avoid sensitivity to floating-pt errors //
        Tdust_eqm = 2.92 * pow(Tdust_ext*e_IR + T_cmb*e_CMB + T_hiegy*e_HiEgy, 1./5.); // approximate equilibrium temp assuming Q~1/lambda [beta=1 opacity law], assuming background IR temp is a fixed constant [relevant in IR-thin limit, but we don't know T_rad, so this is a guess anyways]
    } else { // IR term is not vanishingly small. we will assume the IR radiation temperature is equal to the local Tdust. lacking any direct evolution of that field, this is a good proxy, and exact in the locally-IR-optically-thick limit. in the locally-IR-thin limit it slightly under-estimates Tdust, but usually in that limit the other terms dominate anyways, so this is pretty safe //
        double T0=2.92, q=pow(T0*e_IR,0.25), y=(T_cmb*e_CMB + T_hiegy*e_HiEgy)/(T0*e_IR*q); if(y<=1) {Tdust_eqm=T0*q*(0.8+sqrt(0.04+0.1*y));} else {double y5=pow(y,0.2), y5_3=y5*y5*y5, y5_4=y5_3*y5; Tdust_eqm=T0*q*(1.+15.*y5_4+sqrt(1.+30.*y5_4+25.*y5_4*y5_4))/(20.*y5_3);} // this gives an extremely accurate and exactly-joined solution to the full quintic equation assuming T_rad_IR=T_dust
    }
    return DMAX(DMIN(Tdust_eqm , 2000.) , 1.); // limit at sublimation temperature or some very low temp //
}



/* this function estimates the free electron fraction from heavy ions, assuming a simple mix of cold molecular gas, Mg, and dust, with the ions from singly-ionized Mg, to prevent artificially low free electron fractions */
double return_electron_fraction_from_heavy_ions_adm(int target, double temperature, double density_cgs, double n_elec_HHe)
{
    printf("ADM Alert! cooling_adm, return_electron_fraction_from_heavy_ions_adm\n");
    if(All.ComovingIntegrationOn) {double rhofac=density_cgs/(1000.*COSMIC_BARYON_DENSITY_CGS); if(rhofac<0.2) {return 0;}} // ignore these reactions in the IGM
    double zeta_cr=1.0e-17, f_dustgas=0.01, n_ion_max=4.1533e-5, XH=HYDROGEN_MASSFRAC; // cosmic ray ionization rate (fixed as constant for non-CR runs) and dust-to-gas ratio
#ifdef COSMIC_RAYS
    if(target>=0) {double u_cr=0; int k; for(k=0;k<N_CR_PARTICLE_BINS;k++) {u_cr += SphP[target].CosmicRayEnergyPred[k];}
        zeta_cr = u_cr * 2.2e-6 * ((SphP[target].Density*All.cf_a3inv / P[target].Mass) * (UNIT_PRESSURE_IN_CGS));} // convert to ionization rate
#endif
#ifdef METALS
    if(target>=0) {f_dustgas=0.5*P[target].Metallicity[0];} // constant dust-to-metals ratio
#ifdef COOL_METAL_LINES_BY_SPECIES
    if(target>=0) {n_ion_max = (All.SolarAbundances[6]/24.3)/XH;} // limit, to avoid over-ionization at low metallicities
#endif
#endif
    /* Regime I: highly/photo-ionized, any contributions here would be negligible -- no need to continue */
    if(n_elec_HHe > 0.01) {return n_ion_max;} // contribute something negligible, doesn't matter here //
    double a_grain_micron=0.1, m_ion=24.305*PROTONMASS, mu_eff=2.38, m_neutrals=mu_eff*PROTONMASS, m_grain=4.189e-12*(2.4)*a_grain_micron*a_grain_micron*a_grain_micron, ngrain_ngas=(m_neutrals/m_grain)*f_dustgas; // effective size of grains that matter at these densities, and ions [here Mg] that dominate
    double k_ei=9.77e-8, y=sqrt(m_ion/ELECTRONMASS), ln_y=log(y), psi_0=-ln_y+(1.+ln_y)*log(1.+ln_y)/(2.+ln_y), k_eg_00=0.0195*a_grain_micron*a_grain_micron*sqrt(temperature), k_eg_0=k_eg_00*exp(psi_0), n_crit=k_ei*zeta_cr/(k_eg_0*ngrain_ngas*k_eg_0*ngrain_ngas), n_eff=density_cgs/m_neutrals;

    /* Regime II: CR-ionized with high enough ionization fraction that gas-phase recombinations dominate */
    if(n_eff < 0.01*n_crit) {return DMIN(n_ion_max , sqrt(zeta_cr / (k_ei * n_eff)) / (XH*mu_eff) );} // CR-dominated off gas with recombination -- put into appropriate units here
    if(n_eff < 100.*n_crit) {return DMIN(n_ion_max , 0.5 * (sqrt(4.*n_crit/n_eff + 1.) - 1.) * (k_eg_0*ngrain_ngas)/(k_ei*XH*mu_eff));} // interpolates between gas-recombination and dust-dominated regimes

    /* Regime III: recombination dominated by dust, but dust has a 'normal' efficiency of absorbing grains */
    double psi_fac=16.71/(a_grain_micron*temperature), alpha=zeta_cr*psi_fac/(k_eg_00 * ngrain_ngas*ngrain_ngas * n_eff), alpha_min=0.02, alpha_max=10.; /* Z*psi_fac = Z*e^2/(a_grain*kB*T) to dimensionless [psi] units; alpha=prefactor in charge equation: psi = alpha*(exp[-psi] - y/(1-psi)) */
    if(alpha > alpha_max) {return DMIN(n_ion_max , zeta_cr / (k_eg_0*ngrain_ngas * XH*mu_eff*n_eff ) );}

    /* Regime IV: recombination dominated by dust and grains dominate the charge, strongly suppressing the free charges */
    if(alpha < 1.e-4) {return DMIN(n_ion_max , zeta_cr / (k_eg_00*ngrain_ngas * XH*mu_eff*n_eff) );} // psi->small, negligible correction here
    if(alpha < alpha_min) {double psi=0.5*(1.-sqrt(1.+4.*y*alpha)); return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff) );} // small-psi-limit
    double psi_xmin=0.5*(1.-sqrt(1.+4.*y*alpha_min)), psi=psi_0 + (psi_xmin-psi_0)*2./(1.+alpha_min/alpha); // this interpolates between the asymptotic limmits at low and high alpha, where we can obtain high-accuracy solutions here
    return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff));
}




/* this function evaluates Compton heating+cooling rates and synchrotron cooling for thermal gas populations, accounting for the
    explicitly-evolved radiation field if it is evolved (otherwise assuming a standard background), and B-fields if they
    are evolved, as well as the proper relativistic or non-relativistic effects and two-temperature plasma effects. */
double evaluate_Compton_heating_cooling_rate_adm(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg)
{
    printf("ADM Alert! cooling_adm, evaluate_Compton_heating_cooling_rate_adm\n");
    double Lambda = 0;
    //double compton_prefac_eV = 2.16e-35 / nHcgs; // multiply field in eV/cm^3 by this and temperature difference to obtain rate

    //double e_CMB_eV=0.262*All.cf_a3inv/All.cf_atime, T_cmb = 2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    //Lambda += compton_prefac_eV * n_elec * e_CMB_eV * (T-T_cmb);

    //double e_UVB_eV = shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    //Lambda += compton_prefac_eV * n_elec * e_UVB_eV * (T-2.e4); // assume very crude approx Compton temp ~2e4 for UVB

    double compton_prefac_eV = 1.9e-37/nHcgs; // multiply field in eV/cm^3 by this and temperature difference to obtain rate

    double e_CMB_eV=pow(ELECTRONMASS/All.ADM_ElectronMass,3.0)*pow(All.ADM_FineStructure/0.01,2), T_cmb = 1.35/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    Lambda += compton_prefac_eV * n_elec * e_CMB_eV * (T-T_cmb) * pow(T_cmb, 4.0);

/////////////////////////////////////////////
// COMMENTED OUT ALL RT_USE_GRAVTREE TERMS //
// //////////////////////////////////////////

//    double e_UVB_eV = shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
//    Lambda += compton_prefac_eV * n_elec * e_UVB_eV * (T-2.e4); // assume very crude approx Compton temp ~2e4 for UVB

//#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
/*    if(target >= 0)
    {
        int k; double E_tot_to_evol_eVcgs = (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            double e_tmp = SphP[target].Rad_E_gamma_Pred[k] * E_tot_to_evol_eVcgs, Teff = 0;
*/
//#if defined(GALSF_FB_FIRE_RT_LONGRANGE) /* three-band (UV, OPTICAL, IR) approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
//            if(k==RT_FREQ_BIN_FIRE_IR) {Teff=30.;}
//            if(k==RT_FREQ_BIN_FIRE_OPT) {Teff=4000.;}
//            if(k==RT_FREQ_BIN_FIRE_UV) {Teff=15000.;}
//#endif
//#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
//            if(k==RT_FREQ_BIN_INFRARED) {Teff=SphP[target].Dust_Temperature;}
//#endif
//#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
//            if(k==RT_FREQ_BIN_OPTICAL_NIR) {Teff=2800.;}
//#endif
//#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
//            if(k==RT_FREQ_BIN_NUV) {Teff=12000.;}
//#endif
//#if defined(RT_PHOTOELECTRIC) /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
//            if(k==RT_FREQ_BIN_PHOTOELECTRIC) {Teff=24400.;}
//#endif
//#if defined(RT_LYMAN_WERNER) /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
//            if(k==RT_FREQ_BIN_LYMAN_WERNER) {Teff=28800.;}
//#endif
//#if defined(RT_CHEM_PHOTOION) /* Hydrogen and Helium ionizing bands: H0 here */
//            if(k==RT_FREQ_BIN_H0) {Teff=2340.*rt_nu_eff_eV[k];}
//#endif
//#if defined(RT_PHOTOION_MULTIFREQUENCY) /* Hydrogen and Helium ionizing bands: He bands */
//            if(k==RT_FREQ_BIN_He0 || k==RT_FREQ_BIN_He1 || k==RT_FREQ_BIN_He2) {Teff=2340.*rt_nu_eff_eV[k];}
//#endif
//#if defined(RT_SOFT_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
//            if(k==RT_FREQ_BIN_SOFT_XRAY) {Teff=3.6e6;}
//#endif
//#if defined(RT_HARD_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
//            if(k==RT_FREQ_BIN_HARD_XRAY) {Teff=1.7e7;}
//#endif
//            if(Teff < 3.e4) {e_tmp *= n_elec;} // low-energy radiation acts inefficiently on neutrals here
///            Lambda += compton_prefac_eV * e_tmp * (T - Teff); // add to compton heating/cooling terms
//        }
//    }
//#else // no explicit RHD terms evolved, so assume a MW-like ISRF instead
//    double e_IR_eV=0.31, T_IR=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
//    double e_OUV_eV=0.66, T_OUV=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
//    Lambda += compton_prefac_eV * n_elec * (e_IR_eV*(T-T_IR) + e_OUV_eV*(T-T_OUV));
//#endif

//#ifdef BH_COMPTON_HEATING /* custom band to represent (non)relativistic X-ray compton cooling from an AGN source without full RHD */
//    if(target >= 0)
//    {
//        double e_agn = (SphP[target].Rad_Flux_AGN * UNIT_FLUX_IN_CGS) / (C_LIGHT * ELECTRONVOLT_IN_ERGS), T_agn=2.e7; /* approximate from Sazonov et al. */
//        Lambda += compton_prefac_eV * e_agn * (T-T_agn); // since the heating here is primarily hard X-rays, and the cooling only relevant for very high temps, do not have an n_elec here
//    }
//#endif

//#ifdef MAGNETIC /* include sychrotron losses as well as long as we're here, since these scale more or less identically just using the magnetic instead of radiation energy */
//    if(target >= 0)
//    {
//        double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG;
//        Lambda += compton_prefac_eV * U_mag_ev * T; // synchrotron losses proportional to temperature (non-relativistic here), as inverse compton, just here without needing to worry about "T-T_eff", as if T_eff->0
//    }
//#endif

//    double T_eff_for_relativistic_corr = T; /* used below, but can be corrected */
//    if(Lambda > 0) /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
//    {
//        double Lambda_limiter_var = 1.483e34 * Lambda*Lambda*T; /* = (Lambda/2.6e-22)^2 * (T/1e9): if this >> 1, follow CAFG and cap at cooling rate assuming equilibrium e- temp from Coulomb exchange balancing compton */
//        if(Lambda_limiter_var > 1) {Lambda_limiter_var = 1./pow(Lambda_limiter_var,0.2); Lambda*=Lambda_limiter_var; T_eff_for_relativistic_corr*=Lambda_limiter_var;}
//    }
//    if(T_eff_for_relativistic_corr > 3.e7) {Lambda *= (T_eff_for_relativistic_corr/1.5e9) / (1-exp(-T_eff_for_relativistic_corr/1.5e9));} /* relativistic correction term, becomes important at >1e9 K, enhancing rate */

    return Lambda;
}







/*  this function computes the self-consistent temperature and electron fraction */
double ThermalProperties_adm(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess)
{

printf("ADM Alert! cooling_adm, thermal properties adm\n");
#if defined(CHIMES)
    int i = target; *ne_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_elec]]; *nH0_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HI]];
    *nHp_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HII]]; *nHe0_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeI]];
    *nHep_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeII]]; *nHepp_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeIII]];
    double temp = ChimesGasVars[target].temperature;
    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp, rho, nH0_guess, ne_guess, 0, target);
    return temp;
#else
    if(target >= 0) {*ne_guess=SphP[target].Ne; *nH0_guess = DMAX(0,DMIN(1,1.-( *ne_guess / 1.2 )));} else {*ne_guess=1.; *nH0_guess=0.;}
    rho *= UNIT_DENSITY_IN_CGS; u *= UNIT_SPECEGY_IN_CGS;   /* convert to physical cgs units */
    double temp = convert_u_to_temp_adm(u, rho, target, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess);
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING) && !defined(CHIMES_HII_REGIONS)
//    if(target >= 0) {if(SphP[target].DelayTimeHII > 0) {SphP[target].Ne = 1.0 + 2.0*yhelium(target); *nH0_guess=0; nHe0_guess=0;}} /* fully ionized [if using older model] */
//    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp, rho, nH0_guess, ne_guess, 0, target);
//#endif
    return temp;
#endif
}



/* function to return the local multiplier relative to the UVB model to account in some local RHD models for local ionizing sources */
double return_local_gammamultiplier_adm(int target)
{
printf("ADM Alert! cooling_adm, local_gammamultiplier\n");
// Added this first line to automatically turn off any local UVB multiplier for ADM. Edit this for future, feedback-related runs!
return 1;

#if defined(GALSF_FB_FIRE_RT_UVHEATING) && !defined(CHIMES)
    if((target >= 0) && (gJH0_adm > 0))
    {
        double local_gammamultiplier = SphP[target].Rad_Flux_EUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
        local_gammamultiplier = 1.0 + local_gammamultiplier / gJH0_adm; // this needs to live here in cooling.c where gJH0_adm is declared as a global shared variable!
        if(!isfinite(local_gammamultiplier)) {local_gammamultiplier=1;}
        return DMAX(1., DMIN(1.e20, local_gammamultiplier));
    }
#endif
    return 1;
}


/* function to attenuate the UVB to model self-shielding in optically-thin simulations */
double return_uvb_shieldfac_adm(int target, double gamma_12, double nHcgs, double logT)
{
    printf("ADM Alert! cooling_adm, return_uvb_shieldfac\n");
#ifdef GALSF_EFFECTIVE_EQS
    return 1; // self-shielding is implicit in the sub-grid model already //
#endif
//#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION > 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
//    if(target>=0) {if(SphP[target].DelayTimeHII > 0) {return 1;}} // newer HII region model irradiates and removes shielding for regions, but allows cooling function to evolve //
//#endif
    double NH_SS_z, NH_SS = 0.0123; /* CAFG: H number density above which we assume no ionizing bkg (proper cm^-3) */
    if(gamma_12>0) {NH_SS_z = NH_SS*pow(gamma_12,0.66)*pow(10.,0.173*(logT-4.));} else {NH_SS_z = NH_SS*pow(10.,0.173*(logT-4.));}
    double q_SS = nHcgs/NH_SS_z;
#ifdef COOL_UVB_SELFSHIELD_RAHMATI
    return 0.98 / pow(1 + pow(q_SS,1.64), 2.28) + 0.02 / pow(1 + q_SS*(1.+1.e-4*nHcgs*nHcgs*nHcgs*nHcgs), 0.84); // from Rahmati et al. 2012: gives gentler cutoff at high densities. but we need to modify it with the extra 1+(nHcgs/10)^4 denominator term since at very high nH, this cuts off much too-slowly (as nH^-0.84), which means UVB heating can be stronger than molecular cooling even at densities >> 1e4
#else
    return 1./(1.+q_SS*(1.+q_SS/2.*(1.+q_SS/3.*(1.+q_SS/4.*(1.+q_SS/5.*(1.+q_SS/6.*q_SS)))))); // this is exp(-q) down to ~1e-5, then a bit shallower, but more than sufficient approximation here //
#endif
}



#ifdef CHIMES
/* This routine updates the ChimesGasVars structure for particle target. */
void chimes_update_gas_vars(int target)
{
  double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);
  double u_old_cgs = DMAX(All.MinEgySpec, SphP[target].InternalEnergy) * UNIT_SPECEGY_IN_CGS;
  double rho_cgs = SphP[target].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

#ifdef COOL_METAL_LINES_BY_SPECIES
  double H_mass_fraction = 1.0 - (P[target].Metallicity[0] + P[target].Metallicity[1]);
#else
  double H_mass_fraction = XH;
#endif

  ChimesGasVars[target].temperature = (ChimesFloat) chimes_convert_u_to_temp_adm(u_old_cgs, rho_cgs, target);
  ChimesGasVars[target].nH_tot = (ChimesFloat) (H_mass_fraction * rho_cgs / PROTONMASS);
  ChimesGasVars[target].ThermEvolOn = All.ChimesThermEvolOn;

  // If there is an EoS, need to set TempFloor to that instead.
  ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);
//#if (GALSF_FB_FIRE_STELLAREVOLUTION <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
//    if (SphP[target].DelayTimeHII > 0) {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(HIIRegion_Temp, 10.1);}
//        else {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);}
//#endif

  // Flag to control how the temperature
  // floor is implemented in CHIMES.
  ChimesGasVars[target].temp_floor_mode = 1;

  // Extragalactic UV background
  ChimesGasVars[target].isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0];
  ChimesGasVars[target].isotropic_photon_density[0] *= chimes_rad_field_norm_factor;

  ChimesGasVars[target].G0_parameter[0] = chimes_table_spectra.G0_parameter[0];
  ChimesGasVars[target].H2_dissocJ[0] =  chimes_table_spectra.H2_dissocJ[0];

#ifdef CHIMES_STELLAR_FLUXES
  int kc;
  for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++)
    {
      ChimesGasVars[target].isotropic_photon_density[kc + 1] = (ChimesFloat) (SphP[target].Chimes_fluxPhotIon[kc] / C_LIGHT);

#ifdef CHIMES_HII_REGIONS
    if(SphP[target].DelayTimeHII > 0) {
        ChimesGasVars[target].isotropic_photon_density[kc + 1] += (ChimesFloat) (SphP[target].Chimes_fluxPhotIon_HII[kc] / C_LIGHT);
        ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) ((SphP[target].Chimes_G0[kc] + SphP[target].Chimes_G0_HII[kc]) / DMAX((SphP[target].Chimes_fluxPhotIon[kc] + SphP[target].Chimes_fluxPhotIon_HII[kc]), 1.0e-300));
	} else {ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (SphP[target].Chimes_G0[kc] / DMAX(SphP[target].Chimes_fluxPhotIon[kc], 1.0e-300));}
    ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#else
      ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (SphP[target].Chimes_G0[kc] / DMAX(SphP[target].Chimes_fluxPhotIon[kc], 1.0e-300));
      ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#endif
    }
#endif

  ChimesGasVars[target].cr_rate = (ChimesFloat) cr_rate;  // For now, assume a constant cr_rate.
  ChimesGasVars[target].hydro_timestep = (ChimesFloat) (dt * UNIT_TIME_IN_CGS);

  ChimesGasVars[target].ForceEqOn = ChimesEqmMode;
  ChimesGasVars[target].divVel = (ChimesFloat) P[target].Particle_DivVel / UNIT_TIME_IN_CGS;
  if (All.ComovingIntegrationOn)
    {
      ChimesGasVars[target].divVel *= (ChimesFloat) All.cf_a2inv;
      ChimesGasVars[target].divVel += (ChimesFloat) (3 * All.cf_hubble_a / UNIT_TIME_IN_CGS);  /* Term due to Hubble expansion */
    }
  ChimesGasVars[target].divVel = (ChimesFloat) fabs(ChimesGasVars[target].divVel);

#ifndef COOLING_OPERATOR_SPLIT
  ChimesGasVars[target].constant_heating_rate = ChimesGasVars[target].nH_tot * ((ChimesFloat) SphP[target].DtInternalEnergy);
#else
  ChimesGasVars[target].constant_heating_rate = 0.0;
#endif

#ifdef CHIMES_SOBOLEV_SHIELDING
  double surface_density;
  surface_density = evaluate_NH_from_GradRho(P[target].GradRho,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs
  ChimesGasVars[target].cell_size = (ChimesFloat) (shielding_length_factor * surface_density / rho_cgs);
#else
  ChimesGasVars[target].cell_size = 1.0;
#endif

  ChimesGasVars[target].doppler_broad = 7.1;  // km/s. For now, just set this constant. Thermal broadening is also added within CHIMES.

  ChimesGasVars[target].InitIonState = ChimesInitIonState;

#ifdef CHIMES_HII_REGIONS
  if (SphP[target].DelayTimeHII > 0.0) {ChimesGasVars[target].cell_size = 1.0;} // switch of shielding in HII regions
#endif

#ifdef CHIMES_TURB_DIFF_IONS
  chimes_update_turbulent_abundances(target, 0);
#endif

#if defined(COOL_METAL_LINES_BY_SPECIES) && !defined(GALSF_FB_NOENRICHMENT)
  chimes_update_element_abundances(target);
#endif

  return;
}

#ifdef COOL_METAL_LINES_BY_SPECIES
/* This routine re-computes the element abundances from
 * the metallicity array and updates the individual ion
 * abundances accordingly. */
void chimes_update_element_abundances(int i)
{
  double H_mass_fraction = 1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1]);

  /* Update the element abundances in ChimesGasVars. */
  ChimesGasVars[i].element_abundances[0] = (ChimesFloat) (P[i].Metallicity[1] / (4.0 * H_mass_fraction));   // He
  ChimesGasVars[i].element_abundances[1] = (ChimesFloat) (P[i].Metallicity[2] / (12.0 * H_mass_fraction));  // C
  ChimesGasVars[i].element_abundances[2] = (ChimesFloat) (P[i].Metallicity[3] / (14.0 * H_mass_fraction));  // N
  ChimesGasVars[i].element_abundances[3] = (ChimesFloat) (P[i].Metallicity[4] / (16.0 * H_mass_fraction));  // O
  ChimesGasVars[i].element_abundances[4] = (ChimesFloat) (P[i].Metallicity[5] / (20.0 * H_mass_fraction));  // Ne
  ChimesGasVars[i].element_abundances[5] = (ChimesFloat) (P[i].Metallicity[6] / (24.0 * H_mass_fraction));  // Mg
  ChimesGasVars[i].element_abundances[6] = (ChimesFloat) (P[i].Metallicity[7] / (28.0 * H_mass_fraction));  // Si
  ChimesGasVars[i].element_abundances[7] = (ChimesFloat) (P[i].Metallicity[8] / (32.0 * H_mass_fraction));  // S
  ChimesGasVars[i].element_abundances[8] = (ChimesFloat) (P[i].Metallicity[9] / (40.0 * H_mass_fraction));  // Ca
  ChimesGasVars[i].element_abundances[9] = (ChimesFloat) (P[i].Metallicity[10] / (56.0 * H_mass_fraction)); // Fe

  ChimesGasVars[i].metallicity = (ChimesFloat) (P[i].Metallicity[0] / 0.0129);  // In Zsol. CHIMES uses Zsol = 0.0129.
  ChimesGasVars[i].dust_ratio = ChimesGasVars[i].metallicity;

#ifdef CHIMES_METAL_DEPLETION
#ifdef _OPENMP
  int ThisThread = omp_get_thread_num();
#else
  int ThisThread = 0;
#endif
  chimes_compute_depletions(ChimesGasVars[i].nH_tot, ChimesGasVars[i].temperature, ThisThread);

  ChimesGasVars[i].element_abundances[1] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[0]; // C
  ChimesGasVars[i].element_abundances[2] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[1]; // N
  ChimesGasVars[i].element_abundances[3] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[2]; // O
  ChimesGasVars[i].element_abundances[5] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[3]; // Mg
  ChimesGasVars[i].element_abundances[6] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[4]; // Si
  ChimesGasVars[i].element_abundances[7] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[5]; // S
  ChimesGasVars[i].element_abundances[9] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[6]; // Fe

  ChimesGasVars[i].dust_ratio *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDustRatio;
#endif // CHIMES_METAL_DEPLETION

  /* The element abundances may have changed, so use the check_constraint_equations()
   * routine to update the abundance arrays accordingly. If metals have been injected
   * into a particle, it is distributed across all of the metal's atomic/ionic/molecular
   * species, preserving the ion and molecule fractions of that element. */

  check_constraint_equations(&(ChimesGasVars[i]), &ChimesGlobalVars);
}
#endif // COOL_METAL_LINES_BY_SPECIES

#ifdef CHIMES_TURB_DIFF_IONS
/* mode == 0: re-compute the CHIMES abundance array from the ChimesNIons and ChimesHtot
 *            arrays that are used to track turbulent diffusion of ions and molecules.
 * mode == 1: update the ChimeSNIons array from the current CHIMES abundance array.
 */
void chimes_update_turbulent_abundances(int i, int mode)
{
  int k_species;
  double NHtot = (1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1])) * (P[i].Mass * UNIT_MASS_IN_CGS) / PROTONMASS;

  if (mode == 0)
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	ChimesGasVars[i].abundances[k_species] = (ChimesFloat) (SphP[i].ChimesNIons[k_species] / NHtot);
    }
  else
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	SphP[i].ChimesNIons[k_species] = ((double) ChimesGasVars[i].abundances[k_species]) * NHtot;
    }
}
#endif // CHIMES_TURB_DIFF_IONS

#ifdef CHIMES_METAL_DEPLETION
void chimes_init_depletion_data(void)
{
  int i;

#ifdef _OPENMP
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(maxThreads * sizeof(struct Chimes_depletion_data_structure));
#else
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(sizeof(struct Chimes_depletion_data_structure));
#endif

  // Elements in Jenkins (2009) in the order
  // C, N, O, Mg, Si, P, S, Cl, Ti, Cr, Mn, Fe,
  // Ni, Cu, Zn, Ge, Kr
  // Solar abundances are as mass fractions,
  // taken from the Cloudy default values, as
  // used in CHIMES.
  double SolarAbund[DEPL_N_ELEM] = {2.07e-3, 8.36e-4, 5.49e-3, 5.91e-4, 6.83e-4, 7.01e-6, 4.09e-4, 4.72e-6, 3.56e-6, 1.72e-5, 1.12e-5, 1.1e-3, 7.42e-5, 7.32e-7, 1.83e-6, 2.58e-7, 1.36e-7};

  // Fit parameters, using equation 10 of Jenkins (2009).
  // Where possible, we take the updated fit parameters
  // A2 and B2 from De Cia et al. (2016), which we convert
  // to the Ax and Bx of J09 (with zx = 0). Otherwise, we
  // use the original J09 parameters. We list these in the
  // order Ax, Bx, zx.
  double DeplPars[DEPL_N_ELEM][3] = {{-0.101, -0.193, 0.803}, // C
				     {0.0, -0.109, 0.55},     // N
				     {-0.101, -0.172, 0.0},     // O
				     {-0.412, -0.648, 0.0},     // Mg
				     {-0.426, -0.669, 0.0},     // Si
				     {-0.068, -0.091, 0.0},      // P
				     {-0.189, -0.324, 0.0},     // S
				     {-1.242, -0.314, 0.609}, // Cl
				     {-2.048, -1.957, 0.43},  // Ti
				     {-0.892, -1.188, 0.0},      // Cr
				     {-0.642, -0.923, 0.0},      // Mn
				     {-0.851, -1.287, 0.0},     // Fe
				     {-1.490, -1.829, 0.599}, // Ni
				     {-0.710, -1.102, 0.711}, // Cu
				     {-0.182, -0.274, 0.0},       // Zn
				     {-0.615, -0.725, 0.69},  // Ge
				     {-0.166, -0.332, 0.684}}; // Kr

  int n_thread;
#ifdef _OPENMP
  n_thread = maxThreads;
#else
  n_thread = 1;
#endif

  for (i = 0; i < n_thread; i++)
    {
      memcpy(ChimesDepletionData[i].SolarAbund, SolarAbund, DEPL_N_ELEM * sizeof(double));
      memcpy(ChimesDepletionData[i].DeplPars, DeplPars, DEPL_N_ELEM * 3 * sizeof(double));
      // DustToGasSaturated is the dust to gas ratio when F_star = 1.0, i.e. at maximum depletion onto grains.
      ChimesDepletionData[i].DustToGasSaturated = 5.9688e-03;
    }
}

/* Computes the linear fits as in Jenkins (2009)
 * or De Cia et al. (2016). Note that this returns
 * log10(fraction in the gas phase) */
double chimes_depletion_linear_fit(double nH, double T, double Ax, double Bx, double zx)
{
  // First, compute the parameter F_star, using the
  // best-fit relation from Fig. 16 of Jenkins (2009).
  double F_star = 0.772 + (0.461 * log10(nH));
  double depletion;

  // Limit F_star to be no greater than unity
  if (F_star > 1.0)
    F_star = 1.0;

  // All metals are in the gas phase if T > 10^6 K
  if (T > 1.0e6)
    return 0.0;
  else
    {
      depletion = Bx + (Ax * (F_star - zx));

      // Limit depletion to no greater than 0.0 (remember: it is a log)
      if (depletion > 0.0)
	return 0.0;
      else
	return depletion;
    }
}

void chimes_compute_depletions(double nH, double T, int thread_id)
{
  int i;
  double pars[DEPL_N_ELEM][3];
  memcpy(pars, ChimesDepletionData[thread_id].DeplPars, DEPL_N_ELEM * 3 * sizeof(double));

  // ChimesDepletionFactors array is for the metals
  // in the order C, N, O, Mg, Si, S, Fe. The other
  // metals in CHIMES are not depleted.
  ChimesDepletionData[thread_id].ChimesDepletionFactors[0] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[0][0], pars[0][1], pars[0][2])); // C
  ChimesDepletionData[thread_id].ChimesDepletionFactors[1] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[1][0], pars[1][1], pars[1][2])); // N
  ChimesDepletionData[thread_id].ChimesDepletionFactors[2] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[2][0], pars[2][1], pars[2][2])); // O
  ChimesDepletionData[thread_id].ChimesDepletionFactors[3] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[3][0], pars[3][1], pars[3][2])); // Mg
  ChimesDepletionData[thread_id].ChimesDepletionFactors[4] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[4][0], pars[4][1], pars[4][2])); // Si
  ChimesDepletionData[thread_id].ChimesDepletionFactors[5] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[6][0], pars[6][1], pars[6][2])); // S
  ChimesDepletionData[thread_id].ChimesDepletionFactors[6] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[11][0], pars[11][1], pars[11][2])); // Fe

  // The dust abundance as used in CHIMES will be the
  // metallicity in solar units multiplied by ChimesDustRatio
  ChimesDepletionData[thread_id].ChimesDustRatio = 0.0;
  for (i = 0; i < DEPL_N_ELEM; i++)
    ChimesDepletionData[thread_id].ChimesDustRatio += ChimesDepletionData[thread_id].SolarAbund[i] * (1.0 - pow(10.0, chimes_depletion_linear_fit(nH, T, pars[i][0], pars[i][1], pars[i][2])));

  // The above sum gives the dust to gas mass ratio.
  // Now normalise it by the saturated value.
  ChimesDepletionData[thread_id].ChimesDustRatio /= ChimesDepletionData[thread_id].DustToGasSaturated;
}
#endif // CHIMES_METAL_DEPLETION


void chimes_gizmo_exit(void)
{
  endrun(56275362);
}
#endif // CHIMES

#endif // ADM
#endif // COOLING