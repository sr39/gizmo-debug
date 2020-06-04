#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif

/* Routines for models that require stellar evolution: luminosities, mass loss, SNe rates, etc. 
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifdef GALSF


/* return the light-to-mass ratio [in units of Lsun/Msun] of a star or stellar population with a given age; used throughout the code below */
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i)
{
#ifdef SINGLE_STAR_SINK_DYNAMICS // SINGLE-STAR VERSION: calculate single-star luminosity (and convert to solar luminosity-to-mass ratio, which this output assumes)
    double m0=P[i].Mass; if(P[i].Type == 5) {m0=P[i].BH_Mass;}
    return calculate_individual_stellar_luminosity(0, m0, i) / m0 * (UNIT_LUM_IN_SOLAR) / (UNIT_MASS_IN_SOLAR);

#else // STELLAR-POPULATION VERSION: compute integrated mass-to-light ratio of an SSP
    double lum=1; if(stellar_age_in_gyr < 0.01) {lum=1000;} // default to a dumb imf-averaged 'young/high-mass' vs 'old/low-mass' distinction
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION // fit to updated SB99 tracks: including rotation, new mass-loss tracks, etc.
    if(stellar_age_in_gyr < 0.0035) {lum=1136.59;} else {double log_age=log10(stellar_age_in_gyr/0.0035); lum=1500.*pow(10.,-1.8*log_age+0.3*log_age*log_age-0.025*log_age*log_age*log_age);}
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
    double t1=0.0012, t2=0.0037, f1=800., f2=1100.*pow(Z_for_stellar_evol(i),-0.1), tx=log10(stellar_age_in_gyr/t2), t_g=log10(stellar_age_in_gyr/1.2)/0.05;
    if(stellar_age_in_gyr<=t1) {lum=f1;} else if(stellar_age_in_gyr<=t2) {lum=f1*pow(stellar_age_in_gyr/t1,log(f2/f1)/log(t2/t1));} else {lum=f2*pow(10.,-1.82*tx+0.42*tx*tx-0.07*tx*tx*tx)*(1.+1.2*exp(-0.5*t_g*t_g));}
#endif
#endif
    if(stellar_age_in_gyr<0.033) {lum*=calculate_relative_light_to_mass_ratio_from_imf(stellar_age_in_gyr,i);} // account for IMF variation model [if used]
    return lum;
#endif
}


/* subroutine to calculate luminosity of an individual star, according to accretion rate, 
    mass, age, etc. Modify your assumptions about main-sequence evolution here. ONLY relevant for SINGLE-STAR inputs. */
double calculate_individual_stellar_luminosity(double mdot, double mass, long i)
{
#if !defined(SINGLE_STAR_SINK_DYNAMICS)
    return 0; /* not defined */
#endif
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) && (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2) /* this is pre-calculated, simply return it */
    return P[i].StarLuminosity_Solar / UNIT_LUM_IN_SOLAR;
#endif
    /* if above flags not defined, estimate accretion + main-sequence luminosity as simply as possible */
    double lum=0, lum_sol=0, c_code = C_LIGHT_CODE, m_solar = mass * UNIT_MASS_IN_SOLAR;
    double rad_eff_protostar = 5.0e-7; /* if below the deuterium burning limit, just use the potential energy efficiency at the surface of a jupiter-density object */
    if(m_solar < 0.012) {rad_eff_protostar = 5.e-8 * pow(m_solar/0.00095,2./3.);}
    lum = rad_eff_protostar * mdot * c_code*c_code;
    if(m_solar >= 0.012) /* now for pre-main sequence and main sequence, need to also check the mass-luminosity relation */
    {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;}
        else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
        else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);}
        else {lum_sol = 32000. * m_solar;}
    }
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) && (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1) // now, account for pre-main sequence evolution and calculate accretion luminosity using protostellar radius
    if(i > 0) {if(P[i].Type == 5) {
        double eps_protostar=1.0, T4000_4 = pow(m_solar , 0.55), l_kh = 0.2263 * P[i].ProtoStellarRadius_inSolar*P[i].ProtoStellarRadius_inSolar * T4000_4; // protostellar temperature along Hayashi track and luminosity from KH contraction
        lum = DMAX(lum_sol,l_kh) / UNIT_LUM_IN_SOLAR + eps_protostar * (All.G * P[i].Mass / (P[i].ProtoStellarRadius_inSolar / UNIT_LENGTH_IN_SOLAR)) * mdot; // assume GM/r liberated per unit mass. Note we need radius in code units here since everything else in 'lum' is code-units as well. for pre-ms evolution, if Hayashi-temp luminosity exceeds MS luminosity, use it. otherwise use main sequence luminosity, and assume the star is moving along the Henyey track
        P[i].StarLuminosity_Solar = lum * UNIT_LUM_IN_SOLAR; //store total luminosity of the star in solar units
    }}
#endif
    return lum;

}


/* return the light-to-mass ratio, for the IMF of a given particle, relative to the Chabrier/Kroupa IMF.
    ONLY relevant for STELLAR POPULATION integrated inputs. */
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i)
{
#ifdef GALSF_SFR_IMF_VARIATION // fitting function from David Guszejnov's IMF calculations (ok for Mturnover in range 0.01-100) for how mass-to-light ratio varies with IMF shape/effective turnover mass 
    double log_mimf = log10(P[i].IMF_Mturnover);
    return (0.051+0.042*(log_mimf+2)+0.031*(log_mimf+2)*(log_mimf+2)) / 0.31;
#endif
#ifdef GALSF_SFR_IMF_SAMPLING // account for IMF sampling model if not evolving individual stars
    double mu = 0.01 * P[i].Mass * UNIT_MASS_IN_SOLAR; // 1 O-star per 100 Msun
    if(stellar_age_in_gyr > 0.003) {mu *= 0.326 * (0.003 / stellar_age_in_gyr);} // expectation value is declining with time, so 'effective multiplier' is larger
    return P[i].IMF_NumMassiveStars / mu;
#endif
    return 1; // Chabrier or Kroupa IMF //
}


#if defined(GALSF_FB_FIRE_RT_HIIHEATING) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
/* routine to compute the -ionizing- luminosity coming from either individual stars or an SSP */
double particle_ionizing_luminosity_in_cgs(long i)
{
#ifdef SINGLE_STAR_SINK_DYNAMICS /* SINGLE STAR VERSION: use effective temperature as a function of stellar mass and size to get ionizing photon production */
    double l_sol=bh_lum_bol(0,P[i].Mass,i)*(UNIT_LUM_IN_SOLAR), m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR, r_sol=pow(m_sol,0.738); // L/Lsun, M/Msun, R/Rsun
    double T_eff=5780.*pow(l_sol/(r_sol*r_sol),0.25), x0=157800./T_eff, fion=0; // ZAMS effective temperature; x0=h*nu/kT for nu>13.6 eV; fion=fraction of blackbody emitted above x0
    if(x0 < 30.) {double q=18./(x0*x0) + 1./(8. + x0 + 20.*exp(-x0/10.)); fion = exp(-1./q);} // accurate to <10% for a Planck spectrum to x0>30, well into vanishing flux //
    return fion * l_sol * SOLAR_LUM; // return value in cgs, as desired for this routine [l_sol is in L_sun, by definition above] //

#else /* STELLAR POPULATION VERSION: use updated SB99 tracks: including rotation, new mass-loss tracks, etc. */
    
    if(P[i].Type != 5)
    {
        double lm_ssp=0, star_age=evaluate_stellar_age_Gyr(P[i].StellarAge), t0=0.0035, tmax=0.02;
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
        tmax=0.15; lm_ssp=evaluate_light_to_mass_ratio(star_age,i); if(star_age<t0) {lm_ssp*=0.5;} else {lm_ssp*=0.5*pow(star_age/t0,-2.9);} /* slightly revised fit scales simply with Lbol [easier to modify]; see same references for stellar wind mass-loss rates; and extends to later ages (though most comes out at <100 Myr) */
#else
        if(star_age < t0) {lm_ssp=500.;} else {double log_age=log10(star_age/t0); lm_ssp=470.*pow(10.,-2.24*log_age-4.2*log_age*log_age) + 60.*pow(10.,-3.6*log_age);}
        lm_ssp *= calculate_relative_light_to_mass_ratio_from_imf(star_age, i);
#endif
        if(star_age >= tmax) {return 0;} // skip since old stars don't contribute
        return lm_ssp * SOLAR_LUM * (P[i].Mass*UNIT_MASS_IN_SOLAR); // converts to cgs luminosity [lm_ssp is in Lsun/Msun, here]
    } // (P[i].Type != 5)
#ifdef BH_HII_HEATING /* AGN template: light-to-mass ratio L(>13.6ev)/Mparticle in Lsun/Msun, above is dNion/dt = 5.5e54 s^-1 (Lbol/1e45 erg/s) */
    if(P[i].Type == 5) {return 0.18 * bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i) * UNIT_LUM_IN_CGS;}
#endif

#endif
    return 0; // catch
}
#endif



/* this routine tells the feedback algorithms what to 'inject' when a stellar feedback event occurs.
    you must define the mass, velocity (which defines the momentum and energy), and metal content (yields)
    of the ejecta for the event[s] of interest. Mass [Msne] and velocity [SNe_v_ejecta] should
    be in code units. yields[k] should be defined for all metal species [k], and in dimensionless units
    (mass fraction of the ejecta in that species). */
void particle2in_addFB_fromstars(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration)
{
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    if(fb_loop_iteration == 0) {particle2in_addFB_SNe(in,i);}
    if(fb_loop_iteration == 1) {particle2in_addFB_winds(in,i);}
    if(fb_loop_iteration == 2) {particle2in_addFB_Rprocess(in,i);}
    return;
#endif
    if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    // 'dummy' example model assumes all SNe are identical with IMF-averaged properties from the AGORA model (Kim et al., 2016 ApJ, 833, 202)
    in->Msne = P[i].SNe_ThisTimeStep * (14.8/UNIT_MASS_IN_SOLAR); // assume every SNe carries 14.8 solar masses (IMF-average)
    in->SNe_v_ejecta = 2607. / UNIT_VEL_IN_KMS; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun], which is IMF-averaged
#ifdef SINGLE_STAR_SINK_DYNAMICS // if single-star exploding or returning mass, use its actual mass & assumed energy to obtain the velocity
    in->Msne = DMIN(1.,P[i].SNe_ThisTimeStep) * P[i].Mass; // mass fraction of star being returned this timestep
    in->SNe_v_ejecta = sqrt(2.*(1.e51/UNIT_ENERGY_IN_CGS)/P[i].Mass); // for SNe [total return], simple v=sqrt(2E/m)should be fine without relativistic corrections
    if(P[i].SNe_ThisTimeStep<1) {double m_msun=P[i].Mass*UNIT_MASS_IN_SOLAR; in->SNe_v_ejecta = (616. * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / UNIT_VEL_IN_KMS;} // scaling from size-mass relation+eddington factor, assuming line-driven winds //
#endif
#ifdef METALS
    int k; for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
    if(NUM_METAL_SPECIES>=10) {in->yields[1] = 0.4;} // (catch for Helium, which the above scaling would give bad values for)
#endif
#endif
}


/* this routine calculates the event rates for different types of mechanical/thermal feedback
    algorithms. things like SNe rates, which determine when energy/momentum/mass are injected, should go here.
    you can easily modify this to accomodate any form of thermal or mechanical feedback/injection of various
    quantities from stars. */
double mechanical_fb_calculate_eventrates(int i, double dt)
{
#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION) // FIRE-specific stellar population version: separate calculation for SNe, stellar mass loss, R-process injection //
    double R_SNe = mechanical_fb_calculate_eventrates_SNe(i,dt);
    mechanical_fb_calculate_eventrates_Winds(i,dt);
    mechanical_fb_calculate_eventrates_Rprocess(i,dt);
    return R_SNe;
#endif
    
#ifdef SINGLE_STAR_SINK_DYNAMICS /* SINGLE-STAR version: simple implementation of single-star wind mass-loss and SNe rates */
    double m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR, l_sol=bh_lum_bol(0,P[i].Mass,i)*UNIT_LUM_IN_SOLAR;
#ifdef SINGLE_STAR_FB_WINDS
    double gam=DMIN(0.5,3.2e-5*l_sol/m_sol), alpha=0.5+0.4/(1.+16./m_sol), q0=(1.-alpha)*gam/(1.-gam), k0=1./30.; // Eddington factor (~L/Ledd for winds), capped at 1/2 for sanity reasons, approximate scaling for alpha factor with stellar type (weak effect)
    P[i].SNe_ThisTimeStep = DMIN(0.5, (2.338 * alpha * pow(l_sol,7./8.) * pow(m_sol,0.1845) * (1./q0) * pow(q0*k0,1./alpha) / m_sol) * (dt*UNIT_TIME_IN_GYR)); // Castor, Abbot, & Klein scaling
#endif
#ifdef SINGLE_STAR_FB_SNE
    double t_lifetime_Gyr = 10.*(m_sol/l_sol) + 0.003; /* crude estimate of main-sequence lifetime, capped at 3 Myr*/
    if(evaluate_stellar_age_Gyr(P[i].StellarAge) >= t_lifetime_Gyr) {P[i].SNe_ThisTimeStep=1;}
#endif
    return 1;
#endif
    
#ifdef GALSF_FB_THERMAL /* STELLAR-POPULATION version: pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly */
    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
    if(evaluate_stellar_age_Gyr(P[i].StellarAge) < 0.005) {return 0;} // enforce age limit of 5 Myr
    P[i].SNe_ThisTimeStep = P[i].Mass*UNIT_MASS_IN_SOLAR / 91.; // 1 event per 91 solar masses
    return 1;
#endif
    
#ifdef GALSF_FB_MECHANICAL /* STELLAR-POPULATION version: mechanical feedback: 'dummy' example model below assumes a constant SNe rate for t < 30 Myr, then nothing. experiment! */
    double RSNe=0, star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age < 0.03)
    {
        RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
        double p = RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
    }
    return RSNe;
#endif

    return 0;
}



#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION)
/* functions below contain pre-calculation of event rates and energetics, masses, etc, for FIRE mechanical feedback modules */

double mechanical_fb_calculate_eventrates_SNe(int i, double dt) 
{
#if defined(SINGLE_STAR_SINK_DYNAMICS) && (!defined(SINGLE_STAR_FB_SNE) || defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION)) /* no single-star module to use here, for these flags its in the spawn routine */
    return 0;
#endif
    if(All.SNe_Energy_Renormalization <= 0) return 0;
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), RSNe=0, agemin=0.003401, agebrk=0.01037, agemax=0.03753; // some ages in Gyr used below
    /* here we are determining an expected SNe rate, so SNe occur stochastically but with an age dependence in the population */
    // calculate: NSNe/Myr *if* each SNe had exactly 10^51 ergs; really from the energy curve [do this so we are gauranteed to get the correct SNe energy] //
#if (GALSF_FB_FIRE_STELLAREVOLUTION == 1) || (GALSF_FB_FIRE_STELLAREVOLUTION == 2)
    if(star_age>agemin) {if(star_age<=agebrk) {RSNe=5.408e-4;} else {if(star_age<=agemax) {RSNe=2.516e-4;}}} // core-collapse rate [super-simple 2-piece constant //
    if(star_age>agemax) {RSNe=5.3e-8 + 1.6e-5*exp(-0.5*((star_age-0.05)/0.01)*((star_age-0.05)/0.01));} // Ia (prompt Gaussian+delay, Manucci+06)
#elif (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
    agemin=0.0037; agebrk=0.7e-2; agemax=0.044; double f1=3.9e-4, f2=5.1e-4, f3=1.8e-4;
    if(star_age<agemin) {RSNe=0;} else if(star_age<=agebrk) {RSNe=f1*pow(star_age/agemin,log(f2/f1)/log(agebrk/agemin));}
        else if(star_age<=agemax) {RSNe=f2*pow(star_age/agebrk,log(f3/f2)/log(agemax/agebrk));} else {RSNe=0;} // core-collapse; updated with same stellar evolution models for wind mass loss [see there for references]. simple 2-part power-law provides extremely-accurate fit. models predict a totally negligible metallicity-dependence.
    double t_Ia_min=agemax, norm_Ia=1.6e-3; // t_Ia_min = delay time to first Ia, in Gyr; norm_Ia = Hubble-time integrated number of Ia's per solar mass
    if(star_age>t_Ia_min) {RSNe += norm_Ia * 7.94e-5 * pow(star_age,-1.1) / fabs(pow(t_Ia_min/0.1,-0.1) - 0.61);} // Ia DTD following Maoz & Graur 2017, ApJ, 848, 25
#endif
    
    double renorm = calculate_relative_light_to_mass_ratio_from_imf(star_age,i); // account for higher # of O-stars with a different IMF
    if(star_age<agemax) {RSNe *= renorm;}
#ifdef GALSF_SFR_IMF_SAMPLING
    if(star_age>agemax && P[i].IMF_NumMassiveStars>0) {RSNe += 2.516e-4*renorm;} // account for un-exploded O-stars
#endif
    double p = RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
    double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
#ifdef GALSF_SFR_IMF_SAMPLING
    if(star_age<agemax && P[i].IMF_NumMassiveStars<n_sn_0) {n_sn_0=P[i].IMF_NumMassiveStars;} // limit to number of O-stars for SNe //
#endif
    P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
    return RSNe;
}


void mechanical_fb_calculate_eventrates_Rprocess(int i, double dt)
{
#ifdef GALSF_FB_FIRE_RPROCESS
    /* we'll use the maximum rate here, then in the -yields- setting, 'cut' each down to its sub-population */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age>=0.003) // rate is zero at <3e6 yr
    {
        double p = 3.0e-5 / (1000.*star_age) * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // Nsne/Myr for a 1 Msun population
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID + 7) < p) {n_sn_0++;} // if > 1, this cuts that part off so we get appropriate n > 1 solution
        P[i].RProcessEvent_ThisTimeStep = n_sn_0; // assign event
    }
#endif
}


void mechanical_fb_calculate_eventrates_Winds(int i, double dt)
{
    double D_RETURN_FRAC=1.e-15, p=0;
    
#if defined(SINGLE_STAR_FB_WINDS) /* SINGLE-STAR VERSION: single-star wind mass-loss rates */
    double fire_wind_rel_mass_res = 1e-4; //relative mass resolution of winds, essentially the wind will get spawned in packets of fire_wind_rel_mass_res*(gas_mass_resolution) mass
    D_RETURN_FRAC = fire_wind_rel_mass_res * (2.0*All.MinMassForParticleMerger)/ P[i].Mass;
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION /* for 'fancy' multi-stage modules, have a separate subroutine to compute this */
    if(P[i].wind_mode != 2) {return;} /* only some eligible particles have winds in this module */
    p = single_star_wind_mdot(i,0) * dt / P[i].Mass; /* actual mdot from its own subroutine, given in code units */
#else /* otherwise use standard scaling from e.g. Castor, Abbot, & Klein */
    double m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR, l_sol=bh_lum_bol(0,P[i].Mass,i)*UNIT_LUM_IN_SOLAR; /* luminosity in solar */
    double gam=DMIN(0.5,3.2e-5*l_sol/m_sol), alpha=0.5+0.4/(1.+16./m_sol), q0=(1.-alpha)*gam/(1.-gam), k0=1./30.; // Eddington factor (~L/Ledd for winds), capped at 1/2 for sanity reasons, approximate scaling for alpha factor with stellar type (weak effect)
    p = (2.338 * alpha * pow(L_sol,7./8.) * pow(M_sol,0.1845) * (1./q0) * pow(q0*k0,1./alpha) / m_sol) * (dt*UNIT_TIME_IN_GYR); // mdot in M_sun/Gyr, times dt
#endif
    p=1.-exp(-p);
    
#else /* STELLAR POPULATION VERSION: now do stellar-population-averaged inputs */
    
    D_RETURN_FRAC = 0.01; // fraction of particle mass to return on a recycling step //
    if(All.StellarMassLoss_Rate_Renormalization <= 0) {return;}
    double star_age=evaluate_stellar_age_Gyr(P[i].StellarAge), ZZ=Z_for_stellar_evol(i);
#if (GALSF_FB_FIRE_STELLAREVOLUTION == 1)
    if(star_age<=0.001){p=4.76317;} else {if(star_age<=0.0035){p=4.76317*pow(10.,1.838*(0.79+log10(ZZ))*(log10(star_age)-(-3.00)));} else {
        if(star_age<=0.1){p=29.4*pow(star_age/0.0035,-3.25)+0.0041987;} else {p=0.41987*pow(star_age,-1.1)/(12.9-log(star_age));}}} // normalized  to give expected return fraction from stellar winds alone (~17%)
#elif (GALSF_FB_FIRE_STELLAREVOLUTION == 2)
    if(star_age<=0.001){p=4.76317*ZZ;} else {if(star_age<=0.0035){p=4.76317*ZZ*pow(10.,1.838*(0.79+log10(ZZ))*(log10(star_age)-(-3.00)));} else {
        if(star_age<=0.1){p=29.4*pow(star_age/0.0035,-3.25)+0.0041987;} else {p=0.41987*pow(star_age,-1.1)/(12.9-log(star_age));}}} // normalized  to give expected return fraction from stellar winds alone (~17%)
#elif (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
    /* updated fit. separates the more robust line-driven winds [massive-star-dominated] component, and -very- uncertain AGB. extremely good fits to updated STARBURST99 result for a 3-part Kroupa IMF (0.3,1.3,2.3 slope, 0.01-0.08-0.5-100 Msun, 8-120 SNe/BH cutoff, wind model evolution, Geneva v40 [rotating, Geneva 2013 updated tracks, at all metallicities available, ~0.1-1 solar], sampling times 1e4-2e10 yr at high resolution */
    double f1=3.*pow(ZZ,0.87), f2=20.*pow(ZZ,0.45), f3=0.6*ZZ, t1=0.0017, t2=0.004, t3=0.02, t=star_age; /* fit parameters for 'massive star' mass-loss */
    if(t<=t1) {p=f1;} else if(t<=t2) {p=f1*pow(t/t1,log(f2/f1)/log(t2/t1));} else if(t<=t3) {p=f2*pow(t/t2,log(f3/f2)/log(t3/t2));} else {p=f3*pow(t/t3,-3.1);} /* piecewise continuous function linking constant early and rapid late decay */
    double f_agb=0.01, t_agb=1.; p += f_agb/((1. + pow(t/t_agb,1.1)) * (1. + 0.01/(t/t_agb))); /* add AGB component. note that essentially no models [any of the SB99 geneva or padova tracks, or NuGrid, or recent other MESA models] predict a significant dependence on metallicity (that shifts slightly when the 'bump' occurs, but not the overall loss rate), so this term is effectively metallicity-independent */
#endif
    if(star_age < 0.1) {p *= calculate_relative_light_to_mass_ratio_from_imf(star_age,i);} // late-time independent of massive stars
    p *= All.StellarMassLoss_Rate_Renormalization * (dt*UNIT_TIME_IN_GYR); // fraction of particle mass expected to return in the timestep //
    p = 1.0 - exp(-p); // need to account for p>1 cases //
#endif
    
    double n_wind_0=(double)floor(p/D_RETURN_FRAC); p-=n_wind_0*D_RETURN_FRAC; // if p >> return frac, should have > 1 event, so we inject the correct wind mass
    P[i].MassReturn_ThisTimeStep += n_wind_0*D_RETURN_FRAC; // add this in, then determine if there is a 'remainder' to be added as well
    if(get_random_number(P[i].ID + 5) < p/D_RETURN_FRAC) {P[i].MassReturn_ThisTimeStep += D_RETURN_FRAC;} // add the 'remainder' stochastically
}






void particle2in_addFB_Rprocess(struct addFB_evaluate_data_in_ *in, int i)
{
#ifdef GALSF_FB_FIRE_RPROCESS
    if(P[i].RProcessEvent_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    int k; double star_age=evaluate_stellar_age_Gyr(P[i].StellarAge), p=get_random_number(P[i].ID + 8), pcrit, tcrit;
    for(k=0;k<NUM_RPROCESS_SPECIES;k++)
    {
        if(k==0) {tcrit=0.03;  pcrit=0.3333333333;}  // model 0: tmin > 3e7 yr, rate = 1e-5
        if(k==1) {tcrit=0.003; pcrit=0.3333333333;}  // model 1: tmin > 3e6 yr, rate = 1e-5
        if(k==2) {tcrit=0.03;  pcrit=0.03333333333;} // model 2: tmin > 3e7 yr, rate = 1e-6
        if(k==3) {tcrit=0.003; pcrit=0.03333333333;} // model 3: tmin > 3e6 yr, rate = 1e-6
        if((star_age>=tcrit)&&(p<=pcrit)) {in->yields[NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES+k]=1;} // units irrelevant//
    }
    in->Msne = 0.01 * (double)P[i].RProcessEvent_ThisTimeStep / ((double)(UNIT_MASS_IN_SOLAR)); // mass ejected ~0.01*M_sun; only here for bookkeeping //
#endif
}



void particle2in_addFB_SNe(struct addFB_evaluate_data_in_ *in, int i)
{
    int k; if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    int SNeIaFlag=0; if(evaluate_stellar_age_Gyr(P[i].StellarAge) > 0.03753) {SNeIaFlag=1;}; /* assume SNe before critical time are core-collapse, later are Ia */
    double Msne=10.5; if(SNeIaFlag) {Msne=1.4;} // average ejecta mass for single event (normalized to give total mass loss correctly)
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
    Msne=8.72; if(SNeIaFlag) {Msne=1.4;} // updated table of SNe rates and energetics, this is the updated mean mass per explosion to give the correct total SNe mass
#endif
    double SNeEgy = All.SNe_Energy_Renormalization*P[i].SNe_ThisTimeStep * 1.0e51/UNIT_ENERGY_IN_CGS; // assume each SNe has 1e51 erg
#ifdef METALS
    double yields[NUM_METAL_SPECIES]={0};
    if(NUM_METAL_SPECIES>=10) {
        // All, then He,C,N,O,Ne,Mg,Si,S,Ca,Fe
        if(SNeIaFlag) {
            /* SNIa */ /* from Iwamoto et al. 1999; 'W7' models */
            yields[0]=1.4;/* total metal mass */
            yields[1]=0.0;/*He*/ yields[2]=0.049;/*C*/ yields[3]=1.2e-6;/*N*/ yields[4]=0.143;/*O*/
            yields[5]=0.0045;/*Ne*/ yields[6]=0.0086;/*Mg*/ yields[7]=0.156;/*Si*/
            yields[8]=0.087;/*S*/ yields[9]=0.012;/*Ca*/ yields[10]=0.743;/*Fe*/
        } else {
            /* SNII (IMF-averaged... may not be the best approx on short timescales..., Nomoto 2006 (arXiv:0605725) */
            yields[0]=2.0;/*Z [total metal mass]*/
            yields[1]=3.87;/*He*/ yields[2]=0.133;/*C*/ yields[3]=0.0479;/*N*/ yields[4]=1.17;/*O*/
            yields[5]=0.30;/*Ne*/ yields[6]=0.0987;/*Mg*/ yields[7]=0.0933;/*Si*/
            yields[8]=0.0397;/*S*/ yields[9]=0.00458;/*Ca*/ yields[10]=0.0741;/*Fe*/
            // metal-dependent yields:
            if(P[i].Metallicity[0]<0.033) {yields[3]*=P[i].Metallicity[0]/All.SolarAbundances[0];} else {yields[3]*=1.65;} // N scaling is strongly dependent on initial metallicity of the star //
            yields[0] += yields[3]-0.0479; // correct total metal mass for this correction //
        }
    }
    if(NUM_METAL_SPECIES==3 || NUM_METAL_SPECIES==4)
    {
        if(SNeIaFlag) {
            yields[0]=1.4; yields[1]=0.0086; yields[2]=0.743; // All Z, Mg, Fe in total mass (SnIa)
        } else {
            yields[0]=2.0; yields[1]=0.12; yields[2]=0.0741; // SnII (per-SNe IMF-weighted averages)
        }
    }
    if(NUM_METAL_SPECIES==1) {if(SNeIaFlag) {yields[0]=1.4;} else {yields[0]=2.0;}}
    for(k=0;k<NUM_METAL_SPECIES;k++) {yields[k]=yields[k]/Msne;} // normalize to mass fraction //
    /* add a check to allow for larger abundances in the progenitor stars (usually irrelevant) */
    for(k=0;k<NUM_METAL_SPECIES;k++) {yields[k]=yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]);}
    if(SNeIaFlag) {if(NUM_METAL_SPECIES>=10) {yields[1]=0.0;}} // no He yield for Ia SNe //
    for(k=0;k<NUM_METAL_SPECIES;k++) {if(yields[k]<0) {yields[k]=0.0;} if(yields[k]>1) {yields[k]=1;} in->yields[k]=yields[k];}
#endif
    in->Msne = P[i].SNe_ThisTimeStep * (Msne/UNIT_MASS_IN_SOLAR); // total mass in code units
#ifdef SINGLE_STAR_SINK_DYNAMICS
    in->Msne = P[i].Mass; // conserve mass and destroy star completely
#endif
    in->SNe_v_ejecta = sqrt(2.0*SNeEgy/in->Msne); // v_ej in code units
}



void particle2in_addFB_winds(struct addFB_evaluate_data_in_ *in, int i)
{
    int k; if(P[i].MassReturn_ThisTimeStep<=0) {in->Msne=0; return;} // no event
#ifdef METALS
    /* assume track initial metallicity; turn on COOL_METAL_LINES_BY_SPECIES for more detailed tracking of light elements */
    double yields[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) {yields[k]=P[i].Metallicity[k];} // return surface abundances, to leading order //
    if(NUM_METAL_SPECIES>=10)
    {
        /* All, then He,C,N,O,Ne,Mg,Si,S,Ca,Fe ;; follow AGB/O star yields in more detail for the light elements */
        /*   the interesting species are He & CNO: below is based on a compilation of van den Hoek & Groenewegen 1997, Marigo 2001, Izzard 2004 */
        yields[1]=0.36; /*He*/ yields[2]=0.016; /*C*/ yields[3]=0.0041; /*N*/ yields[4]=0.0118; /*O*/
        // metal-dependent yields: O scaling is strongly dependent on initial metallicity of the star //
        if(P[i].Metallicity[0]<0.033) {yields[4] *= P[i].Metallicity[0]/All.SolarAbundances[0];} else {yields[4] *= 1.65;}
        for(k=1;k<=4;k++) {yields[k]=yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]); if(yields[k]<0) {yields[k]=0.0;} if(yields[k]>1) {yields[k]=1;} in->yields[k]=yields[k];} // enforce yields obeying pre-existing surface abundances, and upper/lower limits //
        yields[0]=0.0; for(k=2;k<NUM_METAL_SPECIES;k++) {yields[0]+=yields[k];}
    } else {
        yields[0]=0.032; for(k=1;k<NUM_METAL_SPECIES;k++) {yields[k]=0.0;}
    }
    for(k=0;k<NUM_METAL_SPECIES;k++) in->yields[k]=yields[k];
#endif
    in->Msne = P[i].Mass * P[i].MassReturn_ThisTimeStep; // mass (in code units) returned
    
    /* STELLAR POPULATION-AVERAGED VERSION: calculate wind kinetic luminosity + internal energy (hot winds from O-stars, slow from AGB winds) */
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && (GALSF_FB_FIRE_STELLAREVOLUTION <= 2)
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), E_wind_tscaling=0.0013;
    if(star_age <= 0.1) {E_wind_tscaling=0.0013 + 16.0/(1+pow(star_age/0.0025,1.4)+pow(star_age/0.01,5.0));} // stellar population age dependence of specific wind energy, in units of an effective internal energy/temperature
    in->SNe_v_ejecta = sqrt(2.0 * (All.StellarMassLoss_Energy_Renormalization * E_wind_tscaling * (3.0e7/((5./3.-1.)*U_TO_TEMP_UNITS)))); // get the actual wind velocity (multiply specific energy by units, user-set normalization, and convert)
#elif (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // ??
    double t=evaluate_stellar_age_Gyr(P[i].StellarAge), Z=Z_for_stellar_evol(i), f0=pow(Z,0.12); /* setup: updated fit here uses the same stellar evolution models/tracks as used to compute mass-loss rates. see those for references here. */
    in->SNe_v_ejecta = sqrt(All.StellarMassLoss_Energy_Renormalization) * f0 * (3000./(1.+pow(t/0.003,2.5)) + 600./(1.+pow(sqrt(Z)*t/0.05,6.)+pow(Z/0.2,1.5)) + 30.) / UNIT_VEL_IN_KMS; /* interpolates smoothly from OB winds through AGB, also versus Z */
#endif


#if defined(SINGLE_STAR_FB_WINDS) /* SINGLE-STAR VERSION: instead of a stellar population, this is wind from a single star */
    double m_msun=P[i].Mass*UNIT_MASS_IN_SOLAR;
    in->SNe_v_ejecta = (616. * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / UNIT_VEL_IN_KMS; // scaling from size-mass relation+eddington factor, assuming line-driven winds //
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION)
    in->SNe_v_ejecta = single_star_wind_velocity(i); /* for fancy models, wind velocity in subroutine, based on Leitherer 1992 and stellar evolutions tage, size, etc. */
#endif
#endif

}


double Z_for_stellar_evol(int i)
{
    if(i<0) {return 1;}
    double Z_solar = P[i].Metallicity[0]/All.SolarAbundances[0]; // use total metallicity
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) && defined(COOL_METAL_LINES_BY_SPECIES) // ??
    int i_Fe=10; Z_solar = P[i].Metallicity[i_Fe]/All.SolarAbundances[i_Fe]; // use Fe, specifically, for computing stellar properties, as its most relevant here. MAKE SURE this is set to the correct abundance in the list, to match Fe!!!
#endif
    return DMIN(DMAX(Z_solar,0.01),3.);
}

#endif // GALSF_FB_MECHANICAL+GALSF_FB_FIRE_STELLAREVOLUTION

    


#ifdef SINGLE_STAR_FB_JETS
double single_star_jet_velocity(int n){
    /*Calculates the launch velocity of jets*/
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) /* use the fancy stellar evolution modules to calculate these for stars or protostars */
    return (All.BAL_f_launch_v * sqrt(All.G * P[n].BH_Mass / (P[n].ProtoStellarRadius_inSolar / UNIT_LENGTH_IN_SOLAR)) * All.cf_atime); // we use the flag as a multiplier times the Kepler velocity at the protostellar radius. Really we'd want v_kick = v_kep * m_accreted / m_kicked to get the right momentum
#else
    return (All.BAL_f_launch_v * sqrt(All.G * P[n].BH_Mass / (10. / UNIT_LENGTH_IN_SOLAR)) * All.cf_atime); // we use the flag as a multiplier times the Kepler velocity at the protostellar radius. Really we'd want v_kick = v_kep * m_accreted / m_kicked to get the right momentum; without a better guess, assume fiducial protostellar radius of 10*Rsun, as in Federrath 2014
#endif
}
#endif


#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) /* begins large block of 'fancy' protostar-through-MS stellar evolution models */

/* 'master' function to update the size (and other properties like effective temperature) of accreting protostars along relevant stellar evolution tracks */
double singlestar_subgrid_protostellar_evolution_update_track(int n, double dm, double dt)
{
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1)
    /* this is the simple version written by Phil: intentionally simplified PS evolution tracks, designed to make it easy to understand and model the evolution and reduce un-necessary complications */
    double lum_sol = 0.0, m_initial = DMAX(1.e-37 , (BPP(n).BH_Mass - dm)), mu = DMAX(0, dm/m_initial), m_solar = BPP(n).BH_Mass*UNIT_MASS_IN_SOLAR, T4000_4 = pow(m_solar, 0.55); // m_initial = mass before the accretion, mu = relative mass accreted, m_solar = mass in solar units, T4000_4 = (temperature/4000K)^4 along Hayashi track
    if(m_solar > 0.012) // below this limit, negligible luminosity //
    {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;} else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
          else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);} else {lum_sol = 32000. * m_solar;}
    }
    double R_Hayashi_Henyey = 2.1 * sqrt(lum_sol / T4000_4); // size below which, at the temperature above, contraction must occur along the Henyey track at constant luminosity
    double t_R_evol = 0, contraction_factor = 0; // timescale for contraction
    if(BPP(n).ProtoStellarRadius_inSolar <= R_Hayashi_Henyey)
    {// currently on Henyey track, contracting at constant Luminosity
        t_R_evol = 1.815e7 * m_solar*m_solar / (BPP(n).ProtoStellarRadius_inSolar * lum_sol) / (UNIT_TIME_IN_YR); // contraction timescale
        contraction_factor = 1. / (1. + dt/t_R_evol);
    } else {// currently on Hayashi track, contracting at constant Temperature
        t_R_evol = 8.021e7 * m_solar*m_solar / (BPP(n).ProtoStellarRadius_inSolar*BPP(n).ProtoStellarRadius_inSolar*BPP(n).ProtoStellarRadius_inSolar * T4000_4) / (UNIT_TIME_IN_YR); // contraction timescale
        contraction_factor = 1. / pow(1 + 3.*dt/t_R_evol, 1./3.);
    }
    double r_new = DMAX(10. * m_solar, 5.24 * pow(m_solar, 1./3)), R_main_sequence_ignition; // r_new = size of newly-formed protostar; R_main_sequence_ignition = main sequence radius - where contraction should halt
    if (m_solar < 0.012) {r_new =  5.24 * pow(m_solar, 1./3);} // constant density for Jupiter-type objects
    BPP(n).ProtoStellarRadius_inSolar = (BPP(n).ProtoStellarRadius_inSolar * contraction_factor + r_new * mu) / (1. + mu); // new size (contraction + accretion both accounted for)
    if(m_solar <= 1) {R_main_sequence_ignition = pow(m_solar,0.8);} else {R_main_sequence_ignition = pow(m_solar,0.57);}
    if(BPP(n).ProtoStellarRadius_inSolar <= R_main_sequence_ignition)
    {
        BPP(n).ProtoStellarRadius_inSolar = R_main_sequence_ignition; BPP(n).ProtoStellarStage = 5; //using same notation for MS as SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1
#ifdef SINGLE_STAR_PROMOTION
        P[n].Type = 4; P[n].StellarAge = All.Time; P[n].Mass = DMAX(P[n].Mass , BPP(n).BH_Mass + BPP(n).BH_Mass_AlphaDisk); // convert type, mark the new ZAMS age according to the current time, and accrete remaining mass
#endif
    }

#elif (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2) /* Protostellar evolution model based on the ORION version, see Offner 2009 Appendix B */
    
    double frad = 0.18; //limit for forming radiative barrier, based on Offner+MckKee 2011 source code
    double fk = 0.5; //fraction of kinetic energy that is radiated away in the inner disk before reaching the surface, using default ORION value here as it is not a GIZMO input parameter
    double f_acc = 0.5; //fraction of accretion power that is radiated away instead of being used to drive winds, using default ORION value here as it is not a GIZMO input parameter
#ifdef SINGLE_STAR_FB_JETS //We can convert the GIZMO parameters into their ORION versions, FWIND=(1.0-All.BAL_f_accretion) and .
    f_acc = ( (1.0-All.BAL_f_accretion)*(1.0-All.BAL_f_launch_v*All.BAL_f_launch_v) + All.BAL_f_accretion )/All.BAL_f_accretion; //The energy not used to drive winds must be radiated away (neglecting thermal energy). It is from the leftover from the jet material and the material not launched. But we need to be careful because the f_acc here only sees the modified mdot = All.BAL_f_accretion * mdot_now, so we need to divide with All.BAL_f_accretion. For our nominal All.BAL_f_accretion=0.7 and BAL_f_launch_v=0.3 this yields 1.39 (0.97 for the full mdot)
#endif
    double max_rel_dr = 0.01; //Maximum relative change in radius per step, if the change over a single timestep is larger than this than we subcycle the evolution of the stellar radius
    double mass = BPP(n).BH_Mass; //mass of star/protostar at the end of the timestep
    double mass_D = BPP(n).Mass_D; //amount of D in the protostar
    double mdot = BPP(n).BH_Mdot; //accretion rate, shorter to write it this way
    double mdot_m_solar_per_year = mdot * UNIT_MASS_IN_SOLAR / UNIT_TIME_IN_YR; // accretion rate in msolar/yr
    double m_solar = mass * UNIT_MASS_IN_SOLAR; // mass in units of Msun
    double m_initial = DMAX(1.e-37 , (mass - dm)); // mass before accretion
    int stage = BPP(n).ProtoStellarStage; /*what stage of stellar evolution the particle is in 0: pre collapse, 1: no burning, 2: fixed Tc burnig, 3: variable Tc burning, 4: shell burning, 5: main sequence, see Offner 2009 Appendix B*/
    if (stage == 0) {//set the radius for the pre-collapse phase according to Eq B1 in Offner 2009, this overwrites the original prescription from sfr_eff.c
        BPP(n).ProtoStellarRadius_inSolar = DMAX(2.5 * pow(mdot_m_solar_per_year*1e5,0.33),2.0); //radius always at least 2 R_sun
    }
    double r = BPP(n).ProtoStellarRadius_inSolar / UNIT_LENGTH_IN_SOLAR; // star radius in code units
    int stage_increase = 0;
    double lum_Hayashi = ps_lum_Hayashi_BB(mass, r); //blackbody radiation assuming the star follows the Hayashi track
    double lum_MS = ps_lum_MS(mass); //luminosity of main sequence star of m mass
    double lum_int = DMAX(lum_Hayashi, lum_MS+(f_acc*fk*All.G*mass*mdot/r )); //luminosity from the stellar interior
    double lum_I = ps_lum_I(mdot); //luminosity needed to ionize the accreted material
    if (stage < 5){ //not a main sequence star
        if (stage >= 1){ //We only evolve those that are beyond the pre-collapse phase
            int loop_subcycle=0;
            int n_subcycle=1; //no subcycle by default
            double dm_curr = dm; double dt_curr = dt;
            mass = m_initial; //value at the beginning o the timestep
            double n_ad, ag, rhoc, Pc, Tc, beta, dlogbeta_dlogm, lum_D, dm_D, rel_dr; //need to declare them here or the compiler reduces their scope to the loop
            do{//we use a do-while loop here because we first try the full timestep and if dr is too large we subcycle
            //Note: this subcycling does not update stage, so a protostar could overshoot
                //Evolve mass
                mass += dm_curr; //increase mass first
                double dm_rel = dm_curr/(mass-dm_curr);
                //Get properties for stellar evolution
                lum_Hayashi = ps_lum_Hayashi_BB(mass, r); //blackbody radiation assuming the star follows the Hayashi track
                lum_MS = ps_lum_MS(mass); //luminosity of main sequence star of m mass
                lum_int = DMAX(lum_Hayashi, lum_MS+(f_acc*fk*All.G*mass*mdot/r )); //luminosity from the stellar interior
                n_ad = ps_adiabatic_index(stage, mdot); //get adiabatic index. Note: ORION does not seem to update this, but I think it is worthwhile as mdot can vary over time
                ag = 3.0/(5.0-n_ad); //shorthand
                rhoc = ps_rhoc(mass, n_ad, r); //central density
                Pc = ps_Pc(mass, n_ad, r); //central pressure
                Tc = ps_Tc(rhoc,Pc); //central temperature
                beta = ps_beta(mass, n_ad, rhoc, Pc); //mean ratio of gas pressure to total pressure
                dlogbeta_dlogm = ps_dlogbeta_dlogm(mass, r, n_ad, beta, rhoc, Pc); // d log beta/ d log m
                //Calculate luminosity from D burning
                lum_D = 0; //luminosity from D burning
                dm_D = dm_curr; //by default we burn no D (stage 1)
                if (stage==2){ //burning at fixed Tc, lum_D set to keep the central temperature constant
                    double dlogbetaperbetac_dlogm = ps_dlogbetaperbetac_dlogm(mass, r, n_ad, beta, rhoc, Pc, Tc); // ratio of gas pressure to total pressure at the center
                    lum_D = lum_int + lum_I + (All.G*mass*mdot/r) * ( 1.-fk-0.5*ag*beta * (1.+dlogbetaperbetac_dlogm) ); // Eq B8 of Offner 2009
                    //Change in available deuterium mass
                    dm_D = dm_curr - dt_curr * (lum_D*UNIT_LUM_IN_SOLAR/15.) * (1e-5) / (UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR) ;
                }
                else{ if (stage>2){
                    //burning all accreted D for stages above 2
                    lum_D = (15./UNIT_LUM_IN_SOLAR) * (mdot_m_solar_per_year/(1e-5));
                    dm_D = 0; //all new D is burned
                    mass_D = 0; //no D left in protostar
                    }
                }                
                //Let's evolve the stellar radius
                rel_dr = ( dm_rel * (1.-(1.-fk)/(ag*beta)+0.5*dlogbeta_dlogm) - dt_curr/(ag*beta)*r/(All.G*mass*mass) * (lum_int+lum_I-lum_D) ); //Eq B4 of Offner 2009 divided by r, and corrected by a factor of 2 as per the ORION source used for Offner+McKee 2011
                //Let's check if we need to subcycle
                if (fabs(rel_dr) > max_rel_dr){
                    n_subcycle = (int) DMAX(ceil(rel_dr/max_rel_dr), 2.0*n_subcycle); //number of subcycle steps, at least 2, either double the previous number or estimated from dr
                    //reset protostar properties, restart loop
                    loop_subcycle = 0;
                    mass = m_initial; mass_D = BPP(n).Mass_D; 
                    dm_curr = dm/((double)n_subcycle); dt_curr = dt/((double)n_subcycle);
                    r = BPP(n).ProtoStellarRadius_inSolar / UNIT_LENGTH_IN_SOLAR;
                }
                else{
                    loop_subcycle++;
                    r *= (1.0+rel_dr);
                    mass_D += dm_D;
                }
            } while(loop_subcycle<n_subcycle); //repeat for the number of subcycle steps
            //Update stellar properties
            BPP(n).ProtoStellarRadius_inSolar = r * UNIT_LENGTH_IN_SOLAR;
            BPP(n).Mass_D = mass_D;
            //Debug message
#ifdef PS_EVOL_OUTPUT_MOREINFO
            if (n_subcycle>1){
                dm_D = mass_D - BPP(n).Mass_D; //get the tota change in D mass in the protostar
                rel_dr = r/(BPP(n).ProtoStellarRadius_inSolar/UNIT_LENGTH_IN_SOLAR)-1.0; //get the actual relative change over the whole timestep
            }
            printf("PS evolution t: %g sink ID: %u mass: %g radius_solar: %g stage: %d mdot_m_solar_per_year: %g mD: %g rel_dr: %g dm: %g dm_D: %g Tc: %g Pc: %g rhoc: %g beta: %g dt: %g n_ad: %g lum_int: %g lum_I: %g lum_D: %g age_Myr: %g StarLuminosity_Solar: %g BH_Mass_AlphaDisk: %g SinkRadius: %g dlogbeta_dlogm: %g n_subcycle: %d.ZAMS_Mass %g PS_end\n",All.Time, P[n].ID,m_solar,BPP(n).ProtoStellarRadius_inSolar,stage, mdot_m_solar_per_year, BPP(n).Mass_D*UNIT_MASS_IN_SOLAR,rel_dr,dm*UNIT_MASS_IN_SOLAR, dm_D*UNIT_MASS_IN_SOLAR, Tc, Pc*UNIT_PRESSURE_IN_CGS, rhoc*UNIT_DENSITY_IN_CGS, beta, dt*UNIT_TIME_IN_MYR, n_ad, lum_int*UNIT_LUM_IN_SOLAR, lum_I*UNIT_LUM_IN_SOLAR, lum_D*UNIT_LUM_IN_SOLAR, (All.Time-P[n].ProtoStellarAge)*UNIT_TIME_IN_MYR, BPP(n).StarLuminosity_Solar, BPP(n).BH_Mass_AlphaDisk, BPP(n).SinkRadius, dlogbeta_dlogm, n_subcycle, P[n].ZAMS_Mass );
#endif
            //Check whether the star can progress to the next state
            //Move from "no burn" to "burning at fixed Tc" phase when central temperature gets high enough for D ignition
            if ( (stage==1) && (Tc >= 1.5e6) && ((All.Time-BPP(n).StellarAge) > DMAX(3.*dt, 1e-4/(UNIT_TIME_IN_MYR)) ) ){ //further check that the sink has been promoted at least a couple of timesteps and 100 yr ago, so that we don't start D burning immediately after forming the sink (relevant in low res cases)
                stage_increase = 1;//particle qualifies to the "fixed Tc burn" phase
            }
            //Move from "burning at fixed Tc" to "variable Tc burn" phase when D runs out
            if ( (stage==2) && (BPP(n).Mass_D <= 0) ){
                BPP(n).Mass_D = 0;
                stage_increase = 1;//particle qualifies to the "variable Tc burn" phase
            }
            //Move from "variable Tc burn" to "shell burn" phase when radiation becomes strong enough to form a radiative zone
            if ( (stage==3) && ((lum_D/lum_MS) < frad) ){
                stage_increase = 1;//particle qualifies to the "shell burn" phase
                BPP(n).ProtoStellarRadius_inSolar *= 2.1; //star swells due to the formation of radiative barrier
            }
            //Move from "shell burn" to "main sequence" phase when the radius reaches the main sequence radius
            if ( (stage==4) && (BPP(n).ProtoStellarRadius_inSolar <= ps_radius_MS_in_solar(mass)) ){
                stage_increase = 1;//particle qualifies to become a ZAMS star
                BPP(n).ProtoStellarRadius_inSolar = ps_radius_MS_in_solar(mass);
            }
        }
        else{ //the protostar is in the "pre-collapse" state, no internal evolution, just check if it can be promoted to the next stage
            BPP(n).Mass_D = BPP(n).BH_Mass; //no D burned so far
            if (m_solar >= 0.01){ 
            stage_increase = 1; //particle qualifies to the "no burning stage"
            } 
        }
        if (stage_increase){
            BPP(n).ProtoStellarStage += stage_increase;
            BPP(n).StellarAge = All.Time; //store the time of the last promotion
            if (BPP(n).ProtoStellarStage == 5){
                P[n].ZAMS_Mass = P[n].BH_Mass; //store the mass at which we reached the main sequence
            }
#ifdef PS_EVOL_OUTPUT_MOREINFO
            printf("%u promoted to %d \n",P[n].ID,(stage+stage_increase)); //Debug message
#endif
        } //increase evolutionary stage if the particle satisfies the requirements
    }
    else{ // for main sequence stars
        if (BPP(n).ProtoStellarStage == 5){ //MS stars
            BPP(n).ProtoStellarRadius_inSolar = ps_radius_MS_in_solar(mass); //update the mass if the mass changes (unlikely)
#ifdef SINGLE_STAR_FB_SNE
            //Check age and see if we need to blow up this star
            double age_Gyr = evaluate_stellar_age_Gyr(P[n].StellarAge);
            if ( age_Gyr > stellar_lifetime_in_Gyr(n) ){
                BPP(n).ProtoStellarStage = 6; //time to explode
                P[n].Mass_final = P[n].BH_Mass; //record the final mass the star had
#ifdef BH_ALPHADISK_ACCRETION
                BPP(n).BH_Mass_AlphaDisk = 0; //probably does not matter, but let's make sure these don't cause issues
                P[n].Mass = P[n].BH_Mass; 
#endif
                //Save properties of SN progenitor
                fprintf(FdBhSNDetails, "%g %u %g %g %g %g %g %g %g \n", All.Time, P[n].ID, P[n].BH_Mass, P[n].Pos[0], P[n].Pos[1],P[n].Pos[2],P[n].Vel[0], P[n].Vel[1],P[n].Vel[2]);
                fflush(FdBhSNDetails);
            }
#endif
        }
    }
    //Calculate the luminosity of the star
    /*********************************************/
    /* Power based parametrization from ORION (2 params)*/
    //lum_acc = facc * fk * All.G*mass*mdot/r; //luminosity radiated at the accretion shock
    //lum_disk = (1.-fk) * All.G*mass*mdot/r; //luminosity released by material that traverses the inner disk
    //BPP(n).StarLuminosity_Solar = (lum_acc + lum_disk + lum_int) * UNIT_LUM_IN_SOLAR ; //luminosity of the star
    /*********************************************/
    /* Mass flux based parametrization (1 params) */
    /* For our nominal choice of BAL_f_accretion=0.7 this gives very similar results to the ORION parameters of fk=facc=0.5, which is equivalent to 0.75*/
    double eps_protostar=0.75; //fraction of gas that does not get launched out with a jet, default value, although 1.0 would be energy conserving.
#ifdef SINGLE_STAR_FB_JETS
    eps_protostar=1.0; // since mdot is already modified by All.BAL_f_accretion
#endif
    eps_protostar -= f_acc*fk; //Need to deduct the part that is already accounted for in L_int (ORION uses the convention to add lum_acc to that)
#ifndef SINGLE_STAR_FB_DISABLE_MS_HEATING
    BPP(n).StarLuminosity_Solar = (eps_protostar*All.G*mass*mdot/r + lum_int) * UNIT_LUM_IN_SOLAR; //luminosity of the star in solar units
#else
    BPP(n).StarLuminosity_Solar = (eps_protostar*All.G*mass*mdot/r + lum_Hayashi) * UNIT_LUM_IN_SOLAR; //same as above but we don't count H burning for th emission. Thsi way the radial evolution follows the same track as with the full model, but we don't provide feedback from H burning to the nearby gas
#endif
#endif//end of SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2

#ifdef PS_EVOL_OUTPUT_MOREINFO // print out the basic star info
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2)
    if (BPP(n).ProtoStellarStage >= 5) //only for MS stars, for previous stages we will print out the properties before
#endif
    {printf("PS evolution t: %g sink ID: %u mass: %g radius_solar: %g stage: %d mdot_m_solar_per_year: %g mD: 0 rel_dr: 0 dm: %g dm_D: 0 Tc: 0 Pc: 0 rhoc: 0 beta: 0 dt: %g n_ad: 0 lum_int: 0 lum_I: 0 lum_D: 0 age_Myr: %g StarLuminosity_Solar: %g BH_Mass_AlphaDisk: %g SinkRadius: %g dlogbeta_dlogm: 0 n_subcycle: 0.ZAMS_Mass %g PS_end\n",All.Time, P[n].ID,BPP(n).BH_Mass*UNIT_MASS_IN_SOLAR,BPP(n).ProtoStellarRadius_inSolar,BPP(n).ProtoStellarStage, BPP(n).BH_Mdot*UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR , dm*UNIT_MASS_IN_SOLAR, dt*UNIT_TIME_IN_MYR, (All.Time-P[n].ProtoStellarAge)*UNIT_TIME_IN_MYR, BPP(n).StarLuminosity_Solar, BPP(n).BH_Mass_AlphaDisk*UNIT_MASS_IN_SOLAR, BPP(n).SinkRadius, P[n].ZAMS_Mass );}
#endif
}

#if defined(SINGLE_STAR_FB_WINDS)
/* Let's get the wind mass loss rate for MS stars. n is the index of the particle (P[n]). mode is 1 when called by the wind spawning routine (blackhole.c) and 2 if called by the FIRE wind module (in this file, mechanical_fb_calculate_eventrates_Winds). The function decides which type of wind feedback is appropriate for the current star and will only give a nonzero mdot to one of these */
double single_star_wind_mdot(int n, int set_mode){ //if set_mode is zero then the wind mode is not changed by calling this function
    double minimum_stellarmass_for_winds_solar  = 2.0;  // minimum stellar mass allowed to have winds
    int    model_wolf_rayet_phase_explicitly    = 1;    // assumes that O stars turn into WR stars at the end of their lifetime, increasing their mass loss rate
    double n_particles_for_discrete_wind_spawn  = 1e-2; // parameter for switching between wind spawning and just depositing momentum to nearby gas (FIRE winds) -- particle number required to trigger 'explicit' spawn module. Setting it to 0 ensures that we always spawn winds, while a high value (e.g. 1e6) ensures we always use the FIRE wind module
    double spawning_min_wind_jet_mom_ratio = 10.0; //If winds are much more powerful than jets ( (wind momentum injection/jet momentum injection) > this value) then we can safely spawn the winds and neglect the jets if we want to
    
    double wind_mass_loss_rate=0; //mass loss rate in code units
    if (P[n].Type != 5) {return 0;}
    if(BPP(n).ProtoStellarStage != 5) {return 0;}
    double m_solar = BPP(n).Mass * UNIT_MASS_IN_SOLAR; // mass in units of Msun
    if (m_solar < minimum_stellarmass_for_winds_solar){return 0.0;} //no winds for low mass stars
    //Winds are for MS only: we are assuming that METALS are also on
    double ZZ = BPP(n).Metallicity[0]/All.SolarAbundances[0]; //relative metallicity to solar
    double logmdot_wind; // log10(Mdot / (Msun/yr))

    // phenomenological prescription: "de Jager / 3" model from Smith 2014, with limiter for "weak-wind problem"
    logmdot_wind = -6 + 1.5 * log10(BPP(n).StarLuminosity_Solar / 1e6) + 0.69 * log10(ZZ); // "de Jager / 3"
    logmdot_wind = DMIN(-7.65 + 2.9*log10(BPP(n).StarLuminosity_Solar/ 1e5), logmdot_wind); // weak-wind problem
    wind_mass_loss_rate = pow(10.0,logmdot_wind) / (UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR); //reducing the rate to be more in line with observations, see Nathan Smith 2014, conversion to code units from Msun/yr
    
    if(model_wolf_rayet_phase_explicitly) {if(evaluate_stellar_age_Gyr(P[n].StellarAge) > (stellar_lifetime_in_Gyr(n)-singlestar_WR_lifetime_Gyr(n))){wind_mass_loss_rate*=10;}} //Our star is in the WR phase, for now use the simple prescription of having 10x higher wind loss rates based on Smith 2014
    //Let's deal with the case of undefined wind mode (just promoted to MS or restart from snapshot)
    if ( set_mode && (wind_mass_loss_rate>0) ){
        //Let's calculate N_wind = Mdot_wind * t_wind / dm_wind, where t_wind is solved from: Mdot_wind * t_wind = material swept up = 4/3 pi rho (v_wind*t_wind)^3
        double v_wind = single_star_wind_velocity(n);
        double t_wind =sqrt( wind_mass_loss_rate * (3.0/(4.0*M_PI*P[n].DensAroundStar)) / (v_wind*v_wind*v_wind));
        double N_wind = wind_mass_loss_rate * t_wind / All.BAL_wind_particle_mass;
        
        int old_wind_mode = P[n].wind_mode;
        if (N_wind >= n_particles_for_discrete_wind_spawn){
            P[n].wind_mode = 1; //we can spawn enough particles per wind time
        } else{
            P[n].wind_mode = 2; //we can't spawn enough particles per wind time, switching to FIRE wind module to reduce burstiness
        }
#ifdef SINGLE_STAR_FB_JETS
        if ( (P[n].wind_mode == 1) && (P[n].BH_Mdot>0) ){ //we want to spawn winds but we have jets too
            double jet_mom_inj = single_star_jet_velocity(n) * P[n].BH_Mdot;
            double wind_mom_inj = v_wind * wind_mass_loss_rate;
            if (spawning_min_wind_jet_mom_ratio < (spawning_min_wind_jet_mom_ratio * jet_mom_inj) ){ P[n].wind_mode = 2;} //we switch back to the FIRE wind injection so that we can spawn the jet and have winds at the same time
        }
#endif
        if (old_wind_mode != P[n].wind_mode){
            printf("Wind mode change for star %llu to %d at %g. Mdot_wind %g\n",P[n].ID,P[n].wind_mode,All.Time, wind_mass_loss_rate);
        }
    }
    return wind_mass_loss_rate;
}

double singlestar_WR_lifetime_Gyr(int n){ //Calculate lifetime for star in Wolf-Rayet Phase
    double m_solar = BPP(n).Mass * UNIT_MASS_IN_SOLAR; // mass in units of Msun
    double ZZ = P[n].Metallicity[0]/All.SolarAbundances[0]; //relative metallicity to solar
    if (m_solar<=20.0){return 0.;} //No WR phase below that
    //Using prescription based on Fig 7 from Meynet & Maeder 2005, all >10 Msun star spend the end of their lifetime as WR
    return DMAX(0., 1.5e-3 * DMIN(1., ((m_solar-20.)/80.)) * pow(ZZ,0.5) );
}

double single_star_wind_velocity(int n){
    /* Let's get the wind velocity for MS stars */
    double T_eff = 5814.33 * pow( P[n].StarLuminosity_Solar/(P[n].ProtoStellarRadius_inSolar*P[n].ProtoStellarRadius_inSolar), 0.25 ); //effective temperature in K
    double v_esc = (617.7 / UNIT_VEL_IN_KMS) * sqrt(P[n].Mass*UNIT_MASS_TO_SOLAR / P[n].ProtoStellarRadius_inSolar); // surface escape velocity - wind escape velocity should be O(1) factor times this, factor given below
    if(T_eff < 1.25e4){ return 0.7 * v_esc;}  // Lamers 1995
    else if (T_eff < 2.1e4) {return 1.3 * v_esc;}
    else {return 2.6 * v_esc;}        
}
#endif


double stellar_lifetime_in_Gyr(int n){ //Estimate lifetime of star, using simple MS approximation t/Gyr ~ 9.6 M/L in solar
    double m_solar = BPP(n).Mass * UNIT_MASS_IN_SOLAR; // mass in units of Msun
    return 9.6 * (m_solar / P[n].StarLuminosity_Solar) + 0.003;     // gives ~10Gyr for solar-type stars, ~40Myr for 8msun ZAMS, and asymptotes to 3Myr at high mass
}


#if defined(SINGLE_STAR_FB_SNE)
double single_star_SN_velocity(int n){ // Initial velocity of SNe ejecta: 10^51 erg/SN, distributed evenly among the mass
    return sqrt(SINGLE_STAR_FB_SNE * (2e51/UNIT_ENERGY_IN_CGS)/P[n].Mass_final); //simple v=sqrt(2E/m)should be fine without relativistic corrections
}

void single_star_SN_init_directions(void){
    /* routine to initialize the distribution of spawned wind particles during SNe. This is essentially a copy of the function rt_init_intensity_directions() in rt_utilities.c */
    int n_polar = SINGLE_STAR_FB_SNE_N_EJECTA_QUADRANT;
    if(n_polar < 1) {printf("Number of SN ejecta particles is invalid (<1). Terminating.\n"); endrun(53463431);}
    double mu[n_polar]; int i,j,k,l,n=0,n_oct=n_polar*(n_polar+1)/2,n_tot=8*n_oct;
    double SN_Ejecta_Direction_tmp[n_oct][3];
    for(j=0;j<n_polar;j++) {mu[j] = sqrt( (j + 1./6.) / (n_polar - 1./2.) );}
    for(i=0;i<n_polar;i++)
    {
        for(j=0;j<n_polar-i;j++)
        {
            k=n_polar-1-i-j;
            SN_Ejecta_Direction_tmp[n][0]=mu[i]; SN_Ejecta_Direction_tmp[n][1]=mu[j]; SN_Ejecta_Direction_tmp[n][2]=mu[k];
            n++;
        }
    }
    n=0;
    for(i=0;i<2;i++)
    {
        double sign_x = 1 - 2*i;
        for(j=0;j<2;j++)
        {
            double sign_y = 1 - 2*j;
            for(k=0;k<2;k++)
            {
                double sign_z = 1 - 2*k;
                for(l=0;l<n_oct;l++)
                {
                    All.SN_Ejecta_Direction[n][0] = SN_Ejecta_Direction_tmp[l][0] * sign_x;
                    All.SN_Ejecta_Direction[n][1] = SN_Ejecta_Direction_tmp[l][1] * sign_y;
                    All.SN_Ejecta_Direction[n][2] = SN_Ejecta_Direction_tmp[l][2] * sign_z;
                    n++;
                }
            }
        }
    }
}
#endif


#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2) /* Functions for protosteller evolution model based on Offner 2009 */
/* Calculate the mean ratio of the gas pressure to the gas+radiation pressure, either by solving the Eddington quartic (for n_ad=3, Eq B5) or by using tabulated values, based on Offner 2009, code taken from ORION */
double ps_beta(double m, double n_ad, double rhoc, double Pc) {
    double mass=m*UNIT_MASS_IN_SOLAR;//in units of solar mass
    double Pc_cgs=Pc*All.cf_a3inv*UNIT_PRESSURE_IN_CGS, rhoc_cgs=rhoc*All.cf_a3inv*UNIT_DENSITY_IN_CGS;
    if(n_ad==3.0) {// In this case we solve the Eddington quartic, P_c^3 = (3/a) (k / (mu mH))^4 (1 - beta) / beta^4 rho_c^4 for beta
        int JMAX=40, j; double BETAMIN=1.0e-4, BETAMAX=1.0, TOL=1.0e-7, dx, f, fmid, xmid, rtb, x1=BETAMIN, x2=BETAMAX, coef=3.0/7.56e-15*pow(BOLTZMANN*rhoc_cgs/(0.613*PROTONMASS),4);
        f = pow(Pc_cgs,3) - coef * (1.0-x1)/pow(x1,4); fmid = pow(Pc_cgs,3) - coef * (1.0-x2)/pow(x2,4); rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
        for (j=1;j<=JMAX;j++) {
          xmid=rtb+(dx *= 0.5); fmid = pow(Pc_cgs,3) - coef * (1.0-xmid)/pow(xmid,4);
          if(fmid <= 0.0) rtb=xmid;
          if(fabs(dx) < TOL*fabs(xmid) || fmid == 0.0) return rtb;
        }
        printf("ps_beta: bisection solve failed to converge"); return(-1);
    } else {
        // For n != 3, we use a table lookup. The values of beta have been pre-computed with mathematica. The table goes from M=5 to 50 solar masses in steps of 2.5 M_sun, and from n=1.5 to n=3 in steps of 0.5. We should never call this routine with M > 50 Msun, since by then the star should be fully on the main sequence.
        double MTABMIN=5.0, MTABMAX=50.0, MTABSTEP=2.5, NTABMIN=1.5, NTABMAX=3.0, NTABSTEP=0.5, MBETMIN=0.1 ;
        //if (mass < MBETMIN){return (1.0+ 0.25*log(mass/MBETMIN)/log(0.01/MBETMIN) );}  // Setting from Offner+Mckee2011, not sure why, does not make much sense above 1, probably to fit to previous results. I made it change continously to avoid big drops in R at 0.1 Msun, value adjusted from 1.15
        if (mass < MTABMIN){return (1.0);}  // Set beta = 1 for M < 5 Msun
        if ((mass >= MTABMAX) || (n_ad >= NTABMAX)) {printf("ps_beta: too high protostar mass, m: %g n_ad %g",m, n_ad); return(-1.0);}
        static double betatab[19][4] = {{0.98785, 0.988928, 0.98947, 0.989634}, {0.97438, 0.976428, 0.977462, 0.977774}, {0.957927, 0.960895, 0.962397, 0.962846},
          {0.939787, 0.943497, 0.945369, 0.945922}, {0.92091, 0.925151, 0.927276, 0.927896}, {0.901932, 0.906512, 0.908785, 0.909436}, {0.883254, 0.888017, 0.890353, 0.891013},
          {0.865111, 0.86994, 0.872277, 0.872927}, {0.847635, 0.852445, 0.854739, 0.855367}, {0.830886, 0.835619, 0.837842, 0.838441}, {0.814885, 0.8195, 0.821635, 0.822201},
          {0.799625, 0.804095, 0.806133, 0.806664}, {0.785082, 0.789394, 0.791328, 0.791825}, {0.771226, 0.775371, 0.777202, 0.777665}, {0.758022, 0.761997, 0.763726, 0.764156},
          {0.745433, 0.749238, 0.750869, 0.751268}, {0.733423, 0.73706, 0.738596, 0.738966}, {0.721954, 0.725429, 0.726874, 0.727216}, {0.710993, 0.714311, 0.715671, 0.715987}};
        // Locate ourselves on the table and do a linear interpolation
        int midx = (int) floor((mass-MTABMIN)/MTABSTEP), nidx = (int) floor((n_ad-NTABMIN)/NTABSTEP);
        double mwgt = (mass-(MTABMIN+midx*MTABSTEP)) / MTABSTEP, nwgt = (n_ad-(NTABMIN+nidx*NTABSTEP)) / NTABSTEP;
        return ( betatab[midx][nidx]*(1.0-mwgt)*(1.0-nwgt) + betatab[midx+1][nidx]*mwgt*(1.0-nwgt) + betatab[midx][nidx+1]*(1.0-mwgt)*nwgt + betatab[midx+1][nidx+1]*mwgt*nwgt );
    }
}
/* Sets the adiabatic index for pre burning protostars based on Eq B2 from Offner 2009 */
double ps_adiabatic_index_func(double mdot) {
    double mdot_m_solar_per_year = mdot * UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR; // accretion rate in msolar/yr
    //return ( 5.0 - 3/(1.475+0.07*log10(mdot_m_solar_per_year)) );
    if (mdot_m_solar_per_year < 1.0e-5){
        return (2.5);
    }else{
        return(2.5+0.25*log10(mdot_m_solar_per_year*1.0e5));
    }
}
/* Sets the adiabatic index for protostars based on Appendix B of Offner 2009 */
double ps_adiabatic_index(int stage, double mdot){
    double n_ad;
    switch(stage) {
        case 0: n_ad = ps_adiabatic_index_func(mdot); break;
        case 1: n_ad = ps_adiabatic_index_func(mdot); break;
        case 2: n_ad = 1.75; break;
        case 3: n_ad = 1.75; break;
        case 4: n_ad = 3.0; break;
        case 5: n_ad = 3.0; break;
    }
    if(n_ad > 3.0) {n_ad=3.0;}
    if(n_ad < 1.75) {n_ad=1.75;}
    return n_ad;
}
/* Calculate central temperature for protostar by solving Pc = rho_c*kb*Tc/(mu*mH)+1/3*a*Tc^4 using bisection, based on Offner 2009 Eq B14, code taken from ORION */
double ps_Tc(double rhoc, double Pc) {
    int JMAX=40; double TOL=1.0e-7; // max number of iterations, and error tolerance, respectively
    double Pc_cgs=Pc*All.cf_a3inv*UNIT_PRESSURE_IN_CGS, rhoc_cgs=rhoc*All.cf_a3inv*UNIT_DENSITY_IN_CGS, Tgas=Pc_cgs*0.613*PROTONMASS/(BOLTZMANN*rhoc_cgs), Trad=pow(3*Pc_cgs/7.56e-15, 0.25), dx, f, fmid, xmid, rtb, x1=0, x2; int j;
    x2 = (Trad > Tgas) ? 2*Trad : Tgas;
    f = Pc_cgs - rhoc_cgs*BOLTZMANN*x1/(0.613*PROTONMASS) - 7.56e-15*pow(x1,4)/3.0;
    fmid=Pc_cgs - rhoc_cgs*BOLTZMANN*x2/(0.613*PROTONMASS) - 7.56e-15*pow(x2,4)/3.0;
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=1;j<=JMAX;j++) {
        xmid=rtb+(dx *= 0.5);
        fmid = Pc_cgs - rhoc_cgs*BOLTZMANN*xmid/(0.613*PROTONMASS) - 7.56e-15*pow(xmid,4)/3.0;
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < TOL*fabs(xmid) || fmid == 0.0) return rtb;
    }
    printf("ps_Tc: bisection solve didn't converge, P_c = %e, rho_c = %e, Tgas = %e Trad = %e",Pc_cgs, rhoc_cgs, Tgas, Trad); return(-1);
}
/* Calculate central density for protostar using a pre-computed table for fixed mass, radius and polytropic index, based on Offner 2009, table and code taken from ORION */
double ps_rhoc(double m, double n_ad, double r) { /* Use Tabulated values of rho_c/rho_mean for n=1.5 to 3.1 in intervals of 0.1*/
    static double rhofactab[] = {0.166931, 0.14742, 0.129933, 0.114265, 0.100242, 0.0877, 0.0764968, 0.0665109, 0.0576198, 0.0497216, 0.0427224, 0.0365357, 0.0310837, 0.0262952, 0.0221057, 0.0184553, 0.01529};
    int itab = (int) floor((n_ad-1.5)/0.1); double wgt = (n_ad - (1.5 + 0.1*itab)) / 0.1, rhofac = rhofactab[itab]*(1.0-wgt) + rhofactab[itab+1]*wgt;
    return (m / (4.0/3.0*M_PI*r*r*r) / rhofac);
}
/* Calculate central pressure for protostar using a pre-computed table for fixed mass, radius and polytropic index, based on Offner 2009, table and code taken from ORION */
double ps_Pc(double m, double n_ad, double r) {
    static double pfactab[] = {0.770087, 0.889001, 1.02979, 1.19731, 1.39753, 1.63818, 1.92909, 2.2825, 2.71504, 3.24792, 3.90921, 4.73657, 5.78067, 7.11088, 8.82286, 11.0515, 13.9885};
    int itab = (int) floor((n_ad-1.5)/0.1); double wgt = (n_ad - (1.5 + 0.1*itab)) / 0.1, pfac = pfactab[itab]*(1.0-wgt) + pfactab[itab+1]*wgt;
    return (pfac * All.G * m*m/(r*r*r*r));
}
/* Calculate the mean ratio of the gas pressure to the gas+radiation pressure at the center, based on Offner 2009, code taken from ORION */
double ps_betac(double rhoc, double Pc, double Tc) {
    double Pc_cgs=Pc*All.cf_a3inv*UNIT_PRESSURE_IN_CGS, rhoc_cgs=rhoc*All.cf_a3inv*UNIT_DENSITY_IN_CGS;
    return (rhoc_cgs*BOLTZMANN*Tc/(0.613*PROTONMASS) / Pc_cgs);
}
/* Calculate dlog(beta)/dlog(m) by taking a numerical derivative, based on Offner 2009, code taken from ORION */
double ps_dlogbeta_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc) {
    double eps=0.01, epspone=1.+eps, rhoc2 = ps_rhoc( m*epspone , n_ad, r), Pc2 = ps_Pc( m*epspone , n_ad, r), beta2 = ps_beta( m*epspone, n_ad, rhoc2, Pc2); //slight imprecision here as we do not update the radius
    return (beta2-beta) / (eps*beta);
}
/* Calculate dlog(beta/betac)/dlog(m) by taking a numerical derivative, based on Offner 2009, code taken from ORION */
double ps_dlogbetaperbetac_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc, double Tc) {
    double eps=0.01, epspone=1.+eps, betac = ps_betac(rhoc, Pc, Tc), rhoc2 = ps_rhoc( m*epspone , n_ad, r), Pc2 = ps_Pc( m*epspone , n_ad, r); //slight imprecision here as we do not update the radius
    double Tc2 = ps_Tc(rhoc2, Pc2), beta2 = ps_beta( m*epspone, n_ad, rhoc2, Pc2), betac2 = ps_betac(rhoc2, Pc2, Tc2);
    return (betac*beta2/betac2 - beta) / (eps*beta);
}
/* Calculate the luminosity [in code units] required to ionize the infalling material, based on Offner 2009 */
double ps_lum_I(double mdot) {
    double mdot_m_solar_per_year = mdot * UNIT_MASS_IN_SOLAR/UNIT_TIME_IN_YR; // accretion rate in msolar/yr
    return (2.5/UNIT_LUM_IN_SOLAR) * ( mdot_m_solar_per_year/(1e-5));
}
/* Calculate the blackbody luminosity [in code units] of the star following the Hayashi track */
double ps_lum_Hayashi_BB(double m, double r) {
    double m_solar = m * UNIT_MASS_IN_SOLAR, T4000_4 = pow(m_solar , 0.55), r_solar = r * UNIT_LENGTH_IN_SOLAR; // protostellar temperature along Hayashi track
    T4000_4 = 0.316406; //ORION prescription T=3000K
    return 0.2263 * r_solar * r_solar * T4000_4 / UNIT_LUM_IN_SOLAR; // luminosity from KH contraction
}
/* Calculate the luminosity of a main sequence star in code units: ORION version: fitting formulas of Tout et al (1996) */
double ps_lum_MS(double m) {
    double m_solar = m * UNIT_MASS_IN_SOLAR;
    if(m_solar <= 0.1) {return 0;} // zero luminosity below some threshold for fits, otherwise use formula below //
    return ((0.39704170*pow(m_solar,5.5) + 8.52762600*pow(m_solar,11)) / (0.00025546+pow(m_solar,3)+5.43288900*pow(m_solar,5) + 5.56357900*pow(m_solar,7) + 0.78866060*pow(m_solar,8)+0.00586685*pow(m_solar,9.5))) / UNIT_LUM_IN_SOLAR;
}
/* Calculate the radius of a main sequence star in solar units: ORION version: fitting formulas of Tout et al (1996) */
double ps_radius_MS_in_solar(double m) {
    double m_solar = m * UNIT_MASS_IN_SOLAR;
    return (1.71535900*pow(m_solar,2.5)+6.59778800*pow(m_solar,6.5)+10.08855000*pow(m_solar,11)+1.01249500*pow(m_solar,19)+0.07490166*pow(m_solar,19.5)) /
        (0.01077422+3.08223400*pow(m_solar,2)+17.84778000*pow(m_solar,8.5)+pow(m_solar,18.5)+0.00022582*pow(m_solar,19.5));
}
#endif // (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 2) /* end functions for protosteller evolution model based on Offner 2009 */


#endif //end of protostellar evolution functions







#endif /* GALSF */
