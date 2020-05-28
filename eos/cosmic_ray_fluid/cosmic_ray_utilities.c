#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../../allvars.h"
#include "../../proto.h"

/*! Routines for cosmic ray 'fluid' modules (as opposed to the explicit CR-PIC methods, which are in the grain+particles section of the code)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef COSMIC_RAYS

/* routine which returns the typical absolute value of the rigidity of a given CR [in GV] in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_rigidity_in_GV(int target, int k_CRegy)
{
    double R = 1;
#if (N_CR_PARTICLE_BINS > 1)
    /* insert physics here */
#endif
    return R;
}

/* routine which returns the typical charge of CRs in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_charge_in_e(int target, int k_CRegy)
{
    double Z = 1;
#if (N_CR_PARTICLE_BINS > 1)
    /* insert physics here */
#endif
    return Z;
}

/* routine which determines the fraction of injected CR energy per 'bin' of CR energy. */
double CR_energy_spectrum_injection_fraction(int k_CRegy, int source_PType, double shock_vel)
{
    double f_bin = 1./N_CR_PARTICLE_BINS; /* uniformly distributed */
#if (N_CR_PARTICLE_BINS > 1)
    /* insert physics here */
#endif
    return f_bin;
}

/* routine which gives diffusion coefficient as a function of energy for the 'constant diffusion coefficient' models:
    current default: -extremely- simple power-law, assuming diffusion coefficient increases with CR energy per unit charge as (E/Z)^(1/2) */
double diffusion_coefficient_constant(int target, int k_CRegy)
{
    double dimensionless_kappa_relative_to_GV_protons = 1;
#if (N_CR_PARTICLE_BINS > 1)
    /* insert physics here */
    dimensionless_kappa_relative_to_GV_protons = pow( return_CRbin_CR_rigidity_in_GV(target,k_CRegy) , 0.5 );
#endif
    return All.CosmicRayDiffusionCoeff * dimensionless_kappa_relative_to_GV_protons;
}

/* cosmic ray interactions affecting the -thermal- temperature of the gas are included in the actual cooling/heating functions;
    they are solved implicitly above. however we need to account for energy losses of the actual cosmic ray fluid, here. The
    timescale for this is reasonably long, so we can treat it semi-explicitly, as we do here.
    -- We use the estimate for combined hadronic + Coulomb losses from Volk 1996, Ensslin 1997, as updated in Guo & Oh 2008: */
void CR_cooling_and_losses(int target, double n_elec, double nHcgs, double dtime_cgs)
{
    if(dtime_cgs <= 0) {return;} /* catch */
    int k,k_CRegy; double f_ion=DMAX(DMIN(Get_Gas_Ionized_Fraction(target),1.),0.);
    double a_hadronic = 6.37e-16, b_coulomb_per_GeV = 3.09e-16*(n_elec + 0.57*(1.-f_ion))*HYDROGEN_MASSFRAC; /* some coefficients; a_hadronic is the default coefficient, b_coulomb_per_GeV the default Coulomb+ionization (the two scale nearly-identically) normalization divided by GeV, b/c we need to divide the energy per CR  */
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double CR_coolrate=0, Z=return_CRbin_CR_charge_in_e(target,k_CRegy);
        if(Z > 0) /* protons here [note for now I'm using Z>0 as synonymous with protons, i.e. ignoring positrons, but we could include those if really desired */
        {
#if (N_CR_PARTICLE_BINS > 2) /* note these are currently energy-loss expressions; for truly multi-bin, probably better to work with dp/dt, instead of dE/dt */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), beta=return_CRbin_beta_factor(target,k_CRegy), R_CR_GV=return_CRbin_CR_rigidity_in_GV(target,k_CRegy);
            CR_coolrate += b_coulomb_per_GeV * ((Z*Z)/(beta*E_GeV)) * nHcgs; // all protons Coulomb-interact, can be rapid for low-E
            if(E_GeV>=0.78) {CR_coolrate += a_hadronic * nHcgs;} // only GeV CRs or higher trigger above threshold for collisions
#else
            CR_coolrate = (0.87*a_hadronic + 0.53*b_coulomb_per_GeV) * nHcgs; /* for N<=2, assume a universal spectral shape, the factor here corrects for the fraction above-threshold for hadronic interactions, and 0.53 likewise for averaging  */
#endif
        } else { /* electrons here: note for electrons and positrons, always in the relativistic limit, don't need to worry about beta << 1 limits */
            /* bremsstrahlung [folllowing Blumenthal & Gould, 1970]: dEkin/dt=4*alpha_finestruct*r_classical_elec^2*c * SUM[n_Z,ion * Z * (Z+1) * (ln[2*gamma_elec]-1/3) * E_kin */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), E_rest=0.000511, gamma=E_GeV/E_rest;
            CR_coolrate += n_elec * nHcgs * 1.39e-16 * DMAX(log(2.*gamma)-0.33,0);
            /* synchrotron and inverse compton scale as dE/dt=(4/3)*sigma_Thompson*c*gamma_elec^2*(U_mag+U_rad), where U_mag and U_rad are the magnetic and radiation energy densities, respectively. Ignoring Klein-Nishina corrections here, as they are negligible at <40 GeV and only a ~15% correction up to ~1e5 GeV */
            double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG, U_rad_ev = get_cell_Urad_in_eVcm3(target);
            CR_coolrate += 5.2e-20 * gamma * (U_mag_ev + U_rad_ev); // U_mag_ev=(B^2/8pi)/(eV/cm^(-3)), here; U_rad=U_rad/(eV/cm^-3) //
        }
        
        /* for now, cooling is being treated as energy loss 'within the bin'. with denser bins, should allow for movement -between- bins. needs to be implemented ?? */
        double q_CR_cool = exp(-CR_coolrate * dtime_cgs); if(CR_coolrate * dtime_cgs > 20.) {q_CR_cool = 0;}
        SphP[target].CosmicRayEnergyPred[k_CRegy] *= q_CR_cool; SphP[target].CosmicRayEnergy[k_CRegy] *= q_CR_cool;
#ifdef COSMIC_RAYS_M1
        for(k=0;k<3;k++) {SphP[target].CosmicRayFlux[k_CRegy][k] *= q_CR_cool; SphP[target].CosmicRayFluxPred[k_CRegy][k] *= q_CR_cool;}
#endif
    }
    return;
}



/* routine which gives diffusion coefficient as a function of CR bin for the self-confinement models [in local equilibrium]. mode sets what we assume about the 'sub-grid'
    parameters f_QLT (rescales quasi-linear theory) or f_cas (rescales turbulence strength)
      <=0: fQLT=1 [most naive quasi-linear theory, ruled out by observations],  fcas=1 [standard Goldreich-Shridar cascade]
        1: fQLT=100, fcas=1
        2: fQLT=1, fcas=100
        3: fQLT=1, fcas-K41 from Hopkins et al. 2020 paper, for pure-Kolmogorov isotropic spectrum
        4: fQLT=1, fcas-IK, IK spectrum instead of GS
   if set mode < 0, will also ignore the dust-damping contribution from Squire et al. 2020
 */
#ifndef COSMIC_RAYS_SET_SC_MODEL
#define COSMIC_RAYS_SET_SC_MODEL 1 /* set which mode to return from the SC subroutine here, of the various choices for how to e.g. model fCas, fQLT */
#endif
double diffusion_coefficient_self_confinement(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG,
    double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion)
{
    double vol_inv = SphP[target].Density*All.cf_a3inv / P[target].Mass, fturb_multiplier=1, f_QLT=1, R_CR_GV=return_CRbin_CR_rigidity_in_GV(target,k_CRegy), Z_charge_CR=return_CRbin_CR_charge_in_e(target,k_CRegy);
    double e_CR = SphP[target].CosmicRayEnergyPred[k_CRegy]*vol_inv, n_cgs=rho_cgs/PROTONMASS, cos_Bgrad=0,B2=0,P2=0,EPSILON_SMALL=1.e-50; int k;
#ifdef MAGNETIC
    for(k=0;k<3;k++) {double b0=SphP[target].BPred[k]*vol_inv*All.cf_a2inv, p0=SphP[target].Gradients.CosmicRayPressure[k_CRegy][k]; cos_Bgrad+=b0*p0; B2+=b0*b0; P2+=p0*p0;}
    cos_Bgrad/=sqrt(B2*P2+EPSILON_SMALL);
#else
    B2=e_CR; cos_Bgrad=1; for(k=0;k<3;k++) {double p0=SphP[target].Gradients.CosmicRayPressure[k_CRegy][k]; P2+=p0*p0;} /* this model doesn't really make sense without B-fields, but included for completeness here */
#endif
    double Omega_gyro=0.00898734*b_muG*(All.UnitTime_in_s/All.HubbleParam)/R_CR_GV, r_L=C_LIGHT_CODE/Omega_gyro, kappa_0=r_L*C_LIGHT_CODE;
    double x_LL = DMAX( C_LIGHT_CODE / (Omega_gyro * L_scale), EPSILON_SMALL ), CRPressureGradScaleLength=Get_CosmicRayGradientLength(target,k_CRegy), vA_code=vA_noion, k_turb=1./L_scale, k_L=1./r_L, x_EB_ECR=(0.5*B2+EPSILON_SMALL)/(e_CR+EPSILON_SMALL);
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    if(f_ion>0) {vA_code /= sqrt(f_ion);} // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    if(mode==1) {f_QLT = 100;} // multiplier to account for arbitrary deviation from QLT, applies to all damping mechanisms [100 = favored value in our study; or could use fcas = 100]
    fturb_multiplier = pow(M_A,3./2.); // multiplier to account for different turbulent cascade models (fcas = 1)
    if(mode==2) {fturb_multiplier *= 100.;} // arbitrary multiplier (fcas = 100, here)
    if(mode==3) {fturb_multiplier = pow(M_A,3./2.) * 1./(pow(M_A,1./2.)*pow(x_LL,1./6.));} // pure-Kolmogorov (fcas-K41)
    if(mode==4) {fturb_multiplier = pow(M_A,3./2.) / pow(x_LL,1./10.);} // GS anisotropic but perp cascade is IK (fcas-IK) /

    //if(M_A<1.) {fturb_multiplier*=DMIN(sqrt(M_A),pow(M_A,7./6.)/pow(x_LL,1./6.));} else {fturb_multiplier*=DMIN(1.,1./(pow(M_A,1./2.)*pow(x_LL,1./6.)));} // corrects to Alfven scale, for correct estimate according to Farmer and Goldreich, Lazarian, etc. /* Lazarian+ 2016 multi-part model for where the resolved scales lie on the cascade */
    /* ok now we finally have all the terms needed to calculate the various damping rates that determine the equilibrium diffusivity */
    double U0bar_grain=3., rhograin_int_cgs=1., fac_grain=R_CR_GV*sqrt(n_cgs)*U0bar_grain/(b_muG*rhograin_int_cgs), f_grainsize = DMAX(8.e-4*pow(fac_grain*(temperature/1.e4),0.25), 3.e-3*sqrt(fac_grain)); // b=2, uniform logarithmic grain spectrum over a factor of ~100 in grain size; f_grainsize = 0.07*pow(sqrt(fion*n1)*EcrGeV*T4/BmuG,0.25); // MRN size spectrum
    double G_dust = vA_code*k_L * (P[target].Metallicity[0]/0.014) * f_grainsize; // also can increase by up to a factor of 2 for regimes where charge collisionally saturated, though this is unlikely to be realized
    if(mode<0) {G_dust = 0;} // for this choice, neglect the dust-damping term 
    double G_ion_neutral = 5.77e-11 * n_cgs * (0.97*nh0 + 0.03*nHe0) * sqrt(temperature) * (All.UnitTime_in_s/All.HubbleParam); if(Z_charge_CR > 1) {G_ion_neutral /= sqrt(2.*Z_charge_CR);} // ion-neutral damping: need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]
    double G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI)*cs_thermal/4.) * sqrt(k_turb*k_L) * fturb_multiplier; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)
    double G0 = G_ion_neutral + G_turb_plus_linear_landau + G_dust; // linear terms all add into single G0 term
    double Gamma_effective = G0, phi_0 = (sqrt(M_PI)/6.)*(fabs(cos_Bgrad))*(1./(x_EB_ECR+EPSILON_SMALL))*(cs_thermal*vA_code*k_L/(CRPressureGradScaleLength*G0*G0 + EPSILON_SMALL)); // parameter which determines whether NLL dominates
    if(isfinite(phi_0) && (phi_0>0.01)) {Gamma_effective *= phi_0/(2.*(sqrt(1.+phi_0)-1.));} // this accounts exactly for the steady-state solution for the Thomas+Pfrommer formulation, including both the linear [Landau,turbulent,ion-neutral] + non-linear terms. can estimate (G_nonlinear_landau_effective = Gamma_effective - G0)
    /* with damping rates above, equilibrium transport is equivalent to pure streaming, with v_stream = vA + (diffusive equilibrium part) give by the solution below, proportional to Gamma_effective and valid to O(v^2/c^2) */
    double v_st_eff = vA_code * (1. + f_QLT * 4. * kappa_0 * Gamma_effective * x_EB_ECR * (1. + 2.*vA_code*vA_code/(C_LIGHT_CODE*C_LIGHT_CODE)) / (M_PI*vA_code*vA_code + EPSILON_SMALL)); // effective equilibrium streaming speed for all terms accounted
    return GAMMA_COSMICRAY * v_st_eff * CRPressureGradScaleLength; // convert to effective diffusivity from 'streaming'
}



/* routine which gives diffusion coefficient for extrinsic turbulence models. 'mode' sets whether we assume Alfven modes (mode<0), Fast-mode scattering (mode>0), or both (=0),
     0: 'default' Alfven + Fast modes (both, summing scattering rates linearly)
    -1: 'default' Alfven modes: correctly accounting for an anisotropic Goldreich-Shridar cascade, per Chandran 2000
    -2: Alfven modes in pure Goldreich-Shridhar cascade, ignoring anisotropic effects [*much* higher scattering rate, artificially]
     1: 'default' Fast modes: following Yan & Lazarian 2002, accounting for damping from viscous, ion-neutral, and other effects, and suppression if beta>1
     2: Fast modes following a pure isotropic Kolmogorov cascade down to gyro radius [*much* higher scattering rate, artificially]
 */
double diffusion_coefficient_extrinsic_turbulence(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG,
    double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion)
{
    double f_cas_ET=MAX_REAL_NUMBER, EPSILON_SMALL=1.e-50, h0_kpc=L_scale*(All.UnitLength_in_cm/All.HubbleParam)/3.086e21;
    if(mode <= 0) /* Alfvenic turbulence [default here following Chandran 200, including anisotropy effects] */
    {
        f_cas_ET = 0.007 * C_LIGHT_CODE / vA_noion; /* damping in Alfvenic turbulence following an anisotropic Goldreich-Shridar cascade, per Chandran 2000 */
        if(mode==-2) {f_cas_ET = 1;} /* pure Goldreich-Shridhar cascade, ignoring anisotropic effects per Chandran-00 */ //if(M_A < 1.) {f_cas_ET = pow(1./DMAX(M_A*M_A , x_LL), 1./3.);} /* Lazarian '16 modification for weak cascade in sub-Alfvenic turbulence */
    }
    if(mode >= 0) /* Fast modes [default here following Yan & Lazarian 2002, including damping effects] */
    {
        double R_CR_GV=return_CRbin_CR_rigidity_in_GV(target,k_CRegy);
        double n1=rho_cgs/PROTONMASS, T4=temperature/1.e4, fcasET_colless = 0.04*cs_thermal/vA_noion; /* collisionless [Landau] damping of fast modes */
        double fcasET_viscBrg = 0.03*pow(EPSILON_SMALL + M_A,4./3.)*T4/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*R_CR_GV*T4,1./6.); /* Spitzer/Braginski viscous damping of fast modes */
        double fcasET_viscMol = 0.41*pow(EPSILON_SMALL + M_A,4./3.)*nh0/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*R_CR_GV/(EPSILON_SMALL + T4),1./6.); /* atomic/molecular collisional damping of fast modes */
        double f_cas_ET_fast = fcasET_colless + fcasET_viscBrg + fcasET_viscMol; /* fast modes, accounting for damping, following Yan+Lazarian 2005 */
        double fast_gyrores_dampingsuppression = 1.; // term to account for the fact that small pitch angles become unscattered when neutral fraction is large or beta >~ 1, making kappa blow up rapidly */
        double beta_half = cs_thermal / vA_noion; if(beta_half > 1.) {fast_gyrores_dampingsuppression=0;} else {fast_gyrores_dampingsuppression*=exp(-beta_half*beta_half*beta_half);} // parallel modes strongly damped if beta >~ 1
        double f_neutral_crit = 0.001 * pow(T4,0.25) / (pow(n1*beta_half*beta_half,0.75) * sqrt(h0_kpc)); // neutral fraction above which the parallel modes are strongly damped
        if(nh0 > 2.*f_neutral_crit) {fast_gyrores_dampingsuppression=0;} else {fast_gyrores_dampingsuppression*=exp(-(nh0*nh0*nh0*nh0)/(f_neutral_crit*f_neutral_crit*f_neutral_crit*f_neutral_crit));} // suppression very rapid, as exp(-[fn/f0]^4)
        if(mode==2) {fast_gyrores_dampingsuppression=0; f_cas_ET_fast = 0.0009 * pow(R_CR_GV*h0_kpc*h0_kpc/b_muG,1./3.)/(M_A*M_A);} /* kappa~9e28 * (l_alfven/kpc)^(2/3) * RGV^(1/3) * (B/muG)^(-1/3),  follows Jokipii 1966, with our corrections for spectral shape */
        if(mode==3) {fast_gyrores_dampingsuppression=0; f_cas_ET_fast = 0.003 * (h0_kpc + 0.1*pow(R_CR_GV*h0_kpc*h0_kpc/b_muG,1./3.)/(M_A*M_A) + 2.4e-7*R_CR_GV/b_muG);} /* Snodin et al. 2016 -- different expression for extrinsic MHD-turb diffusivity, using the proper definition of L_Alfven for scaling to the correct limit */
        f_cas_ET = 1./(EPSILON_SMALL + 1./(EPSILON_SMALL+f_cas_ET) + fast_gyrores_dampingsuppression / (EPSILON_SMALL+f_cas_ET_fast)); /* combine fast-mode and Alfvenic scattering */
    }
    return 1.e32 * h0_kpc / (EPSILON_SMALL + M_A*M_A) * f_cas_ET;
}



/*!----------------------------------------------------------------------------------------------------------------------------------------------------
 routines below are more general and/or numerical: they generally do NOT need to be modified even if you are changing the
   physical assumptions, energies, or other properties of the CRs
 ----------------------------------------------------------------------------------------------------------------------------------------------------*/


/* utility to estimate -locally- (without multi-pass filtering) the local Alfven Mach number */
double Get_AlfvenMachNumber_Local(int i, double vA_idealMHD_codeunits, int use_shear_corrected_vturb_flag)
{
    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,x_LL,M_A,h0,EPSILON_SMALL=1.e-50,fturb_multiplier=1; // factor which will represent which cascade model we are going to use
    for(i1=0;i1<3;i1++)
    {
        v2_t += SphP[i].VelPred[i1]*SphP[i].VelPred[i1];
        for(i2=0;i2<3;i2++) {dv2_t += SphP[i].Gradients.Velocity[i1][i2]*SphP[i].Gradients.Velocity[i1][i2];}
#ifdef MAGNETIC
        b2_t += Get_Gas_BField(i,i1) * Get_Gas_BField(i,i1);
        for(i2=0;i2<3;i2++) {db2_t += SphP[i].Gradients.B[i1][i2]*SphP[i].Gradients.B[i1][i2];}
#endif
    }
    v2_t=sqrt(v2_t); b2_t=sqrt(b2_t); dv2_t=sqrt(dv2_t); db2_t=sqrt(db2_t); dv2_t/=All.cf_atime; db2_t/=All.cf_atime; b2_t*=All.cf_a2inv; db2_t*=All.cf_a2inv; v2_t/=All.cf_atime; dv2_t/=All.cf_atime;
    h0=Get_Particle_Size(i)*All.cf_atime; // physical units

    if(use_shear_corrected_vturb_flag == 1)
    {
        dv2_t = sqrt((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) *
            (SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) + (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) *
            (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) + (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) * (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
            (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] + SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) - (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] + SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv;
#ifdef MAGNETIC
        db2_t = sqrt((1./2.)*((SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) * (SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) +
            (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) * (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) +
            (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2]) * (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2])) +
            (2./3.)*((SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[0][0] + SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[1][1] +
            SphP[i].Gradients.B[2][2]*SphP[i].Gradients.B[2][2]) - (SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[2][2] +
            SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[1][1] + SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[2][2]))) * All.cf_a3inv;
#endif
        double db_v_equiv = h0 * db2_t * vA_idealMHD_codeunits / (EPSILON_SMALL + b2_t); /* effective delta-v corresponding to delta-B fluctuations */
        double dv_e = sqrt((h0*dv2_t)*(h0*dv2_t) + db_v_equiv*db_v_equiv + EPSILON_SMALL); /* total effective delta-velocity */
        double vA_eff = sqrt(vA_idealMHD_codeunits*vA_idealMHD_codeunits + db_v_equiv*db_v_equiv + EPSILON_SMALL); /* effective Alfven speed including fluctuation-B */
        M_A = (EPSILON_SMALL + dv_e) / (EPSILON_SMALL + vA_eff);
    } else {
        M_A = h0*(EPSILON_SMALL + dv2_t) / (EPSILON_SMALL + vA_idealMHD_codeunits); /* velocity fluctuation-inferred Mach number */
        M_A = DMAX(M_A , h0*(EPSILON_SMALL + db2_t) / (EPSILON_SMALL + b2_t)); /* B-field fluctuation-inferred Mach number [in incompressible B-turb, this will 'catch' where dv is locally low instanteously] */
    }
    M_A = DMAX( EPSILON_SMALL , M_A ); // proper calculation of the local Alfven Mach number
    return M_A;
}



/* cosmic ray heating of gas, from Guo & Oh 2008, following Mannheim & Schlickeiser 1994.
    We assume all the electron losses go into radiation [ignoring ionization for now], as the electron coulomb losses into gas are lower than protons by factor of energy and me/mp.
    For protons, we assume 1/6 of the hadronic losses (based on branching ratios) and all of the Coulomb losses thermalize.
    Do want to make sure that the rates we assume here are consistent with those used in the CR cooling routine above. */
double CR_gas_heating(int target, double n_elec, double nHcgs)
{
    double e_heat=0, e_CR_units_0=(SphP[target].Density*All.cf_a3inv/P[target].Mass) * (All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam) / nHcgs; int k_CRegy;
    double a_hadronic = 6.37e-16, b_coulomb_per_GeV = 3.09e-16*n_elec*HYDROGEN_MASSFRAC, f_heat_hadronic=1./6.; /* some coefficients; a_hadronic is the default coefficient, b_coulomb_per_GeV the default divided by GeV, b/c we need to divide the energy per CR  */
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double Z = return_CRbin_CR_charge_in_e(target,k_CRegy);
        if(Z > 0) /* protons only here */
        {
            double e_cr_units = SphP[target].CosmicRayEnergyPred[k_CRegy] * e_CR_units_0;
#if (N_CR_PARTICLE_BINS > 2)
            double E_GeV = return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), beta = return_CRbin_beta_factor(target,k_CRegy);
            e_heat += b_coulomb_per_GeV * ((Z*Z)/(beta*E_GeV)) * e_cr_units; // all protons Coulomb-heat, can be rapid for low-E
            if(E_GeV>=0.78) {e_heat += f_heat_hadronic * a_hadronic * e_cr_units;} // only GeV CRs or higher trigger above threshold for collisions
#else
            e_heat += (0.87*f_heat_hadronic*a_hadronic + 0.53*b_coulomb_per_GeV) * e_cr_units; /* for N<=2, assume a universal spectral shape, the factor here corrects for the fraction above-threshold for hadronic interactions, and 0.53 likewise for averaging  */
#endif
        }
    }
    return e_heat;
}



/* master routine to assign diffusion coefficients. for the most relevant physical models, we do a lot of utility here but do the more interesting
    (and uncertain) physical calculation in the relevant sub-routines above, so you don't need to modify all of this in most cases */
void CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(int i)
{
    /* first define some very general variables, and calculate some useful quantities that will be used for any model */
    int k_CRegy; double DiffusionCoeff, CR_kappa_streaming, CRPressureGradScaleLength, v_streaming; v_streaming=Get_CosmicRayStreamingVelocity(i);
#if (COSMIC_RAYS_DIFFUSION_MODEL > 0)
    double cs_thermal,M_A,L_scale,vA_code,vA_noion,codelength_to_kpc,Z_charge_CR,gizmo2gauss,Omega_per_GeV,Bmag,unit_kappa_code,b_muG,E_B,f_ion,temperature,EPSILON_SMALL; int k;
    unit_kappa_code=All.UnitVelocity_in_cm_per_s*All.UnitLength_in_cm/All.HubbleParam; gizmo2gauss=sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam); f_ion=1; temperature=0; EPSILON_SMALL=1.e-50; codelength_to_kpc=(All.UnitLength_in_cm/All.HubbleParam)/(3.086e21);
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv; cs_thermal=sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred)); /* quick thermal pressure properties (we'll assume beta=1 if MHD not enabled) */
#ifdef MAGNETIC /* get actual B-field */
    double B[3]={0}; Bmag=0; for(k=0;k<3;k++) {B[k]=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B[k]*B[k];} // B-field in code units (physical)
#endif
    Bmag=sqrt(DMAX(Bmag,0)); b_muG=Bmag*gizmo2gauss/1.e-6; b_muG=sqrt(b_muG*b_muG + 1.e-6); vA_code=sqrt(Bmag*Bmag/(SphP[i].Density*All.cf_a3inv)); vA_noion=vA_code; E_B=0.5*Bmag*Bmag*(P[i].Mass/(SphP[i].Density*All.cf_a3inv)); Omega_per_GeV=0.00898734*b_muG*(All.UnitTime_in_s/All.HubbleParam); /* B-field in units of physical microGauss; set a floor at nanoGauss level */
#ifdef COOLING
    double ne=SphP[i].Ne, nh0=0, nHe0=0, nHepp, nhp, nHeII, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); rho_cgs = rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; // get thermodynamic properties
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    if(SphP[i].DelayTimeHII>0) {nh0=DMIN(1.e-4,nh0); nHe0=DMIN(1.e-5,nHe0);} // account for our effective ionization model here
#endif
    f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
#endif
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_code /= sqrt(f_ion); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    M_A = Get_AlfvenMachNumber_Local(i,vA_noion,0); /* get turbulent Alfven Mach number estimate. 0 or 1 to turn on shear-correction */
    L_scale = Get_Particle_Size(i)*All.cf_atime; /* define turbulent scales [estimation of M_A defined by reference to this scale */
#endif
    
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        DiffusionCoeff=0; CR_kappa_streaming=0; CRPressureGradScaleLength=Get_CosmicRayGradientLength(i,k_CRegy); /* set these for the bin as we get started */
#ifndef COSMIC_RAYS_DISABLE_STREAMING /* self-consistently calculate the diffusion coefficients for cosmic ray fluids; first the streaming part of this (kappa~v_stream*L_CR_grad) following e.g. Wentzel 1968, Skilling 1971, 1975, Holman 1979, as updated in Kulsrud 2005, Yan & Lazarian 2008, Ensslin 2011 */
        CR_kappa_streaming = GAMMA_COSMICRAY * v_streaming * CRPressureGradScaleLength; /* the diffusivity is now just the product of these two coefficients (all physical units) */
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 0) /* set diffusivity to a universal power-law scaling (constant per-bin)  */
        DiffusionCoeff = diffusion_coefficient_constant(i,k_CRegy); //  this is the input value of the diffusivity, for constant-kappa models
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL < 0) /* disable CR diffusion, specifically */
        DiffusionCoeff = 0; // no diffusion (but -can- allow streaming)
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 3) /* Farber et al. 2018 -- higher coeff in neutral gas, lower in ionized gas */
        DiffusionCoeff = (3.e29/unit_kappa_code) * (1.-f_ion + f_ion/30.); // 30x lower in neutral (note use f_ion directly here, not temperature as they do)
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 4) /* Wiener et al. 2017 style pure-streaming but with larger streaming speeds and limited losses, using their scaling for assumption that turbulent+non-linear Landau only dominate damping */
        double ni_m3=f_ion*(rho_cgs/PROTONMASS)/1.e-3, T6=temperature/1.e6, Lturbkpc=L_scale*codelength_to_kpc, Lgradkpc=CRPressureGradScaleLength*(All.UnitLength_in_cm/All.HubbleParam)/3.086e21, h0_fac=Get_Particle_Size(i)*All.cf_atime*All.cf_a2inv*All.UnitVelocity_in_cm_per_s/1.e6, dv2_10=0; for(k=0;k<3;k++) {int j; for(j=0;j<3;j++) {dv2_10 += SphP[i].Gradients.Velocity[j][k]*SphP[i].Gradients.Velocity[j][k]*h0_fac*h0_fac;}}
        double ecr_14 = SphP[i].CosmicRayEnergyPred[k_CRegy] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * ((All.UnitEnergy_in_cgs/All.HubbleParam)/pow(All.UnitLength_in_cm/All.HubbleParam,3)) / 1.0e-14; // CR energy density in CGS units //
        CR_kappa_streaming = GAMMA_COSMICRAY * CRPressureGradScaleLength * (v_streaming + (1.e5/All.UnitVelocity_in_cm_per_s)*(4.1*pow(MIN_REAL_NUMBER+ni_m3*T6,0.25)/pow(MIN_REAL_NUMBER+ecr_14*Lgradkpc,0.5) + 1.2*pow(MIN_REAL_NUMBER+dv2_10*ni_m3,0.75)/(MIN_REAL_NUMBER+ecr_14*sqrt(Lturbkpc)))); // convert to effective diffusivity
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 5) /* streaming at fast MHD wavespeed [just to see what it does] */
        CR_kappa_streaming = GAMMA_COSMICRAY * sqrt(v_streaming*v_streaming + cs_thermal*cs_thermal) * CRPressureGradScaleLength;
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 1) || (COSMIC_RAYS_DIFFUSION_MODEL == 2) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* textbook extrinsic turbulence model: kappa~v_CR*r_gyro * B_bulk^2/(B_random[scale~r_gyro]^2) v_CR~c, r_gyro~p*c/(Z*e*B)~1e12 cm * RGV *(3 muG/B)  (RGV~1 is the magnetic rigidity). assuming a Kolmogorov spectrum */
        int scatter_modes = 0; /* default to using both Alfven+damped-fast modes */
#if (COSMIC_RAYS_DIFFUSION_MODEL==1)
        scatter_modes = -1; /* Alfven modes only*/
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL==2)
        scatter_modes = 1; /* Fast modes only*/
#endif
        DiffusionCoeff = diffusion_coefficient_extrinsic_turbulence(scatter_modes,i,k_CRegy,M_A,L_scale,b_muG,vA_noion,rho_cgs,temperature,cs_thermal,nh0,nHe0,f_ion) / unit_kappa_code;
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 6) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* self-confinement-based diffusivity */
        double Omega_gyro=0.00898734*b_muG*(All.UnitTime_in_s/All.HubbleParam)/return_CRbin_CR_rigidity_in_GV(i,k_CRegy), r_L=C_LIGHT_CODE/Omega_gyro, kappa_0=r_L*C_LIGHT_CODE; // some handy numbers for limiting extreme-kappa below
        CR_kappa_streaming = diffusion_coefficient_self_confinement(COSMIC_RAYS_SET_SC_MODEL,i,k_CRegy,M_A,L_scale,b_muG,vA_noion,rho_cgs,temperature,cs_thermal,nh0,nHe0,f_ion) / unit_kappa_code;
        if(!isfinite(CR_kappa_streaming)) {CR_kappa_streaming = 1.e30/unit_kappa_code;} /* apply some limiters since its very easy for the routine above to give wildly-large-or-small diffusivity, which wont make a difference compared to just 'small' or 'large', but will mess things up numerically */
        CR_kappa_streaming = DMIN( DMAX( DMIN(DMAX(CR_kappa_streaming,kappa_0) , 1.0e6*GAMMA_COSMICRAY*CRPressureGradScaleLength*COSMIC_RAYS_M1) , 1.e25/unit_kappa_code ) , 1.e32/unit_kappa_code );
#endif

        /* -- ok, we've done what we came to do -- everything below here is pure-numerical, not physics, and should generally not be modified -- */

#if (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* 'combined' extrinsic turbulence + self-confinement model: add scattering rates linearly (plus lots of checks to prevent unphysical bounds) */
        CR_kappa_streaming = 1. / (EPSILON_SMALL +  1./(CR_kappa_streaming+EPSILON_SMALL) + 1./(DiffusionCoeff+EPSILON_SMALL) ); DiffusionCoeff=0; if(!isfinite(CR_kappa_streaming)) {CR_kappa_streaming = 1.e30/unit_kappa_code;} // if scattering rates add linearly, this is a rough approximation to the total transport (essentially, smaller of the two dominates)
        CR_kappa_streaming = DMIN( DMAX( CR_kappa_streaming , kappa_0 ) , 1.0e6*GAMMA_COSMICRAY*CRPressureGradScaleLength*COSMIC_RAYS_M1 ); CR_kappa_streaming = DMIN( DMAX( CR_kappa_streaming , 1.e25/unit_kappa_code ) , 1.e32/unit_kappa_code );
#endif
        DiffusionCoeff = DiffusionCoeff + CR_kappa_streaming; //  add 'diffusion' and 'streaming' terms since enter numerically the same way
            
#ifndef COSMIC_RAYS_M1 /* now we apply a limiter to prevent the coefficient from becoming too large: cosmic rays cannot stream/diffuse with v_diff > c */
        // [all of this only applies if we are using the pure-diffusion description: the M1-type description should -not- use a limiter here, or negative kappa]
        double diffusion_velocity_limit=C_LIGHT_CODE, L_eff=DMAX(Get_Particle_Size(i)*All.cf_atime,CRPressureGradScaleLength); /* maximum diffusion velocity (set <c if desired) */
        double kappa_diff_vel = DiffusionCoeff * GAMMA_COSMICRAY_MINUS1 / L_eff; DiffusionCoeff *= 1 / (1 + kappa_diff_vel/diffusion_velocity_limit); /* caps maximum here */
#ifdef GALSF /* for multi-physics problems, we suppress diffusion [in the FLD-like limit] where it is irrelevant for timestepping-sake */
        DiffusionCoeff *= 1 / (1 + kappa_diff_vel/(0.01*diffusion_velocity_limit)); /* caps maximum here */
        double P_cr_Ratio = Get_Gas_CosmicRayPressure(i,k_CRegy) / (MIN_REAL_NUMBER + SphP[i].Pressure), P_min=1.0e-4; if(P_cr_Ratio < P_min) {DiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
        P_min = 1.0e-6; if(P_cr_Ratio < P_min) {DiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
#endif
        DiffusionCoeff /= GAMMA_COSMICRAY_MINUS1; // ensure correct units for subsequent operations //
#endif
        if((DiffusionCoeff<=0)||(isnan(DiffusionCoeff))) {DiffusionCoeff=0;} /* nan check */
        SphP[i].CosmicRayDiffusionCoeff[k_CRegy] = DiffusionCoeff; /* final assignment! */
    } // end CR bin loop
}



/* utility routine which handles the numerically-necessary parts of the CR 'injection' for you */
double inject_cosmic_rays(double CR_energy_to_inject, double injection_velocity, int source_PType, int target, double *dir)
{
    int k_CRegy,k; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double dEcr = CR_energy_to_inject * CR_energy_spectrum_injection_fraction(k_CRegy,source_PType,injection_velocity);
        SphP[target].CosmicRayEnergy[k_CRegy]+=dEcr; SphP[target].CosmicRayEnergyPred[k_CRegy]+=dEcr;
#ifdef COSMIC_RAYS_M1
        double dir_mag=0, flux_mag=dEcr*COSMIC_RAYS_M1, dir_to_use[3]={0};
#ifdef MAGNETIC
        double B_dot_dir=0, Bdir[3]={0}; for(k=0;k<3;k++) {Bdir[k]=SphP[target].BPred[k]; B_dot_dir+=dir[k]*Bdir[k];} // the 'default' direction is projected onto B
        for(k=0;k<3;k++) {dir_to_use[k]=B_dot_dir*Bdir[k];} // launch -along- B, projected [with sign determined] by the intially-desired direction
#else
        for(k=0;k<3;k++) {dir_to_use[k]=dir[k];} // launch in the 'default' direction
#endif
        for(k=0;k<3;k++) {dir_mag += dir_to_use[k]*dir_to_use[k];}
        if(dir_mag <= 0) {dir_to_use[0]=0; dir_to_use[1]=0; dir_to_use[2]=1; dir_mag=1;}
        for(k=0;k<3;k++) {double dflux=flux_mag*dir_to_use[k]/sqrt(dir_mag); SphP[target].CosmicRayFlux[k_CRegy][k]+=dflux; SphP[target].CosmicRayFluxPred[k_CRegy][k]+=dflux;}
#endif
    }
}



/* return CR pressure within a given bin */
double INLINE_FUNC Get_Gas_CosmicRayPressure(int i, int k_CRegy)
{
    if((P[i].Mass > 0) && (SphP[i].Density>0) && (SphP[i].CosmicRayEnergyPred[k_CRegy] > 0))
    {
        return GAMMA_COSMICRAY_MINUS1 * (SphP[i].CosmicRayEnergyPred[k_CRegy] * SphP[i].Density) / P[i].Mass; // cosmic ray pressure = (4/3-1) * e_cr = 1/3 * (E_cr/Vol) //
    } else {return 0;}
}



/* return CR gradient scale length, with various physical limiters applied: not intended for pure numerical gradient-length calculations (where units dont matter), but for preventing some unphysical situations */
double Get_CosmicRayGradientLength(int i, int k_CRegy)
{
    /* now we need the -parallel- cosmic ray pressure or energy density scale length */
    double CRPressureGradMag = 0.0;
    int k; for(k=0;k<3;k++) {CRPressureGradMag += SphP[i].Gradients.CosmicRayPressure[k_CRegy][k]*SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];}
    CRPressureGradMag = sqrt(1.e-46 + CRPressureGradMag); // sqrt to make absolute value
#ifdef MAGNETIC /* with anisotropic transport, we really want the -parallel- gradient scale-length, so need another factor here */
    double B2_tot=0.0, b=0; CRPressureGradMag=0; for(k=0;k<3;k++) {b=Get_Gas_BField(i,k); B2_tot+=b*b; CRPressureGradMag+=b*SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];} // note, this is signed!
    CRPressureGradMag = sqrt((1.e-40 + CRPressureGradMag*CRPressureGradMag) / (1.e-46 + B2_tot)); // divide B-magnitude to get scalar magnitude, and take sqrt[(G.P)^2] to get absolute value
#endif
    
    /* limit the scale length: if too sharp, need a slope limiter at around the particle size */
    double L_gradient_min = Get_Particle_Size(i) * All.cf_atime;
    /* limit this scale length; if the gradient is too shallow, there is no information beyond a few smoothing lengths, so we can't let streaming go that far */
    double L_gradient_max = DMAX(1000.*L_gradient_min, 500.0*PPP[i].Hsml*All.cf_atime);

    /* also, physically, cosmic rays cannot stream/diffuse with a faster coefficient than ~v_max*L_mean_free_path, where L_mean_free_path ~ 2.e20 * (cm^-3/n) */
    double nH_cgs = SphP[i].Density * All.cf_a3inv * ( All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam ) / PROTONMASS ;
    double L_mean_free_path = (3.e25 / nH_cgs) / (All.UnitLength_in_cm / All.HubbleParam);
    L_gradient_max = DMIN(L_gradient_max, L_mean_free_path);
    
    double CRPressureGradScaleLength = Get_Gas_CosmicRayPressure(i,k_CRegy) / CRPressureGradMag * All.cf_atime;
    if(CRPressureGradScaleLength > 0) {CRPressureGradScaleLength = 1.0/(1.0/CRPressureGradScaleLength + 1.0/L_gradient_max);} else {CRPressureGradScaleLength=0;}
    CRPressureGradScaleLength = sqrt(L_gradient_min*L_gradient_min + CRPressureGradScaleLength*CRPressureGradScaleLength);
    return CRPressureGradScaleLength; /* this is returned in -physical- units */
}



/* return the effective CR 'streaming' velocity for sub-grid [unresolved] models with streaming velocity set by e.g. the Alfven speed along the gradient of the CR pressure */
double Get_CosmicRayStreamingVelocity(int i)
{
    /* in the weak-field (high-beta) case, the streaming velocity is approximately the sound speed */
    double v_streaming = sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred)); // thermal ion sound speed //
#ifdef MAGNETIC
    /* in the strong-field (low-beta) case, it's actually the Alfven velocity: interpolate between these */
    double vA_2 = 0.0; double cs_stream = v_streaming;
    int k; for(k=0;k<3;k++) {vA_2 += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
    vA_2 *= All.cf_afac1 / (All.cf_atime * SphP[i].Density);
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_2 /= Get_Gas_Ionized_Fraction(i); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    v_streaming = DMIN(1.0e6*cs_stream, sqrt(MIN_REAL_NUMBER*cs_stream*cs_stream + vA_2)); // limit to Alfven speed //
#endif
#ifdef COSMIC_RAYS_M1
    v_streaming = DMIN(v_streaming , COSMIC_RAYS_M1); // limit to maximum transport speed //
#endif
    v_streaming *= All.cf_afac3; // converts to physical units and rescales according to chosen coefficient //
    return v_streaming;
}


/* routine to quickly estimate the atomic mass of the CR particles (in proton masses): making assumptions here [e.g. no positrons] -- could always hard-code some choice here for more complicated scenarios */
double return_CRbin_CRmass_in_mp(int target, int k_CRegy)
{
    double Z = return_CRbin_CR_charge_in_e(target,k_CRegy), m_cr_mp=1;
    if(Z<0) {m_cr_mp=0.000544618;} else {if(Z>1.5) {m_cr_mp=2.*Z;}} // mass relative to proton mass [assuming no psoitrons here -- could make this more user-specified, but ignoring here]
    return m_cr_mp;
}


/* routine which returns the dimensionless velocity beta = v/c of the CRs (this is not the RSOL, but true c, for e.g. cooling, energies, etc) */
double return_CRbin_beta_factor(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double q = return_CRbin_CR_rigidity_in_GV(target,k_CRegy)*1.06579*fabs(return_CRbin_CR_charge_in_e(target,k_CRegy))/m_cr_mp; // dimensionless factor to convert from R to beta
    double gamma = sqrt(1.+q*q), beta = q/gamma;
    return beta;
}


/* routine which returns the dimensionless lorentz factor gamma=1/sqrt[1-beta^2] of the CRs (this is not the RSOL, but true c, for e.g. cooling, energies, etc) */
double return_CRbin_gamma_factor(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double q = return_CRbin_CR_rigidity_in_GV(target,k_CRegy)*1.06579*fabs(return_CRbin_CR_charge_in_e(target,k_CRegy))/m_cr_mp; // dimensionless factor to convert from R to beta
    return sqrt(1.+q*q);
}


/* routine which returns the CR kinetic energy in GeV, for a given rigidity, etc. */
double return_CRbin_kinetic_energy_in_GeV(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double R_GV = return_CRbin_CR_rigidity_in_GV(target,k_CRegy), Z = fabs(return_CRbin_CR_charge_in_e(target,k_CRegy));
    double q = R_GV*Z*1.06579/m_cr_mp; // dimensionless factor to convert from R to beta
    double gamma = sqrt(1.+q*q), beta = q/gamma;
    double KE_fac = DMAX(0.,(1.-sqrt(DMAX(0.,1.-beta*beta)))) / DMAX(beta,MIN_REAL_NUMBER); if(beta < 0.01) {KE_fac = 0.5*beta;} // non-relativistic expansion used when gamma very close to one, to prevent numerical errors
    return R_GV * Z * KE_fac;
}


/* routine which returns the CR number density at a given bin in cm^-3 */
double return_CRbin_numberdensity_in_cgs(int target, int k_CRegy)
{
    double e_CR_tot = SphP[target].CosmicRayEnergyPred[k_CRegy] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * (All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam);
    double ECR_per  = return_CRbin_kinetic_energy_in_GeV(target,k_CRegy) * 0.00160218; /* converts to energy in erg */
    return e_CR_tot / ECR_per;
}


/* handy functoin that just returns the B-field magnitude in microGauss, physical units. purely here to save us time re-writing this */
double get_cell_Bfield_in_microGauss(int i)
{
    double Bmag=0, gizmo2mugauss_2=4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam*1.e12;
#ifdef MAGNETIC
    int k; for(k=0;k<3;k++) {double B=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B*B*gizmo2mugauss_2;} // actual B-field in code units
#else
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv*gizmo2mugauss_2; // assume equipartition
#endif
    return sqrt(DMAX(Bmag,0));
}


/* handy functoin that just returns the radiation energy density in eV/cm^-3, physical units. purely here to save us time re-writing this */
double get_cell_Urad_in_eVcm3(int i)
{
    double erad = 0.26*All.cf_a3inv/All.cf_atime; // default with the CMB energy density, which we assume here is a baseline minimum
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    int kfreq; double e_units = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam * 6.24151e11;
#ifdef RT_EVOLVE_ENERGY
    for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad+=SphP[i].Rad_E_gamma_Pred[kfreq]*e_units;}
#else
    for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad+=SphP[i].Rad_E_gamma[kfreq]*e_units;}
#endif
#else
    double uRad_MW = 0.31 + 0.66; /* dust (0.31) and stars (0.66) for Milky way ISRF from Draine (2011); want this to decay as we approach the IGM (where CMB totally dominates) */
    double prefac_rad=1, rho_cgs=SphP[i].Density*All.cf_a3inv*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam; if(All.ComovingIntegrationOn) {double rhofac = rho_cgs / (1000.*All.OmegaBaryon*(All.HubbleParam*HUBBLE_CGS)*(All.HubbleParam*HUBBLE_CGS)*(3./(8.*M_PI*GRAVITY_G))*All.cf_a3inv);
        if(rhofac < 0.2) {prefac_rad=0;} else {if(rhofac > 200.) {prefac_rad=1;} else {prefac_rad=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off stellar background for any gas with density unless it's >1000 times the cosmic mean density
    prefac_rad *= rho_cgs/(0.01*PROTONMASS + rho_cgs); // cut off below low densities, ~1e-2
    erad += uRad_MW * prefac_rad;
#endif
    return erad;
}



/* routine to do the drift/kick operations for CRs: mode=0 is kick, mode=1 is drift */
#if !defined(COSMIC_RAYS_ALFVEN)
double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode)
{
    if(dt_entr <= 0) {return 0;} // no update

    int k_CRegy;
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        int k; double eCR, u0; k=0; if(mode==0) {eCR=SphP[i].CosmicRayEnergy[k_CRegy]; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred[k_CRegy]; u0=SphP[i].InternalEnergyPred;} // initial energy
        if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforced throughout code
        if(eCR < 0) {eCR=0;} // limit to physical values
        
#if defined(COSMIC_RAYS_M1) /* CR FLUX VECTOR UPDATE */
        // this is the exact solution for the CR flux-update equation over a finite timestep dt: it needs to be solved this way [implicitly] as opposed to explicitly for dt because in the limit of dt_cr_dimless being large, the problem exactly approaches the diffusive solution
        double DtCosmicRayFlux[3]={0}, flux[3]={0}, CR_veff[3]={0}, CR_vmag=0, q_cr = 0, cr_speed = COSMIC_RAYS_M1;
        cr_speed = DMAX( All.cf_afac3*SphP[i].MaxSignalVel , DMIN(COSMIC_RAYS_M1 , 10.*fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy])/(Get_Particle_Size(i)*All.cf_atime)));
        for(k=0;k<3;k++) {DtCosmicRayFlux[k] = -fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy]) * (P[i].Mass/SphP[i].Density) * (SphP[i].Gradients.CosmicRayPressure[k_CRegy][k]/GAMMA_COSMICRAY_MINUS1);}
#ifdef MAGNETIC // do projection onto field lines
        double B0[3]={0}, Bmag2=0, DtCRDotBhat=0;
        for(k=0;k<3;k++)
        {
            if(mode==0) {B0[k]=SphP[i].B[k];} else {B0[k]=SphP[i].BPred[k];}
            DtCRDotBhat += DtCosmicRayFlux[k] * B0[k]; Bmag2 += B0[k]*B0[k];
        }
        if(Bmag2 > 0) {for(k=0;k<3;k++) {DtCosmicRayFlux[k] = DtCRDotBhat * B0[k] / Bmag2;}}
#endif
        double dt_cr_dimless = dt_entr * cr_speed*cr_speed * GAMMA_COSMICRAY_MINUS1 / (MIN_REAL_NUMBER + fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy]));
        if((dt_cr_dimless > 0)&&(dt_cr_dimless < 20.)) {q_cr = exp(-dt_cr_dimless);} // factor for CR interpolation
        if(mode==0) {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFlux[k_CRegy][k];}} else {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFluxPred[k_CRegy][k];}}
#ifdef MAGNETIC // do projection onto field lines
        double fluxmag=0, fluxdot=0; for(k=0;k<3;k++) {fluxmag+=flux[k]*flux[k]; fluxdot+=flux[k]*B0[k];}
        if(fluxmag>0) {fluxmag=sqrt(fluxmag);} else {fluxmag=0;}
        if(fluxdot<0) {fluxmag*=-1;} // points down-field
        if(Bmag2>0) {for(k=0;k<3;k++) {flux[k] = fluxmag * B0[k] / sqrt(Bmag2);}} // re-assign to be along field
#endif
        for(k=0;k<3;k++) {flux[k] = q_cr*flux[k] + (1.-q_cr)*DtCosmicRayFlux[k];} // updated flux
        for(k=0;k<3;k++) {CR_veff[k]=flux[k]/(eCR+MIN_REAL_NUMBER); CR_vmag+=CR_veff[k]*CR_veff[k];} // effective streaming speed
        if((CR_vmag <= 0) || (isnan(CR_vmag))) // check for valid numbers
        {
            for(k=0;k<3;k++) {flux[k]=0; for(k=0;k<3;k++) {CR_veff[k]=0;}} // zero if invalid
        } else {
            CR_vmag = sqrt(CR_vmag); if(CR_vmag > cr_speed) {for(k=0;k<3;k++) {flux[k]*=cr_speed/CR_vmag; CR_veff[k]*=cr_speed/CR_vmag;}} // limit flux to free-streaming speed [as with RT]
        }
        if(mode==0) {for(k=0;k<3;k++) {SphP[i].CosmicRayFlux[k_CRegy][k]=flux[k];}} else {for(k=0;k<3;k++) {SphP[i].CosmicRayFluxPred[k_CRegy][k]=flux[k];}}
#endif
    
        /* update scalar CR energy. first update the CR energies from fluxes. since this is positive-definite, some additional care is needed */
        double dCR_dt = SphP[i].DtCosmicRayEnergy[k_CRegy], gamma_eff = GAMMA_COSMICRAY, eCR_tmp = eCR;
        double dCR = dCR_dt*dt_entr, dCRmax = 1.e10*(eCR_tmp+MIN_REAL_NUMBER);
#if defined(GALSF)
        dCRmax = DMAX(2.0*eCR_tmp , 0.1*u0*P[i].Mass);
#endif
        if(dCR > dCRmax) {dCR=dCRmax;} // don't allow excessively large values
        if(dCR < -eCR_tmp) {dCR=-eCR_tmp;} // don't allow it to go negative
        double eCR_00 = eCR_tmp; eCR_tmp += dCR; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy]=eCR_tmp;} else {SphP[i].CosmicRayEnergyPred[k_CRegy]=eCR_tmp;} // updated energy
        double eCR_0 = eCR_tmp; // save this value for below
        
        /* now need to account for the adiabatic heating/cooling of the 'fluid', here, with gamma=gamma_eff */
        double d_div = (-(gamma_eff-1.) * P[i].Particle_DivVel*All.cf_a2inv) * dt_entr;
        if(All.ComovingIntegrationOn) {d_div += (-3.*(gamma_eff-1.) * All.cf_hubble_a) * dt_entr;} /* adiabatic term from Hubble expansion (needed for cosmological integrations */
        double dCR_div = DMIN(eCR_tmp*d_div , 0.5*u0*P[i].Mass); // limit so don't take away all the gas internal energy [to negative values]
        if(dCR_div + eCR_tmp < 0) {dCR_div = -eCR_tmp;} // check against energy going negative
        eCR_tmp += dCR_div; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        dCR_div = eCR_tmp - eCR_0; // actual change that is going to be applied
        if(dCR_div < -0.5*P[i].Mass*u0) {dCR_div=-0.5*P[i].Mass*u0;} // before re-coupling, ensure this will not cause negative energies
        if(dCR_div < -0.9*eCR_00) {dCR_div=-0.9*eCR_00;} // before re-coupling, ensure this will not cause negative energies
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy] += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;} else {SphP[i].CosmicRayEnergyPred[k_CRegy] += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;}
            
    } // loop over CR bins complete
    return 1;
}
#endif



#if 0
/* this routine does the CR cooling/losses and "heating"/re-acceleration for multi-bin CR spectra: i.e. exchanging CR number
    between bins in the multi-bin approximation and modifying the spectral slope within each bin */
void CR_cooling_and_losses_multibin(int target, double n_elec, double nHcgs, double dtime_cgs)
{
    /*! LOSS TYPES:
     - hadronic+catastrophic: simply remove energy from the bin (N and E decrease together, preserving spectral slope)
     - adiabatic+bremstrahhlung: pure multiplicative: Edot ~ E (so instantaneously conserves slope, shifts pmin,p0,pmax, need to calculate flux)
     - coulomb+ionization: scale identically, for low-E protons important, messy dependence, use fitting function from Girichidis, dp/dt ~ (1 + p^(-1.9))
     - inverse compton+synchrotron: Edot ~ E^2, also pdot~p^2 [ultra-rel b/c e-]: modifies slope
     */
    return;
}
#endif




#endif // closes block for entire file


