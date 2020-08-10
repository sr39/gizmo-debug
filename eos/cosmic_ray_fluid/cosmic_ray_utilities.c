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

#if defined(COSMIC_RAYS_EVOLVE_SPECTRUM)
/* routine which defines the actual bin list for the multi-bin spectral CR models. note the number of entries MUST match the hard-coded N_CR_PARTICLE_BINS defined in allvars.h */
void CR_spectrum_define_bins(void)
{
    /* note, we by default don't go to extremely low-energy proton bins since those have incredibly rapid Coulomb loss times, which without continuous injection grind the code down and just give zero energy in the bins */
    double R[N_CR_PARTICLE_BINS] = {3.16227766e-03, 1.00000000e-02, 3.16227766e-02, 1.00000000e-01, 3.16227766e-01, 1.00000000e+00, 3.16227766e+00, 1.00000000e+01, 3.16227766e+01, 1.00000000e+02, 3.16227766e+02,
                                                                                    1.00000000e-01, 3.16227766e-01, 1.00000000e+00, 3.16227766e+00, 1.00000000e+01, 3.16227766e+01, 1.00000000e+02, 3.16227766e+02}; // bin-centered rigidity in GV, defined here for each bin
    double Z[N_CR_PARTICLE_BINS] = {            -1,             -1,             -1,             -1,             -1,             -1,             -1,             -1,             -1,             -1,             -1,
                                                                                                +1,             +1,             +1,             +1,             +1,             +1,             +1,             +1}; // charge: -1 for electrons, +1 for protons, defined here for each bin
    /* for ease-of-use purposes rather than adding a bunch of extra flags, we will use charge = -2 to signify positrons. the code will understand what to do with them and treat them with the correct charge, but its just a special
        flag. otherwise the charge should correspond to -1 for electrons, or to the fully-ionized charge of a given nucleus (e.g. its atomic number, 1=p/H, 2=He, etc. [some can have hard-coded behaviors/cross-sections below for various processes designed for their use as tracers, including e.g. 4=Be, 5=B, 6/7/8=C/N/O, etc.)
     can use some slightly-offset values if desired to hard-code special flags, or code a shared weight here [e.g. weight in amu] to represent different isotopes, if desired (e.g. 10Be vs 7Be) */
    int k; for(k=0;k<N_CR_PARTICLE_BINS;k++) {CR_global_rigidity_at_bin_center[k]=R[k]; CR_global_charge_in_bin[k]=Z[k];}
}
#endif

/* routine which returns the typical absolute value of the rigidity of a given CR [in GV] in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_rigidity_in_GV(int target, int k_CRegy)
{
    double R = 1;
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    double Rv[2]={1.8 , 0.6}; R=Rv[k_CRegy]; // approximate peak energies of each from Cummings et al. 2016 Fig 15
#endif
#if (N_CR_PARTICLE_BINS > 2) /* arbitrary numbers of bins for CR spectra, assumes binning same in e- and p */
    if(target >= 0) {R=CR_return_mean_rigidity_in_bin_in_GV(target,k_CRegy);} else {R=CR_global_rigidity_at_bin_center[k_CRegy];} // this is pre-defined globally for this bin list
#endif
#endif
    return R;
}

/* routine which returns the typical charge of CRs in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_charge_in_e(int target, int k_CRegy)
{
    double Z = 1;
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    double Zv[2]={1 , -1}; Z=Zv[k_CRegy]; // proton, then e- bin
#endif
#if (N_CR_PARTICLE_BINS > 2)
    Z = CR_global_charge_in_bin[k_CRegy]; // this is pre-defined globally for this bin list
#endif
#endif
    return Z;
}

/* routine which determines the fraction of injected CR energy per 'bin' of CR energy. */
double CR_energy_spectrum_injection_fraction(int k_CRegy, int source_PType, double shock_vel, int return_index_in_bin)
{
    double f_bin = 1./N_CR_PARTICLE_BINS; /* uniformly distributed */
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    double f_bin_v[2]={0.95 , 0.05}; f_bin=f_bin_v[k_CRegy]; // 5% of injection into e-, roughly motivated by observed spectra and nearby SNRs
#endif
#if (N_CR_PARTICLE_BINS > 2) /* multi-bin spectrum for p and e-: inset assumptions about injection spectrum here! */
    double f_elec = 0.02; // fraction of the energy to put into e- as opposed to p+ at injection [early experiments with 'observed'  fraction ~ 1% give lower e-/p+ actually observed in the end, so tentative favoring closer to equal at injection? but not run to z=0, so U_rad high from CMB; still experimenting here]
    double inj_slope = 4.5; // injection slope with j(p) ~ p^(-inj_slope), so dN/dp ~ p^(2-inj_slope)
    double R_break_e = 0.4; // location of spectral break for injection e- spectrum, in GV
    double inj_slope_lowE_e = 4.1; // injection slope with j(p) ~ p^(-inj_slope), so dN/dp ~ p^(2-inj_slope), for electrons below R_break_e
    double Z=return_CRbin_CR_charge_in_e(-1,k_CRegy), R=return_CRbin_CR_rigidity_in_GV(-1,k_CRegy); // get bin-centered Z, R
    if(Z < 0 && R < R_break_e) {inj_slope = inj_slope_lowE_e;} // follow model injection spectra favored in Strong et al. 2011 (A+A, 534, A54), who argue the low-energy e- injection spectrum must break to a lower slope by ~1 independent of propagation and re-acceleration model
    if(return_index_in_bin) {return 2.-inj_slope;} // this is the index corresponding to our dN/dp ~ p^gamma
    double EGeV = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_CRegy); // get bin-centered E_GeV for normalizing total energy in bin
    f_bin = EGeV * pow(R/R_break_e , 3.-inj_slope) * log(CR_global_max_rigidity_in_bin[k_CRegy] / CR_global_min_rigidity_in_bin[k_CRegy]); // normalize accounting for slope, isotropic spectrum, logarithmic bin width [which can vary], and energy per N
    if(Z < 0) {f_bin *= f_elec;} else {f_bin *= 1.-f_elec;} // normalize depending on e- or p+
#endif
#endif
    if(return_index_in_bin) {return 0;}
    return f_bin;
}


/* routine which gives diffusion coefficient as a function of energy for the 'constant diffusion coefficient' models:
    current default: -extremely- simple power-law, assuming diffusion coefficient increases with CR energy per unit charge as (E/Z)^(1/2) */
double diffusion_coefficient_constant(int target, int k_CRegy)
{
    double dimensionless_kappa_relative_to_GV_protons = 1;
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
    dimensionless_kappa_relative_to_GV_protons = pow( return_CRbin_CR_rigidity_in_GV(-1,k_CRegy) , 0.5 ); // assume a quasi-empirical scaling here //
#endif
    return All.CosmicRayDiffusionCoeff * dimensionless_kappa_relative_to_GV_protons;
}



/* routine which gives diffusion coefficient as a function of CR bin for the self-confinement models [in local equilibrium]. mode sets what we assume about the 'sub-grid'
    parameters f_QLT (rescales quasi-linear theory) or f_cas (rescales turbulence strength)
      <=0: fQLT=1 [most naive quasi-linear theory, ruled out by observations],  fcas=1 [standard Goldreich-Shridar cascade]
        1: fQLT=100, fcas=1
        2: fQLT=1, fcas=100
        3: fQLT=1, fcas-K41 from Hopkins et al. 2020 paper, for pure-Kolmogorov isotropic spectrum
        4: fQLT=1, fcas-IK, IK spectrum instead of GS
   if set mode < 0, will also ignore the dust-damping contribution from Squire et al. 2020.
   coefficient is returned in cgs units
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
    double Omega_gyro=(0.00898734*b_muG/R_CR_GV) * UNIT_TIME_IN_CGS, r_L=C_LIGHT_CODE/Omega_gyro, kappa_0=r_L*C_LIGHT_CODE; /* all in physical -code- units */
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
    double U0bar_grain=3., rhograin_int_cgs=1., fac_grain=R_CR_GV*sqrt(n_cgs)*U0bar_grain/(b_muG*rhograin_int_cgs), f_grainsize = DMAX(8.e-4*pow(fac_grain*(temperature/1.e4),0.25), 3.e-3*sqrt(fac_grain)), Z_sol=1.; // b=2, uniform logarithmic grain spectrum over a factor of ~100 in grain size; f_grainsize = 0.07*pow(sqrt(fion*n1)*EcrGeV*T4/BmuG,0.25); // MRN size spectrum
#ifdef METALS
    Z_sol = P[target].Metallicity[0]/0.014;
#endif
    double G_dust = vA_code*k_L * Z_sol * f_grainsize; // also can increase by up to a factor of 2 for regimes where charge collisionally saturated, though this is unlikely to be realized
    if(mode<0) {G_dust = 0;} // for this choice, neglect the dust-damping term 
    double G_ion_neutral = (5.77e-11 * n_cgs * (0.97*nh0 + 0.03*nHe0) * sqrt(temperature)) * UNIT_TIME_IN_CGS; if(Z_charge_CR > 1) {G_ion_neutral /= sqrt(2.*Z_charge_CR);} // ion-neutral damping: need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]
    double G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI)*cs_thermal/4.) * sqrt(k_turb*k_L) * fturb_multiplier; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)
    double G0 = G_ion_neutral + G_turb_plus_linear_landau + G_dust; // linear terms all add into single G0 term
    double Gamma_effective = G0, phi_0 = (sqrt(M_PI)/6.)*(fabs(cos_Bgrad))*(1./(x_EB_ECR+EPSILON_SMALL))*(cs_thermal*vA_code*k_L/(CRPressureGradScaleLength*G0*G0 + EPSILON_SMALL)); // parameter which determines whether NLL dominates
    if(isfinite(phi_0) && (phi_0>0.01)) {Gamma_effective *= phi_0/(2.*(sqrt(1.+phi_0)-1.));} // this accounts exactly for the steady-state solution for the Thomas+Pfrommer formulation, including both the linear [Landau,turbulent,ion-neutral] + non-linear terms. can estimate (G_nonlinear_landau_effective = Gamma_effective - G0)
    /* with damping rates above, equilibrium transport is equivalent to pure streaming, with v_stream = vA + (diffusive equilibrium part) give by the solution below, proportional to Gamma_effective and valid to O(v^2/c^2) */
    double v_st_eff = vA_code * (1. + f_QLT * 4. * kappa_0 * Gamma_effective * x_EB_ECR * (1. + 2.*vA_code*vA_code/(C_LIGHT_CODE*C_LIGHT_CODE)) / (M_PI*vA_code*vA_code + EPSILON_SMALL)); // effective equilibrium streaming speed for all terms accounted
    return (GAMMA_COSMICRAY * v_st_eff * CRPressureGradScaleLength) * (UNIT_VEL_IN_CGS*UNIT_LENGTH_IN_CGS); // convert to effective diffusivity from 'streaming speed' [this introduces the gamma and gradient length], and convert to CGS from code units
}



/* routine which gives diffusion coefficient [in cgs] for extrinsic turbulence models. 'mode' sets whether we assume Alfven modes (mode<0), Fast-mode scattering (mode>0), or both (=0),
     0: 'default' Alfven + Fast modes (both, summing scattering rates linearly)
    -1: 'default' Alfven modes: correctly accounting for an anisotropic Goldreich-Shridar cascade, per Chandran 2000
    -2: Alfven modes in pure Goldreich-Shridhar cascade, ignoring anisotropic effects [*much* higher scattering rate, artificially]
     1: 'default' Fast modes: following Yan & Lazarian 2002, accounting for damping from viscous, ion-neutral, and other effects, and suppression if beta>1
     2: Fast modes following a pure isotropic Kolmogorov cascade down to gyro radius [*much* higher scattering rate, artificially]
 */
double diffusion_coefficient_extrinsic_turbulence(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG,
    double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion)
{
    double f_cas_ET=MAX_REAL_NUMBER, EPSILON_SMALL=1.e-50, h0_kpc=L_scale*UNIT_LENGTH_IN_KPC;
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


/* cosmic ray interactions affecting the -thermal- temperature of the gas are included in the actual cooling/heating functions;
    they are solved implicitly above. however we need to account for energy losses of the actual cosmic ray fluid, here. The
    timescale for this is reasonably long, so we can treat it semi-explicitly, as we do here.
    -- We use the estimate for combined hadronic + Coulomb losses from Volk 1996, Ensslin 1997, as updated in Guo & Oh 2008: */
void CR_cooling_and_losses(int target, double n_elec, double nHcgs, double dtime_cgs)
{
#if defined(COSMIC_RAYS_EVOLVE_SPECTRUM)
    CR_cooling_and_losses_multibin(target, n_elec, nHcgs, dtime_cgs, 0); // call the multi-binned version of this routine, passing the relevant parameters, and exit
    return;
#endif

    if(dtime_cgs <= 0) {return;} /* catch */
    int k_CRegy; double f_ion=DMAX(DMIN(Get_Gas_Ionized_Fraction(target),1.),0.);
    double a_hadronic = 6.37e-16, b_coulomb_per_GeV = 3.09e-16*(n_elec + 0.57*(1.-f_ion))*HYDROGEN_MASSFRAC; /* some coefficients; a_hadronic is the default coefficient, b_coulomb_per_GeV the default Coulomb+ionization (the two scale nearly-identically) normalization divided by GeV, b/c we need to divide the energy per CR  */
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double CR_coolrate=0, Z=return_CRbin_CR_charge_in_e(target,k_CRegy);
        if(Z > 0) /* protons here [note for now I'm using Z>0 as synonymous with protons, i.e. ignoring positrons, but we could include those if really desired */
        {
#if (N_CR_PARTICLE_BINS > 2) /* note these are currently energy-loss expressions; for truly multi-bin, probably better to work with dp/dt, instead of dE/dt */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), beta=return_CRbin_beta_factor(target,k_CRegy);
            CR_coolrate += b_coulomb_per_GeV * ((Z*Z)/(beta*E_GeV)) * nHcgs; // all protons Coulomb-interact, can be rapid for low-E
            if(E_GeV>=0.28) {CR_coolrate += a_hadronic * nHcgs;} // only GeV CRs or higher trigger above threshold for collisions
#else
            CR_coolrate = (0.87*a_hadronic + 0.53*b_coulomb_per_GeV) * nHcgs; /* for N<=2, assume a universal spectral shape, the factor here corrects for the fraction above-threshold for hadronic interactions, and 0.53 likewise for averaging  */
#endif
        } else { /* electrons here: note for electrons and positrons, always in the relativistic limit, don't need to worry about beta << 1 limits */
            /* bremsstrahlung [following Blumenthal & Gould, 1970]: dEkin/dt=4*alpha_finestruct*r_classical_elec^2*c * SUM[n_Z,ion * Z * (Z+1) * (ln[2*gamma_elec]-1/3) * E_kin */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), E_rest=0.000511, gamma=1.+E_GeV/E_rest;
            CR_coolrate += n_elec * nHcgs * 1.39e-16 * DMAX(log(2.*gamma)-0.33,0);
            /* synchrotron and inverse compton scale as dE/dt=(4/3)*sigma_Thompson*c*gamma_elec^2*(U_mag+U_rad), where U_mag and U_rad are the magnetic and radiation energy densities, respectively. Ignoring Klein-Nishina corrections here, as they are negligible at <40 GeV and only a ~15% correction up to ~1e5 GeV */
            double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG, U_rad_ev = get_cell_Urad_in_eVcm3(target);
            CR_coolrate += 5.2e-20 * gamma * (U_mag_ev + U_rad_ev); // U_mag_ev=(B^2/8pi)/(eV/cm^(-3)), here; U_rad=U_rad/(eV/cm^-3) //
        }
        
        /* here, cooling is being treated as energy loss 'within the bin'. with denser bins, should allow for movement -between- bins, using the spectral method below */
        double q_CR_cool = exp(-CR_coolrate * dtime_cgs); if(CR_coolrate * dtime_cgs > 20.) {q_CR_cool = 0;}
        SphP[target].CosmicRayEnergyPred[k_CRegy] *= q_CR_cool; SphP[target].CosmicRayEnergy[k_CRegy] *= q_CR_cool;
#ifdef COSMIC_RAYS_M1
        int k; for(k=0;k<3;k++) {SphP[target].CosmicRayFlux[k_CRegy][k] *= q_CR_cool; SphP[target].CosmicRayFluxPred[k_CRegy][k] *= q_CR_cool;}
#endif
    }
    return;
}



/* utility to estimate -locally- (without multi-pass filtering) the local Alfven Mach number */
double Get_AlfvenMachNumber_Local(int i, double vA_idealMHD_codeunits, int use_shear_corrected_vturb_flag)
{
    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,M_A,h0,EPSILON_SMALL=1.e-50; // factor which will represent which cascade model we are going to use
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
    double e_heat=0, e_CR_units_0=(SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / nHcgs; int k_CRegy;
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
            if(E_GeV>=0.28) {e_heat += f_heat_hadronic * a_hadronic * e_cr_units;} // only GeV CRs or higher trigger above threshold for collisions
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
    double cs_thermal,M_A,L_scale,vA_code,vA_noion,gizmo2gauss,Omega_per_GeV,Bmag,unit_kappa_code,b_muG,E_B,f_ion,temperature,EPSILON_SMALL; int k; k=0;
    unit_kappa_code=UNIT_VEL_IN_CGS*UNIT_LENGTH_IN_CGS; gizmo2gauss=UNIT_B_IN_GAUSS; f_ion=1; temperature=0; EPSILON_SMALL=1.e-50;
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv; cs_thermal=sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred)); /* quick thermal pressure properties (we'll assume beta=1 if MHD not enabled) */
#ifdef MAGNETIC /* get actual B-field */
    double B[3]={0}; Bmag=0; for(k=0;k<3;k++) {B[k]=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B[k]*B[k];} // B-field in code units (physical)
#endif
    Bmag=sqrt(DMAX(Bmag,0)); b_muG=Bmag*gizmo2gauss/1.e-6; b_muG=sqrt(b_muG*b_muG + 1.e-6); vA_code=sqrt(Bmag*Bmag/(SphP[i].Density*All.cf_a3inv)); vA_noion=vA_code; E_B=0.5*Bmag*Bmag*(P[i].Mass/(SphP[i].Density*All.cf_a3inv));
    Omega_per_GeV=(0.00898734*b_muG) * UNIT_TIME_IN_CGS; /* B-field in units of physical microGauss; set a floor at nanoGauss level. convert to physical code units */
#ifdef COOLING
    double ne=1, nh0=0, nHe0=0, nHepp, nhp, nHeII, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); rho_cgs=rho*UNIT_DENSITY_IN_CGS; // get thermodynamic properties
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
        double ni_m3=f_ion*(rho_cgs/PROTONMASS)/1.e-3, T6=temperature/1.e6, Lturbkpc=L_scale*UNIT_LENGTH_IN_KPC, Lgradkpc=CRPressureGradScaleLength*UNIT_LENGTH_IN_KPC, h0_fac=0.1*Get_Particle_Size(i)*All.cf_atime*All.cf_a2inv*UNIT_VEL_IN_KMS, dv2_10=0; for(k=0;k<3;k++) {int j; for(j=0;j<3;j++) {dv2_10 += SphP[i].Gradients.Velocity[j][k]*SphP[i].Gradients.Velocity[j][k]*h0_fac*h0_fac;}}
        double ecr_14 = SphP[i].CosmicRayEnergyPred[k_CRegy] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS / 1.0e-14; // CR energy density in CGS units //
        CR_kappa_streaming = GAMMA_COSMICRAY * CRPressureGradScaleLength * (v_streaming + (1./UNIT_VEL_IN_KMS)*(4.1*pow(MIN_REAL_NUMBER+ni_m3*T6,0.25)/pow(MIN_REAL_NUMBER+ecr_14*Lgradkpc,0.5) + 1.2*pow(MIN_REAL_NUMBER+dv2_10*ni_m3,0.75)/(MIN_REAL_NUMBER+ecr_14*sqrt(Lturbkpc)))); // convert to effective diffusivity
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
#if defined(COSMIC_RAYS_SET_ET_MODEL)
        scatter_modes = COSMIC_RAYS_SET_ET_MODEL; /* set to user-defined value */
#endif
        DiffusionCoeff = diffusion_coefficient_extrinsic_turbulence(scatter_modes,i,k_CRegy,M_A,L_scale,b_muG,vA_noion,rho_cgs,temperature,cs_thermal,nh0,nHe0,f_ion) / unit_kappa_code;
#endif
#if (COSMIC_RAYS_DIFFUSION_MODEL == 6) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* self-confinement-based diffusivity */
        double Omega_gyro=(0.00898734*b_muG/return_CRbin_CR_rigidity_in_GV(i,k_CRegy)) * UNIT_TIME_IN_CGS, r_L=C_LIGHT_CODE/Omega_gyro, kappa_0=r_L*C_LIGHT_CODE; // some handy numbers for limiting extreme-kappa below. all in -physical- code units //
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
void inject_cosmic_rays(double CR_energy_to_inject, double injection_velocity, int source_PType, int target, double *dir)
{
    if(CR_energy_to_inject <= 0) {return;}
    double f_injected[N_CR_PARTICLE_BINS]; f_injected[0]=1; int k_CRegy;;
#if (N_CR_PARTICLE_BINS > 1) /* add a couple steps to make sure injected energy is always normalized properly! */
    double sum_in=0.0; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]=CR_energy_spectrum_injection_fraction(k_CRegy,source_PType,injection_velocity,0); sum_in+=f_injected[k_CRegy];}
    if(sum_in>0.0) {for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]/=sum_in;}} else {for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]=1./N_CR_PARTICLE_BINS;}}
#endif
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double dEcr = CR_energy_to_inject * f_injected[k_CRegy]; // normalized properly to sum to unity
        if(dEcr <= 0) {continue;}
#if defined(COSMIC_RAYS_EVOLVE_SPECTRUM) // update the evolved slopes with the injection spectrum slope: do a simple energy-weighted mean for the updated/mixed slope here
        double E_GeV = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_CRegy), egy_slopemode = 1, xm = CR_global_min_rigidity_in_bin[k_CRegy] / CR_global_rigidity_at_bin_center[k_CRegy], xp = CR_global_max_rigidity_in_bin[k_CRegy] / CR_global_rigidity_at_bin_center[k_CRegy], xm_e=xm, xp_e=xp; // values needed for bin injection parameters
        if(CR_check_if_bin_is_nonrelativistic(k_CRegy)) {egy_slopemode=2; xm_e=xm*xm; xp_e=xp*xp;} // values needed to scale from slope injected to number and back
        double slope_inj = CR_energy_spectrum_injection_fraction(k_CRegy,source_PType,injection_velocity,1); // spectral slope of injected CRs
        double gamma_one = slope_inj + 1., xm_gamma_one = pow(xm, gamma_one), xp_gamma_one = pow(xp, gamma_one); // variables below
        double ntot_inj = (dEcr / E_GeV) * ((gamma_one + egy_slopemode) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp_e - xm_gamma_one*xm_e); // injected number in bin
        /* // below for old method where we explicitly evolved slope instead of number
        double egy_prev = SphP[target].CosmicRayEnergy[k_CRegy], etot_new = egy_prev + dEcr,
        double slope_prev = SphP[target].CosmicRay_PwrLaw_Slopes_in_Bin[k_CRegy]; // spectral slope of initial CRs
        gamma_one = slope_prev + 1.; xm_gamma_one = pow(xm, gamma_one); xp_gamma_one = pow(xp, gamma_one); // variables below
        double ntot_prev = (egy_prev / E_GeV) * ((gamma_one + egy_slopemode) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp_e - xm_gamma_one*xm_e); // CR number before injection
        SphP[target].CosmicRay_PwrLaw_Slopes_in_Bin[k_CRegy] = CR_return_slope_from_number_and_energy_in_bin(etot_new, ntot_prev + ntot_inj, E_GeV, k_CRegy);
        */
        #pragma omp atomic
        SphP[target].CosmicRay_Number_in_Bin[k_CRegy] += ntot_inj; // simply update injected number. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
#endif
        #pragma omp atomic
        SphP[target].CosmicRayEnergy[k_CRegy] += dEcr; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
        #pragma omp atomic
        SphP[target].CosmicRayEnergyPred[k_CRegy] += dEcr; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
#ifdef COSMIC_RAYS_M1
        double dir_mag=0, flux_mag=dEcr*COSMIC_RAYS_M1, dir_to_use[3]={0}; int k;
#ifdef MAGNETIC
        double B_dot_dir=0, Bdir[3]={0}; for(k=0;k<3;k++) {Bdir[k]=SphP[target].BPred[k]; B_dot_dir+=dir[k]*Bdir[k];} // the 'default' direction is projected onto B
        for(k=0;k<3;k++) {dir_to_use[k]=B_dot_dir*Bdir[k];} // launch -along- B, projected [with sign determined] by the intially-desired direction
#else
        for(k=0;k<3;k++) {dir_to_use[k]=dir[k];} // launch in the 'default' direction
#endif
        for(k=0;k<3;k++) {dir_mag += dir_to_use[k]*dir_to_use[k];}
        if(dir_mag <= 0) {dir_to_use[0]=0; dir_to_use[1]=0; dir_to_use[2]=1; dir_mag=1;}
        for(k=0;k<3;k++) {
            double dflux=flux_mag*dir_to_use[k]/sqrt(dir_mag);
            #pragma omp atomic
            SphP[target].CosmicRayFlux[k_CRegy][k]+=dflux; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
            #pragma omp atomic
            SphP[target].CosmicRayFluxPred[k_CRegy][k]+=dflux; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
        }
#endif
    }
    return;
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
    double nH_cgs = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS;
    double L_mean_free_path = (3.e25 / nH_cgs) / UNIT_LENGTH_IN_CGS;
    L_gradient_max = DMIN(L_gradient_max, L_mean_free_path);
    
    double CRPressureGradScaleLength = Get_Gas_CosmicRayPressure(i,k_CRegy) / CRPressureGradMag * All.cf_atime;
    if(CRPressureGradScaleLength > 0) {CRPressureGradScaleLength = 1.0/(1.0/CRPressureGradScaleLength + 1.0/L_gradient_max);} else {CRPressureGradScaleLength=0;}
    CRPressureGradScaleLength = sqrt(L_gradient_min*L_gradient_min + CRPressureGradScaleLength*CRPressureGradScaleLength);
    return CRPressureGradScaleLength; /* this is returned in -physical- units */
}



/* return the effective CR 'streaming' velocity for sub-grid [unresolved] models with streaming velocity set by e.g. the Alfven speed along the gradient of the CR pressure */
double Get_CosmicRayStreamingVelocity(int i)
{
    /* if we don't evolve magnetic fields, we'll assume the streaming velocity is approximately the sound speed, i.e. assume beta~1 */
    double v_streaming = sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred)); // thermal ion sound speed //
#ifdef MAGNETIC  /* since we actually evolve B-fields, it's actually the Alfven velocity: use this */
    double vA_2 = 0.0; double cs_stream = v_streaming;
    int k; for(k=0;k<3;k++) {vA_2 += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);}
    vA_2 *= All.cf_afac1 / (All.cf_atime * SphP[i].Density);
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_2 /= Get_Gas_Ionized_Fraction(i); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    v_streaming = DMIN(1.0e6*cs_stream, sqrt(MIN_REAL_NUMBER*cs_stream*cs_stream + vA_2)); // limit to Alfven speed, but put some limiters for extreme cases //
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
    double e_CR_tot = SphP[target].CosmicRayEnergyPred[k_CRegy] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * (UNIT_PRESSURE_IN_CGS);
    double ECR_per  = return_CRbin_kinetic_energy_in_GeV(target,k_CRegy) * 0.00160218; /* converts to energy in erg */
    return e_CR_tot / ECR_per;
}


/* handy functoin that just returns the B-field magnitude in microGauss, physical units. purely here to save us time re-writing this */
double get_cell_Bfield_in_microGauss(int i)
{
    double Bmag=0;
#ifdef MAGNETIC
    int k; for(k=0;k<3;k++) {double B=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B*B;} // actual B-field in code units
#else
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv; // assume equipartition
#endif
    return UNIT_B_IN_GAUSS * sqrt(DMAX(Bmag,0)) * 1.e6; // return B in microGauss
}


/* handy function that just returns the radiation energy density in eV/cm^-3, physical units. purely here to save us time re-writing this */
double get_cell_Urad_in_eVcm3(int i)
{
    double erad = 0.26*All.cf_a3inv/All.cf_atime; // default with the CMB energy density, which we assume here is a baseline minimum
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    int kfreq; double e_units = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_EV;
    for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad+=SphP[i].Rad_E_gamma_Pred[kfreq]*e_units;}
#else
    double uRad_MW = 0.31 + 0.66; /* dust (0.31) and stars (0.66) for Milky way ISRF from Draine (2011); want this to decay as we approach the IGM (where CMB totally dominates) */
    double prefac_rad=1, rho_cgs=SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; if(All.ComovingIntegrationOn) {double rhofac = rho_cgs / (1000.*COSMIC_BARYON_DENSITY_CGS);
        if(rhofac < 0.2) {prefac_rad=0;} else {if(rhofac > 200.) {prefac_rad=1;} else {prefac_rad=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off stellar background for any gas with density unless it's >1000 times the cosmic mean density
    prefac_rad *= rho_cgs/(0.01*PROTONMASS + rho_cgs); // cut off below low densities, ~1e-2
    erad += uRad_MW * prefac_rad;
#endif
    return erad;
}



/* return pre-factor for CR streaming losses, such that loss rate dE/dt = -E * streamfac */
double CR_get_streaming_loss_rate_coefficient(int target, int k_CRegy)
{
    double streamfac = 0;
#if !defined(COSMIC_RAYS_DISABLE_STREAMING) && !defined(COSMIC_RAYS_ALFVEN)
    double vstream_0 = Get_CosmicRayStreamingVelocity(target), vA=Get_Gas_Alfven_speed_i(target); /* define naive streaming and Alfven speeds */
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA /= Get_Gas_Ionized_Fraction(target); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    if(vA>0) {vstream_0 = DMIN(vA, vstream_0);} /* account for the fact that the loss term is always [or below] the Alfven speed, regardless of the bulk streaming speed */
    double L_cr = Get_CosmicRayGradientLength(target,k_CRegy), v_st_eff = SphP[target].CosmicRayDiffusionCoeff[k_CRegy] / (GAMMA_COSMICRAY * L_cr + MIN_REAL_NUMBER); // maximum possible streaming speed from combined diffusivity
    streamfac = fabs(GAMMA_COSMICRAY_MINUS1 * DMIN(v_st_eff, vstream_0) / L_cr); // if upper-limit to streaming is less than nominal 'default' v_stream/loss term, this should be lower too
#endif
    return streamfac;
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
        double eCR_0, eCR_00; eCR_00 = eCR_tmp; eCR_tmp += dCR; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy]=eCR_tmp;} else {SphP[i].CosmicRayEnergyPred[k_CRegy]=eCR_tmp;} // updated energy
        eCR_0 = eCR_tmp; // save this value for below
        
#if defined(COSMIC_RAYS_EVOLVE_SPECTRUM)
        // add update for CR number if evolved explicitly //
        if(mode==0) // only update on kicks, since we worth with a drift-conserved slope determining the ratio of N and E
        {
            double dN = SphP[i].DtCosmicRay_Number_in_Bin[k_CRegy]*dt_entr, n0 = SphP[i].CosmicRay_Number_in_Bin[k_CRegy], n_new = n0+dN;
            double E_GeV = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_CRegy), xm = CR_global_min_rigidity_in_bin[k] / CR_global_rigidity_at_bin_center[k], xp = CR_global_max_rigidity_in_bin[k] / CR_global_rigidity_at_bin_center[k], xm_e=xm, xp_e=xp;
            if(CR_check_if_bin_is_nonrelativistic(k_CRegy)) {xm_e = xm*xm; xp_e = xp*xp;} // extra power of p in energy equation accounted for here, all that's needed
            double N_min = eCR_tmp / (E_GeV * xp_e * (1.-1.e-4)); // even with arbitrarily large slopes we cannot exceed this limit: all CRs 'piled up' at highest energy
            double N_max = eCR_tmp / (E_GeV * xm_e * (1.+1.e-4)); // even with arbitrarily large slopes we cannot exceed this limit: all CRs 'piled up' at lowest energy
            n_new = DMIN(DMAX(n_new,N_min),N_max); if((n_new<0) || (isnan(n_new))) {n_new=0;}
            SphP[i].CosmicRay_Number_in_Bin[k_CRegy] = n_new; // alright, updated CR number for evolution equations
        }
#endif
        
        /* now need to account for the adiabatic heating/cooling of the 'fluid', here, with gamma=gamma_eff */
        double dCR_div;
#if 1
        double mdivvdt = -dt_entr * P[i].Particle_DivVel*All.cf_a2inv; // get locally-estimated gas velocity divergence for cells - if using non-Lagrangian method, need to modify. take negative of this [for sign of change to energy] and multiply by timestep
        if(All.ComovingIntegrationOn) {mdivvdt += -dt_entr * 3.*All.cf_hubble_a;} // include hubble-flow terms
        mdivvdt = DMAX(-1.5, DMIN(1.5, mdivvdt)); // our timestep limiter should ensure this, but problem is it responds to the -previous- timestep, so we need to impose an additional check here to prevent a numerical divergence when/if the gradients are inaccurate
        
        double fmult, Ui = u0 * P[i].Mass; // factor for multiplication below, and initial thermal energy
        double d_CR = mdivvdt * (gamma_eff-1.) * eCR_tmp; // expected CR change - this is the 'PdV' part of the work, applicable in smooth flows
        double d_Egas = mdivvdt * (GAMMA(i)-1.) * Ui; // expected gas change - in a smooth flow likewise should be entropy-conserving
        double min_IEgy = P[i].Mass * All.MinEgySpec; // minimum internal energy - in total units -
        
        double dtI_hydro = SphP[i].DtInternalEnergy * P[i].Mass * dt_entr; // change given by hydro-step computed delta_InternalEnergy
        if(mdivvdt * dtI_hydro > 0 && mode == 0) // same sign from hydro and from smooth-flow-estimator, suggests we are in a smooth flow, so we'll use stronger assumptions about the effective 'entropy' here
        {
            double abs_limit = fabs(dtI_hydro); // in smooth flow, sum of CR+gas PdV work terms should not exceed this
            double d_sum = fabs(d_CR + d_Egas); // sum of the terms estimated as above
            if(d_sum > abs_limit) {fmult = abs_limit/d_sum; d_CR *= fmult; d_Egas *= fmult;} // limit to not exceed
            
            double limforU = -0.9*Ui; // maximum change in internal energy we will allow, to prevent negative values
            if(Ui <= min_IEgy) {limforU=0;} else {limforU = DMAX(limforU,min_IEgy-Ui);} // actually more restrictive: prevent crossing the minimum enforced temperature in cooling
            if(-d_CR < limforU) {fmult = limforU/d_Egas; d_CR *= fmult; d_Egas *= fmult;} // limit
            
            double limforC = -0.9*eCR_tmp; // maximum change in CR energy we will allow, to prevent negative values
            if(d_CR < limforC) {fmult = limforC/d_CR; d_CR *= fmult; d_Egas *= fmult;} // limit
        } else { // opposite sign, suggests numerical-diffusion dominated; in this case, we will be more conservative about the limits
            double limforU = -0.5*Ui; // maximum change in internal energy we will allow, to prevent negative values
            if(Ui <= min_IEgy) {limforU=0;} else {limforU = DMAX(limforU, 0.5*(min_IEgy-Ui));} // actually more restrictive: prevent crossing the minimum enforced temperature in cooling
            if(-d_CR < limforU) {fmult = limforU/d_Egas; d_CR *= fmult; d_Egas *= fmult;} // limit
            
            double limforC = -0.5*eCR_tmp; // maximum change in CR energy we will allow, to prevent negative values
            if(d_CR < limforC) {fmult = limforC/d_CR; d_CR *= fmult; d_Egas *= fmult;} // limit
        }
        dCR_div = d_CR; // set final value
#else
        /* old treatment here */
        double d_div = (-(gamma_eff-1.) * P[i].Particle_DivVel*All.cf_a2inv) * dt_entr;
        if(All.ComovingIntegrationOn) {d_div += (-3.*(gamma_eff-1.) * All.cf_hubble_a) * dt_entr;} /* adiabatic term from Hubble expansion (needed for cosmological integrations */
        double dCR_div = DMIN(eCR_tmp*d_div , 0.5*u0*P[i].Mass); // limit so don't take away all the gas internal energy [to negative values]
        if(dCR_div + eCR_tmp < 0) {dCR_div = -eCR_tmp;} // check against energy going negative
        eCR_tmp += dCR_div; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        dCR_div = eCR_tmp - eCR_0; // actual change that is going to be applied
        if(dCR_div < -0.5*P[i].Mass*u0) {dCR_div=-0.5*P[i].Mass*u0;} // before re-coupling, ensure this will not cause negative energies
        if(dCR_div < -0.9*eCR_00) {dCR_div=-0.9*eCR_00;} // before re-coupling, ensure this will not cause negative energies
#endif
        double uf = DMAX(u0 - dCR_div/P[i].Mass , All.MinEgySpec); // final updated value of internal energy per above
        //if(mode==0) {SphP[i].InternalEnergy = uf;} else {SphP[i].InternalEnergyPred = uf;} // update gas
        if(mode==0) {SphP[i].DtInternalEnergy += (uf-u0)/(dt_entr+MIN_REAL_NUMBER);} else {SphP[i].InternalEnergyPred = uf;} // update gas
#if !defined(COSMIC_RAYS_EVOLVE_SPECTRUM)
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy] += dCR_div;} else {SphP[i].CosmicRayEnergyPred[k_CRegy] += dCR_div;} // update CRs: note if explicitly evolving spectrum, this is done separately below //
#endif
    } // loop over CR bins complete
    return 1;
}
#endif



#if defined(COSMIC_RAYS_EVOLVE_SPECTRUM)

/* this routine does the CR cooling/losses and "heating"/re-acceleration for multi-bin CR spectra: i.e. exchanging CR number
    between bins in the multi-bin approximation and modifying the spectral slope within each bin */
void CR_cooling_and_losses_multibin(int target, double n_elec, double nHcgs, double dtime_cgs, int mode_driftkick)
{
    /*! loss_mode:
     0 - hadronic+catastrophic: simply remove energy from the bin (N and E decrease together, preserving spectral slope)
     1 - adiabatic+bremstrahhlung: pure multiplicative: Edot ~ E (so instantaneously conserves slope, shifts pmin,p0,pmax, need to calculate flux)
     2 - coulomb+ionization: scale identically, for low-E protons important, messy dependence, use fitting function from Girichidis, dp/dt ~ (1 + p^(-1.9))
     3 - inverse compton+synchrotron: Edot ~ E^2, also pdot~p^2 [ultra-rel b/c e-]: modifies slope
     */
    if(dtime_cgs <= 0 || nHcgs <= 0) {return;} /* catch */
    
    /* define a bunch of general-use variables below */
    //double *bin_slopes; bin_slopes = SphP[target].CosmicRay_PwrLaw_Slopes_in_Bin; // if directly evolving slope
    double *ntot_evolved; ntot_evolved = SphP[target].CosmicRay_Number_in_Bin; double bin_slopes[N_CR_PARTICLE_BINS]; // if directly evolving number
    double *Ucr, Ucr_tot=0; if(mode_driftkick==0) {Ucr=SphP[target].CosmicRayEnergy;} else {Ucr=SphP[target].CosmicRayEnergyPred;}
    int k; for(k=0;k<N_CR_PARTICLE_BINS;k++) {Ucr_tot+=Ucr[k];} // check total energy since some fluid cells can have no CRs
    if(Ucr_tot < MIN_REAL_NUMBER) {return;} // catch - nothing to do here //
    double t=0, dt=0, E_rest_e_GeV=0.000511, f_ion=DMAX(DMIN(Get_Gas_Ionized_Fraction(target),1.),0.), b_muG=get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG, U_rad_ev=get_cell_Urad_in_eVcm3(target);

#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
    double vA=Get_Gas_Alfven_speed_i(target), vA_ion, vA_touse, kappa_max=0, kappa_i[N_CR_PARTICLE_BINS], delta_diffcoeff[N_CR_PARTICLE_BINS], reaccel_coeff[N_CR_PARTICLE_BINS]; vA_ion=vA/sqrt(f_ion); vA_touse=vA;
    for(k=0;k<N_CR_PARTICLE_BINS;k++) {kappa_i[k] = SphP[target].CosmicRayDiffusionCoeff[k]; if(kappa_i[k]>kappa_max) {kappa_max=kappa_i[k];}} // will use diffusion coefficient below
    double reaccel_coeff_0 = 2. * vA_touse*vA_touse / UNIT_TIME_IN_CGS; // below we'll adopt the standard ansatz that the momentum-space diffusion coefficient is related to the spatial coefficient by the usual heuristic expression Dxx*Dpp = <dv^2>
    if(All.Time <= All.TimeBegin || kappa_max*UNIT_LENGTH_IN_CGS*UNIT_VEL_IN_CGS < 1.e25) {reaccel_coeff_0 = 0;} // make sure kappa is initialized and has a physically meaningful value where the prescription here makes sense, or else this will give nonsense values
    reaccel_coeff_0 = DMAX(-0.9/dtime_cgs, DMIN(0.9/dtime_cgs, reaccel_coeff_0)); // our courant-type condition on the timestep should ensure this as well, but just in case, we need to enforce a condition here //
#endif
    
    double gamma_ad_eff=4./3., adiabatic_coeff = (gamma_ad_eff-1.) * (P[target].Particle_DivVel*All.cf_a2inv) / UNIT_TIME_IN_CGS ; // coefficient for adiabatic work [compression/expansion terms]. convert to physical units [a2inv], and then cgs for units here. SIGN is flipped from usual convention since we assume convention where positive coefficients = losses, for convenience with everything else below.
    if(All.ComovingIntegrationOn) {adiabatic_coeff += (gamma_ad_eff-1.) * (3.*All.cf_hubble_a) / UNIT_TIME_IN_CGS;} // adiabatic term from Hubble expansion (needed for cosmological integrations. also converted to physical, cgs, and sign convention we use here.
    double adiabatic_min = -0.5*P[target].Mass*DMAX(DMIN(SphP[target].InternalEnergyPred,SphP[target].InternalEnergy)-All.MinEgySpec,0.) / (Ucr_tot*dtime_cgs + MIN_REAL_NUMBER); if(adiabatic_coeff < adiabatic_min) {adiabatic_coeff = adiabatic_min;} // limit adiabatic -gains- of CRs (careful about sign convention here, negative means gain!) as this leads to too-large thermal losses, prevented by limiters in our step computing the exchange between CRs and gas in adiabatic calc above //
    adiabatic_coeff = DMAX(-0.9/dtime_cgs, DMIN(0.9/dtime_cgs , adiabatic_coeff)); // our timestep limiter should ensure this, but problem is it responds to the -previous- timestep, so we need to impose an additional check here to prevent a numerical divergence when/if the gradients are inaccurate

    double Ucr_i[N_CR_PARTICLE_BINS], Z[N_CR_PARTICLE_BINS], x_m[N_CR_PARTICLE_BINS], x_p[N_CR_PARTICLE_BINS], R0[N_CR_PARTICLE_BINS], E_GeV[N_CR_PARTICLE_BINS], bin_centered_rate_coeff[N_CR_PARTICLE_BINS], streaming_coeff[N_CR_PARTICLE_BINS], brems_coeff[N_CR_PARTICLE_BINS]; int NR_key[N_CR_PARTICLE_BINS];
    double hadronic_coeff = 6.37e-16 * nHcgs; // coefficient for hadronic/catastrophic interactions: dEtot/dt = -(coeff) * Etot, or dPtot/dt = -(coeff) * Ptot (since all p effected are in rel limit, and works by deleting N not by lowering individual E
    double coulomb_coeff = 3.09e-16 * nHcgs * ((n_elec + 0.57*(1.-f_ion))*HYDROGEN_MASSFRAC); // default Coulomb+ionization (the two scale nearly-identically) normalization divided by GeV, b/c we need to divide the energy per CR. needs to be multiplied by ((Z*Z)/(beta*E_GeV))
    double brems_coeff_0 = 1.39e-16  * n_elec * nHcgs; // coefficient for Bremsstrahlung [following Blumenthal & Gould, 1970]: dEkin/dt=4*alpha_finestruct*r_classical_elec^2*c * SUM[n_Z,ion * Z * (Z+1) * (ln[2*gamma_elec]-1/3) * E_kin . this needs to be multiplied by [DMAX(log(2.*gamma)-0.33,0)]; becomes dE/dt = -(coeff) * E, or dP/dt = -(coeff) * P [since all e- in rel limit]
    double synchIC_coeff_0 = 5.2e-20 * (U_mag_ev + U_rad_ev); // synchrotron and inverse compton scale as dE/dt=(4/3)*sigma_Thompson*c*gamma_elec^2*(U_mag+U_rad), where U_mag and U_rad are the magnetic and radiation energy densities, respectively. Ignoring Klein-Nishina corrections here, as they are negligible at <40 GeV and only a ~15% correction up to ~1e5 GeV. U_mag_ev=(B^2/8pi)/(eV/cm^(-3)), here; U_rad=U_rad/(eV/cm^-3). needs to be multiplied by gamma
    double e_ion_coeff = 3.60e-16 * nHcgs * (1.-f_ion), e_coulomb_coeff = 6.40e-16 * nHcgs * n_elec; // electron coulomb + ionization terms - note very similar to proton ionization (slightly different normalization b/c of log terms, and always in relativistic limit. see e.g. Ginzburg and Syrovatskii, 1964; Gould and Burbidge, 1965, Ramaty and Lingenfelter, 1966. note coulomb coeff has a very weak ~0.01*log[E_cr] scaling, but this scaling is basically offset entirely by the beta-dependence of the scaling at energies of interest, and is very weak, so we ignore it. e.g. Gould 72: Edot = (3/2)*sigma_T*ne*me*c^3/(2*beta^2) * (log[me*c^2*beta*sqrt[gamma-1]/(hbar*sqrt(4pi*e^2*ne/me))] + log[2]*(beta^2 / 2 + 1/gamma) + 1/2 + (gamma-1)^2/(16*gamma^2)
    
    double dt_min=dtime_cgs, dt_min_e=dt_min, dt_min_p=dt_min, dt_tmp, CourFac=0.4; // courant-like factor for use in subcycling here //
    int sign_flip_adiabatic_terms = 0, sign_key_for_adiabatic_loop = 1; // key that tells us if the adiabatic+brems+streaming terms have a strong sign flip, in which case we need to do 2 loops instead of 1
    double min_R_e=MAX_REAL_NUMBER, min_R_p=MAX_REAL_NUMBER, max_R_e=MIN_REAL_NUMBER, max_R_p=MIN_REAL_NUMBER; // record the minimum/max bin for e and p, for use in timestepping and limiting below (need to know if we're in the smallest bin)
    for(k=0;k<N_CR_PARTICLE_BINS;k++) /* initialize a bunch of variables for the different bins that we'll refer to below */
    {
        Ucr_i[k] = Ucr[k]; // save initial energy for reference at the end of this loop
        Z[k] = CR_global_charge_in_bin[k]; // want bin-centered values, so give index = -1
        R0[k] = CR_global_rigidity_at_bin_center[k]; // want bin-centered values, so give index = -1
        E_GeV[k] = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k); // want bin-centered values, so give index = -1
        NR_key[k] = CR_check_if_bin_is_nonrelativistic(k); // key to decide whether to use relativistic or non-relativistic scalings
        x_m[k] = CR_global_min_rigidity_in_bin[k] / R0[k]; // ratio of min-to-mid-bin CR rigidity or momentum, used for scaling everything below
        x_p[k] = CR_global_max_rigidity_in_bin[k] / R0[k]; // ratio of max-to-mid-bin CR rigidity or momentum, used for scaling everything below
        bin_slopes[k] = CR_return_slope_from_number_and_energy_in_bin(Ucr[k], ntot_evolved[k], E_GeV[k], k); // initialize slopes to use below from LUT, if not directly evolving them
        if(Z[k] < 0) {if(R0[k]<min_R_e) {min_R_e=R0[k];}} else {if(R0[k]<min_R_p) {min_R_p=R0[k];}} // record minimum rigidity to use in limits applied below
        if(Z[k] < 0) {if(R0[k]>max_R_e) {max_R_e=R0[k];}} else {if(R0[k]>max_R_p) {max_R_p=R0[k];}} // record minimum rigidity to use in limits applied below

        if(Z[k] < 0) {brems_coeff[k] = brems_coeff_0 * DMAX(log(2.*(1.+E_GeV[k]/E_rest_e_GeV))-0.33,0);} else {brems_coeff[k]=0;}
        streaming_coeff[k] = CR_get_streaming_loss_rate_coefficient(target,k) / UNIT_TIME_IN_CGS;
        
        // calculate the timestep limit from all possible bins, maximum step. note no constraint from hadronic here b/c cannot 'cross the bin'
        bin_centered_rate_coeff[k] = 0;
        double adiab_brem_coeff = adiabatic_coeff + streaming_coeff[k] + brems_coeff[k]; bin_centered_rate_coeff[k]+=adiab_brem_coeff; // do constraint from adiabatic + Bremsstrahlung
        if(adiab_brem_coeff < 0) {sign_key_for_adiabatic_loop=-1;} // note the net sign of this term for use below
        if(k>0) {if(adiab_brem_coeff * (adiabatic_coeff + streaming_coeff[k-1] + brems_coeff[k-1]) < 0) {sign_flip_adiabatic_terms = 1;}} // have a sign flip, and its not negligible in magnitude //
        if(Z[k] < 0) // e-: do constraint from synchrotron + IC
        {
            double IC_sync_coeff = (E_GeV[k]/E_rest_e_GeV) * synchIC_coeff_0; bin_centered_rate_coeff[k]+=IC_sync_coeff;
            double ion_coeff = (e_coulomb_coeff + (1.+0.07*log(E_GeV[k])) * e_ion_coeff) / R0[k]; bin_centered_rate_coeff[k]+=ion_coeff; // relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc.
        } else { // p: do constraint from Coulomb + ionization
            if(NR_key[k]==1) // bin is in non-relativistic limit, use those expressions
            {
                double A=1.; if(Z[k]>1) {A=2.*Z[k];}
                double Coul_coeff = ((0.88 * A*A) / (R0[k]*R0[k]*R0[k] * Z[k])) * coulomb_coeff; bin_centered_rate_coeff[k]+=Coul_coeff; // non-relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc. should be multiplied by atomic weight A^2 as well.
            } else { // bin is in relativistic limit, use those expressions
                double Coul_coeff = (Z[k]/R0[k]) * coulomb_coeff; bin_centered_rate_coeff[k]+=Coul_coeff; // relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc.
            }
        }
    }
    if(sign_flip_adiabatic_terms==1) {sign_key_for_adiabatic_loop=-1;} // have sign-flips, so adiabatic term -must- have the oppose sign. otherwise -no- sign flips, so just follow the last sign recorded above

#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
    for(k=0;k<N_CR_PARTICLE_BINS;k++) // additional variables need to be initialized here, after the previous loop definitions
    {
        int minbin_flag=0, maxbin_flag=0;
        if(k>0) {if(Z[k-1] != Z[k]) {minbin_flag=1;}} else {minbin_flag=1;} // previous bin doesn't exist or has opposite sign: lowest-E bin for type
        if(k<N_CR_PARTICLE_BINS-1) {if(Z[k+1] != Z[k]) {maxbin_flag=1;}} else {maxbin_flag=1;} // next bin doesn't exist or has opposite sign: highest-E bin for type
        /* need to estimate the log-slope of the spatial diffusion coefficient vs momentum: delta = dln[kappa] / dln[R], for each species type */
        if(minbin_flag) {
            delta_diffcoeff[k] = log((kappa_i[k+1] + MIN_REAL_NUMBER) / (kappa_i[k] + MIN_REAL_NUMBER)) / log((R0[k+1] + MIN_REAL_NUMBER) / (R0[k] + MIN_REAL_NUMBER)); // linear estimate from next bin
        } else if(maxbin_flag) {
            delta_diffcoeff[k] = log((kappa_i[k] + MIN_REAL_NUMBER) / (kappa_i[k-1] + MIN_REAL_NUMBER)) / log((R0[k] + MIN_REAL_NUMBER) / (R0[k-1] + MIN_REAL_NUMBER)); // linear estimate from previous bin
        } else {
            double y_m = log((kappa_i[k-1] + MIN_REAL_NUMBER) / (kappa_i[k] + MIN_REAL_NUMBER)), y_p = log((kappa_i[k+1] + MIN_REAL_NUMBER) / (kappa_i[k] + MIN_REAL_NUMBER));
            double x_m = log((R0[k-1] + MIN_REAL_NUMBER) / (R0[k] + MIN_REAL_NUMBER)), x_p = log((R0[k+1] + MIN_REAL_NUMBER) / (R0[k] + MIN_REAL_NUMBER));
            delta_diffcoeff[k] = (y_p * x_m*x_m - y_m * x_p*x_p) / (x_m * x_p * (x_m - x_p)); // second-order log-slope estimate
        }
        double delta_limit=0.0315059; delta_diffcoeff[k] = DMIN(DMAX(delta_diffcoeff[k], delta_limit),2.-delta_limit); // need to restrict delta value for the quasi-linear theory model below to make any sense (otherwise just gives unphysical answers because relevant integrals all diverge)
        double alpha_denom = delta_diffcoeff[k] * (4.-delta_diffcoeff[k]) * (4.-delta_diffcoeff[k]*delta_diffcoeff[k]); // quasi-linear isotropic turbulent-type assumption for coefficient relating Dpp and Dxx
        reaccel_coeff[k] = -delta_diffcoeff[k] * reaccel_coeff_0 / ((kappa_i[k] + MIN_REAL_NUMBER) * (alpha_denom + MIN_REAL_NUMBER)); // up to a factor of (R/R0)^(-delta) * (1 - slope_gamma/2), this gives the rate in 1/seconds of the effective momentum-space diffusion. note sign here: generally implies energy gain.
        bin_centered_rate_coeff[k] += reaccel_coeff[k] * (1. - 0.5*DMIN(bin_slopes[k], 2.)); // limit this to ensure we don't miss if this should be larger
    }
#endif

    
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        if((Ucr[k] < 1.e-10 * Ucr_tot) || (bin_centered_rate_coeff[k]==0)) {continue;} // don't bother with timestep limits if the bin contains totally negligible fraction of CR energy
        //if(Z[k] < 0) {if(R0[k] <= 1.001*min_R_e) {continue;}} else {if(R0[k] <= 1.001*min_R_p) {continue;}} // don't need to apply the same timestep in the smallest bin, so long as regulated by the bin above it, since we'll just cool to 0 here in finite time and that's ok
        if(All.ComovingIntegrationOn) {if(All.Time < 0.12) {if(Z[k] < 0) {if(R0[k] >= 0.999*max_R_e) {continue;}}}} // don't need to apply the same timestep in the smallest bin, so long as regulated by the bin above it, since we'll just cool to 0 here in finite time and that's ok
        double abs_bin_coeff_limit = 0.05 * fabs(bin_centered_rate_coeff[k]); // set threshold for fraction of rate where we need to worry about detailed subcycling: sub-dominant processes not important here. find few percent works well here.
        double adiab_brem_coeff = fabs(adiabatic_coeff + streaming_coeff[k] + brems_coeff[k]); // do constraint from adiabatic + Bremsstrahlung
        if(adiab_brem_coeff > abs_bin_coeff_limit) {dt_tmp = CourFac * log(x_p[k]/x_m[k]) / adiab_brem_coeff; if(Z[k]<0) {dt_min_e=DMIN(dt_min_e, dt_tmp);} else {dt_min_p=DMIN(dt_min_p, dt_tmp);}}
        if(Z[k] < 0) // e-: do constraint from synchrotron + IC
        {
           double IC_sync_coeff = (E_GeV[k]/E_rest_e_GeV) * synchIC_coeff_0;
           if(IC_sync_coeff > abs_bin_coeff_limit) {dt_tmp = CourFac * (1./x_m[k] - 1./x_p[k]) / IC_sync_coeff; dt_min_e=DMIN(dt_min_e, dt_tmp);}
           double ion_coeff = (e_coulomb_coeff + (1.+0.07*log(E_GeV[k])) * e_ion_coeff) / R0[k]; // relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc.
           if(ion_coeff > abs_bin_coeff_limit) {dt_tmp = CourFac * DMIN(x_p[k]-x_m[k], x_m[k]) / ion_coeff; dt_min_e=DMIN(dt_min_e, dt_tmp);}
        } else { // p: do constraint from Coulomb + ionization
           if(NR_key[k]==1) // bin is in non-relativistic limit, use those expressions
           {
               double A=1.; if(Z[k]>1) {A=2.*Z[k];}
               double Coul_coeff = ((0.88 * A*A) / (R0[k]*R0[k]*R0[k] * Z[k])) * coulomb_coeff; // non-relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc. should be multiplied by atomic weight A^2 as well.
               if(Coul_coeff > abs_bin_coeff_limit) {dt_tmp = CourFac * (x_p[k]*x_p[k]*x_p[k] - x_m[k]*x_m[k]*x_m[k]) / (3.*Coul_coeff); dt_min_p=DMIN(dt_min_p, dt_tmp);}
           } else { // bin is in relativistic limit, use those expressions
               double Coul_coeff = (Z[k]/R0[k]) * coulomb_coeff; // relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc.
               if(Coul_coeff > abs_bin_coeff_limit) {dt_tmp = CourFac * DMIN(x_p[k]-x_m[k], x_m[k]) / Coul_coeff; dt_min_p=DMIN(dt_min_p, dt_tmp);}
           }
        }
#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
        double rcoeff_bin = fabs(reaccel_coeff[k]) * (1.-0.5*DMIN(bin_slopes[k],1.)); // limit b/c this might change in step, then estimate rate coefficient at bin center
        if(rcoeff_bin > abs_bin_coeff_limit) {dt_tmp = CourFac * (pow(x_p[k],delta_diffcoeff[k]) - pow(x_m[k],delta_diffcoeff[k])) / rcoeff_bin;} // set timestep limit
        if(Z[k]<0) {dt_min_e=DMIN(dt_min_e, dt_tmp);} else {dt_min_p=DMIN(dt_min_p, dt_tmp);} // applies to both e and p
#endif
    }
    if(All.ComovingIntegrationOn) {if(Ucr_tot < 1.e-6*P[target].Mass*SphP[target].InternalEnergy) {dt_min_e*=10.; dt_min_p*=10.; dt_min_e=DMAX(dt_min_e,0.01*dtime_cgs); dt_min_p=DMAX(dt_min_p,0.01*dtime_cgs);}} // allow larger slope errors when the CR energy is a negligible fraction of total
    
    double dt_target = DMIN(dtime_cgs, DMIN(dt_min_e,dt_min_p)); // timescale for subcycling, using constraint above. limit for extreme cases where we might run into expense problems
    if(dt_target < 1.e-4*dtime_cgs)
    {
        printf("WARNING: timestep for subcycling wants to exceed limit: dt_min_e/p=%g/%g dt_tot=%g \n",dt_min_e,dt_min_p,dtime_cgs);
#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
        printf(" ID=%llu mode=%d nH=%g ne=%g vA=%g vA_ion=%g dt=%g Utot=%g hadronic_coeff=%g coulomb_coeff=%g brems_coeff_0=%g synchIC_coeff_0=%g (U_mag_eV=%g U_rad_eV=%g) adiabatic_coeff=%g dtmin=%g signflip=%d signkey=%d \n",P[target].ID,mode_driftkick,nHcgs,n_elec,vA,vA_ion,dtime_cgs,Ucr_tot,hadronic_coeff,coulomb_coeff,brems_coeff_0,synchIC_coeff_0,U_mag_ev,U_rad_ev,adiabatic_coeff,DMIN(dt_min_e,dt_min_p),sign_flip_adiabatic_terms,sign_key_for_adiabatic_loop);
        for(k=0;k<N_CR_PARTICLE_BINS;k++) {printf("  k=%d U=%g Z=%g R0=%g E=%g NR=%d xm=%g xp=%g slope=%g kappa=%g brem_c=%g stream_c=%g reacc_c=%g delta_DxxSlope=%g IC_sync_c=%g el_ion_coul_c=%g p_ion_coul_c_R/NR=%g/%g binratec=%g \n",k,Ucr[k],Z[k],R0[k],E_GeV[k],NR_key[k],x_m[k],x_p[k],bin_slopes[k],kappa_i[k],brems_coeff[k],streaming_coeff[k],reaccel_coeff[k],delta_diffcoeff[k],(E_GeV[k]/E_rest_e_GeV) * synchIC_coeff_0,(e_coulomb_coeff + (1.+0.07*log(E_GeV[k])) * e_ion_coeff) / R0[k], (Z[k]/R0[k]) * coulomb_coeff,0.88/ (R0[k]*R0[k]*R0[k] * Z[k]) * coulomb_coeff,bin_centered_rate_coeff[k]);}
#else
        printf(" ID=%llu mode=%d nH=%g ne=%g dt=%g Utot=%g hadronic_coeff=%g coulomb_coeff=%g brems_coeff_0=%g synchIC_coeff_0=%g (U_mag_eV=%g U_rad_eV=%g) adiabatic_coeff=%g dtmin=%g signflip=%d signkey=%d \n",P[target].ID,mode_driftkick,nHcgs,n_elec,dtime_cgs,Ucr_tot,hadronic_coeff,coulomb_coeff,brems_coeff_0,synchIC_coeff_0,U_mag_ev,U_rad_ev,adiabatic_coeff,DMIN(dt_min_e,dt_min_p),sign_flip_adiabatic_terms,sign_key_for_adiabatic_loop);
        for(k=0;k<N_CR_PARTICLE_BINS;k++) {printf("  k=%d U=%g Z=%g R0=%g E=%g NR=%d xm=%g xp=%g slope=%g bremc=%g strmc=%g binratec=%g \n",k,Ucr[k],Z[k],R0[k],E_GeV[k],NR_key[k],x_m[k],x_p[k],bin_slopes[k],brems_coeff[k],streaming_coeff[k],bin_centered_rate_coeff[k]);}
#endif
        double dt_min_tmp = 1.e-4 * dtime_cgs; // enforce this since our integrators can deal with it, just at some loss of accuracy
        if(dt_target < dt_min_tmp) {dt_min_e=DMAX(dt_min_e, dt_min_tmp); dt_min_p=DMAX(dt_min_p, dt_min_tmp);}
    }

    int proton_key; for(proton_key=0;proton_key<=1;proton_key++) // loop over whether we consider nuclei or electrons, first //
    {
        int n_active=0, bins_sorted[N_CR_PARTICLE_BINS]; double R0_bins[N_CR_PARTICLE_BINS];
        for(k=0;k<N_CR_PARTICLE_BINS;k++) {if((Z[k]>0 && proton_key==1) || (Z[k]<0 && proton_key==0)) {bins_sorted[n_active]=k; R0_bins[n_active]=R0[k]; n_active++;}} // bin is valid: charge matches that desired
        if(n_active<=0) {continue;} // nothing to do here
        if(n_active<N_CR_PARTICLE_BINS) {for(k=n_active;k<N_CR_PARTICLE_BINS;k++) {R0_bins[k]=MAX_REAL_NUMBER;}} // set a dummy value here for sorting purposes below
        //qsort(bins_sorted, n_active, sizeof(int), compare_CR_rigidity_for_sort); // sort on energies from smallest-to-largest [this is hard-coded by requiring the list go in monotonic increasing order for e and p, regardless of how the e and p are themselves ordered //

        t = 0; // reset this before we enter the time integration loop below! //
        if(proton_key==0) {dt_target = DMIN(dtime_cgs, dt_min_e);} else {dt_target = DMIN(dtime_cgs, dt_min_p);} // allow separate timesteps for e and p, since there is no direct exchange between them in the subcycle below
        while(t < dtime_cgs)
        {
            dt = DMIN(dt_target , dtime_cgs - t); // set the subcycle step size
            if(dt <= 0) {break;} // we have reached the end of the timestep - exit loop

            int loss_mode; // this will determine which type[s] of losses [differentiated by their qualitative scalings] we are considering //
            for(loss_mode=0;loss_mode<=5;loss_mode++)
            {
                if(loss_mode==0) {if(proton_key==0) {continue;}} // loss-mode=0 [Hadronic] uses only protons, skip to next in loop
                if(loss_mode==3) {if(proton_key==1) {continue;}} // loss-mode=3 [Compton+Synchrotron] uses only electrons
                if(loss_mode==4) {if(sign_flip_adiabatic_terms == 0) {continue;}} // adiabatic+streaming+brems term walked in same order, so we don't need to do an additional loop here
#if !defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
                if(loss_mode==5) {continue;}
#endif
                
                int order = -1; // default to -descending- energy order [for energy-loss]. but for adiabatic terms may need to use -ascending- order if energy increases
                if(loss_mode==1) {if(sign_key_for_adiabatic_loop < 0) {order = 1;}} // adiabatic: if net energy -gain- in this step [only step where its possible], then switch to ascending order
                if(loss_mode==5) {order = 1;} // reacceleration always produces increasing energies (in the parameter space of interest)
                
                double dn_flux=0, de_flux=0;
                for(k=0;k<n_active;k++)
                {
                    int j = bins_sorted[k]; // target bin for flux, will move 'up the ladder' passing fluxes to higher energies
                    if(order == -1) {j = bins_sorted[n_active-1-k];} // will move 'down the ladder' passing fluxes to lower energies
                    
                    if(loss_mode==0) // hadronic+catastrophic losses.
                    {
                        if(E_GeV[j] > 0.28) {double fac=exp(-DMIN(hadronic_coeff*dt, 60.)); Ucr[j]*=fac; ntot_evolved[j]*=fac;} // only >~GeV trigger threshold for collisions //
                        dn_flux=0; de_flux=0; // make sure these are zero'd for next step
                        continue; // these -destroy- CRs, decreasing N and E identically [ignoring secondary production for now], leaving slope un-modified and -no- bin-to-bin fluxes: easy to solve, just modify total energy+number identically in the bin //
                    }
                    
                    double etot = Ucr[j]; // total energy in bin, one of our key evolved variables
                    double E_bin_NtoE = E_GeV[j]; // use this to normalize the number in a convenient unit (just to keep easier to track, units are arbitrary). save it so we don't forget below.
                    double slope_gamma = bin_slopes[j]; // get the initial slope if we evolve it and can pass it as a conserved quantity
                    double ntot = ntot_evolved[j]; // total number in bin [in our effective units] //

                    double xm = x_m[j], xp = x_p[j], gamma_one = slope_gamma+1., xe_gamma_one = 0, xm_gamma_one = 0, xp_gamma_one = 0; // just for convenience below
                    /* // old mode below, where we directly evolve slope and derived number from it
                    double ntot = 0; // total number in bin. we use the slope we obtain to convert between number and energy, given the bin-centered values
                    xm_gamma_one = pow(xm, gamma_one); xp_gamma_one = pow(xp, gamma_one); // variables needed for conversion below
                    if(NR_key[j]) {ntot = ((gamma_one + 2.) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp*xp - xm_gamma_one*xm*xm);} // non-relativistic map between E and N
                        else {ntot = ((gamma_one + 1.) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp - xm_gamma_one*xm);} // relativistic map between E and N
                    ntot *= (etot / E_bin_NtoE); // dimensional normalization. note the units are arbitrary for the E_bin_NtoE term, as long as we are consistent [we can evolve constant x N, instead of N, for convenience]
                    */
                    
                    double dn_flux_fromprevbin = dn_flux, de_flux_fromprevbin = de_flux; // fluxes coming from the last bin. note these have -already- been cooled, so don't want to double-count that operation, hence we add them at the end of the sub-step
                    double extra_var_topass_for_xe = 0; // dummy variable to pass to subroutine below
                    dn_flux = 0; de_flux = 0; // reset these, they will be re-defined below but in case we skip the relevant loop they need to be zeroed
                    if((etot <= 0 && de_flux_fromprevbin <= 0) || (ntot <= 0 && dn_flux_fromprevbin <= 0)) {continue;} // no CRs to actually work with here [nothing injected yet]!
                    
                    double rate_prefac = 0; // coefficient for the bin-centered loss rate from this mode [-negative- loss rate = energy-gain]
                    if(loss_mode==1 || loss_mode==4) // adiabatic + brems + streaming
                    {
                        if(loss_mode==1) {rate_prefac = adiabatic_coeff;} // adiabatic always in mode=1
                        if((sign_flip_adiabatic_terms==0) || (loss_mode==4)) {rate_prefac += brems_coeff[j] + streaming_coeff[j];} // no sign-flip, or in mode=4 so we do add the strictly-loss terms
                    }
                    if(loss_mode==2) // coulomb + ion
                    {
                        if(Z[j] < 0) // electron ionization losses here
                        {
                            rate_prefac = (e_coulomb_coeff + (1.+0.07*log(E_GeV[j])) * e_ion_coeff) / R0[j]; // always in relativistic limit here. 1/R0 is b/c this coefficient is defined normalized to the bin center in GeV, to make the equations dimensionless
                        } else { // proton ionization + Coulomb losses here
                            rate_prefac = (Z[j]/R0[j]) * coulomb_coeff; // relativistic expression. note this is equation for p evolution, where R is rigidity, so need to be careful with Z factors, etc.
                            if(NR_key[j]) {double A=1.; rate_prefac *= 0.88 *A*A / (R0[j]*R0[j]*Z[j]*Z[j]);} // non-relativistic expression (times constant (A/(R0*Z))^2 //
                        }
                    }
                    if(loss_mode==3) {rate_prefac = (E_GeV[j]/E_rest_e_GeV) * synchIC_coeff_0;} // IC + synch
#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
                    if(loss_mode==5) {if(slope_gamma>2.) {rate_prefac=0;} else {rate_prefac=reaccel_coeff[j]*(1.-0.5*slope_gamma); extra_var_topass_for_xe=delta_diffcoeff[j];}} // we will treat gamma as constant over this integration sub-step
#endif
                    
                    int do_cooling_ops = 1; // option to skip the cooling subcycle part of the step here:
                    if(fabs(rate_prefac) < 0.01*fabs(bin_centered_rate_coeff[j])+MIN_REAL_NUMBER) {do_cooling_ops=0;} // this loss mode is negligible here, skip it [pure optimization]
                    if(fabs(bin_centered_rate_coeff[j])*dt < 1.e-16 * Ucr[j]) {do_cooling_ops=0;} // total loss rate is so small we won't be able to track it to floating-point accuracy, skip it
                    if(loss_mode==1 || loss_mode==4) {if((rate_prefac < 0 && order == -1) || (rate_prefac > 0 && order == 1)) {do_cooling_ops=0;}} // adiabatic: make sure we are operating in correct order (will ignore bin if we have a sign switch, which is ok, since it must be nearly-null in that case
                    
                    if(do_cooling_ops) // enter the 'cooling' (bin-to-bin flux calculations) portion of the subcycle
                    {
                        double rate_dt = rate_prefac * dt; // dimensionless step size with rate prefactor times timestep
                        double x_to_edge = CR_return_new_bin_edge_from_rate(rate_dt, xm, xp, loss_mode, NR_key[j], extra_var_topass_for_xe); // returns the dimensionless 'x' representing the particles furthest from bin edge that 'reach' the edge by end of this sub-step
                        gamma_one = slope_gamma+1.; xe_gamma_one = pow(x_to_edge, gamma_one); xm_gamma_one = pow(xm, gamma_one); xp_gamma_one = pow(xp, gamma_one);
                        
                        double dn_lostfrombin = 0; // calculate the number flux integrating out to that new bin edge
                        if(rate_prefac < 0) {dn_lostfrombin = ntot * (1 - xe_gamma_one / xp_gamma_one) / (1 - xm_gamma_one / xp_gamma_one);}
                            else {dn_lostfrombin = ntot * (xe_gamma_one / xm_gamma_one - 1.) / (xp_gamma_one / xm_gamma_one - 1.);}
                        
                        double etot_final_inbin = etot; // calculate the final energy of all particles still inside bin bounds at end of the sub-step
                        double etot_final_allparticlesfrombin = etot; // calculate the final energy of all particles that began sub-step within bin bounds [difference is net flux of energy out-of-bin]
                        
                        if(loss_mode==1 || loss_mode==4) // adiabatic + brems
                        {
                            double norm_fac = exp(-rate_dt), x1, x0, x1_gamma_one, x0_gamma_one; // now make sure we flip the sign convention back to normal for rate
                            if(rate_dt < 0) {x1=x_to_edge; x1_gamma_one=xe_gamma_one; x0=xm; x0_gamma_one=xm_gamma_one;} else {x1=xp; x1_gamma_one=xp_gamma_one; x0=x_to_edge; x0_gamma_one=xe_gamma_one;} // this is everything that remains in bin. so if decaying, runs from x_m' to xp; if growing [rate < 0] runs from xm to x_p'
                            if(NR_key[j])
                            {
                                norm_fac *= norm_fac; // exp(2*rate*dt) because of KE going as p^2
                                etot_final_inbin = etot * norm_fac * (x1_gamma_one*x1*x1 - x0_gamma_one*x0*x0) / (xp_gamma_one*xp*xp - xm_gamma_one*xm*xm); // extra power of x from extra power of p in KE
                            } else {
                                etot_final_inbin = etot * norm_fac * (x1_gamma_one*x1 - x0_gamma_one*x0) / (xp_gamma_one*xp - xm_gamma_one*xm); // relativistic limit expression
                            }
                            etot_final_allparticlesfrombin = etot * norm_fac; // easy because all particles lose the same fractional energy
                        }
                        if(loss_mode==2) // coulomb + ion
                        {
                            if(NR_key[j])
                            {
                                double norm_fac = etot * (gamma_one+2.) / (xp_gamma_one*xp*xp - xm_gamma_one*xm*xm); // get the pre-factor for the integral from numbers we already have
                                double xmin_remains = pow(3.*rate_dt, 1./3.); // sets absolute minimum value of x for which any energy remains at all, needs to be integrated only from this bound upwards
                                int n_int=10; double xmin=DMAX(xm,xmin_remains), u_int=0, u_int_xe=0, dlnx=log(xp/xm)/((double)n_int), exp_fac=exp(dlnx/2.), x=xmin*exp_fac, f0=CR_coulomb_energy_integrand(xmin,rate_dt,gamma_one), f1, exp_fac_2=exp_fac*exp_fac, x0=xmin, x1;
                                while(x < xp)
                                {
                                    x1 = DMIN(x*exp_fac, xp);
                                    f1 = CR_coulomb_energy_integrand( x1 , rate_dt, gamma_one);
                                    double u_int_fac = log(x1/x0) * (f1 - f0) / (log((f1+MIN_REAL_NUMBER)/(f0+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);
                                    u_int += u_int_fac;
                                    if(x0 < x_to_edge) {if(x1 > x_to_edge) {double f1_e=CR_coulomb_energy_integrand(x_to_edge,rate_dt,gamma_one); u_int_xe+=log(x_to_edge/x0)*(f1_e-f0)/(log((f1_e+MIN_REAL_NUMBER)/(f0+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);} else {u_int_xe+=u_int_fac;}}
                                    x *= exp_fac_2; f0 = f1; x0 = x1;
                                }
                                etot_final_inbin = (u_int - u_int_xe) * norm_fac; etot_final_allparticlesfrombin = u_int * norm_fac; // normalize these appropriately
                            } else {
                                double d_egy_fac = rate_dt * ntot * E_bin_NtoE; // re-normalizes ntot correctly back in our chosen units to the -same- units as etot, for use below
                                etot_final_inbin = etot * (xp_gamma_one*xp - xe_gamma_one*x_to_edge) / (xp_gamma_one*xp - xm_gamma_one*xm)
                                    - d_egy_fac * (xp_gamma_one - xe_gamma_one) / (xp_gamma_one - xm_gamma_one); // relativistic limit expression: accounts both for particles leaving the bin [first term] and total energy lost by all particles that stay in the bin [second term]
                                etot_final_allparticlesfrombin = etot - d_egy_fac; // easy because all particles lose the same absolute energy
                            }
                        }
                        if(loss_mode==3) // IC + synch
                        {
                            double norm_fac = etot * (gamma_one+1.) / (xp_gamma_one*xp - xm_gamma_one*xm); // get the pre-factor for the integral from numbers we already have
                            int n_int=10; double u_int=0, u_int_xe=0, dlnx=log(xp/xm)/((double)n_int), exp_fac=exp(dlnx/2.), x=xm*exp_fac, f0=CR_compton_energy_integrand(xm,rate_dt,gamma_one), f1, exp_fac_2=exp_fac*exp_fac, x0=xm, x1;
                            while(x < xp)
                            {
                                x1 = DMIN(x*exp_fac, xp);
                                f1 = CR_compton_energy_integrand( x1 , rate_dt, gamma_one);
                                double u_int_fac = log(x1/x0) * (f1 - f0) / (log((f1+MIN_REAL_NUMBER)/(f0+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);
                                u_int += u_int_fac;
                                if(x0 < x_to_edge) {if(x1 > x_to_edge) {double f1_e=CR_compton_energy_integrand(x_to_edge,rate_dt,gamma_one); u_int_xe+=log(x_to_edge/x0)*(f1_e-f0)/(log((f1_e+MIN_REAL_NUMBER)/(f0+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);} else {u_int_xe+=u_int_fac;}}
                                x *= exp_fac_2; f0 = f1; x0 = x1;
                            }
                            etot_final_inbin = (u_int - u_int_xe) * norm_fac; etot_final_allparticlesfrombin = u_int * norm_fac; // normalize these appropriately
                        }
                        if(loss_mode==5) // diffusive re-acceleration
                        {
#if defined(COSMIC_RAYS_DIFFUSIVE_REACCELERATION)
                            double norm_fac = etot * (gamma_one+1.) / (xp_gamma_one*xp - xm_gamma_one*xm); // get the pre-factor for the integral from numbers we already have
                            if(NR_key[j]) {norm_fac = etot * (gamma_one+2.) / (xp_gamma_one*xp*xp - xm_gamma_one*xm*xm);} // use correct NR form
                            double delta_j = delta_diffcoeff[j], rate_dt_abs = fabs(rate_dt); // value of delta-slope needed below
                            int n_int=10; double xmin=xm, u_int=0, u_int_xe=0, dlnx=log(xp/xm)/((double)n_int), exp_fac=exp(dlnx/2.), x=xmin*exp_fac, f0=CR_reaccel_energy_integrand(xmin,rate_dt_abs,gamma_one,delta_j,NR_key[j]), f1, exp_fac_2=exp_fac*exp_fac, x0=xmin, x1;
                            while(x < xp)
                            {
                                x1 = DMIN(x*exp_fac, xp);
                                f1 = CR_reaccel_energy_integrand( x1 , rate_dt_abs , gamma_one , delta_j , NR_key[j]);
                                double u_int_fac = log(x1/x0) * (f1 - f0) / (log((f1+MIN_REAL_NUMBER)/(f0+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);
                                u_int += u_int_fac;
                                if(x1 > x_to_edge) {if(x0 < x_to_edge) {double f0_e=CR_reaccel_energy_integrand(x_to_edge, rate_dt_abs, gamma_one, delta_j , NR_key[j]); u_int_xe += log(x1/x_to_edge)*(f1-f0_e)/(log((f1+MIN_REAL_NUMBER)/(f0_e+MIN_REAL_NUMBER)) + MIN_REAL_NUMBER);} else {u_int_xe+=u_int_fac;}} // // integrate from x_edge to f1, or x0+x1 both > x_to_edge, so whole bin is outside range
                                x *= exp_fac_2; f0 = f1; x0 = x1;
                            }
                            etot_final_inbin = (u_int - u_int_xe) * norm_fac; etot_final_allparticlesfrombin = u_int * norm_fac; // normalize these appropriately
#endif
                        }
                        
                        dn_lostfrombin = DMIN(dn_lostfrombin , (1.0-1.e-8)*ntot); // limit to prevent 0's which will give nan's
                        etot_final_inbin = DMAX(etot_final_inbin , 1.e-8*etot); // limit to prevent 0's which will give nan's
                        dn_flux = dn_lostfrombin; // set total outgoing flux of number of particles to next bin
                        de_flux = DMAX(etot_final_allparticlesfrombin - etot_final_inbin, 0); // determine total outgoing flux of energy to next bin

                        ntot -= dn_lostfrombin; // update number
                        etot = etot_final_inbin; // update energy
                    } // close the bin-to-bin flux calculation portion of the subcycle

                    ntot += dn_flux_fromprevbin; // update the bin with the -already cooled- CRs from the previous bin
                    etot += de_flux_fromprevbin; // update the bin with the -already cooled- CRs from the previous bin
                    if(ntot > 0 && etot > 0) {slope_gamma = CR_return_slope_from_number_and_energy_in_bin(etot, ntot, E_bin_NtoE, j);} else {slope_gamma=0;} // get the updated slope for this bin

                    ntot_evolved[j] = ntot; // set finalized updated number in bin
                    Ucr[j] = etot; // set finalized updated energy in bin
                    bin_slopes[j] = slope_gamma; // set finalized updated slope [or N, if that's the conserved quantity we evolve]
                } // loop over active bins
            } // loop over loss mode
            t += dt; // update time
        } // loop over charge
    } // loop over time
    
    if(mode_driftkick==0) {for(k=0;k<N_CR_PARTICLE_BINS;k++) {double fac=1; if(Ucr_i[k]>0 && Ucr[k]>0) {fac=DMAX(1.e-10,DMIN(1.e10,Ucr[k]/Ucr_i[k]));} SphP[target].CosmicRayEnergyPred[k] *= fac;}}
#ifdef COSMIC_RAYS_M1 /* update fluxes. these will self-adapt quickly but should feel losses as well, so we do that here */
    int kdir; for(k=0;k<N_CR_PARTICLE_BINS;k++) {double fac=1; if(Ucr_i[k]>0 && Ucr[k]>0) {fac=DMAX(1.e-2,DMIN(1.e2,Ucr[k]/Ucr_i[k]));}
        for(kdir=0;kdir<3;kdir++) {if(mode_driftkick==0) {SphP[target].CosmicRayFlux[k][kdir] *= fac; SphP[target].CosmicRayFluxPred[k][kdir] *= fac;} else {SphP[target].CosmicRayFluxPred[k][kdir] *= fac;}}}
#endif
    return; // all done!
}




/* initialize CR quantities needed specifically for our multi-bin spectral methods */
void CR_initialize_multibin_quantities(void)
{
    if(ThisTask==0) {printf("Initializing global cosmic ray spectral variables:\n"); fflush(stdout);}
    int k; CR_spectrum_define_bins(); // call this to define the actual list of bins, needs to be done before basically anything else!
    double R0[N_CR_PARTICLE_BINS], Z[N_CR_PARTICLE_BINS]; for(k=0;k<N_CR_PARTICLE_BINS;k++) {R0[k]=CR_global_rigidity_at_bin_center[k]; Z[k]=CR_global_charge_in_bin[k];}
    double R0_bins[N_CR_PARTICLE_BINS], R0_bin_m[N_CR_PARTICLE_BINS], R0_bin_p[N_CR_PARTICLE_BINS];
    int proton_key; for(proton_key=0;proton_key<=1;proton_key++) // loop over whether we consider nuclei or electrons, first //
    {
       int n_active=0, bins_sorted[N_CR_PARTICLE_BINS];
       for(k=0;k<N_CR_PARTICLE_BINS;k++) {if((Z[k]>0 && proton_key==1) || (Z[k]<0 && proton_key==0)) {bins_sorted[n_active]=k; R0_bins[n_active]=R0[k]; n_active++;}} // bin is valid: charge matches that desired
       if(n_active<=0) {continue;} // nothing to do here
       if(n_active<N_CR_PARTICLE_BINS) {for(k=n_active;k<N_CR_PARTICLE_BINS;k++) {R0_bins[k]=MAX_REAL_NUMBER;}} // set a dummy value here for sorting purposes below
       //qsort(bins_sorted, n_active, sizeof(int), compare_CR_rigidity_for_sort); // sort on energies from smallest-to-largest [this is hard-coded by requiring the list go in monotonic increasing order for e and p, regardless of how the e and p are themselves ordered //
       for(k=0;k<n_active;k++)
       {
           int j = bins_sorted[k], j_m=j, j_p=j; // target bin and bin below/above
           if(k > 0) {j_m = bins_sorted[k-1];} // define previous bin 'down'
           if(k < n_active-1) {j_p = bins_sorted[k+1];} // define next bin 'up'
           R0_bin_m[j] = sqrt(R0[j]*R0[j_m]); R0_bin_p[j] = sqrt(R0[j]*R0[j_p]); // take the bin edges at the geometric means (halfway between midpoints in log-space)
           if(j_m==j) {R0_bin_m[j] = R0[j] * pow(R0[j]/R0_bin_p[j], 2);} // lowest bin gets 'padded' in the small-R direction [extends 2x as far in log-space]
           if(j_p==j) {R0_bin_p[j] = R0[j] * pow(R0[j]/R0_bin_m[j], 2);} // highest bin gets 'padded' in the large-R direction [extends 2x as far in log-space]
       }
    }
    for(k=0;k<N_CR_PARTICLE_BINS;k++) {CR_global_min_rigidity_in_bin[k] = R0_bin_m[k]; CR_global_max_rigidity_in_bin[k] = R0_bin_p[k];} // set the variables we just defined

    /* ok, now we need to build the lookup tables */
    double gamma_limit = 120.;
    int n_gamma_sample = 10000, n_table = N_CR_SPECTRUM_LUT-1;
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        double xm = CR_global_min_rigidity_in_bin[k] / CR_global_rigidity_at_bin_center[k]; // dimensionless bin minimum
        double xp = CR_global_max_rigidity_in_bin[k] / CR_global_rigidity_at_bin_center[k]; // dimensionless bin maximum
        double p_power_in_e = 1., xm_e = xm, xp_e = xp; // define variables used below
        if(CR_check_if_bin_is_nonrelativistic(k)) {p_power_in_e = 2.; xm_e = xm*xm; xp_e = xp*xp;} // extra power of p in momentum equation accounted for here, all that's needed
        double gamma_min = -gamma_limit, gamma_max = gamma_limit, d_gamma = (gamma_max-gamma_min) / ((double)n_gamma_sample); int j;

        double gamma = gamma_min, gamma_prev = gamma_min, R_index_prev=0; j=0; int alldone_key=0; // define variables for use in loop below
        while(gamma <= gamma_max)
        {
            double gamma_one = gamma + 1., xm_g = pow(xm, gamma_one), xp_g = pow(xp, gamma_one); // key variables needed
            double R = (gamma_one / (gamma_one + p_power_in_e)) * (xp_g*xp_e - xm_g*xm_e) / (xp_g - xm_g); // ratio of kinetic energy to number times bin-mean energy: this is the key function we need to invert
            double R_index = (log(R / xm_e) / log(xp_e / xm_e)) * n_table; // index one would obtain from this value of R
            if((int)R_index >= j)
            {
                double slope_gamma = gamma + (gamma-gamma_prev) * (((double)j) - R_index) / (R_index - R_index_prev); // linearly interpolate between this and previous step to 'exact' gamma giving desired R
                if(j==1 && slope_gamma < gamma_min) {CR_global_slope_lut[k][j-1] = slope_gamma + gamma_min;}
                CR_global_slope_lut[k][j] = slope_gamma; // set the look-up-table value for use later
                j++; // move to look for next target bin
            }
            if(j >= n_table) {alldone_key=1; break;} // dont need to continue loop, we've filled in all values desired here //
            R_index_prev = R_index; gamma_prev = gamma; // save for next loop
            gamma += d_gamma; // augment
            // check for bad values of gamma that we wish to avoid because they can cause divergences, and just slightly dodge around them //
            double tol = 0.1*d_gamma, badval; // tolerance around the bad values
            badval=-1; if(fabs(gamma - badval) < tol) {if(gamma<badval) {gamma=badval-tol;} else {gamma=badval+tol;}}
            badval=-2; if(fabs(gamma - badval) < tol) {if(gamma<badval) {gamma=badval-tol;} else {gamma=badval+tol;}}
            badval=-3; if(fabs(gamma - badval) < tol) {if(gamma<badval) {gamma=badval-tol;} else {gamma=badval+tol;}}
        }
        if(alldone_key==0 && j<n_table) {int j0=j-1; for(j=j0+1;j<N_CR_SPECTRUM_LUT;j++) {CR_global_slope_lut[k][j]=CR_global_slope_lut[k][j0]+(j-j0)*DMIN(gamma_limit,fabs(gamma_limit-CR_global_slope_lut[k][j0]));}}
            else {CR_global_slope_lut[k][N_CR_SPECTRUM_LUT-1]=gamma_limit;} // set upper end of table value to unity
    }
    
    if(ThisTask==0) {for(k=0;k<N_CR_PARTICLE_BINS;k++) { // print outputs for users
        printf("\n .. bin=%d, charge=%g e, mass=%g mp, rigidity Rmin=%g R0=%g Rmax=%g GV, energy=%g GeV, relativistic?=%d [1=Y/0=N] beta=%g gamma=%g \n",
           k,CR_global_charge_in_bin[k],return_CRbin_CRmass_in_mp(-1,k),CR_global_min_rigidity_in_bin[k],CR_global_rigidity_at_bin_center[k],
           CR_global_max_rigidity_in_bin[k],return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k),1-CR_check_if_bin_is_nonrelativistic(k),return_CRbin_beta_factor(-1,k),
           return_CRbin_gamma_factor(-1,k)); fflush(stdout);
        printf(" .. LUT for CR slopes in this bin: \n"); printf("  .. j  .. R_egy/num .. gamma \n");
        int j; for(j=0;j<N_CR_SPECTRUM_LUT;j++) {printf("  .. %4d  %5.4g %10.3g \n",j,((double)j)/((double)n_table),CR_global_slope_lut[k][j]); fflush(stdout);}
    }}
    
    return;
}



/* input value 'R' = ratio of total CR energy in the bin ('e_tot') to the total CR number times the energy of the CRs with the bin-centered rigidity ('n_tot' x 'E_cr_bin_center_list'), which is a dimensionless function of the slope, whether the bin is relativistic or not, and the bin edges relative to the bin center */
double CR_return_slope_from_number_and_energy_in_bin(double energy_in_code_units, double number_effective_in_code_units, double bin_centered_energy_in_GeV, int k_bin)
{
    if((energy_in_code_units <= 0.) || isnan(energy_in_code_units) || (number_effective_in_code_units <= 0.) || isnan(number_effective_in_code_units)) {return CR_global_slope_lut[k_bin][0];}
    double R = energy_in_code_units / (number_effective_in_code_units * bin_centered_energy_in_GeV + MIN_REAL_NUMBER);
    int n_table = N_CR_SPECTRUM_LUT; // table size
    double xm = CR_global_min_rigidity_in_bin[k_bin] / CR_global_rigidity_at_bin_center[k_bin], xp = CR_global_max_rigidity_in_bin[k_bin] / CR_global_rigidity_at_bin_center[k_bin], xm_e = xm, xp_e = xp;
    if(CR_check_if_bin_is_nonrelativistic(k_bin)) {xm_e=xm*xm; xp_e=xp*xp;} // sets bounds that this value can possibly obtain
    if(R >= xp_e) {return CR_global_slope_lut[k_bin][n_table-1];} // set to maximum
    if(R <= xm_e) {return CR_global_slope_lut[k_bin][0];} // set to minimum
    double n_interp = (log(R / xm_e) / log(xp_e / xm_e)) * (n_table-1); // fraction of the way between min and max in our log-spaced table, returns 0-1
    int n0 = (int) floor(n_interp), n1=n0+1;
    if(n1 > n_table-1) {n1=n_table-1;}
    double slope_gamma = CR_global_slope_lut[k_bin][n0] + (n_interp-(double)n0) * (CR_global_slope_lut[k_bin][n1]-CR_global_slope_lut[k_bin][n0]);
    // check for specific bad values that will nan out and avoid them - this introduces really minimal errors
    double tol = 0.0001, badval; // tolerance around the bad values
    badval=-1; if(fabs(slope_gamma - badval) < tol) {if(slope_gamma<badval) {slope_gamma=badval-tol;} else {slope_gamma=badval+tol;}}
    badval=-2; if(fabs(slope_gamma - badval) < tol) {if(slope_gamma<badval) {slope_gamma=badval-tol;} else {slope_gamma=badval+tol;}}
    badval=-3; if(fabs(slope_gamma - badval) < tol) {if(slope_gamma<badval) {slope_gamma=badval-tol;} else {slope_gamma=badval+tol;}}
    return slope_gamma; // ok, safe to return
}


/* return solution for the distance in dimensionless terms from the bin edge for CRs which can propagate to the edge with the dimensionless rate factor given */
double CR_return_new_bin_edge_from_rate(double rate_dt_dimless, double x_m_bin, double x_p_bin, int loss_mode, int NR_key, double additional_variable_dummy)
{
    if(loss_mode <= 0) {return 0;} // hadronic+catastrophic: should never be called [just should never get to this point since hadronic should be done]
    
    double x_e = 0;
    if(loss_mode == 1 || loss_mode == 4) // adiabatic+brems+streaming
        {if(rate_dt_dimless < 0) {x_e = x_p_bin * exp(rate_dt_dimless);} else {x_e = x_m_bin * exp(rate_dt_dimless);}} // return xp' or xm' depending on sign: solution to dx/dtau=(+/-)x

    if(loss_mode == 2) // coulomb+ionization
    {
        if(NR_key) // non-relativistic Coulomb+ionization terms
            {x_e = pow(x_m_bin*x_m_bin*x_m_bin + 3.*rate_dt_dimless , 1./3.);} // solution to dx/dtau=-1/x^2
        else // relativistic Coulomb+ionization terms
            {x_e = x_m_bin + rate_dt_dimless;} // solution to dx/dtau=-1
    }

    if(loss_mode == 3) // IC+synchrotron
        {if(rate_dt_dimless*x_m_bin < 1.) {x_e = x_m_bin / (1. - rate_dt_dimless * x_m_bin);} else {x_e = x_p_bin;}} // solution for IC+synchrotron: dx/dtau=-x^2
    
    if(loss_mode == 5) // re-acceleration
        {x_e = pow( DMAX(pow(x_p_bin, additional_variable_dummy) + rate_dt_dimless, 0.), 1./additional_variable_dummy);} // solution to dlnx/dtau = x^-delta; here using additional_variable_dummy = delta_diffcoeff[k]
    
    if(x_e < x_m_bin) {x_e=x_m_bin;}
    if(x_e > x_p_bin) {x_e=x_p_bin;}
    return x_e; // catch
}



/* integrand needed for numerical evaluation of non-relativistic Coulomb loss terms; this should be G in int_x0^x1 [G] dlnx_old = E_new, so G = x_old * (dN/dx)_old * E_new, in dimensionless bin-centered units: it's crucial here that 'slope' input is ~df/dlogx = "gamma_one" in the notation above, as opposed to "slope_gamma" */
double CR_coulomb_energy_integrand(double x, double tau, double slope)
{
    return pow(x,slope) * pow(x*x*x - 3.*tau , 2./3.);
}


/* integrand needed for numerical evaluation of re-acceleration energy gain terms; this should be G in int_x0^x1 [G] dlnx_old = E_new, so G = x_old * (dN/dx)_old * E_new, in dimensionless bin-centered units: it's crucial here that 'slope' input is ~df/dlogx = "gamma_one" in the notation above, as opposed to "slope_gamma" */
double CR_reaccel_energy_integrand(double x, double tau, double slope, double delta_slope, int NR_key)
{
    double R_slope = 1./delta_slope;
    if(NR_key) {R_slope *= 2.;}
    return pow(x,slope) * pow( pow(x, delta_slope) + tau , R_slope );
}


/* integrand needed for numerical evaluation of Compton loss terms; this should be G in int_x0^x1 [G] dlnx_old = E_new, so G = x_old * (dN/dx)_old * E_new, in dimensionless bin-centered units: it's crucial here that 'slope' input is ~df/dlogx = "gamma_one" in the notation above, as opposed to "slope_gamma" */
double CR_compton_energy_integrand(double x, double tau, double slope)
{
    return pow(x,slope+1.) / (1. + tau*x);
}


/* quick boolean to determine if we should use relativistic or non-relativistic scalings for a given bin, in our cooling subroutines where we need this to be true across the bin */
int CR_check_if_bin_is_nonrelativistic(int k_bin) // relativistic binflag: all e-, protons/nuclei with R_GV > 1.87655 * A/Z = 1.87655 for protons : need this globally
{
    if((CR_global_charge_in_bin[k_bin] > 0) && (CR_global_rigidity_at_bin_center[k_bin] < 1.87655)) {return 1;} // for now, only protons with division at p=2m0*c, so NR and R expressions equate, qualify as non-relativistic
    return 0; // default to relativistic otherwise
}


/* return number of CRs in bin, in our strange code units, where the 'number' has units of E_code / GeV. to convert to real number need to multiply by this but that can cause lots of overflow issues so we work with this unit */
double CR_return_effective_number_in_bin_in_codeunits(int target, int k_bin)
{
    return SphP[target].CosmicRay_Number_in_Bin[k_bin];
}


/* return the CR bin spectral slope, depending on what we use as our evolved variable */
double CR_return_spectral_slope_target(int target, int k_bin)
{
    //return SphP[target].CosmicRay_PwrLaw_Slopes_in_Bin[k_bin]; // evolving slopes directly
    return CR_return_slope_from_number_and_energy_in_bin(SphP[target].CosmicRayEnergy[k_bin], SphP[target].CosmicRay_Number_in_Bin[k_bin], return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_bin), k_bin); // calculate if not evolving directly
}


/* return true number of CRs in bin, in actual units */
double CR_return_true_number_in_bin(int target, int k_bin)
{
    return CR_return_effective_number_in_bin_in_codeunits(target,k_bin) * UNIT_ENERGY_IN_CGS / (1.0e9*ELECTRONVOLT_IN_ERGS); // does the unit conversion from Ecode/GeV to dimensionless number
}


/* subroutine to return the effective CR number from the bin slope and energy */
double CR_get_number_in_bin_from_slope(int target, int k_bin, double energy, double slope)
{
    double etot = energy;
    double gamma_one = 1. + slope;
    double E_bin_center = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_bin);
    double xm = CR_global_min_rigidity_in_bin[k_bin] / CR_global_rigidity_at_bin_center[k_bin];
    double xp = CR_global_max_rigidity_in_bin[k_bin] / CR_global_rigidity_at_bin_center[k_bin];
    double xm_gamma_one = pow(xm, gamma_one), xp_gamma_one = pow(xp, gamma_one);
    int NR_key = CR_check_if_bin_is_nonrelativistic(k_bin);
    double gamma_fac = 0; // factor to get total number in bin. we use the slope we obtain to convert between number and energy, given the bin-centered values
    if(NR_key) {gamma_fac = ((gamma_one + 2.) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp*xp - xm_gamma_one*xm*xm);} // non-relativistic map between E and N
        else {gamma_fac = ((gamma_one + 1.) / (gamma_one)) * (xp_gamma_one - xm_gamma_one) / (xp_gamma_one*xp - xm_gamma_one*xm);} // relativistic map between E and N
    return (etot / E_bin_center) * gamma_fac; // dimensional normalization. note the units are arbitrary for the E_bin_NtoE term, as long as we are consistent [we can evolve constant x N, instead of N, for convenience]
 }


/* return mean kinetic energy per CR for CRs in bin */
double CR_return_mean_energy_in_bin_in_GeV(int target, int k_bin)
{
    double etot = SphP[target].CosmicRayEnergyPred[k_bin]; // total energy in bin
    double ntot_GeV = CR_return_effective_number_in_bin_in_codeunits(target, k_bin); // total number in units of GeV
    return etot / ntot_GeV; // this returns the mean weighted over the CR spectrum
}


/* return mean rigidity in GV per CR for CRs in bin */
double CR_return_mean_rigidity_in_bin_in_GV(int target, int k_bin)
{
    double slope = CR_return_spectral_slope_target(target, k_bin);
    double gamma_one = 1. + slope;
    double R0 = CR_global_rigidity_at_bin_center[k_bin], xm = CR_global_min_rigidity_in_bin[k_bin] / R0, xp = CR_global_max_rigidity_in_bin[k_bin] / R0, xm_gamma_one = pow(xm, gamma_one), xp_gamma_one = pow(xp, gamma_one);
    double gamma_fac = (gamma_one/(gamma_one+1.)) * (xp_gamma_one*xp - xm_gamma_one*xm) / (xp_gamma_one - xm_gamma_one); // dimensionless term from weighted integral over the bin
    return R0 * gamma_fac; // returns value appropriately weighted over the bin
}


/* sort function if we need to sort CR bins by rigidity from lowest to highest for future compatibility */
int compare_CR_rigidity_for_sort(const void *a, const void *b)
{
    double x = CR_global_rigidity_at_bin_center[ *(int *) a];
    double y = CR_global_rigidity_at_bin_center[ *(int *) b];
    if (x < y) {return -1;} else if(x > y) {return 1;}
    return 0;
}


/* routine to return bin-centered energy using the slope approximation used in the fast cooling sub-cycling, of each bin being in the relativistic or non-relativistic regime */
double return_CRbin_kinetic_energy_in_GeV_binvalsNRR(int k_CRegy)
{
    double R_GV = CR_global_rigidity_at_bin_center[k_CRegy];
    double Zabs = fabs(CR_global_charge_in_bin[k_CRegy]);
    if(CR_check_if_bin_is_nonrelativistic(k_CRegy))
    {
        double fac = 0.5328928536929374; // converts from R_GV to E in GeV for E = p^2/(2m), assuming Z=1, m=mp
        if(Zabs>=1.01) {fac *= Zabs/2.;} // assume A = 2.*Z
        return fac * R_GV*R_GV; // E in GeV in non-relativistic limit
    }
    return R_GV / Zabs; // E in GeV in relativistic limit
}



#endif






#endif // closes block for entire file


