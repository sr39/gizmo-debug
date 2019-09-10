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


/* return the (solar-scaled) light-to-mass ratio of an SSP with a given age; used throughout the code */
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i)
{
    double lum=1; if(stellar_age_in_gyr < 0.01) {lum=1000;} // default to a dumb imf-averaged 'young/high-mass' vs 'old/low-mass' distinction
#ifdef SINGLE_STAR_SINK_DYNAMICS // calculate single-star luminosity (and convert to solar luminosity-to-mass ratio, which this output assumes)
    lum=calculate_individual_stellar_luminosity(0, P[i].BH_Mass, i) / P[i].BH_Mass * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)) / (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS));
#endif
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION // fit to updated SB99 tracks: including rotation, new mass-loss tracks, etc.
    if(stellar_age_in_gyr < 0.0035) {lum=1136.59;} else {double log_age=log10(stellar_age_in_gyr/0.0035); lum=1500.*pow(10.,-1.8*log_age+0.3*log_age*log_age-0.025*log_age*log_age*log_age);}
#endif
    if(stellar_age_in_gyr<0.033) {lum*=calculate_relative_light_to_mass_ratio_from_imf(stellar_age_in_gyr,i);} // account for IMF variation model [if used]
    return lum;
}


/* subroutine to calculate luminosity of an individual star, according to accretion rate,
    mass, age, etc. Modify your assumptions about main-sequence evolution here. */
double calculate_individual_stellar_luminosity(double mdot, double mass, long i)
{
    double lum = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double c_code = C / All.UnitVelocity_in_cm_per_s;
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    /* if below the deuterium burning limit, just use the potential energy efficiency at the surface of a jupiter-density object */
    double rad_eff_protostar = 5.0e-7;
    if(m_solar < 0.012) {rad_eff_protostar = 5.e-8 * pow(m_solar/0.00095,2./3.);}
    lum = rad_eff_protostar * mdot * c_code*c_code;
    /* now for pre-main sequence, need to also check the mass-luminosity relation */
    double lum_sol = 0;
#ifdef SINGLE_STAR_PROMOTION
    if(m_solar >= 0.012)
    {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;}
        else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
        else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);}
        else {lum_sol = 32000. * m_solar;}
    }
#endif
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
    if(i > 0)
    {
        if(P[i].Type == 5) /* account for pre-main sequence evolution */
        {
            double T4000_4 = pow(m_solar , 0.55); // protostellar temperature along Hayashi track
            double l_kh = 0.2263 * P[i].ProtoStellarRadius_inSolar*P[i].ProtoStellarRadius_inSolar * T4000_4; // luminosity from KH contraction
            if(l_kh > lum_sol) {lum_sol = l_kh;} // if Hayashi-temp luminosity exceeds MS luminosity, use it. otherwise use main sequence luminosity, and assume the star is moving along the Henyey track
            // now, calculate accretion luminosity using protostellar radius
#ifdef SINGLE_STAR_FB_JETS
            double eps_protostar=All.BAL_f_accretion; //fraction of gas that does not get launched out with a jet
#else
            double eps_protostar=0.75; //default value, although 1.0 would be energy conserving
#endif
            lum = eps_protostar * (All.G * P[i].Mass / (P[i].ProtoStellarRadius_inSolar * 6.957e10 / All.UnitLength_in_cm)) * mdot; // assume GM/r liberated per unit mass. Note we need radius in code units here since everything else in 'lum' is code-units as well.
        }
    }
#endif
    lum_sol *= SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
    lum += lum_sol;
#endif
    return lum;
}


/* return the light-to-mass ratio, for the IMF of a given particle, relative to the Chabrier/Kroupa IMF */
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i)
{
#ifdef GALSF_SFR_IMF_VARIATION // fitting function from David Guszejnov's IMF calculations (ok for Mturnover in range 0.01-100) for how mass-to-light ratio varies with IMF shape/effective turnover mass
    double log_mimf = log10(P[i].IMF_Mturnover);
    return (0.051+0.042*(log_mimf+2)+0.031*(log_mimf+2)*(log_mimf+2)) / 0.31;
#endif
#ifdef GALSF_SFR_IMF_SAMPLING // account for IMF sampling model if not evolving individual stars
    double mu = 0.01 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam / (1.989e33); // 1 O-star per 100 Msun
    if(stellar_age_in_gyr > 0.003) {mu *= 0.326 * (0.003 / stellar_age_in_gyr);} // expectation value is declining with time, so 'effective multiplier' is larger
    return P[i].IMF_NumMassiveStars / mu;
#endif
    return 1; // Chabrier or Kroupa IMF //
}


#if defined(GALSF_FB_FIRE_RT_HIIHEATING) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
double particle_ionizing_luminosity_in_cgs(long i)
{
    double lm_ssp=0;
    if(P[i].Type != 5)
    {
        /* use updated SB99 tracks: including rotation, new mass-loss tracks, etc. */
        double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age >= 0.02) return 0; // skip since old stars don't contribute
        if(star_age < 0.0035) {lm_ssp=500.;} else {double log_age=log10(star_age/0.0035); lm_ssp=470.*pow(10.,-2.24*log_age-4.2*log_age*log_age) + 60.*pow(10.,-3.6*log_age);}
        lm_ssp *= calculate_relative_light_to_mass_ratio_from_imf(star_age, i);
#ifdef SINGLE_STAR_SINK_DYNAMICS /* use effective temperature as a function of stellar mass and size to get ionizing photon production */
        double l_sol = bh_lum_bol(0,P[i].Mass,i) * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)); // L/Lsun
        double m_sol = P[i].Mass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS)); // M/Msun
        double r_sol = pow(m_sol, 0.738); // R/Rsun
        double T_eff = 5780. * pow(l_sol/(r_sol*r_sol), 0.25); // ZAMS effective temperature
        double x0 = 157800./T_eff; // h*nu/kT for nu>13.6 eV
        double fion = 0.0; // fraction of blackbody emitted above x0
        if(x0 < 30.) {double q=18./(x0*x0) + 1./(8. + x0 + 20.*exp(-x0/10.)); fion = exp(-1./q);} // accurate to <10% for a Planck spectrum to x0>30, well into vanishing flux //
        lm_ssp = fion * l_sol / m_sol; // just needs to be multiplied by the actual stellar luminosity to get luminosity to mass ratio
#endif
    } // (P[i].Type != 5)
#ifdef BH_HII_HEATING /* AGN template: light-to-mass ratio L(>13.6ev)/Mparticle in Lsun/Msun, above is dNion/dt = 5.5e54 s^-1 (Lbol/1e45 erg/s) */
    if(P[i].Type == 5) {lm_ssp = 1.741e6 * bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i) / (P[i].Mass*All.UnitTime_in_Megayears/All.HubbleParam * (C / All.UnitVelocity_in_cm_per_s) * (C / All.UnitVelocity_in_cm_per_s));}
#endif
    lm_ssp *= (1.95*P[i].Mass*All.UnitMass_in_g/All.HubbleParam); // convert to luminosity from L/M
    if((lm_ssp <= 0) || (!isfinite(lm_ssp))) {lm_ssp=0;} // trap for negative values and nans (shouldnt happen)
    return lm_ssp;
}
#endif









/* this routine tells the feedback algorithms what to 'inject' when a stellar feedback event occurs.
    you must define the mass, velocity (which defines the momentum and energy), and metal content (yields)
    of the ejecta for the event[s] of interest. Mass [Msne] and velocity [SNe_v_ejecta] should
    be in code units. yields[k] should be defined for all metal species [k], and in dimensionless units
    (mass fraction of the ejecta in that species). */
void particle2in_addFB_fromstars(struct addFBdata_in *in, int i, int fb_loop_iteration)
{
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
    if(fb_loop_iteration == 0) {particle2in_addFB_SNe(in,i);}
    if(fb_loop_iteration == 1) {particle2in_addFB_winds(in,i);}
    if(fb_loop_iteration == 2) {particle2in_addFB_Rprocess(in,i);}
    if(fb_loop_iteration == 3) {particle2in_addFB_ageTracer(in,i);}
    return;
#endif
    if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    // 'dummy' example model assumes all SNe are identical with IMF-averaged properties from the AGORA model (Kim et al., 2016 ApJ, 833, 202)
    in->Msne = P[i].SNe_ThisTimeStep * (14.8*SOLAR_MASS)/(All.UnitMass_in_g/All.HubbleParam); // assume every SNe carries 14.8 solar masses (IMF-average)
    in->SNe_v_ejecta = 2.607e8 / All.UnitVelocity_in_cm_per_s; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun]
#ifdef METALS
    int k; for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
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
    double RSNe = 0;
#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION)
    // FIRE feedback rates: separate calculation for SNe, stellar mass loss, R-process injection //

    RSNe = mechanical_fb_calculate_eventrates_SNe(i,dt);
    mechanical_fb_calculate_eventrates_Winds(i,dt);
    mechanical_fb_calculate_eventrates_Rprocess(i,dt);
    mechanical_fb_calculate_eventrates_Agetracers(i,dt);
    return RSNe;
#endif
#ifdef GALSF_FB_THERMAL
    // pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly //
    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
    if(evaluate_stellar_age_Gyr(P[i].StellarAge) < 0.005) {return 0;} // enforce age limit of 5 Myr
    P[i].SNe_ThisTimeStep = P[i].Mass * (All.UnitMass_in_g / All.HubbleParam) / (91. * SOLAR_MASS); // 1 event per 91 solar masses
    return 1;
#endif
#ifdef GALSF_FB_MECHANICAL
    // mechanical feedback: 'dummy' example model below assumes a constant SNe rate for t < 30 Myr, then nothing. experiment! //
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age < 0.03)
    {
        RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
        double p = dt * RSNe * P[i].Mass * (All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS; // unit conversion factor
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
    }
#endif
    return RSNe;
}



#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION)
/* functions below contain pre-calculation of event rates and energetics, masses, etc, for FIRE mechanical feedback modules */

double mechanical_fb_calculate_eventrates_SNe(int i, double dt)
{
    if(All.SNeIIEnergyFrac <= 0) return 0;
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
#ifdef SINGLE_STAR_SINK_DYNAMICS
    /* here we are determining SNe for individual stars, so it happens deterministically at the end of their lives */
    double m_sol = P[i].Mass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS)); // M/Msun
    if(m_sol > 8.) // minimum mass for SNe
    {
        double l_sol = bh_lum_bol(0,P[i].Mass,i) * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)); // L/Lsun
        double lifetime = 9.6 * (m_sol/l_sol); // standard lifetime (in Gyr): this gives first SNe at 3Myr
        if(star_age >= lifetime) {P[i].SNe_ThisTimeStep = 1;}
    }
    return 0;
#else
    /* here we are determining an expected SNe rate, so SNe occur stochastically but with an age dependence in the population */
    double agemin=0.003401, agebrk=0.01037, agemax=0.03753, RSNe=0; // in Gyr //
    // calculate: NSNe/Myr *if* each SNe had exactly 10^51 ergs; really from the energy curve; below for 1Msun pop //
    if(star_age > agemin)
    {
        if(star_age<=agebrk) {RSNe=5.408e-4;} else {if(star_age<=agemax) {RSNe=2.516e-4;}} // core-collapse
        if(star_age>agemax) {RSNe=5.3e-8 + 1.6e-5*exp(-0.5*((star_age-0.05)/0.01)*((star_age-0.05)/0.01));} // Ia (prompt Gaussian+delay, Manucci+06)

        double renorm = calculate_relative_light_to_mass_ratio_from_imf(star_age,i); // account for higher # of O-stars with a different IMF
        if(star_age<agemax) {RSNe *= renorm;}
#ifdef GALSF_SFR_IMF_SAMPLING
        if(star_age>agemax && P[i].IMF_NumMassiveStars>0) {RSNe += 2.516e-4*renorm;} // account for un-exploded O-stars
#endif
        double p = dt * RSNe * P[i].Mass * (All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS; // unit conversion factor
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
#ifdef GALSF_SFR_IMF_SAMPLING
        if(star_age<agemax && P[i].IMF_NumMassiveStars<n_sn_0) {n_sn_0=P[i].IMF_NumMassiveStars;} // limit to number of O-stars for SNe //
#endif
        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
    }
    return RSNe;
#endif
}


void mechanical_fb_calculate_eventrates_Rprocess(int i, double dt)
{
#ifdef GALSF_FB_FIRE_RPROCESS
    /* we'll use the maximum rate here, then in the -yields- setting, 'cut' each down to its sub-population */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age>=0.003) // rate is zero at <3e6 yr
    {
        double unit_fac = (All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS;
        double p = 3.0e-5 / (1000.*star_age); // Nsne/Myr for a 1 Msun population
        p *= dt * P[i].Mass * unit_fac; // actual probability in the timestep
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID + 7) < p) {n_sn_0++;} // if > 1, this cuts that part off so we get appropriate n > 1 solution
        P[i].RProcessEvent_ThisTimeStep = n_sn_0; // assign event
    }
#endif
}


#ifdef GALSF_FB_FIRE_AGE_TRACERS
int get_age_tracer_bin(const double age){
/* Finds the age tracer bin (passive scalar) corresponding
   to the given star's age (in Myr) */
  int index = -1;
  const int too_old_flag = -9;
#ifndef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
  // Bins are log-spaced do not need to perform a search
  const double binstart = log10(All.AgeTracerBinStart);
  const double binend   = log10(All.AgeTracerBinEnd);
  const double log_bin_dt = (binend - binstart) / (1.0*NUM_AGE_TRACERS);

  if( age <= All.AgeTracerBinStart){
      index = 0;
  } else if (age >= All.AgeTracerBinEnd){
//      Do nothing here and continue. Likely ONLY happens in test problems,
//      otherwise stellar ages should really never be larger than largest bin
      return too_old_flag;
  } else{
      index = floor( (log10(age) - binstart)/ log_bin_dt); // find bin
  };

#else
  // Bins are custom with arbitrary sizes - perform a search
  if (age < All.AgeTracerTimeBins[0]){
    index = 0;
  } else if (age >= All.AgeTracerTimeBins[NUM_AGE_TRACERS]){
    return too_old_flag;
  } else {
    // search for the bin:
    index = binarySearch(All.AgeTracerTimeBins, age, 0, NUM_AGE_TRACERS+1,
                                                      NUM_AGE_TRACERS+1);
    if (age < All.AgeTracerTimeBins[index]){
        printf("Binary search not working %d  %f  %f  %f\n",index, age,
                                                All.AgeTracerTimeBins[index],
                                                All.AgeTracerTimeBins[index+1]);
        endrun(8888);
    }
  }

#endif
  if ( index < 0 ){
    printf("Age tracer binary search not working %d  %f\n",index, age);
    endrun(8888);
  }

  return index;
}
#endif


void mechanical_fb_calculate_eventrates_Agetracers(int i, double dt)
{
#ifdef GALSF_FB_FIRE_AGE_TRACERS
    P[i].AgeDeposition_ThisTimeStep = 0;
    if (P[i].Type != 4) return;

    // return; // this stops things
    if (All.AgeTracerRateLimitThreshold <= 0) {
      P[i].AgeDeposition_ThisTimeStep = 1.0;
      return;
    }

      // find bin - make sure dt is not comparable to width
    const double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge) * 1000.0; // Age in Myr
    dt                   *= All.UnitTime_in_Megayears; // convert to Myr
    const int k           = get_age_tracer_bin(star_age);
#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
    const double bin_dt   = All.AgeTracerTimeBins[k+1] - All.AgeTracerTimeBins[k];
#else
    const double log_bin_dt = (log10(All.AgeTracerBinEnd - log10(All.AgeTracerBinStart)))
                                        / (1.0*NUM_AGE_TRACERS);
    const double bin_dt = pow(10.0,log10(All.AgeTracerBinStart+(k+1)*log_bin_dt)) -
                          pow(10.0,log10(All.AgeTracerBinStart+(k  )*log_bin_dt));
#endif
    // if dt is large compared to bin spacing, make full event and return
    if (dt / bin_dt > All.AgeTracerRateLimitThreshold){
      P[i].AgeDeposition_ThisTimeStep = 1.0; // make event
      return;
    }

    // otherwise, rate limit according to the return fraction
    //    this represents fraction of time star should be depositing
    if(get_random_number(P[i].ID + 3) < All.AgeTracerReturnFraction){
        // deposit age tracer and use this scaling to increase
        // the normalization of the tracer field
        P[i].AgeDeposition_ThisTimeStep = 1.0 / All.AgeTracerReturnFraction;
    } else {
        // don't do anything this timestep
        P[i].AgeDeposition_ThisTimeStep = 0;
    }
#endif
    return;
}

void mechanical_fb_calculate_eventrates_Winds(int i, double dt)
{
    if(All.GasReturnFraction <= 0) return;
#ifdef NO_WIND_YIELDS
    return;
#endif
    double D_RETURN_FRAC = 0.01; // fraction of particle mass to return on a recycling step //
#ifdef SINGLE_STAR_SINK_DYNAMICS
    D_RETURN_FRAC = 1.0e-7; // needs to be much smaller to have quasi-continuous winds on these scales //
    /* use a standard scaling from e.g. Castor, Abbot, & Klein */
    double L_sol = bh_lum_bol(0, P[i].Mass, i) * All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM); // L in solar
    double M_sol = P[i].Mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS); // M in solar
    double gam = DMIN(0.5,3.2e-5*L_sol/M_sol); // Eddington factor (~L/Ledd for winds), capped at 1/2 for sanity reasons
    double alpha = 0.5 + 0.4/(1. + 16./M_sol); // approximate scaling for alpha factor with stellar type (weak effect)
    double q0 = (1.-alpha)*gam / (1.-gam); double k0=1./30.; //k is a normalization factor in the model
    double mdot = 2.338 * alpha * pow(L_sol,7./8.) * pow(M_sol,0.1845) * (1./q0) * pow(q0*k0,1./alpha); // in Msun/Gyr
    double p = mdot / M_sol; // mass fraction returned per Gyr
    p *= All.GasReturnFraction * (dt*0.001*All.UnitTime_in_Megayears/All.HubbleParam); // fraction of particle mass expected to return in the timestep //
    p = 1.0 - exp(-p); // need to account for p>1 cases //
#else
    double p=0, star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), ZZ = P[i].Metallicity[0]/All.SolarAbundances[0];
#ifdef NO_METAL_DEP_YIELDS
    ZZ = 1.0; // AJE: fix to solar
#endif
    if(ZZ>3) {ZZ=3;}
    if(ZZ<0.01) {ZZ=0.01;}
    if(star_age<=0.001){p=11.6846;} else {
        if(star_age<=0.0035){p=11.6846*ZZ*
            pow(10.,1.838*(0.79+log10(ZZ))*(log10(star_age)-(-3.00)));} else {
                if(star_age<=0.1){p=72.1215*pow(star_age/0.0035,-3.25)+0.0103;} else {
                    p=1.03*pow(star_age,-1.1)/(12.9-log(star_age));
                }}}
    if(star_age < 0.1) {p *= calculate_relative_light_to_mass_ratio_from_imf(star_age,i);} // late-time independent of massive stars
    p *= All.GasReturnFraction * (dt*0.001*All.UnitTime_in_Megayears/All.HubbleParam); // fraction of particle mass expected to return in the timestep //
    p = 1.0 - exp(-p); // need to account for p>1 cases //
    p  = DMIN(1.0, p);
    p *= 1.4 * 0.291175; // to give expected return fraction from stellar winds alone (~17%) //

    /* // updated fit from M Grudic. More accurate for early times.
     //     Needs to add the above call for later times (t >~ 0.02-0.1 Gyr) since late-time AGB loss is not strongly
     //     metallicity-dependent (as fit below only includes line-driven winds).
     double f1 = 4.68 * pow(ZZ, 0.87); // fit for fractional mass loss in first 1.5Myr
     double f3 = 0.44 * pow(ZZ, 0.77); // fit fractional mass loss from 20Myr onward
     if(star_age<=0.0015){p = f1;} else {
     if(star_age<=0.004){p = f1 * pow(star_age/0.0015,2.1);} else {
     if(star_age<=0.02){p = f1 * 7.844 * pow(star_age/0.004, 0.621335*log(0.1275*f3/f1));} else {
     p = f3 * pow(star_age/0.02, -1.1);
     }}}
     if(star_age < 0.1) {p *= calculate_relative_light_to_mass_ratio_from_imf(i);} // late-time independent of massive stars
     p *= All.GasReturnFraction * (dt*0.001*All.UnitTime_in_Megayears/All.HubbleParam); // fraction of particle mass expected to return in the timestep //
     p = 1.0 - exp(-p); // need to account for p>1 cases //
     */
#endif
    double n_wind_0=(double)floor(p/D_RETURN_FRAC); p-=n_wind_0*D_RETURN_FRAC; // if p >> return frac, should have > 1 event, so we inject the correct wind mass
    P[i].MassReturn_ThisTimeStep += n_wind_0*D_RETURN_FRAC; // add this in, then determine if there is a 'remainder' to be added as well
    if(get_random_number(P[i].ID + 5) < p/D_RETURN_FRAC) {P[i].MassReturn_ThisTimeStep += D_RETURN_FRAC;} // add the 'remainder' stochastically
}






void particle2in_addFB_Rprocess(struct addFBdata_in *in, int i)
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
        if((star_age>=tcrit)&&(p<=pcrit)) {in->yields[NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES-NUM_AGE_TRACERS+k]=1;} // units irrelevant//
    }
    in->Msne = 0.01 * (double)P[i].RProcessEvent_ThisTimeStep / ((double)((All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS)); // mass ejected ~0.01*M_sun; only here for bookkeeping //
#endif
}




#ifdef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
int binarySearch(const double * arr, const double x, 
                 const int l, const int r, const int total)
{

  if (r<l){ endrun(7777);}

  const int w     = r-l;
  const int mid   = l + w/2;

  if (w  <= 1){
    if (mid < total-1)
      if (x >= arr[mid] && x < arr[mid+1])
        return mid;
    if (mid >0)
      if(x >= arr[mid-1] && x < arr[mid])
        return mid-1;
  } // endif w < 1
  if (arr[mid] > x){
     return binarySearch(arr,x,l,mid-1,total);
  } else if (arr[mid]<x){
    return binarySearch(arr,x,mid+1,r,total);
  } else {
    return mid;
  }

  return -1;
}
#endif


void particle2in_addFB_ageTracer(struct addFBdata_in *in, int i)
{

  in->Msne = 0.0;
/*
  Routine injects passive scalar tracer field for each star particle
  according to the age of that particle for post processing stellar
  abundances. Finds the appropriate time bin(s) to inject into given
  the star's age and the timestep. Normalization of these fields is
  somewhat arbitrary, but normalizing such that each star deposits
  1 unit in each time bin (regardless of size) per solar mass. Thus,
  scaling is simple (just multiply field by mass of metal species
  produced per solar mass of star formation to get an abundance)
*/

// AJE: NOTE - do I need to worry about HubbleParam in the yield normalization?

#ifdef GALSF_FB_FIRE_AGE_TRACERS
    int k; if(P[i].AgeDeposition_ThisTimeStep<=0) {in->Msne=0; return;} // no deposition
    if (P[i].Type != 4) return; // do nothing!
#ifdef METALS
    const int k_age_start = NUM_METAL_SPECIES-NUM_AGE_TRACERS;
    const int k_age_end   = NUM_METAL_SPECIES;
    /* pull surface "abundances" - really... these are surface "ages" */
    double yields[NUM_AGE_TRACERS]; for(k=0;k<NUM_AGE_TRACERS;k++) {yields[k]=0.0;} // AJE: Normalization here might be hard - test without P[i].Metallicity[k+k_age_start];}
    // set default abundances
    for (k=0; k<NUM_AGE_TRACERS; k++){
        in->yields[k+k_age_start] = yields[k];
    }
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge) * 1000.0; // Age in Myr

    const double dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a * All.UnitTime_in_Megayears; // get particle timestep in Myr

    /* AJE: may need to switch to normalizing over logbin spacing if bins are large
      too avoid small number issues  - this requires undoing the normalization in post */
    const double M_norm = P[i].Mass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS))
                                    * P[i].AgeDeposition_ThisTimeStep;

    // get age tracer bin
    k = get_age_tracer_bin(star_age); if (k == -9){printf("age greater than bin\n"); return;} // flag that age > max bin

    // Now deposit the tracer field, with a check to make sure multiple
    // bins are crossed this dt if dt is large or age happens to be near a bin
    // edge
#ifndef GALSF_FB_FIRE_AGE_TRACERS_CUSTOM
    // bin size in linear time (better normalization)
    const double binstart = log10(All.AgeTracerBinStart);
    const double binend   = log10(All.AgeTracerBinEnd);
    const double log_bin_dt = (binend - binstart) / (1.0*NUM_AGE_TRACERS);

    double bin_dt = pow(10.0, binstart+(k+1)*log_bin_dt) - (k==0? 0.0 :pow(10.0, binstart + k*log_bin_dt));

    if (star_age + dt > pow(10.0, binstart+(k+1)*log_bin_dt)){ // goes over multiple bins!!
       double tstart = star_age;

       while ( tstart < star_age + dt){
           // do fractional / full injection for each bin we cross:
           bin_dt      =   pow(10.0, binstart + (k+1)*log_bin_dt) - pow(10.0, binstart+k*log_bin_dt);
           double tend =   DMIN(   pow(10.0,binstart+(k+1)*log_bin_dt),   star_age + dt);
           in->yields[k+k_age_start] += ((tend - tstart) / bin_dt) * M_norm;
           tstart = tend;
           k++;
       }
    } else {
        // normalization is somewhat arbitrary, but choosing "1" unit
        // per bin per solar mass of star for convenience
        in->yields[k + k_age_start] += ( dt / bin_dt ) * M_norm;
    }

#else
   double bin_dt      =   All.AgeTracerTimeBins[k+1] - All.AgeTracerTimeBins[k];
   if (star_age + dt > All.AgeTracerTimeBins[k+1] ){ // goes over multiple bins!!
      double tstart = star_age;
      while ( tstart < star_age + dt){
          // do fractional here:
          bin_dt      =   All.AgeTracerTimeBins[k+1] - All.AgeTracerTimeBins[k];
          double tend =   DMIN(   All.AgeTracerTimeBins[k+1],   star_age + dt);
          in->yields[k+k_age_start] += ((tend - tstart) / bin_dt) * M_norm;
          tstart = tend;
          k++;
      }
   } else {
       // normalization is somewhat arbitrary, but choosing "1" unit per bin per solar mass of star for convenience
       in->yields[k + k_age_start] += ( dt / bin_dt ) * M_norm;   // add to this yield bin
   }
#endif // custom bins

#ifdef GALSF_FB_FIRE_AGE_TRACERS_SURFACE_YIELDS
    /* Now (optionally) compute the return of surface yields. Since accumulation
       of metals in stars, and thus surface abundance return, is non-linear
       we cannot simply return some trivial factor times the current tracer
       abundances. Must assume a mass loss rate model: */
   // To Do: change flag to integer parameter to switch between potential models

   // use native FIRE wind model [CHANGE THIS IF WIND MODEL CHANGES]
   double p=0, ZZ = P[i].Metallicity[0]/All.SolarAbundances[0];
   star_age *= 0.001; // Convert to Gyr for the below
#ifdef NO_METAL_DEP_YIELDS
   ZZ = 1.0; // AJE: fix to solar
#endif
   if(ZZ>3) {ZZ=3;}
   if(ZZ<0.01) {ZZ=0.01;}
   if(star_age<=0.001){p=11.6846;} else {
       if(star_age<=0.0035){p=11.6846*ZZ*
           pow(10.,1.838*(0.79+log10(ZZ))*(log10(star_age)-(-3.00)));} else {
               if(star_age<=0.1){p=72.1215*pow(star_age/0.0035,-3.25)+0.0103;} else {
                   p=1.03*pow(star_age,-1.1)/(12.9-log(star_age));
               }}}
   if(star_age < 0.1) {p *= calculate_relative_light_to_mass_ratio_from_imf(star_age,i);} // late-time independent of massive stars
   p *= All.GasReturnFraction * (dt*0.001);
   //p = 1.0 - exp(-p); // need to account for p>1 cases //
   //p  = DMIN(1.0, p);
   p *= 1.4 * 0.291175 * P[i].AgeDeposition_ThisTimeStep;
   for(k=0; k<NUM_AGE_TRACERS;k++){
     in->yields[k + k_age_start] += (p*P[i].Metallicity[k+k_age_start]);
   }
#endif

#endif // metals
    in->Msne = 1.0E-10 * SOLAR_MASS / All.UnitMass_in_g / All.HubbleParam; // small number just to be nonzero

#endif // age tracer model
    return;
}



void particle2in_addFB_SNe(struct addFBdata_in *in, int i)
{
    int k; if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    int SNeIaFlag=0; if(evaluate_stellar_age_Gyr(P[i].StellarAge) > 0.03753) {SNeIaFlag=1;}; /* assume SNe before critical time are core-collapse, later are Ia */
    double Msne=10.5; if(SNeIaFlag) {Msne=1.4;} // average ejecta mass for single event (normalized to give total mass loss correctly)
    double SNeEgy = All.SNeIIEnergyFrac*P[i].SNe_ThisTimeStep * 1.0e51/(All.UnitEnergy_in_cgs/All.HubbleParam); // assume each SNe has 1e51 erg
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
#ifndef NO_METAL_DEP_YIELDS
            // metal-dependent yields:
            if(P[i].Metallicity[0]<0.033) {yields[3]*=P[i].Metallicity[0]/All.SolarAbundances[0];} else {yields[3]*=1.65;} // N scaling is strongly dependent on initial metallicity of the star //
            yields[0] += yields[3]-0.0479; // correct total metal mass for this correction //
#endif
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

#ifndef AGE_TRACER_TEST_MODE
    /* add a check to allow for larger abundances in the progenitor stars (usually irrelevant) */
    for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {yields[k]=yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]);}
#else
    // direct add in surface abundances without scaling to do better test
#ifndef NO_SURFACE_ABUNDANCES
    for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {yields[k]=yields[k] + P[i].Metallicity[k];}
#endif    
#endif


    if(SNeIaFlag) {if(NUM_METAL_SPECIES>=10) {yields[1]=0.0;}} // no He yield for Ia SNe //
    for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {if(yields[k]<0) {yields[k]=0.0;} if(yields[k]>1) {yields[k]=1;} in->yields[k]=yields[k];}
#endif
    in->Msne = P[i].SNe_ThisTimeStep * (Msne*SOLAR_MASS)/(All.UnitMass_in_g/All.HubbleParam); // total mass in code units
#ifdef SINGLE_STAR_SINK_DYNAMICS
    in->Msne = P[i].Mass; // conserve mass and destroy star completely
#endif
    in->SNe_v_ejecta = sqrt(2.0*SNeEgy/in->Msne); // v_ej in code units
}



void particle2in_addFB_winds(struct addFBdata_in *in, int i)
{
    int k; if(P[i].MassReturn_ThisTimeStep<=0) {in->Msne=0; return;} // no event
#ifdef METALS
    /* assume track initial metallicity; turn on COOL_METAL_LINES_BY_SPECIES for more detailed tracking of light elements */
    double yields[NUM_METAL_SPECIES-NUM_AGE_TRACERS]={0.0};
#ifndef NO_SURFACE_ABUNDANCES
    for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {yields[k]=P[i].Metallicity[k];} // return surface abundances, to leading order //
#endif
    if(NUM_METAL_SPECIES>=10)
    {
        /* All, then He,C,N,O,Ne,Mg,Si,S,Ca,Fe ;; follow AGB/O star yields in more detail for the light elements */
        /*   the interesting species are He & CNO: below is based on a compilation of van den Hoek & Groenewegen 1997, Marigo 2001, Izzard 2004 */
#ifndef NO_WIND_YIELDS
        yields[1]=0.36; /*He*/ yields[2]=0.016; /*C*/ yields[3]=0.0041; /*N*/ yields[4]=0.0118; /*O*/
#ifndef NO_METAL_DEP_YIELDS
        // metal-dependent yields: O scaling is strongly dependent on initial metallicity of the star //
        if(P[i].Metallicity[0]<0.033) {yields[4] *= P[i].Metallicity[0]/All.SolarAbundances[0];} else {yields[4] *= 1.65;}
#endif // metal dep yields

#ifndef AGE_TRACER_TEST_MODE // turn this scaling off when testing tracers!!!
        for(k=1;k<=4;k++) {yields[k]=yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]); if(yields[k]<0) {yields[k]=0.0;} if(yields[k]>1) {yields[k]=1;} in->yields[k]=yields[k];} // enforce yields obeying pre-existing surface abundances, and upper/lower limits //
//  - no #else statement needed since default yield abundances are surface (see above)
#endif // test mode

        yields[0]=0.0; for(k=2;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {yields[0]+=yields[k];}
#endif // wind yields

    } else {
        yields[0]=0.032; for(k=1;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {yields[k]=0.0;}
    }
    for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) in->yields[k]=yields[k];
#endif
    in->Msne = P[i].Mass * P[i].MassReturn_ThisTimeStep; // mass (in code units) returned
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double m_msun = P[i].Mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    in->SNe_v_ejecta = (616.e5 * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / All.UnitVelocity_in_cm_per_s; // scaling from size-mass relation+eddington factor, assuming line-driven winds //
#else
    /* calculate wind kinetic luminosity + internal energy (hot winds from O-stars, slow from AGB winds) */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), E_wind_tscaling=0.0013;
    if(star_age <= 0.1) {E_wind_tscaling=0.0013 + 16.0/(1+pow(star_age/0.0025,1.4)+pow(star_age/0.01,5.0));} // stellar population age dependence of specific wind energy, in units of an effective internal energy/temperature
    in->SNe_v_ejecta = sqrt(2.0 * (All.AGBGasEnergy * E_wind_tscaling * (3.0e7*(1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS) * All.UnitMass_in_g/All.UnitEnergy_in_cgs))); // get the actual wind velocity (multiply specific energy by units, user-set normalization, and convert)
#endif
}

#endif // GALSF_FB_MECHANICAL+GALSF_FB_FIRE_STELLAREVOLUTION


#endif /* GALSF */
