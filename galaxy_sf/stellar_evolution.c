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
    double c_code = C_LIGHT_CODE;
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    /* if below the deuterium burning limit, just use the potential energy efficiency at the surface of a jupiter-density object */
    double rad_eff_protostar = 5.0e-7;
    if(m_solar < 0.012) {rad_eff_protostar = 5.e-8 * pow(m_solar/0.00095,2./3.);}
    lum = rad_eff_protostar * mdot * c_code*c_code;
    /* now for pre-main sequence, need to also check the mass-luminosity relation */
    double lum_sol = 0;
#if (defined(SINGLE_STAR_PROMOTION) && (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 0))  
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
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) && (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 0)
            double T4000_4 = pow(m_solar , 0.55); // protostellar temperature along Hayashi track
            double l_kh = 0.2263 * P[i].ProtoStellarRadius_inSolar*P[i].ProtoStellarRadius_inSolar * T4000_4; // luminosity from KH contraction
            if(l_kh > lum_sol) {lum_sol = l_kh;} // if Hayashi-temp luminosity exceeds MS luminosity, use it. otherwise use main sequence luminosity, and assume the star is moving along the Henyey track
            // now, calculate accretion luminosity using protostellar radius
#ifdef SINGLE_STAR_FB_JETS
            double eps_protostar=1.0; // since mdot is already modified by All.BAL_f_accretion 
#else
            double eps_protostar=0.75; //fraction of gas that does not get launched out with a jet, default value, although 1.0 would be energy conserving
#endif
            lum = eps_protostar * (All.G * P[i].Mass / (P[i].ProtoStellarRadius_inSolar * 6.957e10 / All.UnitLength_in_cm)) * mdot; // assume GM/r liberated per unit mass. Note we need radius in code units here since everything else in 'lum' is code-units as well.

#elif (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1)
            lum = P[i].StarLuminosity_Solar * SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s); //get pre-calculated luminosity of the star
#endif
        }
    }
#endif
#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION) && (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 0)
    lum_sol *= SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
    lum += lum_sol;
    P[i].StarLuminosity_Solar = lum / ( SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s) ); //store total luminosity of the star in solar units
#endif
    
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
    if(P[i].Type == 5) {lm_ssp = 1.741e6 * bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i) / (P[i].Mass*All.UnitTime_in_Megayears/All.HubbleParam*C_LIGHT_CODE*C_LIGHT_CODE);}
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
    in->Msne = P[i].SNe_ThisTimeStep * (14.8*SOLAR_MASS)/(All.UnitMass_in_g/All.HubbleParam); // assume every SNe carries 14.8 solar masses (IMF-average)
    in->SNe_v_ejecta = 2.607e8 / All.UnitVelocity_in_cm_per_s; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun]
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
    double RSNe = 0;
#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION)
    // FIRE feedback rates: separate calculation for SNe, stellar mass loss, R-process injection //
    RSNe = mechanical_fb_calculate_eventrates_SNe(i,dt);
    mechanical_fb_calculate_eventrates_Winds(i,dt);
    mechanical_fb_calculate_eventrates_Rprocess(i,dt);
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

void mechanical_fb_calculate_eventrates_Winds(int i, double dt)
{
    if(All.GasReturnFraction <= 0) return;
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
    in->Msne = 0.01 * (double)P[i].RProcessEvent_ThisTimeStep / ((double)((All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS)); // mass ejected ~0.01*M_sun; only here for bookkeeping //
#endif
}



void particle2in_addFB_SNe(struct addFB_evaluate_data_in_ *in, int i)
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
    in->Msne = P[i].SNe_ThisTimeStep * (Msne*SOLAR_MASS)/(All.UnitMass_in_g/All.HubbleParam); // total mass in code units
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
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double m_msun = P[i].Mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    in->SNe_v_ejecta = (616.e5 * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / All.UnitVelocity_in_cm_per_s; // scaling from size-mass relation+eddington factor, assuming line-driven winds //
#else
    /* calculate wind kinetic luminosity + internal energy (hot winds from O-stars, slow from AGB winds) */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), E_wind_tscaling=0.0013;
    if(star_age <= 0.1) {E_wind_tscaling=0.0013 + 16.0/(1+pow(star_age/0.0025,1.4)+pow(star_age/0.01,5.0));} // stellar population age dependence of specific wind energy, in units of an effective internal energy/temperature
    in->SNe_v_ejecta = sqrt(2.0 * (All.AGBGasEnergy * E_wind_tscaling * (3.0e7/((5./3.-1.)*U_TO_TEMP_UNITS)))); // get the actual wind velocity (multiply specific energy by units, user-set normalization, and convert)
#endif
}

#endif // GALSF_FB_MECHANICAL+GALSF_FB_FIRE_STELLAREVOLUTION

    





#if defined(SINGLE_STAR_PROTOSTELLAR_EVOLUTION)

/* 'master' function to update the size (and other properties like effective temperature) of accreting protostars along relevant stellar evolution tracks */
double singlestar_subgrid_protostellar_evolution_update_track(int n, double dm, double dt)
{
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION==0)
    /* this is the simple version written by Phil: intentionally simplified PS evolution tracks, designed to make it easy to understand and model the evolution and reduce un-necessary complications */
    double lum_sol = 0.0, m_initial = DMAX(1.e-37 , (BPP(n).BH_Mass - dm)), mu = DMAX(0, dm/m_initial), m_solar = BPP(n).BH_Mass * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS)), T4000_4 = pow(m_solar, 0.55); // m_initial = mass before the accretion, mu = relative mass accreted, m_solar = mass in solar units, T4000_4 = (temperature/4000K)^4 along Hayashi track
    if(m_solar > 0.012) // below this limit, negligible luminosity //
    {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;} else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
          else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);} else {lum_sol = 32000. * m_solar;}
    }
    double R_Hayashi_Henyey = 2.1 * sqrt(lum_sol / T4000_4); // size below which, at the temperature above, contraction must occur along the Henyey track at constant luminosity
    double t_R_evol = 0, contraction_factor = 0; // timescale for contraction
    if(BPP(n).ProtoStellarRadius_inSolar <= R_Hayashi_Henyey)
    {// currently on Henyey track, contracting at constant Luminosity
        t_R_evol = 1.815e7 * m_solar*m_solar / (BPP(n).ProtoStellarRadius_inSolar * lum_sol) / (All.UnitTime_in_s/(All.HubbleParam * SEC_PER_YEAR)); // contraction timescale
        contraction_factor = 1. / (1. + dt/t_R_evol);
    } else {// currently on Hayashi track, contracting at constant Temperature
        t_R_evol = 8.021e7 * m_solar*m_solar / (BPP(n).ProtoStellarRadius_inSolar*BPP(n).ProtoStellarRadius_inSolar*BPP(n).ProtoStellarRadius_inSolar * T4000_4) / (All.UnitTime_in_s/(All.HubbleParam * SEC_PER_YEAR)); // contraction timescale
        contraction_factor = 1. / pow(1 + 3.*dt/t_R_evol, 1./3.);
    }
    double r_new = 100. * m_solar, R_main_sequence_ignition; // r_new = size of newly-formed protostar; R_main_sequence_ignition = main sequence radius - where contraction should halt
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

#elif (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1)
    /* Protostellar evolution model based on the ORION version, see Offner 2009 Appendix B */
    
    const double frad = 0.33; //limit for forming radiative barrier
    const double fk = 0.5; //fraction of kinetic energy that is radiated away in the inner disk before reaching the surface, using default ORION value here as it is not a GIZMO input parameter
    const double max_rel_dr = 0.01; //Maximum relative change in radius per step, if the change over a single timestep is larger than this than we subcycle the evolution of the stellar radius
    double mass = BPP(n).BH_Mass; //mass of star/protostar at the end of the timestep
    double mass_D = BPP(n).Mass_D; //amount of D in the protostar
    double mdot = BPP(n).BH_Mdot; //accretion rate, shorter to write it this way
    double mdot_m_solar_per_year = mdot * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/All.UnitTime_in_s*SEC_PER_YEAR; // accretion rate in msolar/yr
    double m_solar = mass * (All.UnitMass_in_g / SOLAR_MASS); // mass in units of Msun
    double m_initial = DMAX(1.e-37 , (mass - dm)); // mass before accretion
    int stage = BPP(n).ProtoStellarStage; /*what stage of stellar evolution the particle is in 0: pre collapse, 1: no burning, 2: fixed Tc burnig, 3: variable Tc burning, 4: shell burning, 5: main sequence, see Offner 2009 Appendix B*/
    if (stage == 0){
        //set the radius for the pre-collapse phase according to Eq B1 in Offner 2009, this overwrites the original prescription from sfr_eff.c
        BPP(n).ProtoStellarRadius_inSolar = DMAX(2.5 * pow(mdot_m_solar_per_year*1e5,0.2),2.0); //radius always at least 2 R_sun
    }
    double r = BPP(n).ProtoStellarRadius_inSolar * SOLAR_RADIUS/All.UnitLength_in_cm; // star radius in code units
    int stage_increase = 0;
    double lum_Hayashi = ps_lum_Hayashi_BB(mass, r); //blackbody radiation assuming the star follows the Hayashi track
    double lum_MS = ps_lum_MS(mass); //luminosity of main sequence star of m mass
    double lum_int = DMAX(lum_Hayashi, lum_MS); //luminosity from the stellar interior
    double lum_I = ps_lum_I(mdot); //luminosity needed to ionize the accreted material
    if (stage < 5){ //not a main sequence star
        if (stage >= 1){ //We only evolve those that are beyond the pre-collapse phase
            int loop_subcycle=0;
            int n_subcycle=1; //no subcycle by default
            double dm_curr = dm; double dt_curr = dt;
            mass = m_initial; //value at the beginning o the timestep
            double n_ad, ag, rhoc, Pc, Tc, beta, dlogbeta_dlogm, lum_D, dm_D, rel_dr; //need to declare them here or the compiler reduces their scope to the loop???
            do{//we use a do-while loop here because we first try the full timestep and if dr is too large we subcycle
            //Note: this subcycling does not update stage, so a protostar could overshoot
                //Evolve mass
                mass += dm_curr; //increase mass first
                double dm_rel = dm_curr/(mass-dm_curr);
                //Get properties for stellar evolution
                lum_Hayashi = ps_lum_Hayashi_BB(mass, r); //blackbody radiation assuming the star follows the Hayashi track
                lum_MS = ps_lum_MS(mass); //luminosity of main sequence star of m mass
                lum_int = DMAX(lum_Hayashi, lum_MS); //luminosity from the stellar interior
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
                    dm_D = dm_curr - dt_curr * lum_D / (15.*SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)) * (1e-5) / ((All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/All.UnitTime_in_s*SEC_PER_YEAR) ;
                }
                else{ if (stage>2){
                    //burning all accreted D for stages above 2
                    lum_D = 15.*SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s) * (mdot_m_solar_per_year/(1e-5));
                    dm_D = 0; //all new D is burned
                    mass_D = 0; //no D left in protostar
                    }
                }                
                //Let's evolve the stellar radius
                rel_dr = 2. * ( dm_rel * (1.-(1.-fk)/(ag*beta)+0.5*dlogbeta_dlogm) - dt_curr/(ag*beta)*r/(All.G*mass*mass) * (lum_int+lum_I-lum_D) ); //Eq B4 of Offner 2009 divided by r
                //Let's check if we need to subcycle
                if (rel_dr > max_rel_dr){
                    n_subcycle = (int) DMAX(ceil(rel_dr/max_rel_dr), 2.0*n_subcycle); //number of subcycle steps, at least 2, either double the previous number or estimated from dr
                    //reset protostar properties, restart loop
                    loop_subcycle = 0;
                    mass = m_initial; mass_D = BPP(n).Mass_D; 
                    dm_curr = dm/((double)n_subcycle); dt_curr = dt/((double)n_subcycle);
                    r = BPP(n).ProtoStellarRadius_inSolar * SOLAR_RADIUS/All.UnitLength_in_cm;
                }
                else{
                    loop_subcycle++;
                    r *= rel_dr;
                    mass_D += dm_D;
                }
            } while(loop_subcycle<n_subcycle); //repeat for the number of subcycle steps
            //Update stellar properties
            BPP(n).ProtoStellarRadius_inSolar *= (1.0+rel_dr);
            dm_D = mass_D - BPP(n).Mass_D; //get the tota change in D mass in the protostar
            BPP(n).Mass_D = mass_D;
            //Debug message
            printf("PS evolution t: %g sink ID: %u mass: %g radius_solar: %g stage: %d mdot_m_solar_per_year: %g mD: %g rel_dr: %g dm: %g dm_D: %g Tc: %g beta: %g dt: %g n_ad: %g lum_int: %g lum_I: %g lum_D: %g age_Myr: %g StarLuminosity_Solar: %g BH_Mass_AlphaDisk: %g SinkRadius: %g dlogbeta_dlogm: %g n_subcycle: %d PS_end\n",All.Time, P[n].ID,m_solar,BPP(n).ProtoStellarRadius_inSolar,stage, mdot_m_solar_per_year, BPP(n).Mass_D*(All.UnitMass_in_g / SOLAR_MASS),rel_dr,dm* (All.UnitMass_in_g / SOLAR_MASS), dm_D* (All.UnitMass_in_g / SOLAR_MASS), Tc, beta, dt*All.UnitTime_in_Megayears, n_ad, lum_int / (SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)), lum_I/ (SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)), lum_D/ (SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)), (All.Time-P[n].ProtoStellarAge)*All.UnitTime_in_Megayears, BPP(n).StarLuminosity_Solar, BPP(n).BH_Mass_AlphaDisk, BPP(n).SinkRadius, dlogbeta_dlogm, n_subcycle );
            //Check whether the star can progress to the next state
            //Move from "no burn" to "burning at fixed Tc" phase when central temperature gets high enough for D ignition
            if ( (stage==1) && (Tc >= 1.5e6) && ((All.Time-BPP(n).StellarAge) > DMAX(3.*dt, 1e-4/All.UnitTime_in_Megayears) ) ){ //further check that the sink has been promoted at least a couple of timesteps and 100 yr ago, so that we don't start D burning immediately after forming the sink (relevant in low res cases)
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
            //Debug message
            printf("%u promoted to %d \n",P[n].ID,(stage+stage_increase));
        } //increase evolutionary stage if the particle satisfies the requirements
    }
    else{ // for main sequence stars
        BPP(n).ProtoStellarRadius_inSolar = ps_radius_MS_in_solar(mass); //update the mass if the mass changes (unlikely)
    }
    //Calculate the luminosity of the star
    /*********************************************/
    /* Power based parametrization from ORION (2 params)*/
    //lum_acc = facc * fk * All.G*mass*mdot/r; //luminosity radiated at the accretion shock
    //lum_disk = (1.-fk) * All.G*mass*mdot/r; //luminosity released by material that traverses the inner disk
    //BPP(n).StarLuminosity_Solar = (lum_acc + lum_disk + lum_int) / (SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)) ; //luminosity of the star
    /*********************************************/
    /* Mass flux based parametrization (1 params) */
    /* For our nominal choice of BAL_f_accretion=0.7 this gives very similar results to the ORION parameters of fk=facc=0.5, which is equivalent to 0.75*/
    double eps_protostar=0.75; //fraction of gas that does not get launched out with a jet, default value, although 1.0 would be energy conserving
#ifdef SINGLE_STAR_FB_JETS
    eps_protostar=1.0; // since mdot is already modified by All.BAL_f_accretion
#endif
    BPP(n).StarLuminosity_Solar = (eps_protostar*All.G*mass*mdot/r + lum_int)/ (SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)); //luminosity of the star in solar units
#endif
}


#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1) /* Functions for protosteller evolution model based on Offner 2009 */
/* Calculate the mean ratio of the gas pressure to the gas+radiation pressure, either by solving the Eddington quartic (for n_ad=3, Eq B5) or by using tabulated values, based on Offner 2009, code taken from ORION */
double ps_beta(double m, double n_ad, double rhoc, double Pc) {
    double mass=m*All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS);//in units of solar mass
    if(n_ad==3.0) {// In this case we solve the Eddington quartic, P_c^3 = (3/a) (k / (mu mH))^4 (1 - beta) / beta^4 rho_c^4 for beta
        int JMAX=40, j; double BETAMIN=1.0e-4, BETAMAX=1.0, TOL=1.0e-7, dx, f, fmid, xmid, rtb, x1=BETAMIN, x2=BETAMAX, coef=3.0/7.56e-15*pow(BOLTZMANN*rhoc/(0.613*PROTONMASS),4);
        f = pow(Pc,3) - coef * (1.0-x1)/pow(x1,4); fmid = pow(Pc,3) - coef * (1.0-x2)/pow(x2,4); rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
        for (j=1;j<=JMAX;j++) {
          xmid=rtb+(dx *= 0.5); fmid = pow(Pc,3) - coef * (1.0-xmid)/pow(xmid,4);
          if(fmid <= 0.0) rtb=xmid;
          if(fabs(dx) < TOL*fabs(xmid) || fmid == 0.0) return rtb;
        }
        printf("ps_beta: bisection solve failed to converge"); return(-1);
    } else {
        // For n != 3, we use a table lookup. The values of beta have been pre-computed with mathematica. The table goes from M=5 to 50 solar masses in steps of 2.5 M_sun, and from n=1.5 to n=3 in steps of 0.5. We should never call this routine with M > 50 Msun, since by then the star should be fully on the main sequence.
        double MTABMIN=5.0, MTABMAX=50.0, MTABSTEP=2.5, NTABMIN=1.5, NTABMAX=3.0, NTABSTEP=0.5;
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
    double mdot_m_solar_per_year = mdot * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/All.UnitTime_in_s*SEC_PER_YEAR; // accretion rate in msolar/yr
    return ( 5.0 - 3/(1.475+0.07*log10(mdot_m_solar_per_year)) );
}
/* Sets the adiabatic index for protostars based on Appendix B of Offner 2009 */
double ps_adiabatic_index(int stage, double mdot){
    double n_ad;
    switch(stage) {
        case 0: n_ad = ps_adiabatic_index_func(mdot); break;
        case 1: n_ad = ps_adiabatic_index_func(mdot); break;
        case 2: n_ad = 1.5; break;
        case 3: n_ad = 1.5; break;
        case 4: n_ad = 3.0; break;
        case 5: n_ad = 3.0; break;
    }
    if(n_ad > 3.0) {n_ad=3.0;}
    if(n_ad < 1.5) {n_ad=1.5;}
    return n_ad;
}
/* Calculate central temperature for protostar by solving Pc = rho_c*kb*Tc/(mu*mH)+1/3*a*Tc^4 using bisection, based on Offner 2009 Eq B14, code taken from ORION */
double ps_Tc(double rhoc, double Pc) {
    int JMAX=40; double TOL=1.0e-7; // max number of iterations, and error tolerance, respectively
    double Pc_cgs=Pc*All.UnitPressure_in_cgs, rhoc_cgs=rhoc * All.UnitDensity_in_cgs, Tgas=Pc_cgs*0.613*PROTONMASS/(BOLTZMANN*rhoc_cgs), Trad=pow(3*Pc_cgs/7.56e-15, 0.25), dx, f, fmid, xmid, rtb, x1=0, x2; int j; char errstr[256];
    x2 = (Trad > Tgas) ? 2*Trad : 2*Tgas;
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
    return m / (4.0/3.0*M_PI*r*r*r) / rhofac;
}
/* Calculate central pressure for protostar using a pre-computed table for fixed mass, radius and polytropic index, based on Offner 2009, table and code taken from ORION */
double ps_Pc(double m, double n_ad, double r) {
    static double pfactab[] = {0.770087, 0.889001, 1.02979, 1.19731, 1.39753, 1.63818, 1.92909, 2.2825, 2.71504, 3.24792, 3.90921, 4.73657, 5.78067, 7.11088, 8.82286, 11.0515, 13.9885};
    int itab = (int) floor((n_ad-1.5)/0.1); double wgt = (n_ad - (1.5 + 0.1*itab)) / 0.1, pfac = pfactab[itab]*(1.0-wgt) + pfactab[itab+1]*wgt;
    return pfac * All.G * m*m/(r*r*r*r);
}
/* Calculate the mean ratio of the gas pressure to the gas+radiation pressure at the center, based on Offner 2009, code taken from ORION */
double ps_betac(double rhoc, double Pc, double Tc) {
    return rhoc*BOLTZMANN*Tc/(0.613*PROTONMASS) / Pc;
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
    double mdot_m_solar_per_year = mdot * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/All.UnitTime_in_s*SEC_PER_YEAR; // accretion rate in msolar/yr
    return 2.5*SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s) * ( mdot_m_solar_per_year/(1e-5));
}
/* Calculate the blackbody luminosity [in code units] of the star following the Hayashi track */
double ps_lum_Hayashi_BB(double m, double r) {
    double m_solar = m * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS), T4000_4 = pow(m_solar , 0.55); // protostellar temperature along Hayashi track
    return 0.2263 * r * r * T4000_4 * SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s); // luminosity from KH contraction
}
/* Calculate the luminosity of a main sequence star in code units: ORION version: fitting formulas of Tout et al (1996) */
double ps_lum_MS(double m) {
    double m_solar = m * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS), unit_L = SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
    if(m_solar <= 0.1) {return 0;} // zero luminosity below some threshold for fits, otherwise use formulae below //
    return ((0.39704170*pow(m_solar,5.5) + 8.52762600*pow(m_solar,11)) / (0.00025546+pow(m_solar,3)+5.43288900*pow(m_solar,5) + 5.56357900*pow(m_solar,7) + 0.78866060*pow(m_solar,8)+0.00586685*pow(m_solar,9.5))) * unit_L;
}
/* Calculate the radius of a main sequence star in solar units: ORION version: fitting formulas of Tout et al (1996) */
double ps_radius_MS_in_solar(double m) {
    double m_solar = m * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    return (1.71535900*pow(m_solar,2.5)+6.59778800*pow(m_solar,6.5)+10.08855000*pow(m_solar,11)+1.01249500*pow(m_solar,19)+0.07490166*pow(m_solar,19.5)) /
        (0.01077422+3.08223400*pow(m_solar,2)+17.84778000*pow(m_solar,8.5)+pow(m_solar,18.5)+0.00022582*pow(m_solar,19.5)); //*SOLAR_RADIUS/All.UnitLength_in_cm);
}
#endif // (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1) /* end functions for protosteller evolution model based on Offner 2009 */


#endif //end of protostellar evolution functions







#endif /* GALSF */
