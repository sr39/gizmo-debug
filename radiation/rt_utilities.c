#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif
#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/*! \file rt_utilities.c
 *  \brief useful functions for radiation modules
 *
 *  This file contains a variety of useful functions having to do with radiation in different modules
 *    A number of the radiative transfer subroutines and more general mass-to-light ratio calculations
 *    will refer to these routines.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#define tiny 1e-10


#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
/* routine which returns the luminosity for the desired source particles, as a function of whatever the user desires, in the relevant bands */
int rt_get_source_luminosity(MyIDType i, double sigma_0, double *lum)
{
    int active_check = 0;
    
#ifdef GALSF_FB_RT_PHOTONMOMENTUM
    if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
    {
        if(sigma_0<0) return 1;
        active_check = 1;
        double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        double L = P[i].Mass * evaluate_l_over_m_ssp(star_age) * calculate_relative_light_to_mass_ratio_from_imf(i);
        double f_uv, f_op;
#ifndef RT_FIRE_FIX_SPECTRAL_SHAPE
        double sigma_eff = sigma_0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,P[i].DensAroundStar,PPP[i].NumNgb,0);
        if(star_age <= 0.0025) {f_op=0.09;} else {
            if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));
            } else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
        
        double tau_uv = sigma_eff*KAPPA_UV; double tau_op = sigma_eff*KAPPA_OP;
        f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)/(1+0.8*tau_uv+0.85*tau_uv*tau_uv));
        f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)/(1+0.8*tau_op+0.85*tau_op*tau_op);
        /*
         // above is a fitting function for tau_disp~0.22 'tail' w. exp(-tau) 'core', removes expensive functions
         f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)*exp(-tau_uv));
         f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)*exp(-tau_op);
         // leakage for P(tau) ~ exp(-|logtau/tau0|/sig):
         f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)/ (1 + pow(tau_uv,1./(4.*tau_disp))/(3.*tau_disp) + pow(2.*tau_disp*tau_uv,1./tau_disp)));
         f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)/ (1 + pow(tau_op,1./(4.*tau_disp))/(3.*tau_disp) + pow(2.*tau_disp*tau_op,1./tau_disp));
         */
#else
        f_uv = All.PhotonMomentum_fUV;
        f_op = All.PhotonMomentum_fOPT;
#endif
        lum[0] = L * f_uv;
        lum[1] = L * f_op;
        lum[2] = L * (1-f_uv-f_op);
    }
#endif
    
    
#if defined(RT_CHEM_PHOTOION)
    /* ionizing flux (crudely estimated) from old (Petkova & Springel) implementation */
    if((1 << P[i].Type) & (RT_PHOTOION_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k]=0;} // begin from zero //
        double fac = 0;
#if defined(BLACK_HOLES)
        if(P[i].Type==5) {if(sigma_0<0) return 1; active_check=1; fac += 0.0001 * P[i].BH_Mdot * pow(C / All.UnitVelocity_in_cm_per_s, 2);}
#endif
#if defined(GALSF)
        if(P[i].Type==0) {if(sigma_0<0) return 1; active_check=1; fac += SphP[i].Sfr * All.IonizingLumPerSFR * All.UnitTime_in_s / All.HubbleParam;} // flux from gas according to SFR
#endif
        if(P[i].Type==4) {if(sigma_0<0) return 1; active_check=1; fac += (P[i].Mass * All.UnitMass_in_g / SOLAR_MASS) * All.IonizingLumPerSolarMass * All.UnitTime_in_s / All.HubbleParam;} // flux from star particles according to mass
#if defined(RT_PHOTOION_MULTIFREQUENCY)
        // we should have pre-tabulated how much luminosity gets assigned to each different waveband according to the following function //
        for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] += fac * precalc_stellar_luminosity_fraction[k];}
#else
        lum[0] += fac;
#endif
        if(P[i].Type == 0) {if(sigma_0<0) return 1; active_check=1; rt_get_lum_gas(i,lum);}     /* re-distributes flux as a blackbody */
    }
#endif
    
    return active_check;
}




/* calculate the opacity for use in radiation transport operations [in physical code units = L^2/M] */
double rt_kappa(MyIDType i, int k_freq)
{
#ifdef GALSF_FB_RT_PHOTONMOMENTUM
    double fac = All.UnitMass_in_g * All.HubbleParam / (All.UnitLength_in_cm * All.UnitLength_in_cm);
    if(k_freq==0) return KAPPA_UV * fac;
    if(k_freq==1) return KAPPA_OP * fac;
    if(k_freq==2) return KAPPA_IR * fac;
#endif
    
#ifdef RT_CHEM_PHOTOION
    double nH_over_Density = HYDROGEN_MASSFRAC / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
    double kappa = nH_over_Density * (SphP[i].HI + tiny) * rt_sigma_HI[k_freq];
#if defined(RT_CHEM_PHOTOION_HE) && defined(RT_PHOTOION_MULTIFREQUENCY)
    kappa += nH_over_Density * ((SphP[i].HeI + tiny) * rt_sigma_HeI[k_freq] + (SphP[i].HeII + tiny) * rt_sigma_HeII[k_freq]);
#endif
    return kappa;
#endif
    
    return 0;
}


/* rate of photon absorption [absorptions per unit time per photon]: this, times the timestep dt, times the photon energy density E, 
    gives the change in the energy density from absorptions (the sink term) */
double rt_absorption_rate(MyIDType i, int k_freq)
{
    /* should be equal to (C * Kappa_opacity * rho) */
    return, (C/All.UnitVelocity_in_cm_per_s) * rt_kappa(i, k_freq) * SphP[i].Density*All.a3inv;
}

#endif








#ifdef RT_CHEM_PHOTOION
/* this tabulates the fraction of the ionizing luminosity (normalized to some total)
 from stars with a given effective temperature T_eff [in Kelvin], across the bands specified by the intervals nu[*] */
void rt_get_lum_for_spectral_bin_stars(double T_eff, double luminosity_fraction[N_RT_FREQ_BINS])
{
    int i, j, integral = 10000;
    double I_nu, sum, d_nu, e, e_start, e_end, hc = C * PLANCK, R_eff = 7.0e11;
    
    for(i = 0, sum = 0; i < N_RT_FREQ_BINS; i++)
    {
        e_start = nu[i]; if(i==N_RT_FREQ_BINS-1) {e_end = nu[N_RT_FREQ_BINS-1] + 500.;} else {e_end = nu[i+1];} // nu defines intervals to integrate
        d_nu = (e_end - e_start) / (float)(integral - 1);
        for(j=0, luminosity_fraction[i]=0.0; j<integral; j++)
        {
            e = e_start + j*d_nu; // integrate the Planck function to get the total photon number in the interval //
            I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc) / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);
            luminosity_fraction[i] += 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / e * d_nu / PLANCK; // number/s
        }
        sum += luminosity_fraction[i];
    }
    /* normalize to unity; note that if we wanted a fast version of this function, the normalization means all the constants above are
     redundant and can trivially be thrown out, to [slightly] save time. but right now only called once, so no big deal */
    for(i = 0; i < N_RT_FREQ_BINS; i++) {luminosity_fraction[i] /= sum;}
    if(ThisTask == 0)
    {
        printf("Calc Stellar spectra, for T_eff=%g, \n",T_eff);
        for(i = 0; i < N_RT_FREQ_BINS; i++) {printf(" -- spectra: bin=%d nu=%g luminosity_fraction=%g \n",i,nu[i],luminosity_fraction[i]);}
        fflush(stdout);
    }
}


void rt_get_lum_gas(int target, double *je)
{
#ifdef RT_COOLING_PHOTOHEATING // routine below only makes sense (currently) with its coupled cooling routines //
    int j;
    double temp, entropy, molecular_weight, kT, BB_l, BB_r, next, sigma_SB=5.6704e-5, u_cooling, u_BB, d_nu;
 
    /* ??? needs fixing with fixed cooling routines */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[target].Ne);
    temp =  SphP[target].Pressure * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);
    kT = temp * BOLTZMANN;
    entropy = SphP[target].Pressure * pow(SphP[target].Density * All.cf_a3inv, -1 * GAMMA);
    u_cooling = rt_get_cooling_rate(target, entropy) / (SphP[target].Density * All.cf_a3inv);
    
    u_BB = sigma_SB * pow(temp, 4) ; // erg/cm^2
    u_BB /= P[target].Mass * All.UnitEnergy_in_cgs / All.HubbleParam; //energy/cm^2/mass
    
    double prefac = -u_cooling/u_BB * M_PI / (C*C*PLANCK*PLANCK*PLANCK) * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
    for(j = 0; j < N_RT_FREQ_BINS; j++)
    {
        if(j==N_RT_FREQ_BINS-1) {next = nu[N_RT_FREQ_BINS-1] + 500.;} else {next = nu[j+1];}
        d_nu = next - nu[j];
        
        BB_r = pow(next * ELECTRONVOLT_IN_ERGS, 3) / (exp(next * ELECTRONVOLT_IN_ERGS / (kT)) - 1.0);
        BB_l = pow(nu[j] * ELECTRONVOLT_IN_ERGS, 3) / (exp(nu[j] * ELECTRONVOLT_IN_ERGS / (kT)) - 1.0);
        je[j] += prefac * (BB_l / nu[j] + BB_r / next) * d_nu * rt_nu_eff_eV[j]; // gives energy per unit time, as it should //
    }
#endif
}
#endif




#ifdef RADTRANSFER
/* this function initializes some of the variables we need */
void rt_set_simple_inits(void)
{
    int i; for(i = 0; i < N_gas; i++)
    {
        if(P[i].Type == 0)
        {
            //int j; for(j = 0; j < N_RT_FREQ_BINS; j++) {SphP[i].E_gamma[j] = tiny;}
#ifdef RT_CHEM_PHOTOION
            SphP[i].HII = tiny;
            SphP[i].HI = 1.0 - SphP[i].HII;
            SphP[i].Ne = SphP[i].HII;
#ifdef RT_CHEM_PHOTOION_HE
            double fac = (1-HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
            SphP[i].HeIII = tiny * fac;
            SphP[i].HeII = tiny * fac;
            SphP[i].HeI = (1.0 - SphP[i].HeII - SphP[i].HeIII) * fac;
            SphP[i].Ne += SphP[i].HeII + 2.0 * SphP[i].HeIII;
#endif
#endif
        }
    }
}
#endif



