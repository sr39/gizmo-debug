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

/*Functions for protosteller evolution model based on Offner 2009*/
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 1)
/*Physical constants used in formulas, undefined later*/
#define A_SB      7.56e-15 //4 sigma_SB/c
#define MU     0.613  /* Mean molecular weight of a fully ionized gas of solar composition*/

/* Sets the adiabatic index for pre burning protostars based on Eq B2 from Offner 2009 */
double ps_adiabatic_index_func(double mdot){
    double mdot_m_solar_per_year = mdot * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/UnitTime_in_s*SEC_PER_YEAR; // accretion rate in msolar/yr
    return ( 5.0 - 3/(1.475+0.07*log10(mdot_m_solar_per_year)) )
}

/*Sets the adiabatic index for protostars based on Appendix B of Offner 2009*/
double ps_adiabatic_index(int stage, double mdot){
    double n_ad;
    switch(stage){
        case 0: n_ad = ps_adiabatic_index_func(mdot); break;
        case 1: n_ad = ps_adiabatic_index_func(mdot); break;
        case 2: n_ad = 1.5; break;
        case 3: n_ad = 1.5; break;
        case 4: n_ad = 3.0; break;
        case 5: n_ad = 3.0; break;
    }
    if (n_ad > 3.0){n_ad=3.0;}
    if (n_ad < 1.5){n_ad=1.5;}
    return n_ad;
}

/*Calculate central density for protostar using a pre-computed table for fixed mass, radius and polytropic index, based on Offner 2009, table and code taken from ORION*/
double ps_rhoc(double m, double n_ad, double r){
    /*Tabulated values of rho_c/rho_mean for n=1.5 to 3.1 in intervals of 0.1*/
    static double rhofactab[] = {
    0.166931, 0.14742, 0.129933, 0.114265, 0.100242,
    0.0877, 0.0764968, 0.0665109, 0.0576198, 0.0497216,
    0.0427224, 0.0365357, 0.0310837, 0.0262952, 0.0221057,
    0.0184553, 0.01529
    };
    int itab = (int) floor((n_ad-1.5)/0.1);
    double wgt = (n_ad - (1.5 + 0.1*itab)) / 0.1;
    double rhofac = rhofactab[itab]*(1.0-wgt) + rhofactab[itab+1]*wgt;
    return( mass / (4.0/3.0*M_PI*r*r*r) / rhofac );
}

/*Calculate central pressure for protostar using a pre-computed table for fixed mass, radius and polytropic index, based on Offner 2009, table and code taken from ORION*/
double ps_Pc(double m, double n_ad, double r){
    static double pfactab[] = {
    0.770087, 0.889001, 1.02979, 1.19731, 1.39753,
    1.63818, 1.92909, 2.2825, 2.71504, 3.24792, 3.90921,
    4.73657, 5.78067, 7.11088, 8.82286, 11.0515, 13.9885
    };
    int itab = (int) floor((n_ad-1.5)/0.1);
    double wgt = (n_ad - (1.5 + 0.1*itab)) / 0.1;
    double pfac = pfactab[itab]*(1.0-wgt) + pfactab[itab+1]*wgt;
    return( pfac * G * m*m/(r*r*r*r) );
}

/*Calculate central temperature for protostar by solving Pc = rho_c*kb*Tc/(mu*mH)+1/3*a*Tc^4 using bisection, based on Offner 2009 Eq B14, code taken from ORION*/
double ps_Tc(double rhoc, double Pc){
#define JMAX 40 //max number of iterations
#define TOL 1.0e-7 //error tolerance
    double Pc_cgs = Pc * UnitPressure_in_cgs;
    double rhoc_cgs = rhoc * UnitDensity_in_cgs;
    double Tgas, Trad;
    int j;
    double dx, f, fmid, xmid, rtb;
    double x1, x2;
    char errstr[256];
    x1 = 0.0;
    Tgas = Pc_cgs*MU*PROTONMASS/(BOLTZMANN*rhoc_cgs);
    Trad = pow(3*Pc_cgs/A_SB, 0.25);
    x2 = (Trad > Tgas) ? 2*Trad : 2*Tgas;
    f = Pc_cgs - rhoc_cgs*BOLTZMANN*x1/(MU*PROTONMASS) - A_SB*pow(x1,4)/3.0;
    fmid=Pc_cgs - rhoc_cgs*BOLTZMANN*x2/(MU*PROTONMASS) - A_SB*pow(x2,4)/3.0;
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=1;j<=JMAX;j++) {
        xmid=rtb+(dx *= 0.5);
        fmid = Pc_cgs - rhoc_cgs*BOLTZMANN*xmid/(MU*PROTONMASS) - A_SB*pow(xmid,4)/3.0;
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < TOL*fabs(xmid) || fmid == 0.0) return rtb;
    }
    printf("ps_Tc: bisection solve didn't converge, P_c = %e, rho_c = %e, Tgas = %e Trad = %e",Pc_cgs, rhoc_cgs, Tgas, Trad);
    return(-1);
#undef JMAX
#undef TOL
}


/*Calculate the mean ratio of the gas pressure to the gas+radiation pressure, either by solving the Eddington quartic (for n_ad=3, Eq B5) or by using tabulated values, based on Offner 2009, code taken from ORION*/
double ps_beta(double m, double n_ad, double rhoc, double Pc){
    double mass=m*All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS);//in units of solar mass
if (n_ad==3.0) {
    // In this case we solve the Eddington quartic, P_c^3 = (3/a) (k / (mu mH))^4 (1 - beta) / beta^4 rho_c^4 for beta
#define JMAX 40
#define BETAMIN 1.0e-4
#define BETAMAX 1.0
#define TOL 1.0e-7
    int j;
    double dx, f, fmid, xmid, rtb;
    double x1, x2;
    double coef;

    coef = 3.0/A_SB*pow(BOLTZMANN*rhoc/(MU*PROTONMASS),4);
    x1=BETAMIN; x2=BETAMAX;
    f = pow(Pc,3) - coef * (1.0-x1)/pow(x1,4);
    fmid = pow(Pc,3) - coef * (1.0-x2)/pow(x2,4);
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=1;j<=JMAX;j++) {
      xmid=rtb+(dx *= 0.5);
      fmid = pow(Pc,3) - coef * (1.0-xmid)/pow(xmid,4);
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < TOL*fabs(xmid) || fmid == 0.0) return rtb;
    }
    printf("ps_beta: bisection solve failed to converge");
    return(-1);
#undef JMAX
#undef BETAMIN
#undef BETAMAX
#undef TOL
  } else {
    // For n != 3, we use a table lookup. The values of beta have been pre-computed with mathematica. The table goes from M=5 to 50 solar masses in steps of 2.5 M_sun, and from n=1.5 to n=3 in steps of 0.5. We should never call this routine with M > 50 Msun, since by then the star should be fully on the main sequence.
#define MTABMIN  5.0)
#define MTABMAX  50.0
#define MTABSTEP 2.5
#define NTABMIN  1.5
#define NTABMAX  3.0
#define NTABSTEP 0.5
    if (mass < MTABMIN) return(1.0);  // Set beta = 1 for M < 5 Msun
    if ((mass >= MTABMAX) || (n >= NTABMAX)) {
        printf("ps_beta: too high protostar mass, m: %g n_ad %g",mass, n_ad)
        return(-1.0);
    }
    static double betatab[19][4] = {
      {0.98785, 0.988928, 0.98947, 0.989634}, 
      {0.97438, 0.976428, 0.977462, 0.977774}, 
      {0.957927, 0.960895, 0.962397, 0.962846}, 
      {0.939787, 0.943497, 0.945369, 0.945922}, 
      {0.92091, 0.925151, 0.927276, 0.927896}, 
      {0.901932, 0.906512, 0.908785, 0.909436}, 
      {0.883254, 0.888017, 0.890353, 0.891013}, 
      {0.865111, 0.86994, 0.872277, 0.872927}, 
      {0.847635, 0.852445, 0.854739, 0.855367}, 
      {0.830886, 0.835619, 0.837842, 0.838441}, 
      {0.814885, 0.8195, 0.821635, 0.822201}, 
      {0.799625, 0.804095, 0.806133, 0.806664}, 
      {0.785082, 0.789394, 0.791328, 0.791825}, 
      {0.771226, 0.775371, 0.777202, 0.777665}, 
      {0.758022, 0.761997, 0.763726, 0.764156}, 
      {0.745433, 0.749238, 0.750869, 0.751268}, 
      {0.733423, 0.73706, 0.738596, 0.738966}, 
      {0.721954, 0.725429, 0.726874, 0.727216}, 
      {0.710993, 0.714311, 0.715671, 0.715987}
    };

    // Locate ourselves on the table and do a linear interpolation
    int midx = (int) floor((mass-MTABMIN)/MTABSTEP);
    double mwgt = (mass-(MTABMIN+midx*MTABSTEP)) / MTABSTEP;
    int nidx = (int) floor((n_ad-NTABMIN)/NTABSTEP);
    double nwgt = (n_ad-(NTABMIN+nidx*NTABSTEP)) / NTABSTEP;
    return ( betatab[midx][nidx]*(1.0-mwgt)*(1.0-nwgt) +
	     betatab[midx+1][nidx]*mwgt*(1.0-nwgt) +
	     betatab[midx][nidx+1]*(1.0-mwgt)*nwgt +
	     betatab[midx+1][nidx+1]*mwgt*nwgt );
  }
#undef MTABMIN
#undef MTABMAX
#undef MTABSTEP
#undef NTABMIN
#undef NTABMAX
#undef NTABSTEP
    
}

/*Calculate the mean ratio of the gas pressure to the gas+radiation pressure at the center, based on Offner 2009, code taken from ORION*/
double inline ps_betac(double rhoc, double Pc, double Tc){
    return( rhoc*BOLTZMANN*Tc/(MU*PROTONMASS) / Pc );
}


    
#define DM (0.01*m)
/*Calculate dlog beta/d logm by taking a numerical derivative, based on Offner 2009, code taken from ORION*/
double ps_dlogbeta_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc){
    double rhoc2 = ps_rhoc( (m+DM) , n_ad, r) //slight imprecision here as we do not update the radius
    double Pc2 = ps_Pc( (m+DM) , n_ad, r) //slight imprecision here as we do not update the radius
    double beta2 = ps_beta( (m+dm), n_ad, rhoc2, Pc2);
    return( m/beta * (beta2-beta) / DM );
}
/*Calculate dlog (beta/betac)/d logm by taking a numerical derivative, based on Offner 2009, code taken from ORION*/
double ps_dlogbetaperbetac_dlogm(double m, double r, double n_ad, double beta, double rhoc, double Pc, double Tc){
    double betac = ps_betac(rhoc, Pc, Tc)
    double rhoc2 = ps_rhoc( (m+DM) , n_ad, r) //slight imprecision here as we do not update the radius
    double Pc2 = ps_Pc( (m+DM) , n_ad, r) //slight imprecision here as we do not update the radius
    double Tc2 = ps_Tc(rhoc2, Pc2)
    double beta2 = ps_beta( (m+dm), n_ad, rhoc2, Pc2);
    double betac2 = ps_betac(rhoc2, Pc2, Tc2)
    return( m/(beta/betac) * ((beta2/betac2) - (beta/betac)) / DM );

}
#undef DM

/*Calculate the luminosity required to ionize the infalling material, based on Offner 2009*/
double ps_lum_I(double mdot){
    double mdot_m_solar_per_year = mdot * (All.UnitMass_in_g/(All.HubbleParam * SOLAR_MASS))/UnitTime_in_s*SEC_PER_YEAR; // accretion rate in msolar/yr
    return (2.5*SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s) * ( mdot_m_solar_per_year/(1e-5));
}
/*Calculate the blackbody luminosity of the star following the Hayashi track*/
double ps_lum_Hayashi_KH(double m, double r){
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    double T4000_4 = pow(m_solar , 0.55); // protostellar temperature along Hayashi track
    return (0.2263 * r * r * T4000_4 * SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)); // luminosity from KH contraction
}
/*Calculate the luminosity of a main sequence star*/
double ps_lum_MS(double m){
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    double lum_sol = 0;
    
/***************************************************/
/*  Original GIZMO version  */
/*     if(m_solar > 0.012)
        {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;}
        else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
        else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);}
        else {lum_sol = 32000. * m_solar;}
    } */
    
/***************************************************/
/*  ORION version  */
/*Calculate the luminosity of a main sequence star using the fitting formulas of Tout et al (1996), code taken from ORION*/
// Parameters for the main sequence luminosity and radius fitting formulae
#define MS_ALPHA    0.39704170
#define MS_BETA     8.52762600
#define MS_GAMMA    0.00025546
#define MS_DELTA    5.43288900
#define MS_EPSILON  5.56357900
#define MS_ZETA     0.78866060
#define MS_ETA      0.00586685 
    if(m_solar > 0.1){ //minimum mass for the Tout fitting function
        lum_sol = (MS_ALPHA*pow(m_solar,5.5) + MS_BETA*pow(m_solar,11)) / (MS_GAMMA+pow(m_solar,3)+MS_DELTA*pow(m_solar,5) 
        + MS_EPSILON*pow(m_solar,7) + MS_ZETA*pow(m_solar,8)+MS_ETA*pow(m_solar,9.5));
    }
    return(lum_sol*SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s));
#undef MS_ALPHA
#undef MS_BETA
#undef MS_GAMMA
#undef MS_DELTA
#undef MS_EPSILON
#undef MS_ZETA
#undef MS_ETA
    
    
}

/*Calculate the radius of a main sequence star*/
double ps_radius_MS_in_solar(double m){
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
/***************************************************/
/*  ORION version  */
/*Calculate the luminosity of a main sequence star using the fitting formulas of Tout et al (1996), code taken from ORION*/
// Parameters for the main sequence luminosity and radius fitting formulae
#define MS_THETA    1.71535900
#define MS_IOTA     6.59778800
#define MS_KAPPA   10.08855000
#define MS_LAMBDA   1.01249500
#define MS_MU       0.07490166
#define MS_NU       0.01077422
#define MS_XI       3.08223400
#define MS_UPSILON 17.84778000
#define MS_PI       0.00022582
    double rsol = (MS_THETA*pow(m_solar,2.5)+MS_IOTA*pow(m_solar,6.5)+MS_KAPPA*pow(m_solar,11)+
	       MS_LAMBDA*pow(m_solar,19)+MS_MU*pow(m_solar,19.5)) /
    (MS_NU+MS_XI*pow(m_solar,2)+UPSILON*pow(m_solar,8.5)+pow(m_solar,18.5)+
     MS_PI*pow(m_solar,19.5));
  return(rsol); //*SOLAR_RADIUS/All.UnitLength_in_cm);
#undef MS_THETA
#undef MS_IOTA
#undef MS_KAPPA
#undef MS_LAMBDA
#undef MS_MU
#undef MS_NU
#undef MS_XI
#undef MS_UPSILON
#undef MS_PI
}

//remove definitions we don't use
#undef A_SB
#undef MU
#endif //end of protostellar evolution functions


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
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 0)
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
            lum = P[i].StarLuminosity; //get pre-calculated luminosity of the star
#endif
        }
    }
#endif
#if (SINGLE_STAR_PROTOSTELLAR_EVOLUTION == 0)
    lum_sol *= SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
    lum += lum_sol;
    P[i].StarLuminosity = lum; //store total luminosity of the star
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

    
#endif /* GALSF */
