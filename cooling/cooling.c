#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#include "./cooling.h"


/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM, 
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including 
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these 
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel (volker.springel@h-its.org). The code has been modified heavily by 
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; everything except the original metal-free free-free and 
 *   photo-ionization heating physics has been added (or re-written), and the iteration routine to converge to 
 *   temperatures has been significantly modified.
 */


#ifdef COOLING

#define NCOOLTAB  2000

#define SMALLNUM 1.0e-60
#define COOLLIM  0.1
#define HEATLIM	 20.0


static double XH = HYDROGEN_MASSFRAC;	/* hydrogen abundance by mass */
static double yhelium;

#define eV_to_K   11606.0
#define eV_to_erg 1.60184e-12

/* CAFG: H number density above which we assume no ionizing bkg (proper cm^-3) */
#define NH_SS 0.0123

static double mhboltz;		/* hydrogen mass over Boltzmann constant */
static double ethmin;		/* minimum internal energy for neutral gas */

static double Tmin = 0.0;	/* in log10 */
static double Tmax = 9.0;
static double deltaT;

static double *BetaH0, *BetaHep, *Betaff;
static double *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
static double *GammaeH0, *GammaeHe0, *GammaeHep;
#ifdef COOL_METAL_LINES_BY_SPECIES
/* if this is enabled, the cooling table files should be in a folder named 'spcool_tables' in the run directory.
 cooling tables can be downloaded at: https://dl.dropbox.com/u/16659252/spcool_tables.tgz */
static float *SpCoolTable0;
static float *SpCoolTable1;
#endif

static double J_UV = 0, gJH0 = 0, gJHep = 0, gJHe0 = 0, epsH0 = 0, epsHep = 0, epsHe0 = 0;

static double ne, necgs, nHcgs;
static double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
static double gJH0ne, gJHe0ne, gJHepne;
static double nH0, nHp, nHep, nHe0, nHepp;

static double DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;




/* this is just a simple loop if all we're doing is cooling (no star formation) */
void cooling_only(void)
{
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            do_the_cooling_for_particle(i);
        } // if(P[i].Type == 0 && P[i].Mass > 0)
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
} // void cooling_only(void)





/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle(int i)
{
    double unew;
    double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
    double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
    if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
        
        double ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */
        double uold = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);
#ifdef GALSF_FB_HII_HEATING
        double u_to_temp_fac = PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
        double uion = HIIRegion_Temp / u_to_temp_fac;
        if(SphP[i].DelayTimeHII > 0) if(uold<uion) uold=uion; /* u_old should be >= ionized temp if used here */
#endif // GALSF_FB_HII_HEATING
        
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
        double du = SphP[i].DtInternalEnergy * dtime;
        if(du < -0.5*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy = -0.5*SphP[i].InternalEnergy / dtime;}
        if(du >  50.*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy =  50.*SphP[i].InternalEnergy / dtime;}
        /* and convert to cgs before use in the cooling sub-routine */
        SphP[i].DtInternalEnergy *= All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH);
#endif
        
        
#ifndef RT_COOLING_PHOTOHEATING
        unew = DoCooling(uold, SphP[i].Density * All.cf_a3inv, dtime, &ne, i);
#else
        double fac_entr_to_u = pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
        unew = uold + dt * fac_entr_to_u * (rt_DoHeating(i, dt) + rt_DoCooling(i, dt));
#endif // RT_COOLING_PHOTOHEATING
        
        
#ifdef GALSF_FB_HII_HEATING
        /* set internal energy to minimum level if marked as ionized by stars */
        if(SphP[i].DelayTimeHII > 0)
        {
            if(unew<uion)
            {
                unew=uion;
                if(SphP[i].DtInternalEnergy<0) SphP[i].DtInternalEnergy=0;
                //if(SphP[i].dInternalEnergy<0) SphP[i].dInternalEnergy=0; //manifest-indiv-timestep-debug//
            }
            SphP[i].Ne = 1.0 + 2.0*yhelium;
        }
#endif // GALSF_FB_HII_HEATING
        
        
#if defined(BH_THERMALFEEDBACK)
        if(SphP[i].Injected_BH_Energy)
		{
            unew += SphP[i].Injected_BH_Energy / P[i].Mass;
            SphP[i].Injected_BH_Energy = 0;
		}
#endif
        

#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_DISABLE_COOLING)
        /* cosmic ray interactions affecting the -thermal- temperature of the gas are included in the actual cooling/heating functions; 
            they are solved implicitly above. however we need to account for energy losses of the actual cosmic ray fluid, here. The 
            timescale for this is reasonably long, so we can treat it semi-explicitly, as we do here.
            -- We use the estimate for combined hadronic + Coulomb losses from Volk 1996, Ensslin 1997, as updated in Guo & Oh 2008: */
        double ne_cgs = ((0.78 + 0.22*ne*XH) / PROTONMASS) * (SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
        double CR_coolingrate_perunitenergy = -7.51e-16 * ne_cgs * (All.UnitTime_in_s / All.HubbleParam); // converts cgs to code units //
        double CR_Egy_new = SphP[i].CosmicRayEnergyPred * exp(CR_coolingrate_perunitenergy * dtime);
        SphP[i].CosmicRayEnergyPred = SphP[i].CosmicRayEnergy = CR_Egy_new;
#endif
        
        
        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated. 
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        SphP[i].InternalEnergy = unew;
        SphP[i].Ne = ne;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
#ifndef COOLING_OPERATOR_SPLIT
        SphP[i].DtInternalEnergy = 0;
#endif
        
        
#ifdef GALSF_FB_HII_HEATING
        /* count off time which has passed since ionization 'clock' */
        if(SphP[i].DelayTimeHII > 0) SphP[i].DelayTimeHII -= dtime;
        if(SphP[i].DelayTimeHII < 0) SphP[i].DelayTimeHII = 0;
#endif // GALSF_FB_HII_HEATING
        
    } // closes if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}






/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling(double u_old, double rho, double dt, double *ne_guess, int target)
{
  double u, du;
  double u_lower, u_upper;
  double ratefact;
  double LambdaNet;
  int iter=0, iter_upper=0, iter_lower=0;

#ifdef GRACKLE
#ifndef COOLING_OPERATOR_SPLIT
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below
     to split the step */
    du = dt * SphP[target].DtInternalEnergy / (All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH));
    u_old += 0.5*du;
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
    /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
     term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
     the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
     to which we then add the heating term from hydro and compute the solution over a full timestep */
    double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(1-r)/log(r);} r=DMAX(0,DMIN(r,1));
    du*=0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
    return DMAX(u,All.MinEgySpec);
#endif
    
    
  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;


  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;
  u_lower = u;
  u_upper = u;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess, target);

  /* bracketing */

  if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while((iter_upper<MAXITER)&&(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess, target) * dt < 0))
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
        iter_upper++;
	}

    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while((iter_lower<MAXITER)&&(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess, target) * dt > 0))
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
        iter_lower++;
	}
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = CoolingRateFromU(u, rho, ne_guess, target);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	{
	  u_upper = u;
	}
      else
	{
	  u_lower = u;
	}

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %g\n", u);
    }
    while(((fabs(du/u) > 1.0e-3)||((fabs(du/u) > 1.0e-6)&&(iter < 10))) && (iter < MAXITER));
    //while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(10);
    }

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;	/* to internal units */

  return u;
}



/* returns cooling time. 
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess, int target)
{
    double u;
    double ratefact;
    double LambdaNet, coolingtime;
    
#if defined(GRACKLE) && !defined(GALSF_EFFECTIVE_EQS)
    coolingtime = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(coolingtime >= 0) coolingtime = 0.0;
    coolingtime *= All.HubbleParam / All.UnitTime_in_s;
    return coolingtime;
#endif
    
    DoCool_u_old_input = u_old;
    DoCool_rho_input = rho;
    DoCool_ne_guess_input = *ne_guess;
    
    rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
    u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
    
    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;
    u = u_old;
    LambdaNet = CoolingRateFromU(u, rho, ne_guess, target);
    
    /* bracketing */
    
    if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
        return 0;
    
    coolingtime = u_old / (-ratefact * LambdaNet);
    
    coolingtime *= All.HubbleParam / All.UnitTime_in_s;
    
    return coolingtime;
}


/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess, int target)
{
  double m, dm;
  double m_lower, m_upper;
  double ratefact;
  double LambdaNet;
  int iter = 0;

  DoCool_u_old_input = u;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;

  if(fac <= 0)			/* the hot phase is actually colder than the cold reservoir! */
    {
      return 0.01 * m_old;
    }

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;
  fac *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho * fac;

  m = m_old;
  m_lower = m;
  m_upper = m;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess, target);

  /* bracketing */

  if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
      m_upper *= sqrt(1.1);
      m_lower /= sqrt(1.1);
      while(m_upper - m_old -
          m_upper * m_upper / m_old * ratefact * CoolingRateFromU(u, rho * m_upper / m_old,
                                                                  ne_guess, target) * dt < 0)
      {
	m_upper *= 1.1;
	m_lower *= 1.1;
      }
    }

  if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
      m_lower /= sqrt(1.1);
      m_upper *= sqrt(1.1);
      while(m_lower - m_old -
          m_lower * m_lower / m_old * ratefact * CoolingRateFromU(u, rho * m_lower / m_old,
                                                                  ne_guess, target) * dt > 0)
      {
	m_upper /= 1.1;
	m_lower /= 1.1;
      }
    }

  do
    {
      m = 0.5 * (m_lower + m_upper);

        LambdaNet = CoolingRateFromU(u, rho * m / m_old, ne_guess, target);

      if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0)
	{
	  m_upper = m;
	}
      else
	{
	  m_lower = m;
	}

      dm = m_upper - m_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("m= %g\n", m);
    }
  while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      printf("m_old= %g\n", m_old);
      endrun(11);
    }

  return m;
}





void cool_test(void)
{
#if !defined(COOL_METAL_LINES_BY_SPECIES)
    double uin, rhoin, tempin, muin, nein;
    
    uin = 6.01329e+09;
    rhoin = 7.85767e-29;
    tempin = 2034.0025;
    muin = 0.691955;
    nein = (1 + 4 * yhelium) / muin - (1 + yhelium);
    
    double dtin=1.0e-7;
    double uout,uint;
    int i,target;
    for(i=0;i<20;i++) {
        rhoin=SphP[i].Density;
        nein=SphP[i].Ne;
        target=i;
        uin=SphP[i].InternalEnergy;
        uout=DoCooling(uin,rhoin,dtin,&nein,target);
        printf("%d %d : ne: %g %g \n",ThisTask,target,SphP[i].Ne,nein);
        nein=SphP[i].Ne;
        rhoin *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;    /* convert to physical cgs units */
        uint = uin*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
        tempin=convert_u_to_temp(uint, rhoin, &nein, target);
        printf("%d %d : in: : %g %g %g \n",ThisTask,target,uin,rhoin,nein);
        printf("%d %d : out: %g %g %g %g %g \n",ThisTask,target,tempin,
               CoolingRate(log10(tempin),rhoin,&nein,target),
               CoolingRateFromU(uint,rhoin,&nein,target),
               uout,nein);
        fflush(stdout);
    }
#endif
}





/* this function determines the electron fraction, and hence the mean 
 * molecular weight. With it arrives at a self-consistent temperature.
 * Element abundances and the rates for the emission are also computed
 */
double convert_u_to_temp(double u, double rho, double *ne_guess, int target)
{
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu;
  int iter = 0;

  double u_input, rho_input, ne_input;

  u_input = u;
  rho_input = rho;
  ne_input = *ne_guess;

  mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do
    {
      ne_old = *ne_guess;

      find_abundances_and_rates(log10(temp), rho, ne_guess, target);
      temp_old = temp;

      mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max =
	DMAX(max,
	     temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
	  printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
    while(
          ((fabs(temp - temp_old) > 0.1 * temp) ||
           ((fabs(temp - temp_old) > 1.0e-3 * temp) && (temp > 200.))) && iter < MAXITER);

  if(iter >= MAXITER)
      {
	printf("failed to converge in convert_u_to_temp()\n");
	printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
	printf
	  ("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	   DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);

	endrun(12);
      }

    if(temp<=0) temp=pow(10.0,Tmin);
    if(log10(temp)<Tmin) temp=pow(10.0,Tmin);

  return temp;
}



/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess, int target)
{
  double neold, nenew;
  int j, niter;
  double Tlow, Thi, flow, fhi, t;
  double logT_input, rho_input, ne_input;
  double NH_SS_z=NH_SS, shieldfac;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  if(isnan(logT)) logT=Tmin;    /* nan trap (just in case) */
    
  if(logT <= Tmin)		/* everything neutral */
    {
      nH0 = 1.0;
      nHe0 = yhelium;
      nHp = 0;
      nHep = 0;
      nHepp = 0;
      ne = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Tmax)		/* everything is ionized */
    {
      nH0 = 0;
      nHe0 = 0;
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *ne_guess = ne;		/* note: in units of the hydrogen number density */
      return;
    }

  t = (logT - Tmin) / deltaT;
  j = (int) t;
    if(j<0){j=0;}
    if(j>NCOOLTAB){
        printf("warning: j>NCOOLTAB : j=%d t %g Tlow %g Thi %g logT %g Tmin %g deltaT %g \n",j,t,Tmin+deltaT*j,Tmin+deltaT*(j+1),logT,Tmin,deltaT);fflush(stdout);
        j=NCOOLTAB;
    }
  Tlow = Tmin + deltaT * j;
  Thi = Tlow + deltaT;
  fhi = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

    double local_gammamultiplier=1;
#ifdef GALSF_FB_LOCAL_UV_HEATING
    if ((target >= 0) && (gJH0 > 0))
    {
    local_gammamultiplier = SphP[target].RadFluxEUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
    local_gammamultiplier = 1 + local_gammamultiplier/gJH0;
    }
#endif
    
    /* CAFG: this is the density that we should use for UV background threshold */
    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    if(gJH0>0)
        NH_SS_z = NH_SS*pow(local_gammamultiplier*gJH0/1.0e-12,0.66)*pow(10.,0.173*(logT-4.));
    else
        NH_SS_z = NH_SS*pow(10.,0.173*(logT-4.));
    if(nHcgs<100.*NH_SS_z) shieldfac=exp(-nHcgs/NH_SS_z); else shieldfac=0;
#ifdef COOL_LOW_TEMPERATURES
    if(logT < Tmin+1) shieldfac *= (logT-Tmin); // make cutoff towards Tmin more continuous //
    //shieldfac *= 1 - exp(Tmin-logT);
#endif
#ifdef GALSF_EFFECTIVE_EQS
    shieldfac = 1; // self-shielding is implicit in the sub-grid model already //
#endif

  ne = *ne_guess;
  neold = ne;
  niter = 0;
  necgs = ne * nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      aHp = flow * AlphaHp[j] + fhi * AlphaHp[j + 1];
      aHep = flow * AlphaHep[j] + fhi * AlphaHep[j + 1];
      aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
      ad = flow * Alphad[j] + fhi * Alphad[j + 1];
      geH0 = flow * GammaeH0[j] + fhi * GammaeH0[j + 1];
      geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
      geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];
#ifdef COOL_LOW_TEMPERATURES
        // make cutoff towards Tmin more continuous //
        if(logT < Tmin+1) {
            geH0 *= (logT-Tmin);
            geHe0 *= (logT-Tmin);
            geHep *= (logT-Tmin);
        }
#endif

      if(necgs <= 1.e-25 || J_UV == 0)
	{
	  gJH0ne = gJHe0ne = gJHepne = 0;
	}
      else
	{
        /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
        gJH0ne = gJH0 * local_gammamultiplier / necgs * shieldfac;
        gJHe0ne = gJHe0 * local_gammamultiplier / necgs * shieldfac;
        gJHepne = gJHep * local_gammamultiplier / necgs * shieldfac;
	}

      nH0 = aHp / (aHp + geH0 + gJH0ne);	/* eqn (33) */
      nHp = 1.0 - nH0;		/* eqn (34) */

      if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  nHep = 0.0;
	  nHepp = 0.0;
	  nHe0 = yhelium;
	}
      else
	{
	  nHep = yhelium / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
	  nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
	  nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
	}

      neold = ne;

      ne = nHp + nHep + 2 * nHepp;	/* eqn (38) */
      necgs = ne * nHcgs;

      if(J_UV == 0)
	break;

      nenew = 0.5 * (ne + neold);
      ne = nenew;
      necgs = ne * nHcgs;

      if(fabs(ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d\n", ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %g  rho_input= %g  ne_input= %g\n", logT_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(13);
    }

  bH0 = flow * BetaH0[j] + fhi * BetaH0[j + 1];
  bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
  bff = flow * Betaff[j] + fhi * Betaff[j + 1];

  *ne_guess = ne;
}




/*  this function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates 
 *  (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRateFromU(double u, double rho, double *ne_guess, int target)
{
  double temp;
  temp = convert_u_to_temp(u, rho, ne_guess, target);

    return CoolingRate(log10(temp), rho, ne_guess, target);
}


/*  this function computes the self-consistent temperature
 *  and abundance ratios 
 */
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer, int target)
{
  double temp;

  DoCool_u_old_input = u;
  DoCool_rho_input = rho;
  DoCool_ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  temp = convert_u_to_temp(u, rho, ne_guess, target);

  *nH0_pointer = nH0;
  *nHeII_pointer = nHep;

  return temp;
}




extern FILE *fd;





/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRate(double logT, double rho, double *nelec, int target)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;
  double NH_SS_z=NH_SS,shieldfac;
#ifdef COOL_LOW_TEMPERATURES
  double LambdaMol=0;
#endif
#ifdef COOL_METAL_LINES_BY_SPECIES
  double LambdaMetal=0;
  double *Z;
  if(target>=0)
  {
      Z = P[target].Metallicity;
  } else {
      /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
      int k;
      double Zsol[NUM_METAL_SPECIES];
      for(k=0;k<NUM_METAL_SPECIES;k++) Zsol[k]=All.SolarAbundances[k];
      Z = Zsol;
  }
#endif
  double local_gammamultiplier=1;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */


  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

#ifdef GALSF_FB_LOCAL_UV_HEATING
    double LambdaPElec,photoelec=0;
    if((target >= 0) && (gJH0 > 0))
    {
    local_gammamultiplier = SphP[target].RadFluxEUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
    local_gammamultiplier = 1 + local_gammamultiplier/gJH0;
    }
    if(target >= 0) photoelec=SphP[target].RadFluxUV;
#endif
    
    /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
    if(J_UV != 0)
        NH_SS_z=NH_SS*pow(local_gammamultiplier*gJH0/1.0e-12,0.66)*pow(10.,0.173*(logT-4.));
    else
        NH_SS_z=NH_SS*pow(10.,0.173*(logT-4.));
    if(nHcgs<100.*NH_SS_z) shieldfac=exp(-nHcgs/NH_SS_z); else shieldfac=0;
#ifdef GALSF_EFFECTIVE_EQS
    shieldfac = 1; // self-shielding is implicit in the sub-grid model already //
#endif
    
#ifdef BH_COMPTON_HEATING
    double AGN_LambdaPre,AGN_T_Compton;
    AGN_T_Compton = 2.0e7; /* approximate from Sazonov et al. */
    if(target < 0) {
        AGN_LambdaPre = 0;
    } else {
        AGN_LambdaPre = SphP[target].RadFluxAGN * (3.9/2.0) * All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm)*All.HubbleParam*All.cf_a2inv; /* proper units */
        /* now have incident flux, need to convert to relevant pre-factor for heating rate */
        AGN_LambdaPre *= 6.652e-25; /* sigma_T for absorption */
        AGN_LambdaPre *= (4.*1.381e-16)/(9.109e-28*2.998e10*2.998e10); /* times 4*k_B/(me*c^2) */
    }
#endif

    
  T = pow(10.0, logT);
  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec, target);
      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        
      LambdaExcH0 = bH0 * ne * nH0;
      LambdaExcHep = bHep * ne * nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

      LambdaIonH0 = 2.18e-11 * geH0 * ne * nH0;
      LambdaIonHe0 = 3.94e-11 * geHe0 * ne * nHe0;
      LambdaIonHep = 8.72e-11 * geHep * ne * nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

      LambdaRecHp = 1.036e-16 * T * ne * (aHp * nHp);
      LambdaRecHep = 1.036e-16 * T * ne * (aHep * nHep);
      LambdaRecHepp = 1.036e-16 * T * ne * (aHepp * nHepp);
      LambdaRecHepd = 6.526e-11 * ad * ne * nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

      LambdaFF = bff * (nHp + nHep + 4 * nHepp) * ne;

      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

#ifdef COOL_METAL_LINES_BY_SPECIES
        //if((logT > Tmin+0.5*deltaT)&&((logT > 3.87)||(nHcgs<NH_SS_z)))
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
        if((J_UV != 0)&&(logT > Tmin+0.5*deltaT)&&(logT > 4.00))
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies(nHcgs, logT, Z); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= ne;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted;
             ni not used b/c it should vary species-to-species */
            Lambda += LambdaMetal;
        }
#endif
        
#ifdef COOL_LOW_TEMPERATURES
        //if((nHcgs>NH_SS_z)&&(logT <= 5.2)&&(logT > Tmin+0.5*deltaT))
        if((logT <= 5.2)&&(logT > Tmin+0.5*deltaT))
        {
            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
            LambdaMol = //10.0 *
            2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));//*nHcgs*nHcgs;
            LambdaMol *= (1-shieldfac);
            double LambdaDust = 0;
#ifdef COOL_METAL_LINES_BY_SPECIES
            LambdaMol *= (1+Z[0]/All.SolarAbundances[0])*(0.001 + 0.1*nHcgs/(1.0+nHcgs)
                            + 0.09*nHcgs/(1.0+0.1*nHcgs)
                            + (Z[0]/All.SolarAbundances[0])*(Z[0]/All.SolarAbundances[0])/(1.0+nHcgs));
            /* add dust cooling as well */
            double Tdust = 30.;
            if(T > Tdust) {LambdaDust = 0.63e-33 * (T-Tdust) * sqrt(T) * (Z[0]/All.SolarAbundances[0]);}
#endif
            Lambda += LambdaMol + LambdaDust;
            
        }
#endif
        
        
      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;

	  Lambda += LambdaCmptn;
	}
      else
	LambdaCmptn = 0;

#ifdef BH_COMPTON_HEATING
	if(T > AGN_T_Compton)
	{
        LambdaCmptn = AGN_LambdaPre * (T - AGN_T_Compton) * ne/nHcgs;
        if(LambdaCmptn > 2.19e-21/sqrt(T/1.0e8)) LambdaCmptn=2.19e-21/sqrt(T/1.0e8);
        Lambda += LambdaCmptn;
	}
#endif
        
      Heat = 0;
        if(J_UV != 0) {
            /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
            Heat += local_gammamultiplier * (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs * shieldfac;
        }
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_DISABLE_COOLING)
        if(SphP[target].CosmicRayEnergyPred > 0)
        {
            /* cosmic ray heating, from Guo & Oh 2008: this scales proportional to the electron number density and 
                cosmic ray energy density, both of which we quickly evaluate here (make sure we convert to the correct per-atom units) 
                - note that only 1/6 of the hadronic cooling is thermalized, according to their calculation, while all the Coulomb losses heat */
            double Gamma_CR = 1.0e-16 * (0.98 + 1.65*ne*XH) / nHcgs *
                ((SphP[target].CosmicRayEnergyPred / P[target].Mass * SphP[target].Density * All.cf_a3inv) *
                 (All.UnitPressure_in_cgs * All.HubbleParam * All.HubbleParam));
            Heat += Gamma_CR;
        }
#endif
        
#ifdef COOL_LOW_TEMPERATURES
#if !defined(COSMIC_RAYS) || defined(COSMIC_RAYS_DISABLE_COOLING)
        /* if COSMIC_RAYS is not enabled, but low-temperature cooling is on, we account for the CRs as a heating source using
         a more approximate expression (assuming the mean background of the Milky Way clouds) */
        if(logT <= 5.2)
        {
            double Gamma_CR = 1.0e-16 * (0.98 + 1.65*ne*XH) / (1.e-2 + nHcgs) * 9.0e-12;
            // multiplied by background of ~5eV/cm^3 (Goldsmith & Langer (1978),  van Dishoeck & Black (1986) //
            Heat += Gamma_CR;
        }
#endif
#ifdef COOL_METAL_LINES_BY_SPECIES
        /* add dust heating as well */
        double Tdust = 30.;
        if(T < Tdust) {Heat += 0.63e-33 * (Tdust-T) * sqrt(Tdust) * (Z[0]/All.SolarAbundances[0]);}
#endif
#endif
        
        
        
#ifdef BH_COMPTON_HEATING
        if(T < AGN_T_Compton) Heat += AGN_LambdaPre * (AGN_T_Compton - T) / nHcgs;
        /* note this is independent of the free electron fraction */
#endif
#ifdef GALSF_FB_LOCAL_UV_HEATING
        /* photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995) */
        if(T < 1.0e6) {
            LambdaPElec = 1.0e-24*photoelec/nHcgs;
            photoelec *= sqrt(T)/nHcgs;
            LambdaPElec *= 0.049/(1+pow(photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+photoelec/5000.);
            Heat += LambdaPElec;
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.  
         Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
	LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *nelec = ne;		/* note: in units of the hydrogen number density */

      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * ne;

      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  /* add inverse Compton cooling off the microwave background */
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;
	}
      else
	LambdaCmptn = 0;

#ifdef BH_COMPTON_HEATING
        //LambdaCmptn += AGN_LambdaPre * (T - AGN_T_Compton) * ne/nHcgs;
        /* actually at these temperatures want approximation to relativistic compton cooling */
        LambdaCmptn += AGN_LambdaPre * (T - AGN_T_Compton) * (T/1.5e9)/(1-exp(-T/1.5e9)) * ne/nHcgs;
#endif
        
      Lambda = LambdaFF + LambdaCmptn;

      /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be
         the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
      if(Lambda > 2.19e-21/sqrt(T/1.0e8)) Lambda=2.19e-21/sqrt(T/1.0e8);
    }

  /*      
     printf("Lambda= %g\n", Lambda);

     fprintf(fd,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", pow(10, logT),Lambda,
     LambdaExcH0, LambdaExcHep, 
     LambdaIonH0, LambdaIonHe0, LambdaIonHep,
     LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd,
     LambdaFF, LambdaCmptn, Heat,
     ne, nHp, nHep, nHepp);
   */
    
    double Q = Heat - Lambda;
#ifdef COOL_LOW_TEMPERATURES
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
    //if(nHcgs > 0.1) /* don't bother at very low densities, since youre not optically thick */  
    if( (nHcgs > 0.1) && (target >= 0) )  // DAA: protect from target=-1 with GALSF_EFFECTIVE_EQS
    {
        double surface_density = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1);
        surface_density *= All.cf_a2inv * All.UnitDensity_in_cgs * All.HubbleParam * All.UnitLength_in_cm; // converts to cgs
        double effective_area = 2.3 * PROTONMASS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
        if(T < 1500.)
        {
            if(T < 150.) {kappa_eff=0.0027*T*sqrt(T);} else {kappa_eff=5.;}
            kappa_eff *= P[target].Metallicity[0]/All.SolarAbundances[0];
            if(kappa_eff < 0.1) {kappa_eff=0.1;}
        } else {
            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * ne; /* Thompson scattering (non-relativistic) */
            double k_molecular = 0.1 * P[target].Metallicity[0]; /* molecular line opacities */
            double k_Hminus = 1.1e-25 * sqrt(P[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            double k_conductive = 2.6e-7 * ne * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
            kappa_eff = 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif
    
#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for 
        in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) Q += SphP[target].DtInternalEnergy / nHcgs;
#endif

  return Q;
}





double LogTemp(double u, double ne)	/* ne= electron density in terms of hydrogen density */
{
  double T;

  if(u < ethmin)
    u = ethmin;

  T = log10(GAMMA_MINUS1 * u * mhboltz * (1 + 4 * yhelium) / (1 + ne + yhelium));

  return T;
}



void InitCoolMemory(void)
{
  BetaH0 = (double *) mymalloc("BetaH0", (NCOOLTAB + 1) * sizeof(double));
  BetaHep = (double *) mymalloc("BetaHep", (NCOOLTAB + 1) * sizeof(double));
  AlphaHp = (double *) mymalloc("AlphaHp", (NCOOLTAB + 1) * sizeof(double));
  AlphaHep = (double *) mymalloc("AlphaHep", (NCOOLTAB + 1) * sizeof(double));
  Alphad = (double *) mymalloc("Alphad", (NCOOLTAB + 1) * sizeof(double));
  AlphaHepp = (double *) mymalloc("AlphaHepp", (NCOOLTAB + 1) * sizeof(double));
  GammaeH0 = (double *) mymalloc("GammaeH0", (NCOOLTAB + 1) * sizeof(double));
  GammaeHe0 = (double *) mymalloc("GammaeHe0", (NCOOLTAB + 1) * sizeof(double));
  GammaeHep = (double *) mymalloc("GammaeHep", (NCOOLTAB + 1) * sizeof(double));
  Betaff = (double *) mymalloc("Betaff", (NCOOLTAB + 1) * sizeof(double));

#ifdef COOL_METAL_LINES_BY_SPECIES
  long i_nH=41; long i_T=176; long kspecies=(long)NUM_METAL_SPECIES-1;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //kspecies -= 1;
    kspecies -= NUM_RPROCESS_SPECIES;
#endif
  SpCoolTable0 = (float *) mymalloc("SpCoolTable0",(kspecies*i_nH*i_T)*sizeof(float));
  if(All.ComovingIntegrationOn)
    SpCoolTable1 = (float *) mymalloc("SpCoolTable1",(kspecies*i_nH*i_T)*sizeof(float));
#endif
}



void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19 
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i;
    double T,Tfact;
    XH = 0.76;
    yhelium = (1 - XH) / (4 * XH);
    mhboltz = PROTONMASS / BOLTZMANN;
    
    if(All.MinGasTemp > 0.0)
        Tmin = log10(All.MinGasTemp); // Tmin = log10(0.1 * All.MinGasTemp);
    else
        Tmin = 1.0;
    
    deltaT = (Tmax - Tmin) / NCOOLTAB;
    ethmin = pow(10.0, Tmin) * (1. + yhelium) / ((1. + 4. * yhelium) * mhboltz * GAMMA_MINUS1);
    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB; i++)
    {
        BetaH0[i] = BetaHep[i] = Betaff[i] = AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;
        T = pow(10.0, Tmin + deltaT * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
        
        if(118348 / T < 70) BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;
        if(473638 / T < 70) BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;
        
        Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
        //AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
        //AlphaHep[i] = 1.5e-10 * pow(T, -0.6353); /* old Cen92 fit */
        //AlphaHepp[i] = 4. * AlphaHp[i];	/* old Cen92 fit */
        AlphaHp[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHep[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHp[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*alphaHp[1,T/Z^2] */
        
        if(470000 / T < 70) Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));
        if(157809.1 / T < 70) GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
        if(285335.4 / T < 70) GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
        if(631515.0 / T < 70) GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
        
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables(void)
{
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables(0);
    }
}

void ReadMultiSpeciesTables(int iT)
{
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_METAL_SPECIES-1; long i,j,k,r;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //kspecies -= 1;
    kspecies -= NUM_RPROCESS_SPECIES;
#endif
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;
    
    fname=GetMultiSpeciesFilename(iT,0);
    if(ThisTask == 0) printf("Opening Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read species cooling table in file `%s'\n", fname); endrun(456);}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&SpCoolTable0[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached Cooling EOF! \n");fflush(stdout);}
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename(iT,&fname,1);
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
        fname=GetMultiSpeciesFilename(iT+1,0);
        if(ThisTask == 0) printf("Opening (z+) Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read species 1 cooling table in file `%s'\n", fname); endrun(456);}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&SpCoolTable1[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached Cooling EOF! \n");fflush(stdout);}
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename(iT+1,&fname,1);
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

char *GetMultiSpeciesFilename(int i, int hk)
{
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



/* table input (from file TREECOOL) for ionizing parameters */
/* NOTE: we've switched to using the updated TREECOOL from CAFG, june11 version */

#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 250		/* Max # of lines in TREECOOL */

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static int nheattab;		/* length of table */


void ReadIonizeParams(char *fname)
{
  int i;
  FILE *fdcool;

  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      endrun(456);
    }

  for(i = 0; i < TABLESIZE; i++)
    gH0[i] = 0;

  for(i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%g %g %g %g %g %g %g",
	      &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */

  for(i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      nheattab++;
    else
      break;

  if(ThisTask == 0)
    printf("\n\nread ionization table with %d entries in file `%s'.\n\n", nheattab, fname);
}


void IonizeParams(void)
{
  IonizeParamsTable();

  /*
     IonizeParamsFunction();
   */
}



void IonizeParamsTable(void)
{
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;

  if(All.ComovingIntegrationOn)
    redshift = 1 / All.Time - 1;
  else
    {
    /* in non-cosmological mode, still use, but adopt z=0 background */
    redshift = 0;
    /*
         gJHe0 = gJHep = gJH0 = 0;
         epsHe0 = epsHep = epsH0 = 0;
         J_UV = 0;
         return;
    */
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
	ilow = i;
      else
	break;
    }

  dzlow = logz - inlogz[ilow];
  dzhi = inlogz[ilow + 1] - logz;

  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      gJHe0 = gJHep = gJH0 = 0;
      epsHe0 = epsHep = epsH0 = 0;
      J_UV = 0;
      return;
    }
  else
    J_UV = 1.e-21;		/* irrelevant as long as it's not 0 */

  gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
  gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ilow + 1])) / (dzlow + dzhi));
  gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
  epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ilow + 1])) / (dzlow + dzhi));
  epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ilow + 1])) / (dzlow + dzhi));
  epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ilow + 1])) / (dzlow + dzhi));

  return;
}


void SetZeroIonization(void)
{
  gJHe0 = gJHep = gJH0 = 0;
  epsHe0 = epsHep = epsH0 = 0;
  J_UV = 0;
}


void IonizeParamsFunction(void)
{
  int i, nint;
  double a0, planck, ev, e0_H, e0_He, e0_Hep;
  double gint, eint, t, tinv, fac, eps;
  double at, beta, s;
  double pi;

#define UVALPHA         1.0
  double Jold = -1.0;
  double redshift;

  J_UV = 0.;
  gJHe0 = gJHep = gJH0 = 0.;
  epsHe0 = epsHep = epsH0 = 0.;


  if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
      redshift = 1 / All.Time - 1;

      if(redshift >= 6)
	J_UV = 0.;
      else
	{
	  if(redshift >= 3)
	    J_UV = 4e-22 / (1 + redshift);
	  else
	    {
	      if(redshift >= 2)
		J_UV = 1e-22;
	      else
		J_UV = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);
	    }
	}

      if(J_UV == Jold)
	return;


      Jold = J_UV;

      if(J_UV == 0)
	return;


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
	  fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA + 3.);
	  gint += fac * at;
	  eint += fac * (tinv - 1.) * at;
	}

      gJH0 = a0 * gint / planck;
      epsH0 = a0 * eint * (e0_H / planck);
      gJHep = gJH0 * pow(e0_H / e0_Hep, UVALPHA) / 4.0;
      epsHep = epsH0 * pow((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

      at = 7.83e-18;
      beta = 1.66;
      s = 2.05;

      gJHe0 = (at / planck) * pow((e0_H / e0_He), UVALPHA) *
	(beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
      epsHe0 = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA) *
	(beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

      pi = M_PI;
      gJH0 *= 4. * pi * J_UV;
      gJHep *= 4. * pi * J_UV;
      gJHe0 *= 4. * pi * J_UV;
      epsH0 *= 4. * pi * J_UV;
      epsHep *= 4. * pi * J_UV;
      epsHe0 *= 4. * pi * J_UV;
    }
}





void InitCool(void)
{
    if(ThisTask == 0)
        printf("Initializing cooling ...\n");
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();
    
#ifdef GRACKLE
    InitGrackle();
#endif
    
    InitCoolMemory();
    MakeCoolingTable();
    ReadIonizeParams("TREECOOL");
    IonizeParams();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables();
#endif
}




#ifdef COOL_METAL_LINES_BY_SPECIES
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z)
{
    int k;
    double ne_over_nh_tbl=1, Lambda=0;
    int N_species_active = NUM_METAL_SPECIES-1;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //N_species_active -= 1;
    N_species_active -= NUM_RPROCESS_SPECIES;
#endif
    
    ne_over_nh_tbl = GetLambdaSpecies(0,nHcgs,logT,0);
    for (k=1;k<N_species_active;k++)
        Lambda += GetLambdaSpecies(k,nHcgs,logT,0) * Z[k+1]/(All.SolarAbundances[k+1]*0.0127/All.SolarAbundances[0]);
    
    if(ne_over_nh_tbl>0) Lambda /= ne_over_nh_tbl; else Lambda=0;
    return Lambda;
}

double getSpCoolTableVal(long i,long j,long k,long tblK)
{
    long i_nH=41; long i_T=176; long inHT=i_nH*i_T;
    if(tblK==0) {
        //return SpCoolTable0[i][j][k];
        return SpCoolTable0[i*inHT + j*i_T + k];
    }
    if(tblK==1) {
        //return SpCoolTable1[i][j][k];
        return SpCoolTable1[i*inHT + j*i_T + k];
    }
    return 0;
}

double GetLambdaSpecies(int k_species, double nHcgs, double logT, double fHe)
{
    int ix0,iy0,ix1,iy1,ixmax,iymax;/*,ih0,ih1,ihmax;*/
    double i1,i2,j1,j2,w1,w2,u1;
    double v000,v010,v100,v110,v001,v011,v101,v111;
    double dx,dy,dz,mdz;/*,dh,mdh;*/
    dx=dy=dz=0; ix0=ix1=iy0=iy1=0; ixmax=40; iymax=175; /*ihmax=6;*/
    
    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    if(dx<0) dx=0; if(dx>ixmax) dx=ixmax;
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) ix1=ixmax;
    dx=dx-ix0;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    if(dy<0) dy=0; if(dy>iymax) dy=iymax;
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) iy1=iymax;
    dy=dy-iy0;
    /*
     if(k_species<=1) {
     dh = (fHe-0.238)/(0.298-0.238)*ihmax;
     if(dh<0) dh=0; if(dh>ihmax) dh=ihmax;
     ih0=(int)dh; ih1=ih0+1; if(ih1>ihmax) ih1=ihmax;
     dh=dh-ih0;
     } */
    
    v000=getSpCoolTableVal(k_species,ix0,iy0,0);
    v010=getSpCoolTableVal(k_species,ix0,iy1,0);
    v100=getSpCoolTableVal(k_species,ix1,iy0,0);
    v110=getSpCoolTableVal(k_species,ix1,iy1,0);
    /*if(k_species>1) {
     v000=SpCoolTable0[k_species-2][ix0][iy0];
     v010=SpCoolTable0[k_species-2][ix0][iy1];
     v100=SpCoolTable0[k_species-2][ix1][iy0];
     v110=SpCoolTable0[k_species-2][ix1][iy1];
     } else {
     v000=SpCoolTable0_He[k_species][ix0][iy0][ih0]*mdh + SpCoolTable0_He[k_species][ix0][iy0][ih1]*dh;
     v010=SpCoolTable0_He[k_species][ix0][iy1][ih0]*mdh + SpCoolTable0_He[k_species][ix0][iy1][ih1]*dh;
     v100=SpCoolTable0_He[k_species][ix1][iy0][ih0]*mdh + SpCoolTable0_He[k_species][ix1][iy0][ih1]*dh;
     v110=SpCoolTable0_He[k_species][ix1][iy1][ih0]*mdh + SpCoolTable0_He[k_species][ix1][iy1][ih1]*dh;
     }*/
    if(All.ComovingIntegrationOn) {
        if(All.SpeciesTableInUse<48) {
            dz=log10(1/All.Time)*48; dz=dz-(int)dz;
            v001=getSpCoolTableVal(k_species,ix0,iy0,1);
            v011=getSpCoolTableVal(k_species,ix0,iy1,1);
            v101=getSpCoolTableVal(k_species,ix1,iy0,1);
            v111=getSpCoolTableVal(k_species,ix1,iy1,1);
            /*if(k_species>1) {
             v001=SpCoolTable1[k_species-2][ix0][iy0];
             v011=SpCoolTable1[k_species-2][ix0][iy1];
             v101=SpCoolTable1[k_species-2][ix1][iy0];
             v111=SpCoolTable1[k_species-2][ix1][iy1];
             } else {
             v001=SpCoolTable1_He[k_species][ix0][iy0][ih0]*mdh + SpCoolTable1_He[k_species][ix0][iy0][ih1]*dh;
             v011=SpCoolTable1_He[k_species][ix0][iy1][ih0]*mdh + SpCoolTable1_He[k_species][ix0][iy1][ih1]*dh;
             v101=SpCoolTable1_He[k_species][ix1][iy0][ih0]*mdh + SpCoolTable1_He[k_species][ix1][iy0][ih1]*dh;
             v111=SpCoolTable1_He[k_species][ix1][iy1][ih0]*mdh + SpCoolTable1_He[k_species][ix1][iy1][ih1]*dh;
             }*/
        } else {
            v001=v011=v101=v111=dz=0;
        }
    } else {
        v001=v011=v101=v111=dz=0;
    }
    mdz=1-dz;
    i1=v000*mdz + v001*dz;
    i2=v010*mdz + v011*dz;
    j1=v100*mdz + v101*dz;
    j2=v110*mdz + v111*dz;
    w1=i1*(1-dy) + i2*dy;
    w2=j1*(1-dy) + j2*dy;
    u1=w1*(1-dx) + w2*dx;
    return u1;
}

#endif // COOL_METAL_LINES_BY_SPECIES



#ifdef GALSF_FB_LOCAL_UV_HEATING
void selfshield_local_incident_uv_flux(void)
{
    /* include local self-shielding with the following */
    int i;
    double GradRho,sigma_eff_0;
    
    sigma_eff_0 = 0.955*All.UnitMass_in_g*All.HubbleParam / (All.UnitLength_in_cm*All.UnitLength_in_cm);
    if(All.ComovingIntegrationOn) sigma_eff_0 /= All.Time*All.Time;
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type==0)
        {
            if((SphP[i].RadFluxUV>0) && (PPP[i].Hsml>0) && (SphP[i].Density>0) && (P[i].Mass>0) && (All.Time>0))
            {
                GradRho = sigma_eff_0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1);
                SphP[i].RadFluxUV *= 1276.19 * sigma_eff_0 * exp(-KAPPA_UV*GradRho);
                GradRho *= 3.7e6; // 912 angstrom KAPPA_EUV //
                //SphP[i].RadFluxEUV *= 1276.19 * sigma_eff_0 * exp(-GradRho);
                SphP[i].RadFluxEUV *= 1276.19 * sigma_eff_0 * (0.01 + 0.99/(1.0+0.8*GradRho+0.85*GradRho*GradRho));
                // unit conversion and self-shielding (normalized to Habing (local MW) field)
            } else {
                SphP[i].RadFluxUV = 0;
                SphP[i].RadFluxEUV = 0;
            }}}
}
#endif // GALSF_FB_LOCAL_UV_HEATING





#endif
