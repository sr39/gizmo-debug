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

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc */

/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#if defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_GASRETURN)

/* in case you're wondering, here are some conventions that may be useful for solar abundances
    All.SolarAbundances[0]=0.02;        // all metals (by mass); present photospheric abundances from Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
                                        //   also Anders+Grevesse 1989 (older, but hugely-cited compilation; their Z=0.0201, proto-solar=0.0213)
    // note that the 'all metals' above is the only one where solar enters with any direct role; for everything else, 
    //  'solar' is totally arbitrary (the enrichment routines, etc, don't know what solar is 'supposed' to be): so these
    //  are here purely for convenience and to initialize non-zero metallicities
    //
    All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
    All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
    All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
    All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
    All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
    All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
    All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
    All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
    All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
    All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)
*/

#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif


#define SNeIIBW_Radius_Factor 1.0 // (optional) boost cooling radius for resolution-check


struct kernel_addFB
{
    double dp[3];
    double dv[3];
    double r;
    double wk, dwk;
    double hinv, hinv3, hinv4;
};


struct addFBdata_in
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
  MyFloat SNe_v_ejecta;
  MyDouble Msne;
  MyDouble unit_mom_SNe;
  MyFloat area_sum;
#ifdef METALS
  MyDouble yields[NUM_METAL_SPECIES];
#endif
#ifndef DONOTUSENODELIST
  int NodeList[NODELISTLENGTH];
#endif
}
 *AddFBDataIn, *AddFBDataGet;


struct addFBdata_out
{
  MyFloat M_coupled;
}
 *AddFBDataResult, *AddFBDataOut;



void particle2in_addFB(struct addFBdata_in *in, int i, int feedback_type);
void particle2in_addFB_SNe(struct addFBdata_in *in, int i);
void particle2in_addFB_winds(struct addFBdata_in *in, int i);
void particle2in_addFB_Rprocess(struct addFBdata_in *in, int i);
void particle2in_addFB_wt(struct addFBdata_in *in, int i);
void out2particle_addFB(struct addFBdata_out *out, int i, int mode, int feedback_type);

//
// key variable passed throughout is 'feedback type' -- tells the routine which mode we're in:
//     -1 == evaluate weighting factors for feedback
//      0 == SNe (Type I or II)
//      1 == stellar winds
//      2 == rprocess events (NS-NS mergers)
//
void particle2in_addFB(struct addFBdata_in *in, int i, int feedback_type)
{
    if(feedback_type==-1) particle2in_addFB_wt(in,i);
    if(feedback_type==0) particle2in_addFB_SNe(in,i);
    if(feedback_type==1) particle2in_addFB_winds(in,i);
    if(feedback_type==2) particle2in_addFB_Rprocess(in,i);
}


void particle2in_addFB_Rprocess(struct addFBdata_in *in, int i)
{
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    /*
     k=0     no change to defaults (tmin=3e7yr, rate=1e-5, all neighbor particle)
     k=1     as k=0, 1 ngb particles
     k=2     as k=0, 10 ngb particles
     k=3     as k=0, tmin=3e6yr
     k=4     as k=0, tmin=1e7yr
     k=5     as k=0, tmin=1e8yr
     k=6     as k=0, rate=3e-6
     k=7     as k=0, rate=3e-5
    */
    if(P[i].RProcessEvent_ThisTimeStep<=0)
    {
        in->Msne = 0;
        return;
    }
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    for(k=0;k<NUM_METAL_SPECIES;k++) in->yields[k]=0.0;
    /* this is the important operation for setting the yields for each different model */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    double p = get_random_number(P[i].ID + 8);
    double pcrit,tcrit;
    for(k=0;k<NUM_RPROCESS_SPECIES;k++)
    {
        in->yields[NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES+k] = 0.0;
        tcrit=0.03; // default is age > 3e7
        pcrit=0.3333333333;     // rate lower by 1/3 for 'default'
        if(k==3) {tcrit=0.003;} // k=3 requires age > 3e6 yr
        if(k==4) {tcrit=0.01;}  // k=4 requires age > 1e7 yr
        if(k==5) {tcrit=0.1;}   // k=5 requires age > 1e8
        if(k==6) {pcrit=0.1;}   // k=6 has rate lower by 3
        if(k==7) {pcrit=1.0;}   // k=7 has rate higher by 3
        if((star_age>=tcrit)&&(p<=pcrit)&&(P[i].RProcessEvent_ThisTimeStep>0))
        {
            in->yields[NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES+k] = 1.0; // absolute unit is irrelevant, so use 1.0 //
        }
    }
    in->Hsml = PPP[i].Hsml;
    in->Msne = 0.01 * (double)P[i].RProcessEvent_ThisTimeStep / ((double)((All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS)); // mass ejected ~0.01*M_sun; only here for bookkeeping //
    in->unit_mom_SNe = 0;
    in->SNe_v_ejecta = 0.;
    in->area_sum = P[i].Area_weighted_sum;
#endif
}


void particle2in_addFB_wt(struct addFBdata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->Hsml = PPP[i].Hsml;
    in->Msne = P[i].Mass;
    in->unit_mom_SNe = 1;
    in->SNe_v_ejecta = 500.;
    in->area_sum = P[i].Area_weighted_sum;
#ifdef GALSF_TURNOFF_COOLING_WINDS
    /* calculate the 'blast radius' and 'cooling turnoff time' used by this model */
    double n0 = P[i].DensAroundStar*All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
    n0 = SNeIIBW_Radius_Factor * 0.087*pow(n0,-0.36) / (All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime);
    in->Hsml = DMAX(PPP[i].Hsml/2.,DMIN(n0,5.*PPP[i].Hsml));
#endif
}



void particle2in_addFB_SNe(struct addFBdata_in *in, int i)
{
    int k;
    double unitmass_in_msun,agemax,n_sn_0,star_age;
    double Msne,Esne51,unit_egy_SNe,unit_mom_SNe,SNe_v_ejecta;
#ifdef METALS
    double yields[NUM_METAL_SPECIES];
#endif
    
    agemax=0.03753; /* in Gyr */
    unitmass_in_msun=(All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS;
    
    n_sn_0 = P[i].SNe_ThisTimeStep;
    if((n_sn_0<=0)||(P[i].DensAroundStar<=0))
    {
        in->Msne = 0;
        return;
    }
    
    star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    Esne51 = All.SNeIIEnergyFrac*n_sn_0;
    unit_egy_SNe = 1.0e51/(All.UnitEnergy_in_cgs/All.HubbleParam);
    Msne = 10.5; //average ejecta mass for single event (normalized to give total mass loss correctly)
    if(star_age > agemax) Msne=1.4; // SnIa
    
#ifdef METALS
    if(NUM_METAL_SPECIES>=10) {
        // All, then He,C,N,O,Ne,Mg,Si,S,Ca,Fe
        if(star_age > agemax) {
            /* SNIa */ /* from Iwamoto et al. 1999; 'W7' models */
            yields[0]=1.4;/* total metal mass */
            yields[1]=0.0;/*He*/ yields[2]=0.049;/*C*/ yields[3]=1.2e-6;/*N*/ yields[4]=0.143;/*O*/
            yields[5]=0.0045;/*Ne*/ yields[6]=0.0086;/*Mg*/ yields[7]=0.156;/*Si*/
            yields[8]=0.087;/*S*/ yields[9]=0.012;/*Ca*/ yields[10]=0.743;/*Fe*/
        } else {
            /* SNII (IMF-averaged... may not be the best approx on short timescales...) */
            // Woosley & Weaver 1995 'B' models, integrated (log-log interp) over Salpeter IMF, from 10-50 rescaled from 0.1-100 solar //
            //yields[0]=1.70;/*Z*/
            //yields[1]=4.03;/*He*/ yields[2]=0.117;/*C*/ yields[3]=0.0399;/*N*/ yields[4]=1.06;/*O*/
            //yields[5]=0.169;/*Ne*/ yields[6]=0.0596;/*Mg*/ yields[7]=0.0924;/*Si*/
            //yields[8]=0.0408;/*S*/ yields[9]=0.00492;/*Ca*/ yields[10]=0.0842;/*Fe*/
            /* note that the Mg yield here, and perhaps some other products, is systematically low compared to observations (by ~0.4 dex) */
            // Nomoto 2006 (arXiv:0605725) //
            yields[0]=2.0;/*Z*/
            yields[1]=3.87;/*He*/ yields[2]=0.133;/*C*/ yields[3]=0.0479;/*N*/ yields[4]=1.17;/*O*/
            yields[5]=0.30;/*Ne*/ yields[6]=0.0987;/*Mg*/ yields[7]=0.0933;/*Si*/
            yields[8]=0.0397;/*S*/ yields[9]=0.00458;/*Ca*/ yields[10]=0.0741;/*Fe*/
            
            if(P[i].Metallicity[0]<0.033)
            {
                yields[3] *= P[i].Metallicity[0]/All.SolarAbundances[0]; // N scaling is strongly dependent on initial metallicity of the star //
            } else {
                yields[3] *= 1.65;
            }
            yields[0] += yields[3]-0.0479; // correct total metal mass for this correction //
        }
    }
    if(NUM_METAL_SPECIES==3 || NUM_METAL_SPECIES==4)
    {
    if(star_age > agemax) {
        yields[0]=1.4; yields[1]=0.0086; yields[2]=0.743; // All Z, Mg, Fe in total mass (SnIa)
    } else {
        yields[0]=2.0; yields[1]=0.12; yields[2]=0.0741; // SnII (per-SNe IMF-weighted averages)
    }
    }
    if(NUM_METAL_SPECIES==1) {if(star_age > agemax) {yields[0]=1.4;} else {yields[0]=2.0;}}
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    for(k=1;k<=NUM_RPROCESS_SPECIES;k++)
        yields[NUM_METAL_SPECIES-k] = 0.0; // R-process tracker
#endif
    for(k=0;k<NUM_METAL_SPECIES;k++) yields[k]=yields[k]/Msne;
    /* this is almost always irrelevant, but we can add a check to allow for larger abundances in the progenitor stars */
    for(k=0;k<NUM_METAL_SPECIES;k++)
    {
        yields[k]=yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]);
    }
    if(star_age > agemax) {if(NUM_METAL_SPECIES>=10) {yields[1]=0.0;}} // no He yield for Ia SNe // 
    for(k=0;k<NUM_METAL_SPECIES;k++) {if(yields[k]<0) yields[k]=0.0; if(yields[k]>1) yields[k]=1; in->yields[k]=yields[k];}
#endif
    
    Msne *= n_sn_0/unitmass_in_msun;
    SNe_v_ejecta = sqrt(2.0*Esne51*unit_egy_SNe/Msne); // v_ej
    unit_mom_SNe = Msne * SNe_v_ejecta; // total SNe ejecta momentum = M_ej*v_ej
    
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->Hsml = PPP[i].Hsml;
    in->Msne = Msne;
    in->SNe_v_ejecta = SNe_v_ejecta;
    in->unit_mom_SNe = unit_mom_SNe;
    in->area_sum = P[i].Area_weighted_sum;
#ifdef GALSF_TURNOFF_COOLING_WINDS
    /* calculate the 'blast radius' and 'cooling turnoff time' used by this model */
    double n0 = P[i].DensAroundStar*All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
    n0 = SNeIIBW_Radius_Factor * 0.087*pow(n0,-0.36) / (All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime);
    in->Hsml = DMAX(PPP[i].Hsml/2.,DMIN(n0,5.*PPP[i].Hsml));
    in->unit_mom_SNe = 0;
#endif
}



void particle2in_addFB_winds(struct addFBdata_in *in, int i)
{
#ifdef GALSF_FB_GASRETURN
    int k;
    double star_age,T_corr,GasSpecEnergy,M_wind,E_wind,wind_velocity,wind_momentum;
#ifdef METALS
    double yields[NUM_METAL_SPECIES];
#endif
    
    if((P[i].MassReturn_ThisTimeStep<=0)||(P[i].DensAroundStar<=0))
    {
        in->Msne = 0;
        return;
    }
    
    GasSpecEnergy = All.AGBGasEnergy * 3.0e7*(1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS) * All.UnitMass_in_g/All.UnitEnergy_in_cgs;
    
    star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    /* wind internal energy (hot winds from O-stars, slow from AGB winds) */
    if(star_age <= 0.1)
        T_corr=0.0013 + 16.0/(1+pow(star_age/0.0025,1.4)+pow(star_age/0.01,5.0));
    else
        T_corr=0.0013;
    
    /* define the metal yields */
#ifdef METALS
    /* assume track initial metallicity; turn on COOL_METAL_LINES_BY_SPECIES for more detailed tracking of light elements */
    for(k=0;k<NUM_METAL_SPECIES;k++) yields[k]=P[i].Metallicity[k];
    if(NUM_METAL_SPECIES>=10)
    {
        /* All, then He,C,N,O,Ne,Mg,Si,S,Ca,Fe ;; follow AGB/O star yields in more detail for the light elements */
        /*   the interesting species are He & CNO */
        // below is based on a compilation of van den Hoek & Groenewegen 1997, Marigo 2001, Izzard 2004 //
        yields[1]=0.36; /*He*/ yields[2]=0.016; /*C*/ yields[3]=0.0041; /*N*/ yields[4]=0.0118; /*O*/
        if(P[i].Metallicity[0]<0.033)
        {
            yields[4] *= P[i].Metallicity[0]/All.SolarAbundances[0]; // O scaling is strongly dependent on initial metallicity of the star //
        } else {
            yields[4] *= 1.65;
        }
        for(k=1;k<=4;k++)
        {
            yields[k] = yields[k]*(1.-P[i].Metallicity[0]) + (P[i].Metallicity[k]-All.SolarAbundances[k]);
            if(yields[k]<0) yields[k]=0;
            if(yields[k]>1) yields[k]=1;
        }
        yields[0]=0.0; for(k=2;k<=10;k++) yields[0]+=yields[k];
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
        for(k=1;k<=NUM_RPROCESS_SPECIES;k++)
            yields[NUM_METAL_SPECIES-k] = 0.0; // R-process tracker
#endif
    } else {
        for(k=0;k<NUM_METAL_SPECIES;k++) {yields[k]=0.0;}
        yields[0]=0.032;
    }
    for(k=0;k<NUM_METAL_SPECIES;k++) in->yields[k]=yields[k];
#endif
    
    M_wind = P[i].Mass * P[i].MassReturn_ThisTimeStep;
    E_wind = GasSpecEnergy*T_corr*M_wind;
    wind_velocity = sqrt(2.0*E_wind/M_wind);
    wind_momentum = M_wind*wind_velocity;
    
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->Hsml = PPP[i].Hsml;
    in->Msne = M_wind;
    in->SNe_v_ejecta = wind_velocity;
    in->unit_mom_SNe = wind_momentum;
    in->area_sum = P[i].Area_weighted_sum;
#endif // GALSF_FB_GASRETURN //
}



void out2particle_addFB(struct addFBdata_out *out, int i, int mode, int feedback_type)
{
    if(feedback_type==-1)
    {
        ASSIGN_ADD(P[i].Area_weighted_sum, out->M_coupled, mode);
    } else {
        P[i].Mass -= out->M_coupled;
        if(P[i].Mass<0) P[i].Mass=0;
    }
}





void mechanical_fb_calc(int feedback_type)
{
  int j, k, ngrp, ndone, ndone_flag;
  int recvTask, place;
  int save_NextParticle;
  long long n_exported = 0;

  /* allocate buffers to arrange communication */
  long long NTaskTimesNumPart;
  NTaskTimesNumPart = maxThreads * NumPart;
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct addFBdata_in) +
					     sizeof(struct addFBdata_out) +
					     sizemax(sizeof(struct addFBdata_in),
						     sizeof(struct addFBdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  NextParticle = FirstActiveParticle;	/* begin with this index */
  do
    {

      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
#ifdef OMP_NUM_THREADS
      pthread_t mythreads[OMP_NUM_THREADS - 1];
      int threadid[OMP_NUM_THREADS - 1];
      pthread_attr_t attr;

      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      pthread_mutex_init(&mutex_nexport, NULL);
      pthread_mutex_init(&mutex_partnodedrift, NULL);

      TimerFlag = 0;

      for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	{
	  threadid[j] = j + 1;
	  pthread_create(&mythreads[j], &attr, addFB_evaluate_primary, &threadid[j]);
	}
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int mainthreadid = omp_get_thread_num();
#else
	int mainthreadid = 0;
#endif
	addFB_evaluate_primary(&mainthreadid, feedback_type);	/* do local particles and prepare export list */
      }

#ifdef OMP_NUM_THREADS
      for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	pthread_join(mythreads[j], NULL);
#endif


      if(BufferFullFlag)
	{
	  int last_nextparticle = NextParticle;

	  NextParticle = save_NextParticle;

	  while(NextParticle >= 0)
	    {
	      if(NextParticle == last_nextparticle)
		break;

	      if(ProcessedFlag[NextParticle] != 1)
		break;

	      ProcessedFlag[NextParticle] = 2;

	      NextParticle = NextActiveParticle[NextParticle];
	    }

	  if(NextParticle == save_NextParticle)
	    {
	      /* in this case, the buffer is too small to process even a single particle */
	      endrun(116608);
	    }

	  int new_export = 0;

	  for(j = 0, k = 0; j < Nexport; j++)
	    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
	      {
		if(k < j + 1)
		  k = j + 1;

		for(; k < Nexport; k++)
		  if(ProcessedFlag[DataIndexTable[k].Index] == 2)
		    {
		      int old_index = DataIndexTable[j].Index;

		      DataIndexTable[j] = DataIndexTable[k];
		      DataNodeList[j] = DataNodeList[k];
		      DataIndexTable[j].IndexGet = j;
		      new_export++;

		      DataIndexTable[k].Index = old_index;
		      k++;
		      break;
		    }
	      }
	    else
	      new_export++;

	  Nexport = new_export;

	}

      n_exported += Nexport;

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;

      MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      AddFBDataGet = (struct addFBdata_in *) mymalloc("AddFBDataGet", Nimport * sizeof(struct addFBdata_in));
      AddFBDataIn = (struct addFBdata_in *) mymalloc("AddFBDataIn", Nexport * sizeof(struct addFBdata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  particle2in_addFB(&AddFBDataIn[j], place, feedback_type);
#ifndef DONOTUSENODELIST
	  memcpy(AddFBDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif

	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&AddFBDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct addFBdata_in), MPI_BYTE,
			       recvTask, TAG_FBLOOP_A,
			       &AddFBDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct addFBdata_in), MPI_BYTE,
			       recvTask, TAG_FBLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(AddFBDataIn);
      AddFBDataResult =
	(struct addFBdata_out *) mymalloc("AddFBDataResult", Nimport * sizeof(struct addFBdata_out));
      AddFBDataOut =
	(struct addFBdata_out *) mymalloc("AddFBDataOut", Nexport * sizeof(struct addFBdata_out));

      /* now do the particles that were sent to us */
      NextJ = 0;
#ifdef OMP_NUM_THREADS
      for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	pthread_create(&mythreads[j], &attr, addFB_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int mainthreadid = omp_get_thread_num();
#else
	int mainthreadid = 0;
#endif
	addFB_evaluate_secondary(&mainthreadid, feedback_type);
      }

#ifdef OMP_NUM_THREADS
      for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	pthread_join(mythreads[j], NULL);

      pthread_mutex_destroy(&mutex_partnodedrift);
      pthread_mutex_destroy(&mutex_nexport);
      pthread_attr_destroy(&attr);
#endif

      if(NextParticle < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&AddFBDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct addFBdata_out),
			       MPI_BYTE, recvTask, TAG_FBLOOP_B,
			       &AddFBDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct addFBdata_out),
			       MPI_BYTE, recvTask, TAG_FBLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  out2particle_addFB(&AddFBDataOut[j], place, 1, feedback_type);
	}
      myfree(AddFBDataOut);
      myfree(AddFBDataResult);
      myfree(AddFBDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  /* do final operations on results */
  /* (not needed here, since the entire operation is from central particle to the gas) */
/*
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
      }
*/
}






int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		   int *ngblist, int feedback_type)
{
    int startnode, numngb_inbox, listindex = 0;
    int j, k, n;
    double u,r2,h2;
    double v_ejecta_max,kernel_zero,wk,dM,dP,dE;
    double E_coupled,wk_sum,dP_sum,dP_boost_sum;

    struct kernel_addFB kernel;
    struct addFBdata_in local;
    struct addFBdata_out out;
    memset(&out, 0, sizeof(struct addFBdata_out));

    v_ejecta_max = 5000.0 * 1.0e5/ All.UnitVelocity_in_cm_per_s;
    // 'speed limit' to prevent numerically problematic kicks at low resolution //
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1);

    /* Load the data for the particle injecting feedback */
    if(mode == 0)
        particle2in_addFB(&local, target, feedback_type);
    else
        local = AddFBDataGet[target];

    if(local.Msne<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

    // some units (just used below, but handy to define for clarity) //
    double unitlength_in_kpc=All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime;
    double density_to_n=All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
    double unit_egy_SNe = 1.0e51/(All.UnitEnergy_in_cgs/All.HubbleParam);
#ifdef GALSF_TURNOFF_COOLING_WINDS
    double pressure_to_p4 = (1/All.cf_afac1)*density_to_n*(All.UnitEnergy_in_cgs/All.UnitMass_in_g) / 1.0e4;
#endif
    // now define quantities that will be used below //
    double Esne51 = 0.5*local.SNe_v_ejecta*local.SNe_v_ejecta*local.Msne / unit_egy_SNe;
    double r2sne, RsneKPC, RsneKPC_0, RsneMAX;
    r2sne=0; RsneKPC=0.; RsneMAX=local.Hsml;
    RsneKPC_0=(0.0284/unitlength_in_kpc)*pow(Esne51,0.286); //Cioffi: weak external pressure



  /* Now start the actual FB computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = AddFBDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
        numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
        
	  if(numngb_inbox < 0)
	    return -1;

      E_coupled = wk_sum = dP_sum = dP_boost_sum = 0;
	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = ngblist[n];
            if(P[j].Type != 0) continue; // require a gas particle //
            if(P[j].Mass <= 0) continue; // require the particle has mass //
            
            for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef PERIODIC			/* find the closest image in the given box size  */
            NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
            r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
            if(r2<=0) continue; // same particle //
            
            double h2j = PPP[j].Hsml * PPP[j].Hsml;
            if((r2>h2)&&(r2>h2j)) continue; // outside kernel (in both 'directions') //
            
            // calculate kernel quantities //
            kernel.r = sqrt(r2);
            if(kernel.r > DMAX(2.0/unitlength_in_kpc,PPP[j].Hsml)) continue; // no super-long-range effects allowed! (of course this is arbitrary in code units) //
            
            //u = kernel.r * kernel.hinv;
            //kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, -1);
            for(k=0; k<3; k++) kernel.dv[k] = local.Vel[k] - P[j].Vel[k];
            
            /*
            wk = 1./SphP[j].Density; // wt ~ 1 (uniform in SPH terms)
            wk = kernel.wk * P[j].Mass / SphP[j].Density; // psi
            */
            double h_eff_j = Get_Particle_Size(j);
            //wk = h_eff_j * h_eff_j / (r2 + 0.01*h2); // solid-angle weight (actually, because of summation/division below, this really double-downweights further particles)
            wk = h_eff_j * h_eff_j; // area (solid-angle after summation/division below) weight
            //wk = h_eff_j * h_eff_j * h_eff_j; // volume weight
            
            // if feedback_type==-1, this is a pre-calc loop to get the relevant weights for coupling //
            if(feedback_type==-1)
            {
                out.M_coupled += wk;
                continue;
            }
            // NOW do the actual feedback calculation //
                wk /= local.area_sum; // this way wk matches the value summed above for the weighting //
                // need to check to make sure the coupled fraction doesn't exceed the solid angle subtended by the particles //
                //double wkmax = 1.5 * M_PI * h_eff_j * h_eff_j / (4. * M_PI * (0.5625*r2 + 0.005*h2)); if(wk > wkmax) {wk = wkmax;}
            
                dM = wk * local.Msne;
                dP = local.SNe_v_ejecta / kernel.r;
                /* define initial mass and ejecta velocity in this 'cone' */

#ifndef GALSF_TURNOFF_COOLING_WINDS
            //RsneKPC=RsneKPC_0 * pow(SphP[j].Density*density_to_n+1.0e-3,-0.429);
            //RsneKPC=RsneKPC_0 / sqrt(SphP[j].Density*density_to_n+1.0e-3);
            RsneKPC = RsneKPC_0;
            double n0 = SphP[j].Density*density_to_n;
            /* this is tedious, but is a fast approximation (essentially a lookup table) for the -0.429 power above */
            if(n0 < 1.e-3) {RsneKPC *= 19.4;} else {
                if(n0 < 1.e-2) {RsneKPC *= 1.9 + 23./(1.+333.*n0);} else {
                    if(n0 < 1.e-1) {RsneKPC *= 0.7 + 8.4/(1.+33.3*n0);} else {
                        if(n0 < 1) {RsneKPC *= 0.08 + 3.1/(1.+2.5*n0);} else {
                            if(n0 < 10) {RsneKPC *= 0.1 + 1.14/(1.+0.333*n0);} else {
                                if(n0 < 100) {RsneKPC *= 0.035 + 0.43/(1.+0.0333*n0);} else {
                                    if(n0 < 1000) {RsneKPC *= 0.017 + 0.154/(1.+0.00333*n0);} else {
                                        if(n0 < 1.e4) {RsneKPC *= 0.006 + 0.057/(1.+0.000333*n0);} else {
                                            RsneKPC *= pow(n0, -0.429); }}}}}}}}
            
                /*
                if(P[j].Metallicity[0]/All.SolarAbundances[0] < 0.01) {RsneKPC*=2.0;} else {
                    if(P[j].Metallicity[0]<All.SolarAbundances[0]) {RsneKPC*=pow(P[j].Metallicity[0]/All.SolarAbundances[0],-0.15);} else {RsneKPC*=pow(P[j].Metallicity[0]/All.SolarAbundances[0],-0.09);}}
                */
                /* below expression is again just as good a fit to the simulations, and much faster to evaluate */
                double z0 = P[j].Metallicity[0]/All.SolarAbundances[0];
                if(z0 < 0.01)
                {
                    RsneKPC *= 2.0;
                } else {
                    if(z0 < 1)
                    {
                        RsneKPC *= 0.93 + 0.0615 / (0.05 + 0.8*z0);
                    } else {
                        RsneKPC *= 0.8 + 0.4 / (1 + z0);
                    }
                }
                /* calculates cooling radius given density and metallicity in this annulus into which the ejecta propagate */
            
                if(RsneMAX<RsneKPC) RsneKPC=RsneMAX;
                /* limit to Hsml for coupling */

                r2sne = RsneKPC*RsneKPC;
                // if(r2 > r2sne) dP *= pow(r2sne/r2 , 1.625);
                if(r2 > r2sne) dP *= r2sne*RsneKPC / (r2*kernel.r); // just as good a fit, and much faster to evaluate //
                /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium:
                    from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool */
#endif
            
                /* now, add contribution from relative star-gas particle motion to shock energy */
                u = 0.; dE = 0.;
                for(k=0; k<3; k++)
                {
                    // relative outflow-particle velocity = v_ej*x_i/r + v_star - v_gasparticle //
                    u = -dP*kernel.dp[k] + kernel.dv[k]/All.cf_atime;
                    dE += u*u;
                }
                dE *= 0.5 * dM;
                /* now we have the proper energy to couple */
                
                E_coupled += dE;
                wk_sum += 1.;
                out.M_coupled += dM;
            
                /* inject actual mass from mass return */
                if(P[j].Hsml<=0) {if(SphP[j].Density>0){SphP[j].Density*=(1+dM/P[j].Mass);} else {SphP[j].Density=dM*kernel.hinv3;}} else {
                    SphP[j].Density+=kernel_zero*dM/(P[j].Hsml*P[j].Hsml*P[j].Hsml);}
                SphP[j].Density *= 1 + dM/P[j].Mass; // inject mass at constant particle volume //
                P[j].Mass += dM;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[j].MassTrue += dM;
#endif
                /* inject metals */
#ifdef METALS
                u=dM/P[j].Mass;
                if(u>1) u=1;
                for(k=0;k<NUM_METAL_SPECIES;k++)
                {
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
                    u=dM/P[j].Mass;
                    if(u>1) u=1;
                    if(k>=NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES)
                    {
                        if(local.yields[k]>0)
                        {
                            /* for k=1, couple to only 1 neighbor */
                            if(k==NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES+1)
                            {
                                if((mode==0)&&(wk_sum<=1))
                                {
                                    u=local.Msne/P[j].Mass;
                                } else {
                                    u=0;
                                }
                            }
                            /* for k=2, couple to only 10 neighbors */
                            if(k==NUM_METAL_SPECIES-NUM_RPROCESS_SPECIES+2)
                            {
                                if((mode==0)&&(wk_sum<=10))
                                {
                                    u=0.1*local.Msne/P[j].Mass;
                                } else {
                                    u=0;
                                }
                            }
                        }
                    }
#endif
                    P[j].Metallicity[k]=(1-u)*P[j].Metallicity[k]+u*local.yields[k];
                }
#endif
            
#if defined(COSMIC_RAYS) && defined(GALSF_FB_SNE_HEATING)
            if(local.SNe_v_ejecta > 5.0e7 / All.UnitVelocity_in_cm_per_s)
            {
                /* a fraction of the *INITIAL* energy goes into cosmic rays [this is -not- affected by the radiative losses above] */
                double dE_init_coupled = 0.5 * dM * local.SNe_v_ejecta * local.SNe_v_ejecta;
                SphP[j].CosmicRayEnergy += All.CosmicRay_SNeFraction * dE_init_coupled;
                SphP[j].CosmicRayEnergyPred += All.CosmicRay_SNeFraction * dE_init_coupled;
            }
#endif
                /* inject the post-shock energy and momentum (convert to specific units as needed first) */
                dE *= 1 / P[j].Mass;
                SphP[j].InternalEnergy += dE;
                SphP[j].InternalEnergyPred += dE;
#ifdef GALSF_TURNOFF_COOLING_WINDS
                /* if the sub-grid 'cooling turnoff' model is enabled, turn off cooling for the 'blastwave timescale' */
                dP = 7.08 * pow(Esne51*SphP[j].Density*density_to_n,0.34) * pow(SphP[j].Pressure*pressure_to_p4,-0.70)
                        / (All.UnitTime_in_Megayears/All.HubbleParam);
                if(dP>SphP[j].DelayTimeCoolingSNe) SphP[j].DelayTimeCoolingSNe=dP;
#else
                /* inject momentum */
                dP = wk * local.unit_mom_SNe / P[j].Mass;
                dP_sum += dP;
                dP *= sqrt(1. + NORM_COEFF*(SphP[j].Density*RsneKPC*RsneKPC*RsneKPC)/local.Msne);
                /* above is the appropriate factor for the ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                if(dP > v_ejecta_max) dP = v_ejecta_max;
                dP_boost_sum += dP;
                dP *= All.cf_atime / kernel.r;
                for(k=0; k<3; k++)
                {
                    u = wk * local.Msne * kernel.dv[k] / P[j].Mass;
                    if (u > v_ejecta_max*All.cf_atime) u = v_ejecta_max*All.cf_atime;
                    if (u < -v_ejecta_max*All.cf_atime) u = -v_ejecta_max*All.cf_atime;
                    P[j].Vel[k] += -dP*kernel.dp[k] + u;
                    SphP[j].VelPred[k] += -dP*kernel.dp[k] + u;
                }
#endif
            
#ifdef PM_HIRES_REGION_CLIPPING
                dP=0; for(k=0;k<3;k++) dP+=P[j].Vel[k]*P[j].Vel[k]; dP=sqrt(dP);
                if(dP>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) P[j].Mass=0;
                if(dP>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) for(k=0;k<3;k++) P[j].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/dP;
#endif

/* // this is handled for us now in the hydro routine //
#ifdef WAKEUP
                SphP[j].wakeup = 1; // wakeup particle after feedback injection
#endif
*/
	    } // for(n = 0; n < numngb; n++)
        // output some of these; for convenience, we will only output the local mean values, but normalized appropriately //
/*
        if((mode == 0)&&(feedback_type>=0))
        {
            if(dP_sum>0) dP_boost_sum /= dP_sum;
            printf("SNe/Wind/R-Process Feedback: Time=%g FB_Type=%d M_Ejecta=%g v0_Ejecta=%g Hsml_coupling=%g boost_mean=%g \n",
                   All.Time,feedback_type,local.Msne,local.SNe_v_ejecta,local.Hsml,dP_boost_sum); fflush(stdout);
 
            //fprintf(FdGasReturn, "%lg %d %g %g %g %g \n",
            //       All.Time,feedback_type,local.Msne,local.SNe_v_ejecta,local.Hsml,dP_boost_sum); fflush(FdGasReturn);
            // problem: this file can only be accessed from the head node (ThisTask==0); need to MPI share to it to print //
        }
*/
	} // while(startnode >= 0)
        

#ifndef DONOTUSENODELIST
    if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = AddFBDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	} // if(mode == 1)
#endif
    } // while(startnode >= 0)

  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_addFB(&out, target, 0, feedback_type);
  else
    AddFBDataResult[target] = out;

  return 0;
} // int addFB_evaluate 






void *addFB_evaluate_primary(void *p, int feedback_type)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;
  int active_check = 0;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      UNLOCK_NEXPORT;
      if(exitFlag)
	break;
        
    active_check = 0;
    if(PPP[i].NumNgb > 0 && PPP[i].Hsml > 0 && P[i].Mass > 0)
    {
#ifdef GALSF_FB_SNE_HEATING
        if(P[i].SNe_ThisTimeStep>0)
        if(feedback_type==-1 || feedback_type==0)
        active_check = 1;
#endif
#ifdef GALSF_FB_GASRETURN
        if(P[i].MassReturn_ThisTimeStep>0)
        if(feedback_type==-1 || feedback_type==1)
        active_check = 1;
#endif
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
        if(P[i].RProcessEvent_ThisTimeStep>0)
        if(feedback_type==-1 || feedback_type==2)
        active_check = 1;
#endif
    }
        
    if(active_check==1)
    {
        if(addFB_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist, feedback_type) < 0)
            break;		// export buffer has filled up //
    }

      ProcessedFlag[i] = 1; /* particle successfully finished */
    }

  return NULL;
}




void *addFB_evaluate_secondary(void *p, int feedback_type)
{
  int thread_id = *(int *) p;
  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;

  while(1)
    {
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }
      UNLOCK_NEXPORT;

      if(j >= Nimport)
	break;

      addFB_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist, feedback_type);
    }

  return NULL;
}





void determine_where_SNe_occur()
{
    int i;
    double dt,star_age,agemin,agebrk,agemax,RSNe,p,n_sn_0,RSNeFac;
    double npossible,nhosttotal,ntotal,ptotal,dtmean,rmean;
    npossible=nhosttotal=ntotal=ptotal=dtmean=rmean=0;
    double mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean;
    mpi_npossible=mpi_nhosttotal=mpi_ntotal=mpi_ptotal=mpi_dtmean=mpi_rmean=0;
#ifdef GALSF_FB_GASRETURN
    double D_RETURN_FRAC = 0.01; // fraction of particle mass to return on a recycling step //
#endif
    
    if(All.Time<=0) return;
    
    // basic variables we will use //
    agemin=0.003401; agebrk=0.01037; agemax=0.03753; // in Gyr //
    // converts rate to code units //
    RSNeFac=(All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS;
    
    // loop over particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].SNe_ThisTimeStep=0;
#ifdef GALSF_FB_GASRETURN
        P[i].MassReturn_ThisTimeStep=0;
#endif
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
        P[i].RProcessEvent_ThisTimeStep=0;
#endif
        if(All.ComovingIntegrationOn) if(P[i].Type != 4) continue;
        if(All.ComovingIntegrationOn==0) if((P[i].Type<2)||(P[i].Type>4)) continue;
        if(P[i].Mass<=0) continue;
        
        // get particle timestep //
#ifndef WAKEUP
        dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; // dloga to dt_physical
#else
        dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a; //
#endif
        if(dt<=0) continue;
        
        // now use the simple tabulated SNe rates //
        RSNe=0.0;
        star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age<=0) continue;
        npossible++;
        
        // SNe component //
        if(All.SNeIIEnergyFrac>0)
        {
            if(star_age > agemin)
            {
                if((star_age>=agemin)&&(star_age<=agebrk))
                    RSNe = 5.408e-4; // NSNe/Myr *if* each SNe had exactly 10^51 ergs; really from the energy curve //
                if((star_age>=agebrk)&&(star_age<=agemax))
                    RSNe = 2.516e-4; // this is for a 1 Msun population //
                // add contribution from Type-Ia //
                if(star_age>agemax)
                    RSNe = 5.3e-8 + 1.6e-5*exp(-0.5*((star_age-0.05)/0.01)*((star_age-0.05)/0.01));
                // delayed population (constant rate)  +  prompt population (gaussian) //
                p = dt * (RSNe*RSNeFac) * P[i].Mass * calculate_relative_light_to_mass_ratio_from_imf(i);
                ptotal += p;
                n_sn_0 = (float)floor(p);
                p -= n_sn_0;
                if(get_random_number(P[i].ID + 6) < p) n_sn_0++;
                P[i].SNe_ThisTimeStep = n_sn_0;
                
                rmean += RSNe;
                ntotal += n_sn_0;
                if(n_sn_0>0) nhosttotal++;
            } // if(star_age > agemin) //
        }
        dtmean += dt;
        
#ifdef GALSF_FB_GASRETURN
        // Stellar Winds component //
        if(All.GasReturnFraction>0)
        {
            p=0.0;
            if(star_age<=0.001){p=11.6846;} else {
                if(star_age<=0.0035){p=11.6846*(0.01+P[i].Metallicity[0]/All.SolarAbundances[0])*
                    pow(10.,1.838*(0.79+log10(P[i].Metallicity[0]/All.SolarAbundances[0]))*(log10(star_age)-(-3.00)));} else {
                        if(star_age<=0.1){p=72.1215*pow(star_age/0.0035,-3.25)+0.0103;} else {
                            p=1.03*pow(star_age,-1.1)/(12.9-log(star_age));
                        }}}
            p *= 1.4 * 0.291175; // to give expected return fraction from stellar winds alone (~17%) //
            p *= All.GasReturnFraction * (dt*0.001*All.UnitTime_in_Megayears/All.HubbleParam);
            p *= calculate_relative_light_to_mass_ratio_from_imf(i);
            // this is the fraction of the particle mass expected to return in the timestep //
            p = 1.0 - exp(-p); // need to account for p>1 cases //
            if(get_random_number(P[i].ID + 5) < p/D_RETURN_FRAC) // ok, have a mass return event //
                P[i].MassReturn_ThisTimeStep = D_RETURN_FRAC;
        }
#endif
        
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
        p=0.0;
        /* we'll use the maximum rate here, then in the -yields- setting, 'cut' each down to its sub-population */
        if(star_age>=0.003) // rate is zero at <3e6 yr
        {
            p = 3.0e-5 / (1000.*star_age); // Nsne/Myr for a 1 Msun population
            p *= dt * RSNeFac * P[i].Mass * calculate_relative_light_to_mass_ratio_from_imf(i);
            n_sn_0 = (float)floor(p);
            p -= n_sn_0;
            if(get_random_number(P[i].ID + 7) < p) n_sn_0++;
            P[i].RProcessEvent_ThisTimeStep = n_sn_0;
        }
#endif
        
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
    
    
    MPI_Reduce(&dtmean, &mpi_dtmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rmean, &mpi_rmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ptotal, &mpi_ptotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nhosttotal, &mpi_nhosttotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ntotal, &mpi_ntotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&npossible, &mpi_npossible, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(mpi_npossible>0)
        {
            mpi_dtmean /= mpi_npossible;
            mpi_rmean /= mpi_npossible;
            fprintf(FdSneIIHeating, "%lg %g %g %g %g %g %g \n",
                    All.Time,mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean);
            fflush(FdSneIIHeating);
        }
    } // if(ThisTask == 0) //
    
} // void determine_where_SNe_occur() //



#ifdef GALSF_GASOLINE_RADHEATING
/* this routine copies the 'heating' term from young stars included in the Stinson+ 2013 GASOLINE model:
 the integrated luminosity from all stars with age<4 Myr is coupled directly as a heating term to the
 gas, with some efficiency parameter */
void luminosity_heating_gasoline(void)
{
    double Gasoline_LumHeating_Efficiency = 0.1;
    // coupling 'efficiency' (currently hard-coded to 10% to match their model, can change) //
    double dt,star_age,dE,dE_j,wtsum,h;
    int startnode, numngb, i, j, n, dummy=0;
    
    if(All.Time<=0) return;
    Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3))))
        {
#ifndef WAKEUP
            dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
            if((star_age < 0.004)&&(dt>0)&&(P[i].Mass>0))
            {
                dE = evaluate_l_over_m_ssp(star_age) * calculate_relative_light_to_mass_ratio_from_imf(i);
                dE *= Gasoline_LumHeating_Efficiency;
                dE *= (4.0/2.0) * (P[i].Mass*All.UnitMass_in_g/All.HubbleParam); // L in CGS
                dE *= (dt*All.UnitTime_in_s/All.HubbleParam) / (All.UnitEnergy_in_cgs/All.HubbleParam); // dE in code units
                
                h = 1.0*P[i].Hsml;
                numngb=ngb_treefind_variable_threads(P[i].Pos,h,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
                wtsum = 0.0;
                if(numngb>0)
                {
                    for(n = 0; n < numngb; n++)
                    {
                        j = Ngblist[n];
                        if(P[j].Type == 0 && P[j].Mass > 0)
                        {
                            wtsum += 1.0;
                        }
                    }
                    if(wtsum>0)
                    {
                        dE /= wtsum;
                        for(n = 0; n < numngb; n++)
                        {
                            j = Ngblist[n];
                            if(P[j].Type == 0 && P[j].Mass > 0)
                            {
                                dE_j = dE / P[j].Mass;
                                SphP[j].InternalEnergy += dE_j;
                                SphP[j].InternalEnergyPred += dE_j;
                            }
                        }
                    } // if(wtsum>0)
                } // if(numngb>0)
            } // if((star_age < 0.004)&&(dt>0)&&(P[i].Mass>0))
        } // if((P[i].Type == 4)||(P[i].Type == 2)||(P[i].Type == 3))
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    myfree(Ngblist);
} // void luminosity_heating_gasoline(void)
#endif


#endif /* GALSF_FB_SNE_HEATING */

