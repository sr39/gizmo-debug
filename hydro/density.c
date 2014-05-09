#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef COSMIC_RAYS
#include "../cosmic_rays/cosmic_rays.h"
#endif
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

/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first hydro loop", where the gas densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


struct kernel_density
{
  double dx, dy, dz;
  double r;
  double dvx, dvy, dvz;
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
#if defined(AV_CD10_VISCOSITY_SWITCH)
  MyFloat Accel[3];
#endif
  MyFloat Vel[3];
  MyFloat Hsml;
#ifdef GALSF_SUBGRID_WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
  MyFloat BPred[3];
#endif
  int Type;
}
 *DensDataIn, *DensDataGet;

static struct densdata_out
{
    MyLongDouble Ngb;
    MyLongDouble Rho;
    MyLongDouble DhsmlNgb;
    MyLongDouble Particle_DivVel;
    MyFloat NV_T[3][3];
#ifdef HYDRO_SPH
    MyLongDouble DhsmlHydroSumFactor;
#endif
    
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    MyLongDouble EgyRho;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    MyFloat AGS_zeta;
#endif

#if defined(AV_CD10_VISCOSITY_SWITCH)
    MyFloat NV_D[3][3];
    MyFloat NV_A[3][3];
#endif

#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
    MyFloat divB;
#endif

#if defined(BH_POPIII_SEEDS) || defined(GALSF_FB_LOCAL_UV_HEATING) || defined(GALSF_FB_RPWIND_FROMSTARS) || defined(BH_PHOTONMOMENTUM) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN)
    MyFloat GradRho[3];
#endif
    
#if defined(GRAIN_FLUID)
    MyLongDouble RhoGrains;
    MyLongDouble GrainVel[3];
    MyLongDouble SmoothedEntr;
#endif
    
#if defined(BLACK_HOLES)
    int BH_TimeBinGasNeighbor;
#endif

#if defined(VS_TURB) || defined(AB_TURB) || defined(GRAIN_FLUID)
    MyLongDouble GasVel[3];
#endif

}
 *DensDataResult, *DensDataOut;

void particle2in_density(struct densdata_in *in, int i);
void out2particle_density(struct densdata_out *out, int i, int mode);
void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
					struct kernel_density *kernel, int j);


void particle2in_density(struct densdata_in *in, int i)
{
    int k;
    in->Type = P[i].Type;
    in->Hsml = PPP[i].Hsml;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}
    }
    
    if(P[i].Type == 0)
    {
#if defined(AV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            in->Accel[k] = All.cf_a2inv*P[i].GravAccel[k] + SphP[i].HydroAccel[k]; // PHYSICAL units //
#endif
        
#ifdef GALSF_SUBGRID_WINDS
        in->DelayTime = SphP[i].DelayTime;
#endif
        
#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
        for(k = 0; k < 3; k++)
            in->BPred[k] = SphP[i].BPred[k];
#endif
    }
}


void out2particle_density(struct densdata_out *out, int i, int mode)
{
    int j,k;
    ASSIGN_ADD(PPP[i].NumNgb, out->Ngb, mode);
    ASSIGN_ADD(PPPZ[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
    ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel,   mode);
    
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
#endif
    
    if(P[i].Type == 0)
    {
        ASSIGN_ADD(SphP[i].Density, out->Rho, mode);

        for(k = 0; k < 3; k++)
            for(j = 0; j < 3; j++)
                ASSIGN_ADD(SphP[i].NV_T[k][j], out->NV_T[k][j], mode);

#ifdef HYDRO_SPH
        ASSIGN_ADD(SphP[i].DhsmlHydroSumFactor, out->DhsmlHydroSumFactor, mode);
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
        ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
#endif

#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
        ASSIGN_ADD(SphP[i].EgyWtDensity,   out->EgyRho,   mode);
#endif

#if defined(VS_TURB) || defined(AB_TURB)
        for(k = 0; k < 3; k++)
            ASSIGN_ADD(SphP[i].SmoothedVel[k], out->GasVel[k], mode);
#endif

#if defined(AV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            for(j = 0; j < 3; j++)
            {
                ASSIGN_ADD(SphP[i].NV_D[k][j], out->NV_D[k][j], mode);
                ASSIGN_ADD(SphP[i].NV_A[k][j], out->NV_A[k][j], mode);
            }
#endif
        
#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
        ASSIGN_ADD(SphP[i].divB, out->divB, mode);
#endif

    } // P[i].Type == 0 //

#if (defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS))
    if(P[i].Type == 4)
        ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
#endif

#if defined(GRAIN_FLUID)
    if(P[i].Type > 0)
    {
        ASSIGN_ADD(P[i].Gas_Density, out->Rho, mode);
        ASSIGN_ADD(P[i].Grain_Density, out->RhoGrains, mode);
        ASSIGN_ADD(P[i].Gas_InternalEnergy, out->SmoothedEntr, mode);
        for(k = 0; k<3; k++)
        {
            ASSIGN_ADD(P[i].Gas_Velocity[k], out->GasVel[k], mode);
            ASSIGN_ADD(P[i].Grain_Velocity[k], out->GrainVel[k], mode);
        }
    }
#endif

#if defined(GALSF_FB_RPWIND_FROMSTARS) || defined(GALSF_FB_GASRETURN) || defined(GALSF_FB_HII_HEATING) || defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN )
    if((P[i].Type == 4)||(P[i].Type == 2)||(P[i].Type==3))
        ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
#endif

#if defined(BH_POPIII_SEEDS) || defined(GALSF_FB_LOCAL_UV_HEATING) || defined(GALSF_FB_RPWIND_FROMSTARS) || defined(BH_PHOTONMOMENTUM) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN)
    if(P[i].Type != 0)
        for(k = 0; k<3; k++)
            ASSIGN_ADD(P[i].GradRho[k], out->GradRho[k], mode);
#endif
    

#ifdef BLACK_HOLES
    if(P[i].Type == 5)
    {
        ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
        if(mode == 0)
            BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;
        else
        {
            if(BPP(i).BH_TimeBinGasNeighbor > out->BH_TimeBinGasNeighbor)
                BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;
        }
    } /* if(P[i].Type == 5) */
#endif
}



/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  MyFloat *Left, *Right;
  int i, j, k, k1, k2, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac, fac_lim;
  double Tinv[3][3], detT, CNumHolder=0, ConditionNumber=0;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;

#ifdef COSMIC_RAYS
  int CRpop;
#endif
    
  CPU_Step[CPU_DENSMISC] += measure_time();

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(density_isactive(i))
	{
	  Left[i] = Right[i] = 0;

#ifdef BLACK_HOLES
	  P[i].SwallowID = 0;
#endif
	}
    }

  /* allocate buffers to arrange communication */
  size_t MyBufferSize = All.BufferSize;
  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					     sizemax(sizeof(struct densdata_in),
						     sizeof(struct densdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = my_second();

  desnumngb = All.DesNumNgb;
  desnumngbdev = All.MaxNumNgbDeviation;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      NextParticle = FirstActiveParticle;	/* begin with this index */

      do
	{
	  BufferFullFlag = 0;
	  Nexport = 0;
	  save_NextParticle = NextParticle;

	  tstart = my_second();

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
	      pthread_create(&mythreads[j], &attr, density_evaluate_primary, &threadid[j]);
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
	    density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
	  }

#ifdef OMP_NUM_THREADS
	  for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);
#endif

	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);

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
		  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
			 P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
		  if(P[NextParticle].Type == 0)
		    printf("   rho=%g hsml=%g\n",SphP[NextParticle].Density,PPP[NextParticle].Hsml);

		  endrun(112208);
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

	  tstart = my_second();

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

	  tend = my_second();
	  timewait1 += timediff(tstart, tend);

	  for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      Nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", Nimport * sizeof(struct densdata_in));
	  DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", Nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      particle2in_density(&DensDataIn[j], place);

	      memcpy(DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }
	  /* exchange particle data */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = my_second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataIn);
	  DensDataResult =
	    (struct densdata_out *) mymalloc("DensDataResult", Nimport * sizeof(struct densdata_out));
	  DensDataOut =
	    (struct densdata_out *) mymalloc("DensDataOut", Nexport * sizeof(struct densdata_out));

	  report_memory_usage(&HighMark_sphdensity, "SPH_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = my_second();

	  NextJ = 0;

#ifdef OMP_NUM_THREADS
	  for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	    pthread_create(&mythreads[j], &attr, density_evaluate_secondary, &threadid[j]);
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
	    density_evaluate_secondary(&mainthreadid);
	  }

#ifdef OMP_NUM_THREADS
	  for(j = 0; j < OMP_NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);

	  pthread_mutex_destroy(&mutex_partnodedrift);
	  pthread_mutex_destroy(&mutex_nexport);
	  pthread_attr_destroy(&attr);
#endif

	  tend = my_second();
	  timecomp2 += timediff(tstart, tend);

	  if(NextParticle < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = my_second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  tend = my_second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = my_second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = my_second();
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      out2particle_density(&DensDataOut[j], place, 1);
	    }
	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);


        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(density_isactive(i))
            {
                PPPZ[i].DhsmlNgbFactor *= PPP[i].Hsml / (NUMDIMS * PPP[i].NumNgb);
                P[i].Particle_DivVel /= PPP[i].NumNgb;
                /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                PPP[i].NumNgb *= NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS);
                
                // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(PPPZ[i].DhsmlNgbFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
                    PPPZ[i].DhsmlNgbFactor = 1 / (1 + PPPZ[i].DhsmlNgbFactor);
                else
                    PPPZ[i].DhsmlNgbFactor = 1;
                P[i].Particle_DivVel *= PPPZ[i].DhsmlNgbFactor;
            
                if(P[i].Type == 0)
                {
                    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
                    SphP[i].NV_T[1][0]=SphP[i].NV_T[0][1]; SphP[i].NV_T[2][0]=SphP[i].NV_T[0][2]; SphP[i].NV_T[2][1]=SphP[i].NV_T[1][2];
                    /* Now invert the NV_T matrix we just measured */
                    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
                        this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
                    ConditionNumber=CNumHolder=0;
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {ConditionNumber += SphP[i].NV_T[k1][k2]*SphP[i].NV_T[k1][k2];}}
#ifdef ONEDIM
                    /* one-dimensional case */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    detT = SphP[i].NV_T[0][0];
                    if(SphP[i].NV_T[0][0]!=0 && !isnan(SphP[i].NV_T[0][0])) Tinv[0][0] = 1/detT; /* only one non-trivial element in 1D! */
#endif
#ifdef TWODIMS
                    /* two-dimensional case */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    detT = SphP[i].NV_T[0][0]*SphP[i].NV_T[1][1] - SphP[i].NV_T[0][1]*SphP[i].NV_T[1][0];
                    if((detT != 0)&&(!isnan(detT)))
                    {
                        Tinv[0][0] = SphP[i].NV_T[1][1] / detT;
                        Tinv[0][1] = -SphP[i].NV_T[0][1] / detT;
                        Tinv[1][0] = -SphP[i].NV_T[1][0] / detT;
                        Tinv[1][1] = SphP[i].NV_T[0][0] / detT;
                    }
#endif
#if !defined(ONEDIM) && !defined(TWODIMS)
                    /* three-dimensional case */
                    detT = SphP[i].NV_T[0][0] * SphP[i].NV_T[1][1] * SphP[i].NV_T[2][2] +
                        SphP[i].NV_T[0][1] * SphP[i].NV_T[1][2] * SphP[i].NV_T[2][0] +
                        SphP[i].NV_T[0][2] * SphP[i].NV_T[1][0] * SphP[i].NV_T[2][1] -
                        SphP[i].NV_T[0][2] * SphP[i].NV_T[1][1] * SphP[i].NV_T[2][0] -
                        SphP[i].NV_T[0][1] * SphP[i].NV_T[1][0] * SphP[i].NV_T[2][2] -
                        SphP[i].NV_T[0][0] * SphP[i].NV_T[1][2] * SphP[i].NV_T[2][1];
                    /* check for zero determinant */
                    if((detT != 0) && !isnan(detT))
                    {
                        Tinv[0][0] = (SphP[i].NV_T[1][1] * SphP[i].NV_T[2][2] - SphP[i].NV_T[1][2] * SphP[i].NV_T[2][1]) / detT;
                        Tinv[0][1] = (SphP[i].NV_T[0][2] * SphP[i].NV_T[2][1] - SphP[i].NV_T[0][1] * SphP[i].NV_T[2][2]) / detT;
                        Tinv[0][2] = (SphP[i].NV_T[0][1] * SphP[i].NV_T[1][2] - SphP[i].NV_T[0][2] * SphP[i].NV_T[1][1]) / detT;
                        Tinv[1][0] = (SphP[i].NV_T[1][2] * SphP[i].NV_T[2][0] - SphP[i].NV_T[1][0] * SphP[i].NV_T[2][2]) / detT;
                        Tinv[1][1] = (SphP[i].NV_T[0][0] * SphP[i].NV_T[2][2] - SphP[i].NV_T[0][2] * SphP[i].NV_T[2][0]) / detT;
                        Tinv[1][2] = (SphP[i].NV_T[0][2] * SphP[i].NV_T[1][0] - SphP[i].NV_T[0][0] * SphP[i].NV_T[1][2]) / detT;
                        Tinv[2][0] = (SphP[i].NV_T[1][0] * SphP[i].NV_T[2][1] - SphP[i].NV_T[1][1] * SphP[i].NV_T[2][0]) / detT;
                        Tinv[2][1] = (SphP[i].NV_T[0][1] * SphP[i].NV_T[2][0] - SphP[i].NV_T[0][0] * SphP[i].NV_T[2][1]) / detT;
                        Tinv[2][2] = (SphP[i].NV_T[0][0] * SphP[i].NV_T[1][1] - SphP[i].NV_T[0][1] * SphP[i].NV_T[1][0]) / detT;
                    } else {
                        for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    }
#endif
                    
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {CNumHolder += Tinv[k1][k2]*Tinv[k1][k2];}}
                    ConditionNumber = sqrt(ConditionNumber*CNumHolder) / NUMDIMS;
                    if(ConditionNumber<1) ConditionNumber=1;
                    /* this = sqrt( ||NV_T^-1||*||NV_T|| ) :: should be ~1 for a well-conditioned matrix */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {SphP[i].NV_T[k1][k2]=Tinv[k1][k2];}}
                    /* now NV_T holds the inverted matrix elements, for use in hydro */
                } // P[i].Type == 0 //
                
                /* now check whether we had enough neighbours */
                double ncorr_ngb = 1.0;
                if(P[i].Type==0)
                {
                    /* use the previous timestep condition number to correct how many neighbors we should use for stability */
                    if((iter==0)&&(ConditionNumber>SphP[i].ConditionNumber))
                    {
                        /* if we find ourselves with a sudden increase in condition number - check if we have a reasonable 
                            neighbor number for the previous iteration, and if so, use the new (larger) correction */
                        ncorr_ngb = DMIN(2.0,sqrt(1.0 + SphP[i].ConditionNumber/((double)CONDITION_NUMBER_DANGER)));
                        double dn_ngb = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(All.MaxNumNgbDeviation*ncorr_ngb);
                        ncorr_ngb = DMIN(2.0,sqrt(1.0 + ConditionNumber/((double)CONDITION_NUMBER_DANGER)));
                        double dn_ngb_alt = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(All.MaxNumNgbDeviation*ncorr_ngb);
                        dn_ngb = DMIN(dn_ngb,dn_ngb_alt);
                        if(dn_ngb < 10.0) SphP[i].ConditionNumber = ConditionNumber;
                    }
                    ncorr_ngb = DMIN(2.0, sqrt(1.0 + SphP[i].ConditionNumber/((double)CONDITION_NUMBER_DANGER)));
                }
                
                desnumngb = All.DesNumNgb * ncorr_ngb;
                desnumngbdev = All.MaxNumNgbDeviation * ncorr_ngb;
                
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                {
                    desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
                    desnumngbdev = 4 * (All.BlackHoleNgbFactor+1);
                }
#endif
                
#if defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS)
                if(P[i].Type == 4)
                    desnumngb = 64; /* will assign the stellar luminosity to very few (one actually) gas particles */
#endif
                
#if defined(GALSF_FB_RPWIND_FROMSTARS) || defined(GALSF_FB_GASRETURN) || defined(GALSF_FB_HII_HEATING) || defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN)
                /* use a much looser check for N_neighbors when the central point is a star particle,
                 since the accuracy is limited anyways to the coupling efficiency -- the routines use their
                 own estimators+neighbor loops, anyways, so this is just to get some nearby particles */
                if((P[i].Type!=0)&&(P[i].Type!=5))
                {
                    desnumngb = All.DesNumNgb;
                    desnumngbdev = All.DesNumNgb / 4;
                }
#endif
                
                redo_particle = 0;
                
                if(PPP[i].NumNgb < (desnumngb - desnumngbdev) ||
                   (PPP[i].NumNgb > (desnumngb + desnumngbdev) && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
                    redo_particle = 1;
                
                if(PPP[i].NumNgb < (desnumngb - desnumngbdev) && PPP[i].Hsml >= All.MaxHsml)
                {
                    PPP[i].Hsml = All.MaxHsml;
                    redo_particle = 0;
                }
                
                if((redo_particle==0)&&(P[i].Type == 0))
                {
                    /* ok we have reached the desired number of neighbors: save the condition number for next timestep */
                    if(ConditionNumber > 1000.0 * (double)CONDITION_NUMBER_DANGER)
                    {
                        printf("Warning: Condition number=%g CNum_prevtimestep=%g Num_Ngb=%g desnumngb=%g Hsml=%g Hsml_min=%g \n",
                               ConditionNumber,SphP[i].ConditionNumber,PPP[i].NumNgb,desnumngb,PPP[i].Hsml,All.MinGasHsml);
                        fflush(stdout);
                    }
                    SphP[i].ConditionNumber = ConditionNumber;
                }
                
                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        printf("i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                               i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, PPP[i].Hsml, PPPZ[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                        fflush(stdout);
                    }
                    
                    /* need to redo this particle */
                    npleft++;
                    
                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            continue;
                        }
                    
                    if(PPP[i].NumNgb < (desnumngb - desnumngbdev))
                        Left[i] = DMAX(PPP[i].Hsml, Left[i]);
                    else
                    {
                        if(Right[i] != 0)
                        {
                            if(PPP[i].Hsml < Right[i])
                                Right[i] = PPP[i].Hsml;
                        }
                        else
                            Right[i] = PPP[i].Hsml;
                    }
                    
                    // right/left define upper/lower bounds from previous iterations
                    if(Right[i] > 0 && Left[i] > 0)
                    {
                        // geometric interpolation between right/left //
                        PPP[i].Hsml *= exp(PPPZ[i].DhsmlNgbFactor * log( desnumngb / PPP[i].NumNgb ) / NUMDIMS);
                        if((PPP[i].Hsml<=Right[i])||(PPP[i].Hsml>=Left[i]))
                        {
                            if(PPP[i].Hsml>Right[i]) PPP[i].Hsml=Right[i];
                            if(PPP[i].Hsml<Left[i]) PPP[i].Hsml=Left[i];
                            PPP[i].Hsml = pow(PPP[i].Hsml * Left[i] * Right[i] , 1.0/3.0);
                        }
                    }
                    else
                    {
                        if(Right[i] == 0 && Left[i] == 0)
                        {
                            char buf[1000];
                            sprintf(buf, "Right[i] == 0 && Left[i] == 0 && PPP[i].Hsml=%g\n", PPP[i].Hsml);
                            terminate(buf);
                        }
                        
                        if(Right[i] == 0 && Left[i] > 0)
                        {
                            if (PPP[i].NumNgb > 1)
                                fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
                            else
                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                            
                            //if(fabs(PPP[i].NumNgb - desnumngb) < 0.75 * desnumngb)
                            if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                            {
                                fac = fac_lim * PPPZ[i].DhsmlNgbFactor; // account for derivative in making the 'corrected' guess
                                if(iter>=20)
                                    if(PPPZ[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                
                                if(fac < fac_lim+0.231)
                                {
                                    PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                }
                                else
                                {
                                    PPP[i].Hsml *= exp(fac_lim+0.231);
                                    // fac~0.26 leads to expected doubling of number if density is constant,
                                    //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                }
                            }
                            else
                                PPP[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                        }
                        
                        if(Right[i] > 0 && Left[i] == 0)
                        {
                            if (PPP[i].NumNgb > 1)
                                fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                            else
                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                            
                            if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
                            
                            //if(fabs(PPP[i].NumNgb - desnumngb) < 0.75 * desnumngb)
                            if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                            {
                                fac = fac_lim * PPPZ[i].DhsmlNgbFactor; // account for derivative in making the 'corrected' guess
                                if(iter>=20)
                                    if(PPPZ[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                
                                if(fac > fac_lim-0.231)
                                {
                                    PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                }
                                else
                                    PPP[i].Hsml *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
                            }
                            else
                                PPP[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                        }
                    }
                    
                    if(PPP[i].Hsml < All.MinGasHsml)
                        PPP[i].Hsml = All.MinGasHsml;
                    
#ifdef BLACK_HOLES
                    if(P[i].Type == 5)
                        if(Left[i] > All.BlackHoleMaxAccretionRadius)
                        {
                            /* this will stop the search for a new BH smoothing length in the next iteration */
                            PPP[i].Hsml = Left[i] = Right[i] = All.BlackHoleMaxAccretionRadius;
                        }
#endif
                }
                else
                    P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
            }
        }
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0)
            {
                printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
                fflush(stdout);
            }
            if(iter > MAXITER)
            {
                printf("failed to converge in neighbour iteration in density()\n");
                fflush(stdout);
                endrun(1155);
            }
        }
    }
    while(ntot > 0);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Right);
    myfree(Left);
    myfree(Ngblist);
    
    
    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0)
            P[i].TimeBin = -P[i].TimeBin - 1;
    }
    
    
    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results
     ( any quantities that only need to be evaluated once, on the final iteration --
     won't save much b/c the real cost is in the neighbor loop for each particle, but it's something )
     -- also, some results (for example, viscosity suppression below) should not be calculated unless
     the quantities are 'stabilized' at their final values -- */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(density_isactive(i))
        {
            if(P[i].Type == 0 && P[i].Mass > 0)
            {
                if(SphP[i].Density > 0)
                {
#ifdef HYDRO_SPH
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
                    SphP[i].EgyWtDensity /= SphP[i].InternalEnergyPred;
#endif
                    /* need to divide by the sum of x_tilde=1, i.e. numden_ngb */
                    double numden_ngb = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                    SphP[i].DhsmlHydroSumFactor *= PPP[i].Hsml / (NUMDIMS * numden_ngb);
                    SphP[i].DhsmlHydroSumFactor *= -PPPZ[i].DhsmlNgbFactor; /* now this is ready to be called in hydro routine */
#endif
                    
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                    double ndenNGB = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                    PPPZ[i].AGS_zeta *= 0.5 * P[i].Mass * PPP[i].Hsml / (NUMDIMS * ndenNGB) * PPPZ[i].DhsmlNgbFactor;
#endif
                    
#if defined(AV_CD10_VISCOSITY_SWITCH)
                    for(k1 = 0; k1 < 3; k1++)
                        for(k2 = 0; k2 < 3; k2++)
                        {
                            SphP[i].NV_D[k2][k1] *= All.cf_a2inv; // converts to physical velocity/length
                            SphP[i].NV_A[k2][k1] /= All.cf_atime; // converts to physical accel/length
                        }
                    // all quantities below in this block should now be in proper PHYSICAL units, for subsequent operations //
                    double dtDV[3][3], A[3][3], V[3][3], S[3][3];
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            V[k1][k2] = SphP[i].NV_D[k1][0]*Tinv[0][k2] + SphP[i].NV_D[k1][1]*Tinv[1][k2] + SphP[i].NV_D[k1][2]*Tinv[2][k2];
                            A[k1][k2] = SphP[i].NV_A[k1][0]*Tinv[0][k2] + SphP[i].NV_A[k1][1]*Tinv[1][k2] + SphP[i].NV_A[k1][2]*Tinv[2][k2];
                        }
                    SphP[i].NV_DivVel = V[0][0] + V[1][1] + V[2][2];
                    SphP[i].NV_trSSt = 0;
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            dtDV[k1][k2] = A[k1][k2] - (V[k1][0]*V[0][k2] + V[k1][1]*V[1][k2] + V[k1][2]*V[2][k2]);
                            /* S = 0.5*(V+V_transpose) - delta_ij*div_v/3 */
                            S[k1][k2] = 0.5 * (V[k1][k2] + V[k2][k1]);
                            if(k2==k1) S[k1][k2] -= SphP[i].NV_DivVel / NUMDIMS;
                            /* Trace[S*S_transpose] = SSt[0][0]+SSt[1][1]+SSt[2][2] = |S|^2 = sum(Sij^2) */
                            SphP[i].NV_trSSt += S[k1][k2]*S[k1][k2];
                        }
                    SphP[i].NV_dt_DivVel = dtDV[0][0] + dtDV[1][1] + dtDV[2][2];
#endif
                    
                    
#if defined(VS_TURB) || defined(AB_TURB)
                    SphP[i].SmoothedVel[0] /= SphP[i].Density;
                    SphP[i].SmoothedVel[1] /= SphP[i].Density;
                    SphP[i].SmoothedVel[2] /= SphP[i].Density;
#endif
                    
#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
                    SphP[i].divB *= 1 / SphP[i].Density; // add DhsmlNgbFactor per TP correction //
#endif
                }
                
#ifndef HYDRO_SPH
                SphP[i].Density = P[i].Mass * PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
#endif
                SphP[i].Pressure = get_pressure(i);		// should account for density independent pressure

            } // P[i].Type == 0
            
#ifdef PM_HIRES_REGION_CLIPPING
#ifdef BLACK_HOLES
            if (P[i].Type != 5)
            {
#endif
                if(P[i].Type == 0) if ((SphP[i].Density <= 0) || (PPP[i].NumNgb <= 0)) P[i].Mass = 0;
                if ((PPP[i].Hsml <= 0) || (PPP[i].Hsml >= PM_HIRES_REGION_CLIPPING)) P[i].Mass = 0;
                double vmag=0; for(k=0;k<3;k++) vmag+=P[i].Vel[k]*P[i].Vel[k]; vmag = sqrt(vmag);
                if(vmag>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) P[i].Mass=0;
                if(vmag>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) for(k=0;k<3;k++) P[i].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/vmag;
#ifdef BLACK_HOLES
            }
#endif // BLACK_HOLES
#endif // ifdef PM_HIRES_REGION_CLIPPING
            
            
        } // density_isactive(i)
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_DENSCOMPUTE] += timecomp;
    CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}






/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
                     int *ngblist)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
    double r2, h2, u, mass_j, wk;
    struct kernel_density kernel;
    struct densdata_in local;
    struct densdata_out out;
    memset(&out, 0, sizeof(struct densdata_out));
#if defined(BLACK_HOLES)
    out.BH_TimeBinGasNeighbor = TIMEBINS;
#endif
    
    if(mode == 0)
        particle2in_density(&local, target);
    else
        local = DensDataGet[target];
    h2 = local.Hsml * local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox =
            ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
                                          exportnodecount, exportindex, ngblist);
            
            if(numngb_inbox < 0) return -1;
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
#ifdef GALSF_SUBGRID_WINDS
                if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
                    if(!(local.DelayTime > 0))	/* if I'm not wind, then ignore the wind particle */
                        continue;
#endif
                if(P[j].Mass <= 0) continue;
                
                kernel.dx = local.Pos[0] - P[j].Pos[0];
                kernel.dy = local.Pos[1] - P[j].Pos[1];
                kernel.dz = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  find the closest image in the given box size  */
                kernel.dx = NEAREST_X(kernel.dx);
                kernel.dy = NEAREST_Y(kernel.dy);
                kernel.dz = NEAREST_Z(kernel.dz);
#endif
                r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
                
                if(r2 < h2)
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);
                    
                    out.Ngb += kernel.wk;
                    out.Rho += kernel.mj_wk;
                    out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#ifdef HYDRO_SPH
                    double mass_eff = mass_j;
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
                    mass_eff *= SphP[j].InternalEnergyPred;
                    out.EgyRho += kernel.wk * mass_eff;
#endif
                    out.DhsmlHydroSumFactor += -mass_eff * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#endif
                    
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
                    out.AGS_zeta += mass_j * kernel_gravity(u, kernel.hinv, kernel.hinv3, 0);
#endif
                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        if(local.Type == 0)
                        {
                            wk = kernel.wk; /* MAKE SURE THIS MATCHES CHOICE IN GRADIENTS.c!!! */
                            /* the weights for the MLS tensor used for gradient estimation */
                            out.NV_T[0][0] +=  wk * kernel.dx * kernel.dx;
                            out.NV_T[0][1] +=  wk * kernel.dx * kernel.dy;
                            out.NV_T[0][2] +=  wk * kernel.dx * kernel.dz;
                            out.NV_T[1][1] +=  wk * kernel.dy * kernel.dy;
                            out.NV_T[1][2] +=  wk * kernel.dy * kernel.dz;
                            out.NV_T[2][2] +=  wk * kernel.dz * kernel.dz;
                        }
                        kernel.dvx = local.Vel[0] - SphP[j].VelPred[0];
                        kernel.dvy = local.Vel[1] - SphP[j].VelPred[1];
                        kernel.dvz = local.Vel[2] - SphP[j].VelPred[2];
                        out.Particle_DivVel += kernel.dwk * (kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz) / kernel.r;
                        /* not-very accurate SPH div-v estimator: however, it exactly describes the -particle- drift */
                        
                        density_evaluate_extra_physics_gas(&local, &out, &kernel, j);
                    } // kernel.r > 0 //
                }
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        out2particle_density(&out, target, 0);
    else
        DensDataResult[target] = out;
    
    return 0;
}



void *density_evaluate_primary(void *p)
{
    int thread_id = *(int *) p;
    int i, j;
    int *exportflag, *exportnodecount, *exportindex, *ngblist;
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
        
        if(density_isactive(i))
        {
            if(density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
                break;		/* export buffer has filled up */
        }
        ProcessedFlag[i] = 1;	/* particle successfully finished */
    }
    return NULL;
}



void *density_evaluate_secondary(void *p)
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
        
        density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }
    
    return NULL;
    
}




int density_isactive(int n)
{
    /* first check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
    if(P[n].TimeBin < 0)
        return 0;
    
#if defined(GRAIN_FLUID)
    /* all particles can potentially interact with the gas in this mode, if drag > 0 */
    if(P[n].Type >= 0)
        return 1;
#endif
    
#if (defined(RADTRANSFER) && defined(EDDINGTON_TENSOR_STARS))
    if(P[n].Type == 4)
        return 1;
#endif
    
#if defined(GALSF_FB_RPWIND_FROMSTARS) || defined(GALSF_FB_GASRETURN) || defined(GALSF_FB_HII_HEATING) || defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN )
    if(((P[n].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[n].Type == 2)||(P[n].Type==3))))&&(P[n].Mass>0))
    {
#if defined(GALSF_FB_SNE_HEATING)
        /* check if there is going to be a SNe this timestep, in which case, we want the density info! */
        if(P[n].SNe_ThisTimeStep>0) return 1;
#endif
#if defined(GALSF_FB_GASRETURN)
        if(P[n].MassReturn_ThisTimeStep>0) return 1;
#endif
#if defined(GALSF_FB_RPROCESS_ENRICHMENT)
        if(P[n].RProcessEvent_ThisTimeStep>0) return 1;
#endif
        if(P[n].DensAroundStar<=0) return 1;
        // only do stellar age evaluation if we have to //
        float star_age=0;
        star_age = evaluate_stellar_age_Gyr(P[n].StellarAge);
        if(star_age < 0.035) return 1;
    }
#endif
    
#ifdef BLACK_HOLES
    if(P[n].Type == 5)
        return 1;
#endif
    
    if(P[n].Type == 0 && P[n].Mass > 0)
        return 1;
    
    return 0;
}





void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
                                        struct kernel_density *kernel, int j)
{
    kernel->mj_dwk_r = P[j].Mass * kernel->dwk / kernel->r;

    if(local->Type != 0)
    {
        
#if defined(GRAIN_FLUID)
        out->SmoothedEntr += FLT(kernel->mj_wk * SphP[j].InternalEnergy);
#endif
        
#if defined(BLACK_HOLES)
        if(out->BH_TimeBinGasNeighbor > P[j].TimeBin)
            out->BH_TimeBinGasNeighbor = P[j].TimeBin;
#endif
        
#if defined(VS_TURB) || defined(AB_TURB) || defined(GRAIN_FLUID)
        out->GasVel[0] += FLT(kernel->mj_wk * SphP[j].VelPred[0]);
        out->GasVel[1] += FLT(kernel->mj_wk * SphP[j].VelPred[1]);
        out->GasVel[2] += FLT(kernel->mj_wk * SphP[j].VelPred[2]);
#endif
        
#if defined(BH_POPIII_SEEDS) || defined(GALSF_FB_LOCAL_UV_HEATING) || defined(GALSF_FB_RPWIND_FROMSTARS) || defined(BH_PHOTONMOMENTUM) || defined(GALSF_FB_RT_PHOTON_LOCALATTEN)
        /* this is here because for the models of BH growth and self-shielding of stars, we
         just need a quick-and-dirty, single-pass approximation for the gradients (the error from
         using this as opposed to the higher-order gradient estimators is small compared to the
         Sobolev approximation): use only for -non-gas- particles */
        out->GradRho[0] += kernel->mj_dwk_r * kernel->dx;
        out->GradRho[1] += kernel->mj_dwk_r * kernel->dy;
        out->GradRho[2] += kernel->mj_dwk_r * kernel->dz;
#endif
        
    } else { /* local.Type == 0 */
        
#if defined(AV_CD10_VISCOSITY_SWITCH)
        double wk = kernel->wk;
        out->NV_A[0][0] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dx * wk;
        out->NV_A[0][1] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dy * wk;
        out->NV_A[0][2] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dz * wk;
        out->NV_A[1][0] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dx * wk;
        out->NV_A[1][1] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dy * wk;
        out->NV_A[1][2] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dz * wk;
        out->NV_A[2][0] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dx * wk;
        out->NV_A[2][1] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dy * wk;
        out->NV_A[2][2] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dz * wk;
        
        out->NV_D[0][0] += kernel->dvx * kernel->dx * wk;
        out->NV_D[0][1] += kernel->dvx * kernel->dy * wk;
        out->NV_D[0][2] += kernel->dvx * kernel->dz * wk;
        out->NV_D[1][0] += kernel->dvy * kernel->dx * wk;
        out->NV_D[1][1] += kernel->dvy * kernel->dy * wk;
        out->NV_D[1][2] += kernel->dvy * kernel->dz * wk;
        out->NV_D[2][0] += kernel->dvz * kernel->dx * wk;
        out->NV_D[2][1] += kernel->dvz * kernel->dy * wk;
        out->NV_D[2][2] += kernel->dvz * kernel->dz * wk;
#endif
        
        
#if defined(MAGNETIC) && defined(HYDRO_SPH) && (defined(TRACEDIVB) || defined(DIVBCLEANING_DEDNER))
        /* this lives here because the Lagrangian form of the Dedner cleaning scheme depends on
         the -specific- functional form of the div_B estimator */
        out->divB += FLT(-kernel->mj_dwk_r * ((local->BPred[0] - SphP[j].BPred[0]) * kernel->dx +
                                              (local->BPred[1] - SphP[j].BPred[1]) * kernel->dy +
                                              (local->BPred[2] - SphP[j].BPred[2]) * kernel->dz));
#endif
    
    } // Type = 0 check
}




