#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#include "dm_interaction.h"
#include "MT.h"
#define NDEBUG
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
#ifdef DM_BARYON_INTERACTION
struct dm_Conserved_var_Riemann
{
    MyDouble rho;
    MyDouble p;
    MyDouble v[3];
    MyDouble u;
    MyDouble cs;
};


struct dm_kernel_hydra
{
    double dp[3];
    double r, vsig, sound_i, sound_j;
    double dv[3], vrel[3], v_rel, vdotr2;
    double wk_i, wk_j, dwk_i, dwk_j;
    double h_i, h_j, dwk_ij, rho_ij_inv;
    double spec_egy_u_i;
};
struct dm_hydrodata_in
{
    /* basic hydro variables */
    MyDouble Pos[3];
    MyFloat Vel[3];
    MyFloat dm_Hsml;
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat Density;
    MyFloat dm_density;
    MyFloat dm_coll;
    MyFloat vx_dm;
    MyFloat vy_dm;
    MyFloat vz_dm;
    MyFloat vx_baryon;
    MyFloat vy_baryon;
    MyFloat vz_baryon;
    MyFloat dm_count;
    MyFloat baryon_count;
    MyFloat Pressure;
    MyFloat ConditionNumber;
    MyFloat InternalEnergyPred;
    MyFloat SoundSpeed;
    int Timestep;

    
    /* matrix of the conserved variable gradients: rho, u, vx, vy, vz */
    struct
    {
        MyDouble Density[3];
        MyDouble Pressure[3];
        MyDouble Velocity[3][3];
#ifndef DONOTUSENODELIST
    int NodeList[NODELISTLENGTH];
#endif

     }
}
*DM_HydroDataIn, *DM_HydroDataGet;


/* --------------------------------------------------------------------------------- */
/* outputs: this is what the routine needs to return to the particles to set their final values */
/* --------------------------------------------------------------------------------- */

struct dm_hydrodata_out
{
    
    MyLongDouble dm_DtInternalEnergy;
    //MyLongDouble dInternalEnergy; //manifest-indiv-timestep-debug//
    MyDouble baryon_dtVel[3];
    int baryon_count;
     
 }  *DM_HydroDataResult, *DM_HydroDataOut;

static inline void dm_particle2in_hydra(struct dm_hydrodata_in *in, int i);
static inline void dm_out2particle_hydra(struct dm_hydrodata_out *out, int i, int mode);
static inline void dm_particle2in_hydra(struct dm_hydrodata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}
    }
    in->dm_Hsml = PPP[i].dm_Hsml;
    in->Hsml  = PPP[i].Hsml;
    in->Mass = P[i].Mass;
    in->Density = SphP[i].Density;
    in->dm_density = SphP[i].dm_density;
    in->dm_coll = SphP[i].dm_coll;
    in->vx_dm = SphP[i].vx_dm;
    in->vy_dm = SphP[i].vy_dm;
    in->vz_dm = SphP[i].vz_dm;
    in->vx_baryon = SphP[i].vx_baryon;
    in->vy_baryon = SphP[i].vy_baryon;
    in->vz_baryon = SphP[i].vz_baryon;
    in->dm_count = SphP[i].dm_count;
    in->baryon_count = SphP[i].baryon_count;
    in->Pressure = SphP[i].Pressure;
    in->InternalEnergyPred = SphP[i].InternalEnergyPred;
    in->SoundSpeed = Particle_effective_soundspeed_i(i);
    in->Timestep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
    in->ConditionNumber = SphP[i].ConditionNumber;
     

}



static inline void dm_out2particle_hydra(struct dm_hydrodata_out *out, int i, int mode)
{
//    int k;
    
//    for(k = 0; k < 3; k++)
//    {
//      SphP[i].baryon_dtVel[k] += out->baryon_dtVel[k];
//     }
//    SphP[i].dm_DtInternalEnergy += out->dm_DtInternalEnergy;
 
  if(P[i].Type == 0){ 
    SphP[i].count0 += out->baryon_count;   
    }
}


#include "dm_hydro_evaluate.h"

void dm_hydro_final_operations_and_cleanup(void)
{
    int i,k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            double dt;
            dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
        }
        

            
         
/*            for(k=0;k<3;k++)
            {
                SphP[i].DtInternalEnergy -= (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].HydroAccel[k];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[i].DtInternalEnergy += 0.5 * (SphP[i].VelPred[k]/All.cf_atime) * (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].DtMass;
                SphP[i].HydroAccel[k] -= (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].DtMass; 
#endif
                SphP[i].HydroAccel[k] /= P[i].Mass;
            }*/
            
            
            SphP[i].dm_DtInternalEnergy /= P[i].Mass;
            
            
            if(PPP[i].dm_Hsml >= 0.99*All.MaxHsml) {SphP[i].dm_DtInternalEnergy = 0;}
            
               /*         if(All.ComovingIntegrationOn) SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a;*/
            
            
       }
}





void dm_hydro_force(void)
{
    int i, j, k, ngrp, ndone, ndone_flag;
    int recvTask, place;
    double timeall=0, timecomp1=0, timecomp2=0, timecommsumm1=0, timecommsumm2=0, timewait1=0, timewait2=0, timenetwork=0;
    double timecomp, timecomm, timewait, tstart, tend, t0, t1;
    int save_NextParticle;
    long long n_exported = 0;
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type==0)
        {
#ifdef DM_BARYON_INTERACTION
             SphP[i].dm_DtInternalEnergy = 0;
              for(k=0;k<3;k++)
            {
                SphP[i].baryon_dtVel[k] = 0;
            }
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            SphP[i].MaxKineticEnergyNgb = -1.e10;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            SphP[i].DtMass = 0;
            SphP[i].dMass = 0;
            for(k=0;k<3;k++) SphP[i].GravWorkTerm[k] = 0;

#if defined(RT_EVOLVE_NGAMMA_IN_HYDRO)
            for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[i].Dt_E_gamma[k] = 0;}
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Dt_Flux[k][k_dir] = 0;}}
#endif
            
#ifdef MAGNETIC
            SphP[i].divB = 0;
            for(k=0;k<3;k++) {SphP[i].Face_Area[k] = 0;}
#ifdef DIVBCLEANING_DEDNER
            for(k=0;k<3;k++) {SphP[i].DtB_PhiCorr[k] = 0;}
#endif
#ifndef HYDRO_SPH
            for(k=0;k<3;k++) {SphP[i].DtB[k] = 0;}
#ifdef DIVBCLEANING_DEDNER
            SphP[i].DtPhi = 0;
#endif
#endif
#endif
#ifdef WAKEUP
            PPPZ[i].wakeup = 0;
#endif
#endif
        }else{
#ifdef DM_BARYON_INTERACTION
            for(k=0;k<3;k++){P[i].dm_dtVel[k] = 0;}
#endif
}
    
       long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct dm_hydrodata_in) +
                                                             sizeof(struct dm_hydrodata_out) +
                                                             sizemax(sizeof(struct dm_hydrodata_in),
                                                                     sizeof(struct dm_hydrodata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    CPU_Step[CPU_HYDMISC] += measure_time();
    t0 = my_second();
    NextParticle = FirstActiveParticle;
    
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
            pthread_create(&mythreads[j], &attr, dm_hydro_evaluate_primary, &threadid[j]);
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
            dm_hydro_evaluate_primary(&mainthreadid);
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
               
                endrun(115508);
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
        {
            Send_count[j] = 0;
            Recv_count[j] = 0;
        }
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
        DM_HydroDataGet = (struct dm_hydrodata_in *) mymalloc("DM_HydroDataGet", Nimport * sizeof(struct dm_hydrodata_in));
        DM_HydroDataIn = (struct dm_hydrodata_in *) mymalloc("DM_HydroDataIn", Nexport * sizeof(struct dm_hydrodata_in));
        
        
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            dm_particle2in_hydra(&DM_HydroDataIn[j], place);
#ifndef DONOTUSENODELIST
            memcpy(DM_HydroDataIn[j].NodeList,
                   DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif
            
        }
        
               tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                   
                    MPI_Sendrecv(&DM_HydroDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct dm_hydrodata_in), MPI_BYTE,
                                 recvTask, TAG_HYDRO_A,
                                 &DM_HydroDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct dm_hydrodata_in), MPI_BYTE,
                                 recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm1 += timediff(tstart, tend);
        
        myfree(DM_HydroDataIn);
        DM_HydroDataResult = (struct dm_hydrodata_out *) mymalloc("DM_HydroDataResult", Nimport * sizeof(struct dm_hydrodata_out));
        DM_HydroDataOut = (struct dm_hydrodata_out *) mymalloc("DM_HydroDataOut", Nexport * sizeof(struct dm_hydrodata_out));
        report_memory_usage(&HighMark_sphhydro, "SPH_DM_HYDRO");
        
       
        tstart = my_second();
        NextJ = 0;
#ifdef OMP_NUM_THREADS
        for(j = 0; j < OMP_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, dm_hydro_evaluate_secondary, &threadid[j]);
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
            dm_hydro_evaluate_secondary(&mainthreadid);
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
        
               tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  
                    MPI_Sendrecv(&DM_HydroDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct dm_hydrodata_out),
                                 MPI_BYTE, recvTask, TAG_HYDRO_B,
                                 &DM_HydroDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct dm_hydrodata_out),
                                 MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm2 += timediff(tstart, tend);
        
       
        tstart = my_second();
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            dm_out2particle_hydra(&DM_HydroDataOut[j], place, 1);
        }
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        
        myfree(DM_HydroDataOut);
        myfree(DM_HydroDataResult);
        myfree(DM_HydroDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    
      dm_hydro_final_operations_and_cleanup();
 
    
       t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_HYDCOMPUTE] += timecomp;
    CPU_Step[CPU_HYDWAIT] += timewait;
    CPU_Step[CPU_HYDCOMM] += timecomm;
    CPU_Step[CPU_HYDNETWORK] += timenetwork;
    CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}



void *dm_hydro_evaluate_primary(void *p)
{
    int thread_id = *(int *) p;
    int i, j;
    int *exportflag, *exportnodecount, *exportindex, *ngblist;
    
    ngblist = Ngblist + thread_id * NumPart;
    exportflag = Exportflag + thread_id * NTask;
    exportnodecount = Exportnodecount + thread_id * NTask;
    exportindex = Exportindex + thread_id * NTask;
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
        
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
	if(SphP[i].dm_density > 0)
              {	    
                if(dm_hydro_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
                    break;
              }
        }
        ProcessedFlag[i] = 1;
    }
    return NULL;
}



void *dm_hydro_evaluate_secondary(void *p)
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
        
        dm_hydro_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }
    return NULL;
}

#endif
