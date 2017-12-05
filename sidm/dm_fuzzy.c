#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/*! \file dm_fuzzy.c
 *  \brief routines needed for fuzzy-DM implementation
 *         This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef DM_FUZZY



/* --------------------------------------------------------------------------
 Actual evaluation of fluxes from the quantum pressure tensor
 -------------------------------------------------------------------------- */
void do_dm_fuzzy_flux_computation(double *dp, double *GradRho_L, double *GradRho_R,
                                  double rho_L, double rho_R,
                                  double v_L, double v_R, double *Area, double *fluxes)
{
    int m,n;
    double f0=0.5 * 591569. / (All.FuzzyDM_Mass_in_eV * UnitVelocity_in_cm_per_s * All.UnitLength_in_cm/All.HubbleParam); // this encodes the coefficient with the mass of the particle: units vel*L = hbar / particle_mass
    // (0.5=1/2, pre-factor from dimensionless equations; 591569 = hbar/eV in cm^2/s; add mass in eV, and put in code units
    f0 = f0*f0; // squared in terms below //
    double rhoL_i=1./rho_L, rhoR_i=1./rho_R, r2=0, rSi=1./(rho_L+rho_R), flux[3]={0};
    for(m=0;m<3;m++) {r2+=dp[m]*dp[m]; fluxes[m]=0;} /* zero fluxes and calculate separation */
    double r=sqrt(r2), wavespeed=2.*sqrt(f0)*(M_PI/r); // approximate k = 2pi/lambda = 2pi/(2*dr) as the maximum k the code will allow locally */
    /* note that the QPT admits waves parallel to k, with wavespeed omega = pm 2*sqrt[f0]*k, so include these for HLLC solution */
    for(m=0;m<3;m++)
    {
        for(n=0;n<3;n++)
        {
            /* calculate second-derivative of density in pair-wise faction across the face */
            double drho_mn = (GradRho_R[n]-GradRho_L[n]) * dp[m]/ r2;
            /* calculate pressure tensor based on 1st and new 2nd derivative */
            double QPT_L = f0*(rhoL_i*GradRho_L[m]*GradRho_L[n]-drho_mn) * All.cf_a3inv*All.cf_a2inv; // convert grad^2_rho ~ rho/L^2 from code units to physical
            double QPT_R = f0*(rhoR_i*GradRho_R[m]*GradRho_R[n]-drho_mn) * All.cf_a3inv*All.cf_a2inv; // convert grad^2_rho ~ rho/L^2 from code units to physical
            /* calculate 'star' solution (interface moving with contact wave, since we have a Lagrangian code)
             for HLLC reimann problem based on these pressure tensors */
            double P_star = (QPT_L*rho_R + QPT_R*rho_L) * rSi; // if there were no waves and all at rest
            if(v_L > v_R) {P_star += (v_L-v_R)*(v_L-v_R)*rho_L*rho_R*rSi;} // approach velocities lead to up-wind mixing
            ??? (PFH: need to add appropriate dissipation terms from wavespeeds)
            fluxes[m] += -Area[n] * P_star; /* momentum flux into direction 'm' given by Area.Pressure */
        }
    }
    return;
}



/* --------------------------------------------------------------------------
 very quick sub-routine to get the particle densities from their volumes
 -------------------------------------------------------------------------- */
double get_particle_volume_ags(int j)
{
    double L_j = Get_Particle_Size_AGS(j);
#if (NUMDIMS==1)
    return L_j;
#elif (NUMDIMS==2)
    return L_j*L_j;
#else
    return L_j*L_j*L_j;
#endif
}


/* --------------------------------------------------------------------------
 Subroutine here exists to calculate the MFM-like effective faces for purposes of face-interaction evaluation
 -------------------------------------------------------------------------- */

#if !defined(CBE_INTEGRATOR) // if this isn't defined, need to define this function here so can calculate the matrix for gradients and faces //

/* routine to invert the NV_T matrix after neighbor pass */
double do_cbe_nvt_inversion_for_faces(int i)
{
    MyFloat NV_T[3][3]; int j,k;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=P[i].NV_T[j][k];}} // initialize matrix to be inverted //
    double Tinv[3][3], FrobNorm=0, FrobNorm_inv=0, detT=0;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {Tinv[j][k]=0;}}
    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
    NV_T[1][0]=NV_T[0][1]; NV_T[2][0]=NV_T[0][2]; NV_T[2][1]=NV_T[1][2];
    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
     this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm += NV_T[j][k]*NV_T[j][k];}}
#if (NUMDIMS==1) // 1-D case //
    detT = NV_T[0][0];
    if(detT!=0 && !isnan(detT)) {Tinv[0][0] = 1/detT}; /* only one non-trivial element in 1D! */
#endif
#if (NUMDIMS==2) // 2-D case //
    detT = NV_T[0][0]*NV_T[1][1] - NV_T[0][1]*NV_T[1][0];
    if((detT != 0)&&(!isnan(detT)))
    {
        Tinv[0][0] = NV_T[1][1] / detT; Tinv[0][1] = -NV_T[0][1] / detT;
        Tinv[1][0] = -NV_T[1][0] / detT; Tinv[1][1] = NV_T[0][0] / detT;
    }
#endif
#if (NUMDIMS==3) // 3-D case //
    detT = NV_T[0][0] * NV_T[1][1] * NV_T[2][2] + NV_T[0][1] * NV_T[1][2] * NV_T[2][0] +
    NV_T[0][2] * NV_T[1][0] * NV_T[2][1] - NV_T[0][2] * NV_T[1][1] * NV_T[2][0] -
    NV_T[0][1] * NV_T[1][0] * NV_T[2][2] - NV_T[0][0] * NV_T[1][2] * NV_T[2][1];
    /* check for zero determinant */
    if((detT != 0) && !isnan(detT))
    {
        Tinv[0][0] = (NV_T[1][1] * NV_T[2][2] - NV_T[1][2] * NV_T[2][1]) / detT;
        Tinv[0][1] = (NV_T[0][2] * NV_T[2][1] - NV_T[0][1] * NV_T[2][2]) / detT;
        Tinv[0][2] = (NV_T[0][1] * NV_T[1][2] - NV_T[0][2] * NV_T[1][1]) / detT;
        Tinv[1][0] = (NV_T[1][2] * NV_T[2][0] - NV_T[1][0] * NV_T[2][2]) / detT;
        Tinv[1][1] = (NV_T[0][0] * NV_T[2][2] - NV_T[0][2] * NV_T[2][0]) / detT;
        Tinv[1][2] = (NV_T[0][2] * NV_T[1][0] - NV_T[0][0] * NV_T[1][2]) / detT;
        Tinv[2][0] = (NV_T[1][0] * NV_T[2][1] - NV_T[1][1] * NV_T[2][0]) / detT;
        Tinv[2][1] = (NV_T[0][1] * NV_T[2][0] - NV_T[0][0] * NV_T[2][1]) / detT;
        Tinv[2][2] = (NV_T[0][0] * NV_T[1][1] - NV_T[0][1] * NV_T[1][0]) / detT;
    }
#endif
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm_inv += Tinv[j][k]*Tinv[j][k];}}
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=Tinv[j][k];}} // now NV_T holds the inverted matrix elements //
    double ConditionNumber = DMAX(sqrt(FrobNorm * FrobNorm_inv) / NUMDIMS, 1); // = sqrt( ||NV_T^-1||*||NV_T|| ) :: should be ~1 for a well-conditioned matrix //
#ifdef CBE_DEBUG
    if((ThisTask==0)&&(ConditionNumber>100.)) {printf("Condition number == %g (Task=%d i=%d)\n",ConditionNumber,ThisTask,i);}
#endif
    return ConditionNumber;
}

#endif // CBE_INTEGRATOR




/* --------------------------------------------------------------------------
 Everything below here is a giant block to define the sub-routines needed
 to calculate the higher-order matrix gradient estimators for the density
 field, around each DM element (based on its interacting neighbor set,
 within the AGS_Hsml volume). This will give the density gradients
 AGS_Gradients_Density needed to actually compute the quantum pressure tensor
 -------------------------------------------------------------------------- */


/* define a common 'gradients' structure to hold everything we're going to take derivatives of */
struct Quantities_for_Gradients
{
    MyDouble Density;
};

struct DMGraddata_in
{
    MyDouble Pos[3];
    MyDouble AGS_Hsml;
    struct Quantities_for_Gradients GQuant;
}
*DMGradDataIn, *DMGradDataGet;

struct DMGraddata_out
{
    struct Quantities_for_Gradients Gradients[3];
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
}
*DMGradDataResult, *DMGradDataOut;

/* this is a temporary structure for quantities used ONLY in the loop below, for example for computing the slope-limiters (for the Reimann problem) */
static struct temporary_data_topass
{
    struct Quantities_for_Gradients Maxima;
    struct Quantities_for_Gradients Minima;
}
*DMGradDataPasser;

struct kernel_DMGrad
{
    double dp[3],r,wk_i, wk_j, dwk_i, dwk_j,h_i;
};

static inline void particle2in_DMGrad(struct DMGraddata_in *in, int i);

static inline void particle2in_DMGrad(struct DMGraddata_in *in, int i)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->GQuant.Density = P[i].AGS_Density;
}

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

static inline void out2particle_DMGrad(struct DMGraddata_out *out, int i, int mode);

static inline void out2particle_DMGrad(struct DMGraddata_out *out, int i, int mode)
{
    MAX_ADD(DMGradDataPasser[i].Maxima.Density,out->Maxima.Density,mode);
    MIN_ADD(DMGradDataPasser[i].Minima.Density,out->Minima.Density,mode);
    int k; for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(P[i].AGS_Gradients_Density[k],out->Gradients[k].Density,mode);}
}

void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol);

void local_slopelimiter(double *grad, double valmax, double valmin, double alim, double h, double shoot_tol)
{
    double d_abs = 0.0;
    int k; for(k=0;k<3;k++) {d_abs += grad[k]*grad[k];}
    if(d_abs > 0)
    {
        double cfac = 1 / (alim * h * sqrt(d_abs));
        double fabs_max = fabs(valmax);
        double fabs_min = fabs(valmin);
        double abs_min = DMIN(fabs_max,fabs_min);
        if(shoot_tol > 0)
        {
            double abs_max = DMAX(fabs_max,fabs_min);
            cfac *= DMIN(abs_min + shoot_tol*abs_max, abs_max);
        } else {
            cfac *= abs_min;
        }
        if(cfac < 1) {for(k=0;k<3;k++) {grad[k] *= cfac;}}
    }
}

void construct_gradient(double *grad, int i);

void construct_gradient(double *grad, int i)
{
    /* use the NV_T matrix-based gradient estimator */
    int k; double v_tmp[3];
    for(k=0;k<3;k++) {v_tmp[k] = grad[k];}
    for(k=0;k<3;k++) {grad[k] = P[i].NV_T[k][0]*v_tmp[0] + P[i].NV_T[k][1]*v_tmp[1] + P[i].NV_T[k][2]*v_tmp[2];}
}




void DMGrad_gradient_calc(void)
{
    int i, j, k, k1, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
    double timecomp, timecomm, timewait, tstart, tend, t0, t1;
    long long n_exported = 0, NTaskTimesNumPart;
    
    /* allocate buffers to arrange communication */
    DMGradDataPasser = (struct temporary_data_topass *) mymalloc("DMGradDataPasser",N_gas * sizeof(struct temporary_data_topass));
    NTaskTimesNumPart = maxThreads * NumPart;
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct DMGraddata_in) + sizeof(struct DMGraddata_out) + sizemax(sizeof(struct DMGraddata_in),sizeof(struct DMGraddata_out))));
    CPU_Step[CPU_AGSDENSMISC] += measure_time();
    t0 = my_second();
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type==0) {memset(&DMGradDataPasser[i], 0, sizeof(struct temporary_data_topass));}
    
    /* begin the main gradient loop */
    NextParticle = FirstActiveParticle;    /* begin with this index */
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
        tstart = my_second();
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            DMGrad_evaluate_primary(&mainthreadid);    /* do local particles and prepare export list */
        }
        
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            NextParticle = save_NextParticle;
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle) break;
                if(ProcessedFlag[NextParticle] != 1) break;
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle)
            {
                endrun(123708); /* in this case, the buffer is too small to process even a single particle */
            }
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) {k = j + 1;}
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
                else {new_export++;}
            Nexport = new_export;
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
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
        
        /* prepare particle data for export */
        DMGradDataGet = (struct DMGraddata_in *) mymalloc("DMGradDataGet", Nimport * sizeof(struct DMGraddata_in));
        DMGradDataIn = (struct DMGraddata_in *) mymalloc("DMGradDataIn", Nexport * sizeof(struct DMGraddata_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_DMGrad(&DMGradDataIn[j], place);
            memcpy(DMGradDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                    MPI_Sendrecv(&DMGradDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct DMGraddata_in), MPI_BYTE,
                                 recvTask, TAG_GRADLOOP_A,
                                 &DMGradDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct DMGraddata_in), MPI_BYTE,
                                 recvTask, TAG_GRADLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm1 += timediff(tstart, tend);
        myfree(DMGradDataIn);
        DMGradDataResult = (struct DMGraddata_out *) mymalloc("DMGradDataResult", Nimport * sizeof(struct DMGraddata_out));
        DMGradDataOut = (struct DMGraddata_out *) mymalloc("DMGradDataOut", Nexport * sizeof(struct DMGraddata_out));
        report_memory_usage(&HighMark_DMGrad, "GRADIENTS_LOOP");
        
        /* now do the particles that were sent to us */
        tstart = my_second();
        NextJ = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            DMGrad_evaluate_secondary(&mainthreadid);
        }
        
        tend = my_second();
        timecomp2 += timediff(tstart, tend);
        
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
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
                    MPI_Sendrecv(&DMGradDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct DMGraddata_out),
                                 MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                 &DMGradDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct DMGraddata_out),
                                 MPI_BYTE, recvTask, TAG_GRADLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
            out2particle_DMGrad(&DMGradDataOut[j], place, 1);
        }
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        myfree(DMGradDataOut);
        myfree(DMGradDataResult);
        myfree(DMGradDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* now we can properly calculate (second-order accurate) gradients of hydrodynamic quantities from this loop */
        construct_gradient(P[i].AGS_Gradients_Density,i);
        /* finally, we need to apply a sensible slope limiter to the gradients, to prevent overshooting */
        local_slopelimiter(P[i].AGS_Gradients_Density,DMGradDataPasser[i].Maxima.Density,DMGradDataPasser[i].Minima.Density,0.5,PPP[i].AGS_Hsml,0);
    }
    
    myfree(DMGradDataPasser); /* free the temporary structure we created for the MinMax and additional data passing */
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm;
    CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}


int DMGrad_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    /* define variables */
    int startnode, numngb, listindex = 0, j, k, k2, n, swap_to_j;
    double hinv, hinv3, hinv4, r2, u, hinv_j, hinv3_j, hinv4_j;
    struct kernel_DMGrad kernel;
    struct DMGraddata_in local;
    struct DMGraddata_out out;
    /* zero memory and import data for local target */
    memset(&out, 0, sizeof(struct DMGraddata_out));
    memset(&kernel, 0, sizeof(struct kernel_DMGrad));
    if(mode == 0) {particle2in_DMGrad(&local, target);} else {local = DMGradDataGet[target];}
    /* check if we should bother doing a neighbor loop */
    if(local.AGS_Hsml <= 0) return 0;
    if(local.GQuant.Density <= 0) return 0;
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.AGS_Hsml;
    double h2_i = kernel.h_i*kernel.h_i;
    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
    int kernel_mode_i = -1; // only need to calculate wk, by default
    
    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DMGradDataGet[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.AGS_Hsml, target, &startnode, mode, exportflag,
                                                                  exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            if(numngb < 0) {return -1;} /* no neighbors! */
            for(n = 0; n < numngb; n++) /* neighbor loop */
            {
                j = ngblist[n];
                if((P[j].Mass <= 0)||(P[j].AGS_Density <= 0)) {continue;} /* make sure neighbor is valid */
                /* calculate position relative to target */
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0]; kernel.dp[1] = local.Pos[1] - P[j].Pos[1]; kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC            /*  now find the closest image in the given box size  */
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                if((r2 <= 0) || (r2 >= h2_i)) continue;
                /* calculate kernel quantities needed below */
                kernel.r = sqrt(r2); u = kernel.r * hinv;
                kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, kernel_mode_i);
                /* DIFFERENCE & SLOPE LIMITING: need to check maxima and minima of particle values in the kernel, to avoid 'overshoot' with our gradient estimators. this check should be among all interacting pairs */
                double d_rho = P[j].AGS_Density - local.GQuant.Density;
                MINMAX_CHECK(d_rho,out.Minima.Density,out.Maxima.Density);
                /* now output what's actually needed for the gradient */
                for(k=0;k<3;k++)
                {
                    double wk_xyz_i = -kernel.wk_i * kernel.dp[k]; /* sign is important here! */
                    out.Gradients[k].Density += wk_xyz_i * d_rho;
                }
            } // numngb loop
        } // while(startnode)
        /* continue to open leaves if needed */
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DMGradDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
            }
        }
    }
    /* Collect the result at the right place */
    if(mode == 0) {out2particle_DMGrad(&out, target, 0);} else {DMGradDataResult[target] = out;}
    return 0;
}





void *DMGrad_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(P[i].Type==0)
#define EVALUATION_CALL DMGrad_evaluate(i,0,exportflag,exportnodecount,exportindex,ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}

void *DMGrad_evaluate_secondary(void *p)
{
#define EVALUATION_CALL DMGrad_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}


#endif // DM_FUZZY

