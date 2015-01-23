#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"

/*
 
 This module contains the self-contained sub-routines needed for
 grain-specific physics in proto-planetary/proto-stellar/planetary cases.
 It's also potentially use-able for GMC and ISM scales, and terrestrial
 turbulence. Anything where aerodynamic particles are interesting
 
 
 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 
 */



#ifdef GRAIN_FLUID


/* function to apply the drag on the grains from surrounding gas properties */
void apply_grain_dragforce(void)
{
    
    CPU_Step[CPU_MISC] += measure_time();
    if(ThisTask == 0)
    {
        printf("Beginning Grain Drag Force\n");
        fflush(stdout);
    }
    
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 0)
        {
            if(P[i].Gas_Density > 0)
            {
                double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                if(dt > 0)
                {
                    double cs = sqrt( GAMMA * GAMMA_MINUS1 * P[i].Gas_InternalEnergy);
                    double R_grain = P[i].Grain_Size;
                    double rho_gas = P[i].Gas_Density * All.cf_a3inv;
                    double rho_grain = All.Grain_Internal_Density;
                    double vgas_mag = 0.0;
                    for(k=0;k<3;k++) {vgas_mag+=P[i].Gas_Velocity[k]*P[i].Gas_Velocity[k];}

                    if(vgas_mag > 0)
                    {
                        vgas_mag = sqrt(vgas_mag) / All.cf_atime;
                        double x0 = 0.469993 * vgas_mag/cs; // (3/8)*sqrt[pi/2]*|vgas-vgrain|/cs //
                        double tstop_inv = 1.59577 * rho_gas * cs / (R_grain * rho_grain); // 2*sqrt[2/pi] * 1/tstop //
#ifdef GRAIN_EPSTEIN
                        double mu = 2.3 * PROTONMASS;
                        double temperature = mu * (P[i].Gas_InternalEnergy*All.UnitEnergy_in_cgs*All.HubbleParam/All.UnitMass_in_g) / BOLTZMANN;
                        double cross_section = GRAIN_EPSTEIN * 2.0e-15 * (1. + 70./temperature);
                        cross_section /= (All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam*All.HubbleParam));
                        double n_mol = rho_gas / (mu * All.HubbleParam/All.UnitMass_in_g);
                        double mean_free_path = 1 / (n_mol * cross_section); // should be in code units now //
                        double corr_mfp = R_grain / ((9./4.) * mean_free_path);
                        if(corr_mfp > 1) {tstop_inv /= corr_mfp;}
#endif
                        double C1 = (-1-sqrt(1+x0*x0)) / x0;
                        double C2 = C1 * exp( dt * tstop_inv );
                        double xf = -2 * C2 / (C2*C2 -1);
                        double slow_fac = 1 - xf / x0;
                        // note that, with an external (gravitational) acceleration, we can still solve this equation for the relevant update //
                        
                        double delta_egy = 0;
                        double delta_mom[3];
                        for(k=0; k<3; k++)
                        {
                            double vel_new = P[i].Vel[k] + slow_fac * P[i].Gas_Velocity[k];
                            delta_mom[k] = P[i].Mass * (vel_new - P[i].Vel[k]);
                            delta_egy += 0.5*P[i].Mass * (vel_new*vel_new - P[i].Vel[k]*P[i].Vel[k]);
                            P[i].Vel[k] = vel_new;
                        }
                    

                    
#ifdef GRAIN_BACKREACTION
                        int i,k;
                        double dt, dvel, degy, soundspeed, R_grain, t_stop, slow_fac, vel_new, delta_mom[3], delta_egy;
                        int N_MAX_KERNEL,N_MIN_KERNEL,MAXITER_FB,NITER,startnode,dummy,numngb_inbox,jnearest,i,j,k,n;
                        double *pos,h,h2,hinv,hinv3,r2,rho,u,wk;
                        Ngblist = (int *) mymalloc(NumPart * sizeof(int));
                        
                        /* now add in a loop to find particles in same domain, share back the
                         momentum and energy to them (to be properly conservative) */
                        N_MIN_KERNEL=4;N_MAX_KERNEL=20;MAXITER_FB=30;NITER=0;jnearest=0;
                        startnode=All.MaxPart;dummy=0;h=0;numngb_inbox=0;pos=P[i].Pos;
                        h=PPP[i].Hsml; if(h<=0) h=All.SofteningTable[0];
                        do {
                            numngb_inbox = ngb_treefind_variable(pos,h,-1,&startnode,0,&dummy,&dummy);
                            h2=h*h; hinv=1/h; hinv3=hinv*hinv*hinv; rho=0;
                            if((numngb_inbox>=N_MIN_KERNEL)&&(numngb_inbox<=N_MAX_KERNEL))
                            {
                                jnearest=0;r2nearest=1.0e10;
                                for(n=0; n<numngb_inbox; n++)
                                {
                                    j = Ngblist[n];
                                    r2=0;for(k=0;k<3;k++) r2+=(P[i].Pos[k]-P[j].Pos[k])*(P[i].Pos[k]-P[j].Pos[k]);
                                    if((r2<r2nearest)&&(P[j].Mass>0)) {
                                        r2nearest=r2;jnearest=j;
                                    }
                                    if ((r2<=h2)&&(P[j].Mass>0)&&(SphP[j].Density>0)) {
                                        u=sqrt(r2)*hinv;
                                        wk=hinv3*kernel_wk(u);
                                        rho += (P[j].Mass*wk);
                                    }
                                } /* for(n=0; n<numngb_inbox; n++) */
                            } /* if(numngb_inbox>0) */
                            else
                            {
                                startnode=All.MaxPart;
                                if(numngb_inbox<N_MIN_KERNEL)
                                {
                                    if(numngb_inbox<=0) {
                                        h*=2.0;
                                    } else {
                                        if(NITER<=5)
                                            h*=pow((float)numngb_inbox/(float)N_MIN_KERNEL,-1/NUMDIMS);
                                        else
                                            h*=1.26; /* iterate until find appropriate > N_MIN # particles */
                                    }
                                }
                                if(numngb_inbox>N_MAX_KERNEL)
                                {
                                    if(NITER<=5)
                                        h*=pow((float)numngb_inbox/(float)N_MAX_KERNEL,-1/NUMDIMS);
                                    else
                                        h/=1.31; /* iterate until find appropriate < N_MAX # particles */
                                }
                            }
                            NITER++;
                        } while((startnode >= 0)&&(NITER<=MAXITER_FB));
                        if(jnearest != 0) {if((P[jnearest].Mass<=0)||(SphP[jnearest].Density<=0)) jnearest=0;}
                        
                        if((jnearest != 0)&&(numngb_inbox>0))
                        {
                            // share the coupled momentum and energy back to the nearby gas //
                            if(rho>0)
                            {
                                for(n=0; n<numngb_inbox; n++)
                                {
                                    j = Ngblist[n];
                                    r2=0; for(k=0;k<3;k++) r2+=(P[i].Pos[k]-P[j].Pos[k])*(P[i].Pos[k]-P[j].Pos[k]);
                                    if ((r2<=h2)&&(P[j].Mass>0)&&(SphP[j].Density>0))
                                    {
                                        u=sqrt(r2)*hinv;
                                        wk=P[j].Mass * hinv3*kernel_wk(u) / rho;
                                        degy=0;
                                        MyDouble VelPred_j[3];
                                        for(k=0;k<3;k++) {VelPred_j[k]=P[j].Vel[k];}
#ifdef SHEARING_BOX
                                        if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                                        if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                                        for(k=0; k<3; k++)
                                        {
                                            dvel=-wk*delta_mom[k]/P[j].Mass;
                                            degy-=0.5*P[j].Mass*((VelPred_j[k]+dvel)*(VelPred_j[k]+dvel) - VelPred_j[k]*VelPred_j[k]);
                                            P[j].Vel[k] += dvel;
                                            SphP[j].VelPred[k] += dvel;
                                        }
                                        degy += -wk*delta_egy; // energy 'donated' - kinetic energy change = thermal energy change
                                        // degy is now the change in thermal energy, convert to specific energy change //
                                        degy *= 1 / P[j].Mass;
                                        if(degy<-0.9*SphP[j].InternalEnergy) SphP[j].InternalEnergy*=0.1; else SphP[j].InternalEnergy += degy;
                                    }
                                } // for(n=0; n<numngb_inbox; n++)
                            } else {
                                j=jnearest; wk=1;
                                degy=0;
                                MyDouble VelPred_j[3];
                                for(k=0;k<3;k++) {VelPred_j[k]=P[j].Vel[k];}
#ifdef SHEARING_BOX
                                if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                                if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                                for(k=0; k<3; k++)
                                {
                                    dvel=-wk*delta_mom[k]/P[j].Mass;
                                    degy-=0.5*P[j].Mass*((VelPred_j[k]+dvel)*(VelPred_j[k]+dvel) - VelPred_j[k]*VelPred_j[k]);
                                    P[j].Vel[k] += dvel;
                                    SphP[j].VelPred[k] += dvel;
                                }
                                degy += -wk*delta_egy; // energy 'donated' - kinetic energy change = thermal energy change
                                // degy is now the change in thermal energy, convert to specific energy change //
                                degy *= 1 / P[j].Mass;
                                if(degy<-0.9*SphP[j].InternalEnergy) SphP[j].InternalEnergy*=0.1; else SphP[j].InternalEnergy += degy;
                            } // if(rho>0) else
                        } // closes if((jnearest != 0)&&(rho>0)) //
                        
#endif // closes GRAIN_BACKREACTION

                    } // closes check for if(v_mag > 0)
                } // closes check for if(dt > 0)
            } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    
    myfree(Ngblist);
    CPU_Step[CPU_DRAGFORCE] += measure_time();
    
}








#ifdef GRAIN_COLLISIONS

void grain_collisions(void)
{
    int i,k;
    double dt, dvel, degy, soundspeed, R_grain, t_stop, slow_fac, vel_new, delta_mom[3], delta_egy;
    int N_MAX_KERNEL,N_MIN_KERNEL,MAXITER_FB,NITER,startnode,dummy,numngb_inbox,jnearest,i,j,k,n;
    double *pos,h,h2,hinv,hinv3,r2,rho,u,wk;
    Ngblist = (int *) mymalloc(NumPart * sizeof(int));
    
    CPU_Step[CPU_MISC] += measure_time();
    if(ThisTask == 0)
    {
        printf("Beginning Grain Collisions & Interactions \n");
        fflush(stdout);
    }
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 3)
        {
            if(P[i].Grain_Density > 0)
            {
#ifndef WAKEUP
                dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
#else
                dt = P[i].dt_step * All.Timebase_interval;
#endif
                if(dt > 0)
                {
                    soundspeed = Particle_effective_soundspeed_i(i);
                    R_grain = P[i].Grain_Size; // grain size in --code-- units //
                    
                    
                } // closes check for if(dt > 0)
            } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    myfree(Ngblist);
    CPU_Step[CPU_DRAGFORCE] += measure_time();
} /* closes grain_collisions routine */




/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct grain_densdata_in
{
    MyDouble Pos[3];
    MyFloat Vel[3];
    MyFloat Hsml;
    int NodeList[NODELISTLENGTH];
}
*GrnDensDataIn, *GrnDensDataGet;

static struct grain_densdata_out
{
    MyLongDouble RhoGrains;
    MyLongDouble GrainVel[3];
}
*GrnDensDataResult, *GrnDensDataOut;


/*
 calculates density and velocity of surrounding grain particles;
 currently is the full routine allowing for MPI communications, etc;
 but have commented out the iteration to solve for the neighbor number (we just use Hsml
 as defined by the SPH; could easily make this free to iterate itself, at small cost)
 */
void grain_density(void)
{
    MyFloat *Left, *Right;
    int i, j, ndone, ndone_flag, npleft, dt_step, dummy, iter = 0;
    int ngrp, sendTask, recvTask, place, nexport, nimport;
    long long ntot;
    double fac;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
    double timecomp, timecomm, timewait;
    double tstart, tend, t0, t1;
    double desnumngb, valuenorm;
    
    CPU_Step[CPU_DENSMISC] += measure_time();
    
    Left = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(grain_density_isactive(i))
        {
            Left[i] = Right[i] = 0;
        }
    }
    
    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc(NumPart * sizeof(int));
    All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct grain_densdata_in) + sizeof(struct grain_densdata_out) +
                                             sizemax(sizeof(struct grain_densdata_in),
                                                     sizeof(struct grain_densdata_out))));
    DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));
    
    t0 = my_second();
    desnumngb = All.DesNumNgb;
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        i = FirstActiveParticle;	/* begin with this index */
        do
        {
            for(j = 0; j < NTask; j++)
            {
                Send_count[j] = 0;
                Exportflag[j] = -1;
            }
            /* do local particles and prepare export list */
            tstart = my_second();
            for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            {
                if(grain_density_isactive(i))
                {
                    if(grain_density_evaluate(i, 0, &nexport, Send_count) < 0)
                        break;
                }
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
            
            tstart = my_second();
            MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);
            tend = my_second();
            timewait1 += timediff(tstart, tend);
            
            for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
            {
                Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
                nimport += Recv_count[j];
                if(j > 0)
                {
                    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                    Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }
            GrnDensDataGet = (struct grain_densdata_in *) mymalloc(nimport * sizeof(struct grain_densdata_in));
            GrnDensDataIn = (struct grain_densdata_in *) mymalloc(nexport * sizeof(struct grain_densdata_in));
            
            /* prepare particle data for export */
            for(j = 0; j < nexport; j++)
            {
                place = DataIndexTable[j].Index;
                GrnDensDataIn[j].Pos[0] = P[place].Pos[0];
                GrnDensDataIn[j].Pos[1] = P[place].Pos[1];
                GrnDensDataIn[j].Pos[2] = P[place].Pos[2];
                GrnDensDataIn[j].Hsml = PPP[place].Hsml;
                memcpy(GrnDensDataIn[j].NodeList,
                       DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }
            /* exchange particle data */
            tstart = my_second();
            for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;
                
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&GrnDensDataIn[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct grain_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DENS_A,
                                     &GrnDensDataGet[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct grain_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm1 += timediff(tstart, tend);
            myfree(GrnDensDataIn);
            GrnDensDataResult = (struct grain_densdata_out *) mymalloc(nimport * sizeof(struct grain_densdata_out));
            GrnDensDataOut = (struct grain_densdata_out *) mymalloc(nexport * sizeof(struct grain_densdata_out));
            
            /* now do the particles that were sent to us */
            tstart = my_second();
            for(j = 0; j < nimport; j++)
                grain_density_evaluate(j, 1, &dummy, &dummy);
            tend = my_second();
            timecomp2 += timediff(tstart, tend);
            
            if(i < 0)
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
                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&GrnDensDataResult[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct grain_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DENS_B,
                                     &GrnDensDataOut[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct grain_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm2 += timediff(tstart, tend);
            
            /* add the result to the local particles */
            tstart = my_second();
            for(j = 0; j < nexport; j++)
            {
                place = DataIndexTable[j].Index;
                if(P[place].Type == 3)
                {
                    //PPP[place].NumNgb += GrnDensDataOut[j].Ngb;
                    P[place].Grain_Density += GrnDensDataOut[j].RhoGrains;
                    P[place].Grain_Velocity[0] += GrnDensDataOut[j].GrainVel[0];
                    P[place].Grain_Velocity[1] += GrnDensDataOut[j].GrainVel[1];
                    P[place].Grain_Velocity[2] += GrnDensDataOut[j].GrainVel[2];
                }
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            myfree(GrnDensDataOut);
            myfree(GrnDensDataResult);
            myfree(GrnDensDataGet);
        }
        while(ndone < NTask);
        
        
        /* do final operations on results */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(grain_density_isactive(i))
            {
                if(P[i].Grain_Density > 0)
                {
                    P[i].Grain_Velocity[0] /= P[i].Grain_Density;
                    P[i].Grain_Velocity[1] /= P[i].Grain_Density;
                    P[i].Grain_Velocity[2] /= P[i].Grain_Density;
                }
                
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            sumup_large_ints(1, &npleft, &ntot);
            
            /*
             if(ntot > 0)
             {
             iter++;
             
             if(iter > 0 && ThisTask == 0)
             {
             printf("ngb (grain) iteration %d: need to repeat for %d%09d particles.\n", iter,
             (int) (ntot / 1000000000), (int) (ntot % 1000000000));
             fflush(stdout);
             }
             
             if(iter > MAXITER)
             {
             printf("failed to converge in neighbour iteration in grain_density()\n");
             fflush(stdout);
             endrun(1155);
             }
             }
             */
        }
        while(ntot > 0);
        
        
        myfree(DataNodeList);
        myfree(DataIndexTable);
        myfree(Ngblist);
        //myfree(Right);
        //myfree(Left);
        
        /* mark as active again */
        /*
         for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
         if(P[i].TimeBin < 0)
         P[i].TimeBin = -P[i].TimeBin - 1;
         */
        
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
    
    
    
    
    
    
    /*! core of sph density computation, adapted to search for grains now */
    int grain_density_evaluate(int target, int mode, int *nexport, int *nsend_local)
    {
        int j, n;
        int startnode, numngb, numngb_inbox, listindex = 0;
        double h, h2, fac, hinv, hinv3, hinv4, wk, dwk;
        double dx, dy, dz, r, r2, u, mass_j;
        MyLongDouble sum_variable;
        MyLongDouble rho;
        MyLongDouble weighted_numngb;
        MyLongDouble gasvel[3];
        MyDouble *pos;
        MyFloat *vel;
        gasvel[0] = gasvel[1] = gasvel[2] = 0;
        rho = weighted_numngb = 0;
        
        if(mode == 0)
        {
            pos = P[target].Pos;
            h = PPP[target].Hsml;
            vel = P[target].Vel;
        }
        else
        {
            pos = GrnDensDataGet[target].Pos;
            vel = GrnDensDataGet[target].Vel;
            h = GrnDensDataGet[target].Hsml;
        }
        
        h2 = h * h;
        hinv = 1.0 / h;
#ifndef  TWODIMS
        hinv3 = hinv * hinv * hinv;
#else
        hinv3 = hinv * hinv / boxSize_Z;
#endif
        hinv4 = hinv3 * hinv;
        
        if(mode == 0)
        {
            startnode = All.MaxPart;	/* root node */
        }
        else
        {
            startnode = GrnDensDataGet[target].NodeList[0];
            startnode = Nodes[startnode].u.d.nextnode;	/* open it */
        }
        
        numngb = 0;
        while(startnode >= 0)
        {
            while(startnode >= 0)
            {
                numngb_inbox = grain_ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);
                if(numngb_inbox < 0) return -1;
                
                for(n = 0; n < numngb_inbox; n++)
                {
                    j = Ngblist[n];
                    if(P[j].Mass == 0) continue;
                    dx = pos[0] - P[j].Pos[0];
                    dy = pos[1] - P[j].Pos[1];
                    dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
                    dx = NEAREST_X(dx);
                    dy = NEAREST_Y(dy);
                    dz = NEAREST_Z(dz);
#endif
                    r2 = dx*dx + dy*dy + dz*dz;
                    if(r2 < h2)
                    {
                        numngb++;
                        r = sqrt(r2);
                        u = r * hinv;
                        wk = hinv3*kernel_wk(u);
                        dwk = hinv4*kernel_dwk(u);
                        mass_j = P[j].Mass;
                        rho += FLT(mass_j * wk);
                        weighted_numngb += FLT(NORM_COEFF * wk / hinv3);
                        MyDouble VelPred_j[3];
                        for(k=0;k<3;k++) {VelPred_j[k]=SphP[j].VelPred[k];}
#ifdef SHEARING_BOX
                        if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                        if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                        gasvel[0] += FLT(mass_j * wk * VelPred_j[0]);
                        gasvel[1] += FLT(mass_j * wk * VelPred_j[1]);
                        gasvel[2] += FLT(mass_j * wk * VelPred_j[2]);
                    }
                }
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = GrnDensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
    {
        //P[target].Grain_NumNgb = weighted_numngb;
        P[target].Grain_Density = rho;
        P[target].Grain_Velocity[0] = 0;
        P[target].Grain_Velocity[1] = 0;
        P[target].Grain_Velocity[2] = 0;
    }
    else
    {
        //GrnDensDataResult[target].Ngb = weighted_numngb;
        GrnDensDataResult[target].RhoGrains = rho;
        GrnDensDataResult[target].GrainVel[0] = 0;
        GrnDensDataResult[target].GrainVel[1] = 0;
        GrnDensDataResult[target].GrainVel[2] = 0;
    }
    
    return 0;
}



/* code to tell the grain density routine which particles to use */
int grain_density_isactive(int n)
{
    if(P[n].TimeBin < 0) return 0;
    if(P[n].Type == 3) return 1;
    return 0;
}


/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
int grain_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
                                int *nexport, int *nsend_local)
{
    int numngb, no, p, task, nexport_save;
    struct NODE *current;
    MyDouble dx, dy, dz, dist;
    
#ifdef PERIODIC
    MyDouble xtmp;
#endif
    nexport_save = *nexport;
    
    numngb = 0;
    no = *startnode;
    
    while(no >= 0)
    {
        if(no < All.MaxPart)	/* single particle */
        {
            p = no;
            no = Nextnode[no];
            
            /* this is the only line changed from gas case:: use type=3 instead of type=0 */
            if(P[p].Type != 3)
                continue;
            
            if(P[p].Ti_current != All.Ti_Current)
                drift_particle(p, All.Ti_Current);
            
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
            if(dz > dist)
                continue;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
            
            Ngblist[numngb++] = p;
        }
        else
        {
            if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
            {
                if(mode == 1)
                    endrun(12312);
                
                if(target >= 0)	/* if no target is given, export will not occur */
                {
                    if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                    {
                        Exportflag[task] = target;
                        Exportnodecount[task] = NODELISTLENGTH;
                    }
                    
                    if(Exportnodecount[task] == NODELISTLENGTH)
                    {
                        if(*nexport >= All.BunchSize)
                        {
                            *nexport = nexport_save;
                            if(nexport_save == 0)
                                endrun(13004);	/* in this case, the buffer is too small to process even a single particle */
                            for(task = 0; task < NTask; task++)
                                nsend_local[task] = 0;
                            for(no = 0; no < nexport_save; no++)
                                nsend_local[DataIndexTable[no].Task]++;
                            return -1;
                        }
                        Exportnodecount[task] = 0;
                        Exportindex[task] = *nexport;
                        DataIndexTable[*nexport].Task = task;
                        DataIndexTable[*nexport].Index = target;
                        DataIndexTable[*nexport].IndexGet = *nexport;
                        *nexport = *nexport + 1;
                        nsend_local[task]++;
                    }
                    
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                    DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                    
                    if(Exportnodecount[task] < NODELISTLENGTH)
                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }
                
                no = Nextnode[no - MaxNodes];
                continue;
            }
            
            current = &Nodes[no];
            
            if(mode == 1)
            {
                if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                {
                    *startnode = -1;
                    return numngb;
                }
            }
            
            if(current->Ti_current != All.Ti_Current)
                force_drift_node(no, All.Ti_Current);
            
            if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
            {
                if(current->u.d.mass)	/* open cell */
                {
                    no = current->u.d.nextnode;
                    continue;
                }
            }
            
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist = hsml + 0.5 * current->len;;
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
            if(dz > dist)
                continue;
            /* now test against the minimal sphere enclosing everything */
            dist += FACT1 * current->len;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
            
            no = current->u.d.nextnode;	/* ok, we need to open the node */
        }
    }
    
    *startnode = -1;
    return numngb;
}

#endif // closes GRAIN_COLLISIONS



#endif


