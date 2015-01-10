//
//  blackhole_feedback.c
//  BuildingGIZMO
//
//  Created by Paul Torrey on 10/22/14.
//
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#include "blackhole_local.h"

static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_BH_swallowed;

void blackhole_swallow_and_kick_loop(void)
{
    
    int i, j, k;
    int ndone_flag, ndone;
    int ngrp, recvTask, place, nexport, nimport, dummy;
    MPI_Status status;

#ifdef BH_GRAVACCRETION
    double m_tmp_for_bhar;
    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar,mdisk_for_bhar;
    double f0_for_bhar;
#endif
#ifdef BH_SUBGRIDBHVARIABILITY
    long nsubgridvar;
    int jsub;
    double varsg1,varsg2;
    double omega_ri,n0_sgrid_elements,norm_subgrid,time_var_subgridvar;
    gsl_rng *random_generator_forbh;
#endif
#ifdef BH_BONDI
    double norm, soundspeed, bhvel, rho;
#endif
#ifdef KD_FRICTION
    /* add a friction force for the black-holes, accounting for the environment */
    double fac_friction, relvel, accgrv, accfrc;
    double a_erf, lambda;
#endif
    

    int Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed;
    
    
    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),
                                                                     sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    N_gas_swallowed = N_star_swallowed = N_dm_swallowed = N_BH_swallowed = 0;
    Ntot_gas_swallowed = Ntot_star_swallowed = Ntot_dm_swallowed = Ntot_BH_swallowed = 0;
    
    i = FirstActiveParticle;	/* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(P[i].SwallowID == 0)     /* this particle not being swallowed */
                    if(blackhole_swallow_and_kick_evaluate(i, 0, &nexport, Send_count) < 0)
                        break;
        
        qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        BlackholeDataGet = (struct blackholedata_in *) mymalloc("BlackholeDataGet", nimport * sizeof(struct blackholedata_in));
        BlackholeDataIn = (struct blackholedata_in *) mymalloc("BlackholeDataIn", nexport * sizeof(struct blackholedata_in));
        
        
        /* populate the struct to be exported */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
           
            for(k = 0; k < 3; k++)
            {
    #if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].GradRho[k];
    #endif
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
    #ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
    #endif
    #if defined(BH_PHOTONMOMENTUM) 	|| defined(BH_BAL_WINDS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
            BlackholeDataIn[j].BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
    #endif
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
    #ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
    #else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
    #endif
            BlackholeDataIn[j].ID = P[place].ID;
            BlackholeDataIn[j].Mass = P[place].Mass;
            memcpy(BlackholeDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                    MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_DENS_A,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
        
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult", nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_swallow_and_kick_evaluate(j, 1, &dummy, &dummy);  /* set BlackholeDataResult based on BlackholeDataGet */
        
        
        if(i < 0)
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
                    MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_DENS_B,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
                }
            }
        }
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_Mass += BlackholeDataOut[j].Mass;
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_BH_Mass += BlackholeDataOut[j].BH_Mass;
    #ifdef BH_ALPHADISK_ACCRETION
            BPP(place).BH_Mass_AlphaDisk += BlackholeDataOut[j].BH_Mass_AlphaDisk;
    #endif
    #ifdef BH_BUBBLES
            BPP(place).b7.dBH_accreted_BHMass_bubbles += BlackholeDataOut[j].BH_Mass_bubbles;
    #ifdef UNIFIED_FEEDBACK
            BPP(place).b8.dBH_accreted_BHMass_radio += BlackholeDataOut[j].BH_Mass_radio;
    #endif
    #endif
            for(k = 0; k < 3; k++)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_momentum[k] += BlackholeDataOut[j].accreted_momentum[k];
    #ifdef BH_COUNTPROGS
            BPP(place).BH_CountProgs += BlackholeDataOut[j].BH_CountProgs;
    #endif
    #ifdef GALSF
            if(P[place].StellarAge > BlackholeDataOut[j].Accreted_Age)
                P[place].StellarAge = BlackholeDataOut[j].Accreted_Age;
    #endif
        }
        myfree(BlackholeDataOut);
        myfree(BlackholeDataResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);

    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    
    MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_star_swallowed, &Ntot_star_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_dm_swallowed, &Ntot_dm_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(ThisTask == 0)
    {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d BH particles\n",
               Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed);
        fflush(stdout);
    }
    
}




int blackhole_swallow_and_kick_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, bin, listindex = 0;
    MyIDType id;
    MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3];
    MyFloat *pos, h_i, bh_mass;     // hinv
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
    MyFloat *velocity, hinv, hinv3;
#endif

#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat mdot,dt;
#endif

    MyFloat dir[3],norm,mom=0;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat f_accreted;
    f_accreted=0;
    double BH_angle_weighted_kernel_sum, mom_wt;
    MyFloat theta,*Jgas_in_Kernel,BH_disk_hr,kernel_zero,dwk;
    kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1);
#if !defined(BH_ALPHADISK_ACCRETION)
    MyFloat v_kick;
#endif
#endif
#ifdef BH_BUBBLES
    MyLongDouble accreted_BH_mass_bubbles = 0;
    MyLongDouble accreted_BH_mass_radio = 0;
#endif
#ifdef GALSF
    double accreted_age = 1;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat accreted_BH_mass_alphadisk; //bh_mass_alphadisk
#endif
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
    MyFloat m_wind, e_wind;
#endif
    
    int mod_index;
    
    if(mode == 0)
    {
        pos = P[target].Pos;
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
        velocity = P[target].Vel;
#endif
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        bh_mass = BPP(target).BH_Mass;
//#ifdef BH_ALPHADISK_ACCRETION
//        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
//#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#endif

#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = P[target].GradRho;
        BH_disk_hr = P[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
#endif
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
        velocity = BlackholeDataGet[target].Vel;
#endif
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        bh_mass = BlackholeDataGet[target].BH_Mass;
//#ifdef BH_ALPHADISK_ACCRETION
//        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;       /* remove this from comm? */
//#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeDataGet[target].BH_angle_weighted_kernel_sum;
#endif
    }
    
    accreted_mass = 0;
    accreted_BH_mass = 0;
#ifdef BH_ALPHADISK_ACCRETION
    accreted_BH_mass_alphadisk = 0;
#endif
    accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
#ifdef BH_COUNTPROGS
    int accreted_BH_progs = 0;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    mom = All.BlackHoleFeedbackFactor *
    All.BlackHoleRadiativeEfficiency * mdot * dt * (C / All.UnitVelocity_in_cm_per_s);
    mom_wt = 0;
#endif
    
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
    hinv=h_i;
    hinv3=hinv*hinv*hinv;
#endif
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);
            if(numngb < 0) return -1;
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                
                /* we've found a particle to be swallowed.  This could a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == id)
                {
                    /* this is a BH-BH merger */
                    if(P[j].Type == 5)
                    {
                        fprintf(FdBlackHolesDetails,
                                "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
                                ThisTask, All.Time, id, P[j].ID, bh_mass, BPP(j).BH_Mass);
                        
                        accreted_mass    += FLT(P[j].Mass);
                        accreted_BH_mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
                        
#ifdef BH_BUBBLES
                        accreted_BH_mass_bubbles += FLT(BPP(j).BH_Mass_bubbles - BPP(j).BH_Mass_ini);
#ifdef UNIFIED_FEEDBACK
                        accreted_BH_mass_radio += FLT(BPP(j).BH_Mass_radio - BPP(j).BH_Mass_ini);
#endif
#endif
                        
#ifdef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(BPP(j).BH_Mass * P[j].Vel[k]);
#else
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
#endif
                        
                        
#ifdef BH_COUNTPROGS
                        accreted_BH_progs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin;
                        TimeBin_BH_mass[bin] -= BPP(j).BH_Mass;
                        TimeBin_BH_dynamicalmass[bin] -= P[j].Mass;
                        TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0)
                            TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;
                        P[j].Mass = 0;
                        BPP(j).BH_Mass = 0;
                        BPP(j).BH_Mdot = 0;
#ifdef BH_BUBBLES
                        BPP(j).BH_Mass_bubbles = 0;
                        BPP(j).BH_Mass_ini = 0;
#ifdef UNIFIED_FEEDBACK
                        BPP(j).BH_Mass_radio = 0;
#endif
#endif
#ifdef GALSF
                        accreted_age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5)
                    
                    
                        /* this is a DM particle:
                        In this case, no kick, so just zero out the mass and 'get rid of' the
                        particle (preferably by putting it somewhere irrelevant) */
                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {
                        printf("BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n",
                               j,P[j].Type,
                               P[j].Mass,
                               P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],
                               P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2]);
                        fflush(stdout);
                        
                        accreted_mass += FLT(f_accreted*P[j].Mass);
                        accreted_BH_mass += FLT(f_accreted*P[j].Mass);
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(f_accreted*P[j].Mass * P[j].Vel[k]);
                        P[j].Mass = 0;		// zero out particle mass.  it has now been fully swallowed.
                    }

                    
                    
                    /* finally deal with gas or star accretion events */
                    if( (P[j].Type == 0) || (P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {
#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                        printf("BAL kick: j %d Type(j) %d f_acc %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",
                                   j,P[j].Type, All.BAL_f_accretion,P[j].Mass,
                                   P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],
                                   P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2],
                                   All.BAL_v_outflow);
                        fflush(stdout);

                        
#ifdef BH_ALPHADISK_ACCRETION
                        /* all mass goes into the alpha disk, before going into the BH */
                        f_accreted=1.0;    // accrete whole particles to the alpha disk
                        accreted_mass += FLT(f_accreted*P[j].Mass);
                        accreted_BH_mass += 0;
                        accreted_BH_mass_alphadisk += FLT(f_accreted*P[j].Mass);
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(f_accreted*P[j].Mass * P[j].Vel[k]);
                        P[j].Mass = 0;		// zero out particle mass.  it has been fully swallowed.
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue = 0;
#endif  //  HYDRO_MESHLESS_FINITE_VOLUME
#else   //  BH_ALPHADISK_ACCRETION
                        if(P[j].Type==0) v_kick=All.BAL_v_outflow; else v_kick=0.1*All.BAL_v_outflow;
                        if(P[j].Type==0) f_accreted=All.BAL_f_accretion; else f_accreted=1;
                        
                        accreted_mass += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass);
                        accreted_BH_mass += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass);
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass * P[j].Vel[k]);
                        P[j].Mass *= (1-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue *= (1-f_accreted);
#endif  //HYDRO_MESHLESS_FINITE_VOLUME
                        
                        /* BAL kicking operations */
                        dir[0]=dir[1]=dir[2]=0;
                        for(k = 0; k < 3; k++) dir[k]=P[j].Pos[k]-pos[k];
                        for(k = 0, norm = 0; k < 3; k++) norm += dir[k]*dir[k];
                        if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm);}
                        for(k = 0; k < 3; k++)
                        {
                            P[j].Vel[k] += v_kick*All.cf_atime*dir[k]/norm;
                            if(P[j].Type==0) SphP[j].VelPred[k] += v_kick*All.cf_atime*dir[k]/norm;
                        }
#endif  // else BH_ALPHADISK_ACCRETION
        
                        
#else // #if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                        
                        accreted_mass += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass);
#ifdef BH_ALPHADISK_ACCRETION       /* mass goes into the alpha disk, before going into the BH */
                        accreted_BH_mass_alphadisk += FLT(P[j].Mass);
#else                               /* mass goes directly to the BH, not just the parent particle */
                        accreted_BH_mass += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass);
#endif
                        
#ifndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass * P[j].Vel[k]);
#endif
                        
                        P[j].Mass = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue = 0;
#endif
                        
#endif // #else defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                        
                        N_gas_swallowed++;
                    } // if(P[j].Type == 0)
                    
                
                    
                    
                    
                
                } // if(P[j].SwallowID == id)
                
                

                
                
                
                
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0)
                {
                    if(P[j].Type==0)
                    {
                        if((P[j].Mass>0)&&(P[j].SwallowID==0)) // not swallowed!
                        {
                            for(norm=0,k=0;k<3;k++)
                            {
                                dir[k] = (pos[k]-P[j].Pos[k]);
                                norm += dir[k]*dir[k];
                            }
                            if(norm>0)
                            {
                                norm=sqrt(norm); for(k=0;k<3;k++) dir[k]/=norm;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                                /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                for(norm=0,k=0;k<3;k++) norm += dir[k]*Jgas_in_Kernel[k];
                                theta = acos(fabs(norm));
#ifdef BH_PHOTONMOMENTUM
                                /* now we get the weight function based on what we calculated earlier */
                                mom_wt = All.BH_FluxMomentumFactor * bh_angleweight_localcoupling(j,BH_disk_hr,theta) / BH_angle_weighted_kernel_sum;
                                if(BH_angle_weighted_kernel_sum<=0) mom_wt=0;
                                /* add initial L/c optical/UV coupling to the gas at the dust sublimation radius */
                                v_kick = mom_wt * mom / P[j].Mass;
                                for(k = 0; k < 3; k++)
                                {
                                    P[j].Vel[k] += v_kick*All.cf_atime*dir[k];
                                    SphP[j].VelPred[k] += v_kick*All.cf_atime*dir[k];
                                }
#endif // BH_PHOTONMOMENTUM
#if defined(BH_BAL_WINDS) && ( (  (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))  || defined(BH_ALPHADISK_ACCRETION) ) || defined(BH_ALPHADISK_ACCRETION) )
                                /* this is where Alpha Disk feedback is set */
                                
                                mom_wt = bh_angleweight_localcoupling(j,BH_disk_hr,theta) / BH_angle_weighted_kernel_sum;
                                m_wind = mom_wt * All.BlackHoleFeedbackFactor * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * mdot*dt; /* mass to couple */
                                if(BH_angle_weighted_kernel_sum<=0) m_wind=0;
                                /* add wind mass to particle, correcting density as needed */
                                if(P[j].Hsml<=0)
                                {
                                    if(SphP[j].Density>0){SphP[j].Density*=(1+m_wind/P[j].Mass);} else {SphP[j].Density=m_wind*hinv3;}
                                } else {
                                    SphP[j].Density += kernel_zero * m_wind/(P[j].Hsml*P[j].Hsml*P[j].Hsml);
                                }
                                P[j].Mass += m_wind;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                                SphP[j].MassTrue += m_wind;
#endif
                                /* now add wind momentum to particle */
                                for(e_wind=0,k=0;k<3;k++)
                                {
                                    // relative wind-particle velocity (in code units) including BH-particle motion;
                                    norm = All.cf_atime*All.BAL_v_outflow*dir[k] + velocity[k]-P[j].Vel[k];
                                    // momentum conservation gives the following change in velocities
                                    P[j].Vel[k] += norm * m_wind/P[j].Mass;
                                    SphP[j].VelPred[k] += norm * m_wind/P[j].Mass;
                                    // and the shocked wind energy is given by
                                    e_wind += (norm/All.cf_atime)*(norm/All.cf_atime);
                                }
                                e_wind *= 0.5*m_wind;
                                /* now add wind shock energy to particle */
                                e_wind *= 1 / P[j].Mass;
                                SphP[j].InternalEnergy += e_wind;
                                SphP[j].InternalEnergyPred += e_wind;
#endif //defined(BH_BAL_WINDS) && (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))
#endif // defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                            } // norm > 0
                        } // (P[j].Mass>0)&&(P[j].SwallowID==0)
                    } // P[j].Type==0
                } // (mom>0)&&(BH_angle_weighted_kernel_sum>0)
                
                
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        BlackholeTempInfo[mod_index].accreted_Mass = accreted_mass;
        BlackholeTempInfo[mod_index].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BPP(target).BH_Mass_AlphaDisk += accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++)
            BlackholeTempInfo[mod_index].accreted_momentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
        BPP(target).b7.dBH_accreted_BHMass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
        BPP(target).b8.dBH_accreted_BHMass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
        BPP(target).BH_CountProgs += accreted_BH_progs;
#endif
#ifdef KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER
        for(k = 0; k < 3; k++)
            BlackholeTempInfo[mod_index].BH_SwallowPos[k] = bh_swallowpos[k];
#endif
#ifdef GALSF
        if(P[target].StellarAge > accreted_age)
            P[target].StellarAge = accreted_age;
#endif
    }
    else
    {
        BlackholeDataResult[target].Mass = accreted_mass;
        BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeDataResult[target].BH_Mass_AlphaDisk = accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++)
            BlackholeDataResult[target].accreted_momentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
        BlackholeDataResult[target].BH_Mass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
        BlackholeDataResult[target].BH_Mass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
        BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
#ifdef GALSF
        BlackholeDataResult[target].Accreted_Age = accreted_age;
#endif
    }
    
    return 0;
} /* closes bh_evaluate_swallow */






void blackhole_final_loop(void)
{
    int i, k, n, bin;
    double  dt;
#ifdef BH_GRAVACCRETION
    double m_tmp_for_bhar;
    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar,mdisk_for_bhar;
    double f0_for_bhar;
#endif
#ifdef BH_SUBGRIDBHVARIABILITY
    long nsubgridvar;
    int jsub;
    double varsg1,varsg2;
    double omega_ri,n0_sgrid_elements,norm_subgrid,time_var_subgridvar;
    gsl_rng *random_generator_forbh;
#endif
#ifdef BH_BONDI
    double norm, soundspeed, bhvel, rho;
#endif
#ifdef KD_FRICTION
    /* add a friction force for the black-holes, accounting for the environment */
    double fac_friction, relvel, accgrv, accfrc;
    double a_erf, lambda;
#endif


#ifdef BH_REPOSITION_ON_POTMIN
    for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
        if(P[n].Type == 5)
            if(BPP(n).BH_MinPot < 0.5 * BHPOTVALUEINIT)
                for(k = 0; k < 3; k++)
                    P[n].Pos[k] = BPP(n).BH_MinPotPos[k];
#endif
    
    
    for(n = 0; n < TIMEBINS; n++)
    {
        if(TimeBinActive[n])
        {
            TimeBin_BH_mass[n] = 0;
            TimeBin_BH_dynamicalmass[n] = 0;
            TimeBin_BH_Mdot[n] = 0;
            TimeBin_BH_Medd[n] = 0;
        }
    }
    
    
    for(i=0; i<N_active_loc_BHs; i++)
    {
        n = BlackholeTempInfo[i].index;
        if(((BlackholeTempInfo[i].accreted_Mass>0)||(BlackholeTempInfo[i].accreted_BH_Mass>0)) && P[n].Mass > 0)
        {
            for(k = 0; k < 3; k++)
            {
#ifndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                P[n].Vel[k] = (P[n].Vel[k]*P[n].Mass + BlackholeTempInfo[i].accreted_momentum[k]) / (BlackholeTempInfo[i].accreted_Mass + P[n].Mass);
#else
                P[n].Vel[k] = (P[n].Vel[k]*P[n].Mass + BlackholeTempInfo[i].accreted_momentum[k]) / (BlackholeTempInfo[i].accreted_BH_Mass + P[n].Mass);
#endif
            } //for(k = 0; k < 3; k++)
            P[n].Mass += BlackholeTempInfo[i].accreted_Mass;
            BPP(n).BH_Mass += BlackholeTempInfo[i].accreted_BH_Mass;
#ifdef BH_BUBBLES
            BPP(n).BH_Mass_bubbles += BPP(n).b7.BH_accreted_BHMass_bubbles;
#ifdef UNIFIED_FEEDBACK
            BPP(n).BH_Mass_radio += BPP(n).b8.BH_accreted_BHMass_radio;
#endif
#endif
            
        } // if(((BlackholeTempInfo[n].accreted_Mass>0)||(BlackholeTempInfo[n].accreted_BH_Mass>0)) && P[n].Mass > 0)
        
        
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
        /* now that its been added to the BH, use the actual swallowed mass to define the BHAR */
#ifndef WAKEUP
        dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[n].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        
#ifdef BH_ALPHADISK_ACCRETION
        if(dt>0)
            BPP(n).BH_Mdot += BlackholeTempInfo[i].accreted_BH_Mass / dt;
#else
        if(dt>0)
            BPP(n).BH_Mdot = BlackholeTempInfo[i].accreted_BH_Mass / dt;
        else
            BPP(n).BH_Mdot = 0;
#endif // ifdef BH_ALPHADISK_ACCRETION
#endif
        bin = P[n].TimeBin;
        TimeBin_BH_mass[bin] += BPP(n).BH_Mass;
        TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
        TimeBin_BH_Mdot[bin] += BPP(n).BH_Mdot;
        if(BPP(n).BH_Mass > 0)
            TimeBin_BH_Medd[bin] += BPP(n).BH_Mdot / BPP(n).BH_Mass;
#ifdef BH_BUBBLES
        if(BPP(n).BH_Mass_bubbles > 0 && BPP(n).BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * BPP(n).BH_Mass_ini) num_activebh++;
#endif

    } // for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n]) //
    

}

