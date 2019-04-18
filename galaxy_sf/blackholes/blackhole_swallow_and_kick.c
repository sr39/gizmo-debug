/*! \file blackhole_swallow_and_kick.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clarity by parsing the existing code into
 *   smaller files and routines. Some communication and black hole structures were modified
 *   to reduce memory usage. Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

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
    
    int Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed;
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackholedata_out))));
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
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(NEWSINK_J_FEEDBACK)
                BlackholeDataIn[j].Jsink[k] = BPP(place).Jsink[k];
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = BlackholeTempInfo[P[place].IndexMapToTempStruc].Jgas_in_Kernel[k];
#endif
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].ID = P[place].ID;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
	    BlackholeDataIn[j].ProtoStellar_Radius = BPP(place).ProtoStellar_Radius;
#endif	    
#ifdef NEWSINK
            BlackholeDataIn[j].SinkRadius = BPP(place).SinkRadius;
            //Copy info on neighbours
            BlackholeDataIn[j].n_neighbor = BlackholeTempInfo[P[place].IndexMapToTempStruc].n_neighbor;
            memcpy(BlackholeDataIn[j].rgas,BlackholeTempInfo[P[place].IndexMapToTempStruc].rgas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].xgas,BlackholeTempInfo[P[place].IndexMapToTempStruc].xgas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].ygas,BlackholeTempInfo[P[place].IndexMapToTempStruc].ygas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].zgas,BlackholeTempInfo[P[place].IndexMapToTempStruc].zgas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].mgas,BlackholeTempInfo[P[place].IndexMapToTempStruc].mgas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].Hsmlgas,BlackholeTempInfo[P[place].IndexMapToTempStruc].Hsmlgas, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].gasID,BlackholeTempInfo[P[place].IndexMapToTempStruc].gasID, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            memcpy(BlackholeDataIn[j].isbound,BlackholeTempInfo[P[place].IndexMapToTempStruc].isbound, NEWSINK_NEIGHBORMAX * sizeof(int));
            memcpy(BlackholeDataIn[j].f_acc,BlackholeTempInfo[P[place].IndexMapToTempStruc].f_acc, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
#if defined(NEWSINK_J_FEEDBACK)
            memcpy(BlackholeDataIn[j].dv_ang_kick_norm,BlackholeTempInfo[P[place].IndexMapToTempStruc].dv_ang_kick_norm, NEWSINK_NEIGHBORMAX * sizeof(MyFloat));
            BlackholeDataIn[j].t_disc = BPP(place).t_disc;
#endif
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
            BlackholeDataIn[j].BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
#endif
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
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
                                 recvTask, TAG_BH_G,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_G, MPI_COMM_WORLD, &status);
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
                                 MPI_BYTE, recvTask, TAG_BH_H,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_H, MPI_COMM_WORLD, &status);
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
            for(k = 0; k < 3; k++){
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_momentum[k] += BlackholeDataOut[j].accreted_momentum[k];
#ifdef SINGLE_STAR_STRICT_ACCRETION
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_moment[k] += BlackholeDataOut[j].accreted_moment[k];
#endif
#if defined(NEWSINK_J_FEEDBACK)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_J[k] += BlackholeDataOut[j].accreted_J[k];
#endif		
            }
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
    
    if((ThisTask == 0)&&(Ntot_gas_swallowed+Ntot_star_swallowed+Ntot_dm_swallowed+Ntot_BH_swallowed>0))
    {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d BH particles\n",
               Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed);
    }
    
}




int blackhole_swallow_and_kick_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, bin, listindex = 0;
    MyIDType id;
    MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3];
#ifdef SINGLE_STAR_STRICT_ACCRETION
    MyLongDouble accreted_moment[3];
#endif
#ifdef NEWSINK
    MyFloat f_acc_corr=1.0,mdot_avg;
    MyFloat *str_f_acc;
    int n_neighbor;
    MyIDType *str_gasID;
    MyFloat int_zone_radius;
#endif
#if defined(NEWSINK_J_FEEDBACK)
    MyLongDouble accreted_J[3];
    MyFloat dx[3], dv[3], dr;
    MyFloat Jsinktot, dJsinkpred, Jcrossdr[3];
    MyFloat *Jsink, *str_dv_ang_kick_norm;
    MyFloat tdisc;
    MyDouble dv_ang_kick_norm=0; /*Normalization factor for angular momentum feedback kicks*/ 
#endif
#if defined(NEWSINK_STOCHASTIC_ACCRETION)
    double w; int kicked=0;
#endif
#ifdef NEWSINK_JET_OPENING_ANGLE
    double phi_angle, theta_angle;
    double max_theta_angle=NEWSINK_JET_OPENING_ANGLE/2.0*M_PI/180.0; //max of theta is half the opening angles
    MyFloat reldir[3],b_vect1[3],b_vect2[3],b_vect3[3];
#endif
    MyFloat *pos, h_i, bh_mass;
#if (defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)) || defined(NEWSINK_J_FEEDBACK)
    MyFloat *velocity, hinv, hinv3;
#endif
    MyFloat f_accreted=0;
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
    MyFloat protostellar_radius;
#endif    
#if defined(NEWSINK_J_FEEDBACK) || defined(BH_WIND_KICK)
    MyFloat mass;
#ifdef BH_WIND_KICK
    MyFloat v_kick=0;
    MyFloat bh_mass_withdisk;
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk;     // DAA: we need bh_mass_alphadisk for BH_WIND_KICK winds below
#endif
#endif
#endif
#if (defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)) || defined(NEWSINK)
    MyFloat mdot,dt;
#endif
    
    MyFloat dir[3], norm, mom;
    mom=0; norm=0; dir[0]=0;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
    MyFloat *Jgas_in_Kernel;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
    double BH_angle_weighted_kernel_sum, mom_wt;
    MyFloat theta,BH_disk_hr,kernel_zero,dwk;
    kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1);
#endif
#ifdef GALSF
    double accreted_age = 1;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat accreted_BH_mass_alphadisk;   
#endif
    
    int mod_index = 0;
//printf("%d BH swallow line 285\n", ThisTask);
    if(mode == 0)
    {
        pos = P[target].Pos;
#if (defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)) || defined(NEWSINK_J_FEEDBACK)
        velocity = P[target].Vel;
#endif
        h_i = PPP[target].Hsml;
        id = P[target].ID;
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
	protostellar_radius = BPP(target).ProtoStellar_Radius;
#endif	    
#if defined(BH_WIND_KICK) || defined(NEWSINK_J_FEEDBACK)
        mass = P[target].Mass;    
#endif
#if defined(BH_ALPHADISK_ACCRETION) && defined(BH_WIND_KICK)
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        bh_mass = BPP(target).BH_Mass;
#if (defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)) || defined(NEWSINK)
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#endif
#if defined(NEWSINK)
        mdot_avg = BPP(target).BH_Mdot_Avg;
        int_zone_radius = P[target].Hsml * INT_ZONE_TO_HSML;
        n_neighbor = BlackholeTempInfo[P[target].IndexMapToTempStruc].n_neighbor;
        str_f_acc = BlackholeTempInfo[P[target].IndexMapToTempStruc].f_acc;
        str_gasID = BlackholeTempInfo[P[target].IndexMapToTempStruc].gasID;
#if defined(NEWSINK_J_FEEDBACK)
        Jsink = BPP(target).Jsink;
        Jsinktot = sqrt(Jsink[0]*Jsink[0] + Jsink[1]*Jsink[1] +Jsink[2]*Jsink[2]);
        tdisc = BPP(target).t_disc;
        str_dv_ang_kick_norm = BlackholeTempInfo[P[target].IndexMapToTempStruc].dv_ang_kick_norm;
#endif
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
        Jgas_in_Kernel = BlackholeTempInfo[P[target].IndexMapToTempStruc].Jgas_in_Kernel;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
        BH_disk_hr = P[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum;
#endif
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
#if (defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)) || defined(NEWSINK_J_FEEDBACK)
        velocity = BlackholeDataGet[target].Vel;
#endif
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
	protostellar_radius = BlackholeDataGet[target].ProtoStellar_Radius;
#endif	
#if defined(BH_WIND_KICK) || defined(NEWSINK_J_FEEDBACK)
        mass = BlackholeDataGet[target].Mass;
#if defined(BH_ALPHADISK_ACCRETION) && defined(BH_WIND_KICK)
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;      
#endif
#endif
        bh_mass = BlackholeDataGet[target].BH_Mass;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(NEWSINK)
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
#endif
#if defined(NEWSINK)
        mdot_avg = BlackholeDataGet[target].BH_Mdot_Avg;
        int_zone_radius = BlackholeDataGet[target].Hsml * INT_ZONE_TO_HSML;
        n_neighbor = BlackholeDataGet[target].n_neighbor;
        str_f_acc = BlackholeDataGet[target].f_acc;
        str_gasID = BlackholeDataGet[target].gasID;
#if defined(NEWSINK_J_FEEDBACK)
        Jsink = BlackholeDataGet[target].Jsink;
        Jsinktot = sqrt(Jsink[0]*Jsink[0] + Jsink[1]*Jsink[1] +Jsink[2]*Jsink[2]);
        tdisc = BlackholeDataGet[target].t_disc;
        str_dv_ang_kick_norm = BlackholeDataGet[target].dv_ang_kick_norm;
#endif	
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeDataGet[target].BH_angle_weighted_kernel_sum;
#endif
    }
//printf("%d BH swallow line 369\n", ThisTask);
#ifdef BH_WIND_KICK
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif
    
    accreted_mass = 0;
    accreted_BH_mass = 0;
#ifdef BH_ALPHADISK_ACCRETION
    accreted_BH_mass_alphadisk = 0;
#endif
    accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
#ifdef SINGLE_STAR_STRICT_ACCRETION
    accreted_moment[0] = accreted_moment[1] = accreted_moment[2] = 0;
#endif
#if defined(NEWSINK_J_FEEDBACK)
    accreted_J[0] = accreted_J[1] = accreted_J[2] = 0;
    if (Jsinktot>0){
        dJsinkpred = Jsinktot * (1.0 - exp(-dt/tdisc));
        /*Sum up normalization factor for angular momentum feedback*/
//printf("%d BH swallow line 388\n", ThisTask);
        for(k=0;k<n_neighbor;k++){
            if (str_f_acc[k]<1.0){
                dv_ang_kick_norm += str_dv_ang_kick_norm[k]; //we only give feedback to particles we don't swallow completely
            }
        }
    }
//printf("%d BH swallow ang_kick normalization calculated: %g  with %d neighbors \n", ThisTask, dv_ang_kick_norm, n_neighbor );
#endif
    
#ifdef BH_COUNTPROGS
    int accreted_BH_progs = 0;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
#ifdef NEWSINK
    mom = bh_lum_bol(mdot_avg, bh_mass, -1) * dt / (C / All.UnitVelocity_in_cm_per_s);
#else
    mom = bh_lum_bol(mdot, bh_mass, -1) * dt / (C / All.UnitVelocity_in_cm_per_s);
#endif
    mom_wt = 0;
#endif
    
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
    hinv=h_i; hinv3=hinv*hinv*hinv;
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
#if defined(NEWSINK)
            numngb = ngb_treefind_variable_targeted(pos, int_zone_radius, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
#else
            numngb = ngb_treefind_variable_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
#endif
            if(numngb < 0) return -1;
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n]; MyIDType OriginallyMarkedSwallowID = P[j].SwallowID; // record this to help prevent double-counting below
#if defined(NEWSINK_J_FEEDBACK)
                dx[0]=P[j].Pos[0]-pos[0]; dx[1]=P[j].Pos[1]-pos[1]; dx[2]=P[j].Pos[2]-pos[2]; dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
#endif
                /* we've found a particle to be swallowed.  This could be a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == id && P[j].Mass > 0)
                {
#ifndef IO_REDUCED_MODE
                    printf("found particle P[j].ID = %llu with P[j].SwallowID = %llu of type P[j].Type = %d nearby id = %llu with P[j].Mass=%g\n",
                           (unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID, P[j].Type, (unsigned long long) id, P[j].Mass);
#endif
                    if(P[j].Type == 5)  /* this is a BH-BH merger */
                    {
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhMergerDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f\n", All.Time,  id,bh_mass,pos[0],pos[1],pos[2],  P[j].ID,BPP(j).BH_Mass,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]);
#else
#ifndef IO_REDUCED_MODE
                        fprintf(FdBlackHolesDetails,"ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n", ThisTask, All.Time, id, P[j].ID, bh_mass, BPP(j).BH_Mass);
#endif
#endif

#ifdef BH_INCREASE_DYNAMIC_MASS
                        /* the true dynamical mass of the merging BH is P[j].Mass/BH_INCREASE_DYNAMIC_MASS unless exceeded by physical growth
                         - in the limit BPP(j).BH_Mass > BH_INCREASE_DYNAMIC_MASS x m_b, then bh_mass=P[j].Mass on average and we are good as well  */
                        accreted_mass    += FLT( DMAX(BPP(j).BH_Mass, P[j].Mass/BH_INCREASE_DYNAMIC_MASS) );
#else
                        accreted_mass    += FLT(P[j].Mass);
#endif
                        accreted_BH_mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
                        for(k = 0; k < 3; k++){accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);}
#if defined(NEWSINK_J_FEEDBACK)
                        dv[0]=BPP(j).Vel[0]-velocity[0];dv[1]=BPP(j).Vel[1]-velocity[1];dv[2]=BPP(j).Vel[2]-velocity[2];
                        accreted_J[0] += FLT(P[j].Mass *(dx[1]*dv[2] - dx[2]*dv[1]) + BPP(j).Jsink[0]);
                        accreted_J[1] += FLT(P[j].Mass *(dx[2]*dv[0] - dx[0]*dv[2]) + BPP(j).Jsink[1]);
                        accreted_J[2] += FLT(P[j].Mass *(dx[0]*dv[1] - dx[1]*dv[0]) + BPP(j).Jsink[2]);
#endif
#ifdef BH_COUNTPROGS
                        accreted_BH_progs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin; TimeBin_BH_mass[bin] -= BPP(j).BH_Mass; TimeBin_BH_dynamicalmass[bin] -= P[j].Mass; TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0) {TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;}
                        P[j].Mass = 0; BPP(j).BH_Mass = 0; BPP(j).BH_Mdot = 0;
#ifdef GALSF
                        accreted_age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5) -- BH + BH merger


#ifdef BH_GRAVCAPTURE_NONGAS /* DM and star particles can only be accreted ifdef BH_GRAVCAPTURE_NONGAS */
                    /* this is a DM particle: In this case, no kick, so just zero out the mass and 'get rid of' the particle (preferably by putting it somewhere irrelevant) */
                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {
#ifndef IO_REDUCED_MODE
                        printf("BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n", j,P[j].Type,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2]);
#endif
                        accreted_mass += FLT(P[j].Mass); accreted_BH_mass += FLT(P[j].Mass);
                        P[j].Mass = 0;		// zero out particle mass.  it has now been fully swallowed.
                        N_dm_swallowed++;
                    }

                    /* this is a star particle: If there is an alpha-disk, we let them go to the disk. If there is no alpha-disk, stars go to the BH directly and won't affect feedback. (Can be simply modified if we need something different.) */
                    if((P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {
                        accreted_mass += FLT(P[j].Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(P[j].Mass);
#else 
                        accreted_BH_mass += FLT(P[j].Mass);   /* mass goes directly to the BH, not just the parent particle */
#endif
                        P[j].Mass = 0;          // zero out particle mass.  it has now been fully swallowed.
                        N_star_swallowed++;
                    }
#endif // #ifdef BH_GRAVCAPTURE_NONGAS -- BH + DM or Star merger



                    /* this is a gas particle: DAA: we need to see if the gas particle has to be accreted in full or not, depending on BH_WIND_KICK
                     the only difference with BH_ALPHADISK_ACCRETION should be that the mass goes first to the alphadisk */
                    if(P[j].Type == 0)                    
                    {
#ifdef BH_WIND_KICK
                        f_accreted = All.BAL_f_accretion;
#ifndef BH_GRAVCAPTURE_GAS
                        if((All.BlackHoleFeedbackFactor > 0) && (All.BlackHoleFeedbackFactor != 1.)) {f_accreted /= All.BlackHoleFeedbackFactor;} else {if(All.BAL_v_outflow > 0) f_accreted = 1./(1. + fabs(1.*BH_WIND_KICK)*All.BlackHoleRadiativeEfficiency*(C/All.UnitVelocity_in_cm_per_s)/(All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s));}
                        if((bh_mass_withdisk - mass) <= 0) {f_accreted=0;} // DAA: no need to accrete gas particle to enforce mass conservation (we will simply kick),  note that here the particle mass P.Mass is larger than the physical BH mass P.BH_Mass
#endif // #ifdef BH_GRAVCAPTURE_GAS
#else // #ifdef BH_WIND_KICK
                        f_accreted = 1;                           // DAA: no "kick winds" so we need to accrete gas particle in full
#endif

#if defined(NEWSINK)
/*Only take a portion of the mass if we are accreting more than needed*/
                        f_acc_corr = 1.0;
                        for(k=0;k<n_neighbor;k++){ /*Find the accretion factor we prescribed for this particle from list*/
                            if( P[j].ID == str_gasID[k]){
                                f_acc_corr = DMIN( str_f_acc[k], 1.0);
                                if (f_acc_corr < 0) {f_acc_corr=0;}
#if !defined(NEWSINK_STOCHASTIC_ACCRETION)
                                else {if ((1.0-f_acc_corr) < 1e-2) {f_acc_corr=1.0;} //failsafe for weird numerical issues
                                     else {f_accreted *= f_acc_corr;} //change accretion fraction if needed
                                }
#endif
                            }
                        }
#if defined(NEWSINK_STOCHASTIC_ACCRETION) //In this case we stochastically decide whether to accrete the entire particle
                        w = get_random_number(P[j].ID);
                        if(w < f_acc_corr){ f_accreted=1.0;} //this means we fully accrete the particle
                        else{ f_accreted=0.0; } //we don't take this particle
#endif

#ifdef BH_OUTPUT_MOREINFO
                        if ((f_acc_corr != 1.0) && (f_acc_corr != 0.0)) {printf("n=%llu f_acc_corr is: %g for particle with id %llu and mass %g around BH with id %llu\n", (unsigned long long) target, (MyFloat) f_acc_corr,(unsigned long long) P[j].ID, P[j].Mass,(unsigned long long) id);}
#endif
#endif
                        if (f_accreted>0.0){
#if defined(NEWSINK_STOCHASTIC_ACCRETION) && defined(BH_WIND_KICK) //We stochastically determine if this "accreted" particle is really accreted and we take its mass or it gets kicked out
                            w = get_random_number(P[j].ID); kicked=0;
                            if(w > All.BAL_f_accretion){
                                kicked=1;f_accreted=0.0;
                            }
                            else{
#endif
                                accreted_mass += FLT(f_accreted*P[j].Mass);
#ifdef BH_GRAVCAPTURE_GAS
#ifdef BH_ALPHADISK_ACCRETION       /* mass goes into the alpha disk, before going into the BH */
                                accreted_BH_mass_alphadisk += FLT(f_accreted*P[j].Mass);
#else                               /* mass goes directly to the BH, not just the parent particle */
                                accreted_BH_mass += FLT(f_accreted*P[j].Mass);
#endif
                                for(k = 0; k < 3; k++){
                                    accreted_momentum[k] += FLT(f_accreted * P[j].Mass * P[j].Vel[k]);
#ifdef SINGLE_STAR_STRICT_ACCRETION
                                    accreted_moment[k] += FLT(f_accreted * P[j].Mass * P[j].Pos[k]);
#endif
                                }
#if defined(NEWSINK_J_FEEDBACK)
                                dv[0]=P[j].Vel[0]-velocity[0];dv[1]=P[j].Vel[1]-velocity[1];dv[2]=P[j].Vel[2]-velocity[2];
                                accreted_J[0] += FLT(f_accreted * P[j].Mass *(dx[1]*dv[2] - dx[2]*dv[1]) + P[j].Jsink[0]);
                                accreted_J[1] += FLT(f_accreted * P[j].Mass *(dx[2]*dv[0] - dx[0]*dv[2]) + P[j].Jsink[1]);
                                accreted_J[2] += FLT(f_accreted * P[j].Mass *(dx[0]*dv[1] - dx[1]*dv[0]) + P[j].Jsink[2]);
#endif				
/* #ifdef NEWSINK_B_FEEDBACK */
/* 				if(f_accreted == 1.0){ // if the particle is still around after then we leave the flux alone */
/* 				  for(k=0;k<3;k++)   accreted_B[k] += SphP[i].B[k]; */
/* #ifdef DIVBCLEANING_DEDNER */
/* 				  accreted_Phi += SphP[i].Phi; */
/* #endif */
//				}
//#endif				
#endif
                                P[j].Mass *= (1.0-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                                SphP[j].MassTrue *= (1.0-f_accreted);
#endif
#ifdef BH_OUTPUT_MOREINFO
                                if ((1.0-f_accreted)>0) {printf("f_accreted is: %g for particle with id %llu and mass %g around BH with id %llu\n", (MyFloat) f_accreted,(unsigned long long) P[j].ID, P[j].Mass,(unsigned long long) id);}
                                else{printf("Particle with id %llu and mass %g swallowed by BH with id %llu\n", (unsigned long long) P[j].ID, P[j].Mass,(unsigned long long) id);}
#endif
#if defined(NEWSINK_STOCHASTIC_ACCRETION) && defined(BH_WIND_KICK)
                            }//end of else for determining if the particle is kicked
#endif
                            
#if defined(NEWSINK_STOCHASTIC_ACCRETION) //check if we actually kick this particle in the stochastic case
                            if (kicked){
#endif
#ifdef BH_WIND_KICK     /* BAL kicking operations. NOTE: we have two separate BAL wind models, particle kicking and smooth wind model. This is where we do the particle kicking BAL model. This should also work when there is alpha-disk. */
                                v_kick=All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s; //if( !(All.ComovingIntegrationOn) && (All.Time < 0.001)) {v_kick *= All.Time/0.001;}
#ifdef SINGLE_STAR_PROTOSTELLAR_EVOLUTION
				v_kick = sqrt(All.G * mass / (protostellar_radius * 6.957e10 / All.UnitLength_in_cm)); // Kepler velocity at the protostellar radius. Really we'd want v_kick = v_kep * m_accreted / m_kicked to get the right momentum
#endif 
#if defined(NEWSINK) && !defined(NEWSINK_STOCHASTIC_ACCRETION) /*It is possible to accrete only part of the particle so we need to be more careful about our kicks*/
                                if (f_acc_corr<1.0){
                                    v_kick *= f_acc_corr*(1.0-All.BAL_f_accretion)/(1.0-All.BAL_f_accretion*f_acc_corr); /*we wanted to only accrete an f_acc_corr portion, so the imparted momentum is proportional to only f_acc_corr*(1-All.BAL_f_accretion) times the initial mass*/
                                }
#endif
                                dir[0]=dir[1]=dir[2]=0; for(k=0;k<3;k++) {dir[k]=P[j].Pos[k]-pos[k];} // DAA: default direction is radially outwards
#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
                                double dEcr = All.BH_CosmicRay_Injection_Efficiency * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s);
                                SphP[j].CosmicRayEnergy+=dEcr; SphP[j].CosmicRayEnergyPred+=dEcr;
#ifdef COSMIC_RAYS_M1
                                dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]+=dEcr*dir[k]; SphP[j].CosmicRayFluxPred[k]+=dEcr*dir[k];}
#endif
#endif
#if (BH_WIND_KICK < 0)  /* DAA: along polar axis defined by angular momentum within Kernel (we could add finite opening angle) work out the geometry w/r to the plane of the disk */
#if defined(NEWSINK_J_FEEDBACK) /*Use Jsink instead of Jgas_in_Kernel for direction*/
                                if((dir[0]*Jsink[0] + dir[1]*Jsink[1] + dir[2]*Jsink[2]) > 0){for(k=0;k<3;k++) {dir[k]=Jsink[k];}} else {for(k=0;k<3;k++) {dir[k]=-Jsink[k];}}
#else
                                if((dir[0]*Jgas_in_Kernel[0] + dir[1]*Jgas_in_Kernel[1] + dir[2]*Jgas_in_Kernel[2]) > 0){for(k=0;k<3;k++) {dir[k]=Jgas_in_Kernel[k];}} else {for(k=0;k<3;k++) {dir[k]=-Jgas_in_Kernel[k];}}
#endif
#endif
                                for(k=0,norm=0;k<3;k++) {norm+=dir[k]*dir[k];} if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm); dir[0]/=norm;dir[1]/=norm;dir[2]/=norm;}
#if defined(NEWSINK_JET_OPENING_ANGLE) //get the new relative position vector for the particle velocity (from sink)					
                                theta_angle = max_theta_angle * get_random_number(P[j].ID); //uniformly chosen
                                phi_angle=acos(1.0 - 2.0 * get_random_number(P[j].ID)); //chosen in a way to get a uniform distribution on the spherical surface
                                reldir[0]=cos(phi_angle) * sin(theta_angle); reldir[1]=sin(phi_angle) * sin(theta_angle); reldir[2]=cos(theta_angle); //get relative direction from polar axis      
                                //Let's get the other base vectors and get the new velocity direction for the particle. 
                                b_vect3[0]=dir[0];b_vect3[1]=dir[1];b_vect3[2]=dir[2];
                                b_vect1[0] = 0.0; b_vect1[1] = dir[2]; b_vect1[2] = - dir[1]; //We get the first base by taking cross product of dir with +x unit vector
                                for(k=0,norm=0;k<3;k++) {norm+=b_vect1[k]*b_vect1[k];} if(norm<=0) {b_vect1[0]=0;b_vect1[1]=1.0;b_vect1[2]=0;norm=1;} else {norm=sqrt(norm);b_vect1[0]/=norm;b_vect1[1]/=norm;b_vect1[2]/=norm;}
                                //second vector is dir cross b_vect1, and it should be normalized by default as it is the cross product of two orthogonal vectors
                                b_vect2[0] = b_vect3[1] * b_vect1[2] - b_vect3[2] * b_vect1[1]; 
                                b_vect2[1] = b_vect3[0] * b_vect1[2] - b_vect3[2] * b_vect1[0]; 
                                b_vect2[2] = b_vect3[0] * b_vect1[1] - b_vect3[1] * b_vect1[0];
                                //Now we get the new direction
                                for(k=0;k<3;k++) {dir[k]=reldir[0]*b_vect1[k]+reldir[1]*b_vect2[k]+reldir[2]*b_vect3[k];}
#if defined(NEWSINK_RELOCATE_KICKED_PARTICLE)
                                //Let's reposition the particle
                                for(k=0;k<3;k++) {P[j].Pos[k] = pos[k] + dir[k]*int_zone_radius;}//Put the particle at the edge of the interaction zone
#endif
#endif

                                for(k=0;k<3;k++) {P[j].Vel[k]+=v_kick*All.cf_atime*dir[k]; SphP[j].VelPred[k]+=v_kick*All.cf_atime*dir[k];}				
#ifdef NEWSINK
                                for(k=0;k<3;k++) {accreted_momentum[k] -= P[j].Mass * v_kick * All.cf_atime * dir[k]; } // To conserve momentum
#endif				
#ifdef GALSF_SUBGRID_WINDS // if sub-grid galactic winds are decoupled from the hydro, we decouple the BH kick winds as well
                                SphP[j].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
#endif  

#ifndef IO_REDUCED_MODE
                                //printf("BAL kick: P[j].ID %llu BH ID dir: %g %g %g, Jsink: %g %g %g reldir: %g %g %g bvect3 %g %g %g \n", (unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID, dir[0],dir[1],dir[2],Jsink[0],Jsink[1],Jsink[2], reldir[0],reldir[1],reldir[2], b_vect3[0],b_vect3[1],b_vect3[2]);
                                printf("BAL kick: All.BAL_v_outflow %g \t f_acc_corr %g \t v_kick %g\n",(All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s),f_acc_corr,v_kick);
                                printf("BAL kick: P[j].ID %llu BH ID %llu Type(j) %d All.BAL_f_accretion %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",
                                       (unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID,P[j].Type, All.BAL_f_accretion,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2],v_kick);
#endif
#ifdef BH_OUTPUT_MOREINFO
                                fprintf(FdBhWindDetails,"%g  %u %g  %2.7f %2.7f %2.7f  %2.7f %2.7f %2.7f  %g %g %g  %u  %2.7f %2.7f %2.7f\n",
                                        All.Time, P[j].ID, P[j].Mass,  P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],  P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],dir[0],dir[1],dir[2], id, pos[0],pos[1],pos[2]);
#endif
#endif   // #ifdef BH_WIND_KICK
#if defined(NEWSINK_STOCHASTIC_ACCRETION) //continuation of the if (kicked) statement
                            }
                            else{N_gas_swallowed++;} //only count it s swallowed if it actually is
#else
                            N_gas_swallowed++;
#endif //defined(NEWSINK_STOCHASTIC_ACCRETION)
                        } // f_accreted>0.0
                    }  // if(P[j].Type == 0)

                    /* DAA: make sure it is not accreted (or ejected) by the same BH again if inactive in the next timestep */
                    P[j].SwallowID = 0; 
                } // if(P[j].SwallowID == id)  -- particles being entirely or partially swallowed!!!
#if defined(NEWSINK_J_FEEDBACK)
		int n;
                if( Jsinktot > 0 && P[j].Mass > 0 && P[j].Type == 0 ){ /*There is angular mom in the sink and this is gas*/
                /*Let's find if it is on the neighbor list*/
                    for(n=0;n<n_neighbor;n++){
                        if( P[j].ID == str_gasID[n] && str_f_acc[n] < 1.0 ){ /*It should be a particle we don't swallow fully*/
                            Jcrossdr[0] = -Jsink[2]*dx[1] + Jsink[1]*dx[2]; Jcrossdr[1] = Jsink[2]*dx[0] - Jsink[0]*dx[2]; Jcrossdr[2] = -Jsink[1]*dx[0] + Jsink[0]*dx[1]; // L x dx cross product
			    for(k=0; k<3; k++) {
			        dv[k] = dJsinkpred * Jcrossdr[k] / dv_ang_kick_norm; 
                                P[j].Vel[k] += dv[k];  //Eq 22 in Hubber 2013
				SphP[j].VelPred[k] += dv[k]; 
				accreted_momentum[k] -= dv[k]*P[j].Mass; // to conserve momentum
			    }
			    accreted_J[0] -= (dv[2]*dx[1] - dv[1]*dx[2])*P[j].Mass; accreted_J[1] -= (-dv[2]*dx[0] + dv[0]*dx[2])*P[j].Mass; accreted_J[2] -= (dv[1]*dx[0] - dv[0]*dx[1])*P[j].Mass;
                        }
                    }
// #ifdef BH_OUTPUT_MOREINFO
                    // printf("swallow n=%llu last Jsink[2] is: %g \n", (unsigned long long) target, (MyFloat) Jsink[2]);
                    // printf("swallow n=%llu last dv_ang_kick_norm is: %g \n", (unsigned long long) target, (MyFloat) dv_ang_kick_norm);
                    // printf("swallow n=%llu last Jsinktot is: %g \n", (unsigned long long) target, (MyFloat) Jsinktot);
                    // printf("swallow n=%llu last dJsinkpred is: %g \n", (unsigned long long) target, (MyFloat) dJsinkpred);
                    // printf("swallow n=%llu last Jcrossdr[0] is: %g \n", (unsigned long long) target, (MyFloat) Jcrossdr[0]);
                    // printf("swallow n=%llu last accreted_J[0] is: %g \n", (unsigned long long) target, (MyFloat) accreted_J[0]);
                    // printf("swallow n=%llu P[j].Vel[0] is: %g \n", (unsigned long long) target, (MyFloat) P[j].Vel[0]);
                    // printf("swallow n=%llu P[j].Vel[1] is: %g \n", (unsigned long long) target, (MyFloat) P[j].Vel[1]);
                    // printf("swallow n=%llu P[j].Vel[2] is: %g \n", (unsigned long long) target, (MyFloat) P[j].Vel[2]);
                    // printf("swallow n=%llu dv[0] is: %g \n", (unsigned long long) target, (MyFloat) dv[0]);
                    // printf("swallow n=%llu dv[1] is: %g \n", (unsigned long long) target, (MyFloat) dv[1]);
                    // printf("swallow n=%llu dv[2] is: %g \n", (unsigned long long) target, (MyFloat) dv[2]);
                    // printf("swallow n=%llu P[j].Mass is: %g \n", (unsigned long long) target, (MyFloat) P[j].Mass);
                    // printf("swallow n=%llu last sink mass is: %g \n", (unsigned long long) target, (MyFloat) mass);
// #endif
                }
#endif
                
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)                
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0 && mdot>0 && dt>0 && OriginallyMarkedSwallowID==0 && P[j].SwallowID==0 && P[j].Mass>0 && P[j].Type==0) // particles NOT being swallowed!
                {
                    double r=0; for(k=0;k<3;k++) {dir[k]=P[j].Pos[k]-pos[k]; r+=dir[k]*dir[k];} // should be away from BH
                    if(r>0)
                            {
                        r=sqrt(r); for(k=0;k<3;k++) {dir[k]/=r;} /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                        for(norm=0,k=0;k<3;k++) {norm+=dir[k]*Jgas_in_Kernel[k];}
                                theta = acos(fabs(norm));
                                /* now we get the weight function based on what we calculated earlier */
                        mom_wt = bh_angleweight_localcoupling(j,BH_disk_hr,theta,r,h_i) / BH_angle_weighted_kernel_sum;
                                if(BH_angle_weighted_kernel_sum<=0) mom_wt=0;
                                
#ifdef BH_PHOTONMOMENTUM /* inject radiation pressure: add initial L/c optical/UV coupling to the gas at the dust sublimation radius */
                        double v_kick = All.BH_FluxMomentumFactor * mom_wt * mom / P[j].Mass;
                        for(k=0;k<3;k++) {P[j].Vel[k]+=v_kick*All.cf_atime*dir[k]; SphP[j].VelPred[k]+=v_kick*All.cf_atime*dir[k];}
#endif
#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside continuous wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * mom_wt * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s) * mdot*dt;
                                SphP[j].CosmicRayEnergy+=dEcr; SphP[j].CosmicRayEnergyPred+=dEcr;
#ifdef COSMIC_RAYS_M1
                                dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]+=dEcr*dir[k]; SphP[j].CosmicRayFluxPred[k]+=dEcr*dir[k];}
#endif
#endif
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK) /* inject BAL winds, this is the more standard smooth feedback model */
                                double m_wind = mom_wt * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * mdot*dt; /* mass to couple */
                                if(BH_angle_weighted_kernel_sum<=0) m_wind=0;
                                
//1. check if (Vw-V0)*rhat <= 0   [ equivalently, check if   |Vw| <= V0*rhat ]
//2. if (1) is False, the wind will catch the particle, couple mass, momentum, energy, according to the equations above
//3. if (1) is True, the wind will not catch the particle, or will only asymptotically catch it. For the sake of mass conservation in the disk, I think it is easiest to treat this like the 'marginal' case where the wind barely catches the particle. In this case, add the mass normally, but no momentum, and no energy, giving:
                        //dm = m_wind, dV = 0, du = -mu*u0   [decrease the thermal energy slightly to account for adding more 'cold' material to it]
                                
                                double dvr_gas_to_bh, dr_gas_to_bh;
                                for(dvr_gas_to_bh=dr_gas_to_bh=0, k=0;k<3;k++)
                                {
                                    dvr_gas_to_bh += (velocity[k]-P[j].Vel[k]) * (pos[k]-P[j].Pos[k]);
                                    dr_gas_to_bh  += (pos[k]-P[j].Pos[k]) * (pos[k]-P[j].Pos[k]);
                                }
                                dvr_gas_to_bh /= dr_gas_to_bh ;
                                
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
                                if(dvr_gas_to_bh < (All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s))   // gas moving away from BH at v < BAL speed
                                {
                                    double e_wind = 0;
                                    for(k=0;k<3;k++)
                                    {
                                norm = All.cf_atime*(All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s)*dir[k] + velocity[k]-P[j].Vel[k]; // relative wind-particle velocity (in code units) including BH-particle motion;
                                P[j].Vel[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass; // momentum conservation gives updated velocity
                                        SphP[j].VelPred[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass;
                                e_wind += (norm/All.cf_atime)*(norm/All.cf_atime); // -specific- shocked wind energy
                                    }
                            e_wind *= 0.5*m_wind/P[j].Mass; // make total wind energy, add to particle as specific energy of -particle-
                            SphP[j].InternalEnergy += e_wind; SphP[j].InternalEnergyPred += e_wind;
                        } else {    // gas moving away from BH at wind speed (or faster) already.
                                    if(SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass > 0)
                                        SphP[j].InternalEnergy = SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass;
                                }
#endif // if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
                            } // norm > 0
                } // (check if valid gas neighbor of interest)
#endif // defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)                
                

            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}	/* open it */
            }
        }
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        BlackholeTempInfo[mod_index].accreted_Mass = accreted_mass;
        BlackholeTempInfo[mod_index].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        // DAA: could be better to include this in BlackholeTempInfo and update BH_Mass_AlphaDisk only at the end (like Mass and BH_Mass)
        BPP(target).BH_Mass_AlphaDisk += accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++) {
            BlackholeTempInfo[mod_index].accreted_momentum[k] = accreted_momentum[k];
#ifdef SINGLE_STAR_STRICT_ACCRETION
	    BlackholeTempInfo[mod_index].accreted_moment[k] = accreted_moment[k];
#endif
#if defined(NEWSINK_J_FEEDBACK)
            BlackholeTempInfo[mod_index].accreted_J[k] = accreted_J[k];
#endif	    
        }
#ifdef BH_COUNTPROGS
        BPP(target).BH_CountProgs += accreted_BH_progs;
#endif
#ifdef GALSF
        if(P[target].StellarAge > accreted_age) {P[target].StellarAge = accreted_age;}
#endif
    }
    else
    {
        BlackholeDataResult[target].Mass = accreted_mass;
        BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeDataResult[target].BH_Mass_AlphaDisk = accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++) {
            BlackholeDataResult[target].accreted_momentum[k] = accreted_momentum[k];
#ifdef SINGLE_STAR_STRICT_ACCRETION
            BlackholeDataResult[target].accreted_moment[k] = accreted_moment[k];	    
#endif	    
#if defined(NEWSINK_J_FEEDBACK)
            BlackholeDataResult[target].accreted_J[k] = accreted_J[k];
#endif
        }
#ifdef BH_COUNTPROGS
        BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
#ifdef GALSF
        BlackholeDataResult[target].Accreted_Age = accreted_age;
#endif
    }
    
    return 0;
} /* closes bh_evaluate_swallow */



#ifdef BH_WIND_SPAWN
void spawn_bh_wind_feedback(void)
{
    int i, n_particles_split = 0, MPI_n_particles_split, dummy_gas_tag=0;
    for(i = 0; i < NumPart; i++)
        if(P[i].Type==0)
        {
            dummy_gas_tag=i;
            break;
        }
    
    /* don't loop or go forward if there are no gas particles in the domain, or the code will crash */
    if(dummy_gas_tag >= 0)
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type ==5)
            {
#ifndef IO_REDUCED_MODE
                printf("attempting to spawn feedback particles for BH %d on Task %d \n", i, ThisTask);
#endif
                n_particles_split += blackhole_spawn_particle_wind_shell( i , dummy_gas_tag);
            }
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifndef IO_REDUCED_MODE
    if(ThisTask == 0) {printf("Particle BH spawn check: %d particles spawned \n", MPI_n_particles_split);}
#endif
    /* rearrange_particle_sequence -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas   += (long long)MPI_n_particles_split;
    Gas_split       = n_particles_split;                    // specific to the local processor //
    
    rearrange_particle_sequence();
}




/*! this code copies what was used in merge_split.c for the gas particle split case */
int blackhole_spawn_particle_wind_shell( int i, int dummy_sph_i_to_clone )
{
#ifndef IO_REDUCED_MODE
    printf(" splitting BH %d using SphP particle %d\n", i, dummy_sph_i_to_clone);
#endif
    double mass_of_new_particle, total_mass_in_winds, dt;
    int n_particles_split, bin; long j;
    
#ifndef WAKEUP
    dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
    dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
    
    /* here is where the details of the split are coded, the rest is bookkeeping */
    total_mass_in_winds = BPP(i).unspawned_wind_mass;
    n_particles_split   = floor( total_mass_in_winds / All.BAL_wind_particle_mass );
    if( (n_particles_split == 0) || (n_particles_split < 1) ) {return 0;}
    mass_of_new_particle = total_mass_in_winds / n_particles_split;
#ifndef IO_REDUCED_MODE
    printf("want to create %g mass in wind with %d new particles each of mass %g \n", total_mass_in_winds, n_particles_split, mass_of_new_particle);
#endif
    if(NumPart + n_particles_split >= All.MaxPart)
    {
        printf ("On Task=%d with NumPart=%d we tried to split a particle, but there is no space left...(All.MaxPart=%d). Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, All.MaxPart);
        fflush(stdout); endrun(8888);
    }
    
    int k=0;
    double phi = 2.0*M_PI*get_random_number(i+1+ThisTask); // random from 0 to 2pi //
    double cos_theta = 2.0*(get_random_number(i+3+2*ThisTask)-0.5); // random between 1 to -1 //
    double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0]);
#endif
    d_r = DMIN(0.0001, d_r);
    
    for (bin = 0; bin < TIMEBINS; bin++) {if (TimeBinCount[bin] > 0) break;}
    
    /* find the first non-gas particle and move it to the end of the particle list */
    for(j = NumPart; j < NumPart + n_particles_split; j++)
    {
        // i is the BH particle tag
        // j is the new "spawed" particle's location
        // dummy_sph_i_to_clone is a dummy SPH particle's tag to be used to init the wind particle
        k=0;
        phi = 2.0*M_PI*get_random_number(j+1+ThisTask); // random from 0 to 2pi //
        cos_theta = 2.0*(get_random_number(j+3+2*ThisTask)-0.5); // random between 1 to -1 //
        double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // epsilon*Hsml; epsilon<<1, to maintain stability //
        d_r = DMIN(0.0001, d_r);

        /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
        P[j]    = P[dummy_sph_i_to_clone];
        SphP[j] = SphP[dummy_sph_i_to_clone];
        P[j].TimeBin = bin;            // put this particle on the lowest active time bin
        P[j].dt_step = bin ? (((integertime) 1) << bin) : 0;
        
        /* the particle needs to be 'born active' and added to the active set */
        NextActiveParticle[j] = FirstActiveParticle;
        FirstActiveParticle = j;
        NumForceUpdate++;
        
        /* likewise add it to the counters that register how many particles are in each timebin */
        TimeBinCount[bin]++;
        TimeBinCountSph[bin]++;
        PrevInTimeBin[j] = i;
        
        if(FirstInTimeBin[bin] < 0){  // only particle in this time bin on this task
            FirstInTimeBin[bin] = j;
            LastInTimeBin[bin] = j;
            NextInTimeBin[j] = -1;
            PrevInTimeBin[j] = -1;
        } else {                      // there is already at least one particle; add this one "to the front" of the list
            NextInTimeBin[j] = FirstInTimeBin[bin];
            PrevInTimeBin[j] = -1;
            PrevInTimeBin[FirstInTimeBin[bin]] = j;
            FirstInTimeBin[bin] = j;
        }
        
        /* the particle needs an ID: we give it a bit-flip from the original particle to signify the split */
        unsigned int bits;
        int SPLIT_GENERATIONS = 4;
        for(bits = 0; SPLIT_GENERATIONS > (1 << bits); bits++);
        /* correction:  We are using a fixed wind ID, to allow for trivial wind particle identification */
        P[j].ID = All.AGNWindID;
        
        /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
        SphP[j].ConditionNumber *= 10.0;

        SphP[j].Density *= 1e-10; /* will be re-generated anyways */
        SphP[j].Pressure *= 1e-10; /* will be re-generated anyways */
        P[j].Hsml = All.SofteningTable[0]; /* will be re-generated anyways */
        PPP[j].Hsml = All.SofteningTable[0]; /* will be re-generated anyways */
        
        SphP[j].InternalEnergy = All.BAL_internal_temperature / (  PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g  );
        SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;
        
        /* this is a giant pile of variables to zero out. dont need everything here because we cloned a valid particle, but handy anyways */
        for(k=0;k<3;k++) SphP[j].HydroAccel[k] = 0;
        P[i].Particle_DivVel = 0; SphP[j].DtInternalEnergy = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[j].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].dMass = 0; SphP[j].DtMass = 0; SphP[j].MassTrue = P[j].Mass; for(k=0;k<3;k++) SphP[j].GravWorkTerm[k] = 0;
#endif
        
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[j].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        PPP[j].AGS_Hsml = PPP[j].Hsml;
#endif
#endif
#ifdef CONDUCTION
        SphP[j].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        SphP[j].Eta_MHD_OhmicResistivity_Coeff = 0; SphP[j].Eta_MHD_HallEffect_Coeff = 0; SphP[j].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        SphP[j].Eta_ShearViscosity = 0; SphP[j].Zeta_BulkViscosity = 0;
#endif
#ifdef TURB_DIFFUSION
        SphP[j].TD_DiffCoeff = 0;
#endif
#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)
        SphP[j].HostHaloMass = 0;
#endif
#endif
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
        SphP[j].DelayTimeHII = 0;
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        SphP[j].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        SphP[j].Sfr = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        SphP[j].alpha = 0.0;
#endif
#if defined(BH_THERMALFEEDBACK)
        SphP[j].Injected_BH_Energy = 0;
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            SphP[j].E_gamma[k] = 0;
#if defined(RT_EVOLVE_NGAMMA)
            SphP[j].E_gamma_Pred[k] = 0; SphP[j].Dt_E_gamma[k] = 0;
#endif
        }
#endif

        /* note, if you want to use this routine to inject magnetic flux or cosmic rays, do this below */
#ifdef MAGNETIC
        for(k=0;k<3;k++)
        {
            SphP[j].BPred[k] = SphP[j].B[k] = 0; /* add magnetic flux here if desired */
            SphP[j].DtB[k] = 0;
        }
        SphP[j].divB = 0;
#ifdef DIVBCLEANING_DEDNER
        SphP[j].DtPhi = SphP[j].PhiPred = SphP[j].Phi = 0;
#endif
#endif
#ifdef COSMIC_RAYS
        SphP[j].CosmicRayEnergyPred = SphP[j].CosmicRayEnergy = 0; /* add CR energy here if desired */
        SphP[j].DtCosmicRayEnergy = 0;
#endif
        
        
        /* assign masses to both particles (so they sum correctly) */
        P[j].Mass = mass_of_new_particle;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].MassTrue = P[j].Mass;
#endif
        P[i].Mass -= P[j].Mass;
        
        /* shift the particle locations according to the random number we drew above */
        double dx, dy, dz;
        double sin_theta = sqrt(1 - cos_theta*cos_theta);
        dx = d_r * sin_theta * cos(phi);
        dy = d_r * sin_theta * sin(phi);
        dz = d_r * cos_theta;

        P[j].Pos[0] =  P[i].Pos[0] + dx;
        P[j].Pos[1] =  P[i].Pos[1] + dy;
        P[j].Pos[2] =  P[i].Pos[2] + dz;
        
        P[j].Vel[0] =  P[i].Vel[0] + dx / d_r * (All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s) * All.cf_atime;
        P[j].Vel[1] =  P[i].Vel[1] + dy / d_r * (All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s) * All.cf_atime;
        P[j].Vel[2] =  P[i].Vel[2] + dz / d_r * (All.BAL_v_outflow*1e5/All.UnitVelocity_in_cm_per_s) * All.cf_atime;
        SphP[j].VelPred[0] = P[j].Vel[0]; SphP[j].VelPred[1] = P[j].Vel[1]; SphP[j].VelPred[2] = P[j].Vel[2]; 
        
#if defined(BH_COSMIC_RAYS)
        /* inject cosmic rays alongside wind injection */
        double dEcr = All.BH_CosmicRay_Injection_Efficiency * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * (C / All.UnitVelocity_in_cm_per_s)*(C / All.UnitVelocity_in_cm_per_s);
        SphP[j].CosmicRayEnergy+=dEcr; SphP[j].CosmicRayEnergyPred+=dEcr;
#ifdef COSMIC_RAYS_M1
        dEcr*=COSMIC_RAYS_M1; for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k]+=dEcr*dir[k]; SphP[j].CosmicRayFluxPred[k]+=dEcr*dir[k];}
#endif
#endif

        /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
        force_add_star_to_tree(i, j);// (buggy)
        /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
    }
    
    BPP(i).unspawned_wind_mass = 0.0;
    
    return n_particles_split;
}
#endif
