//
//  blackhole_feed.c
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




void blackhole_feed_loop(void)
{
    
    int i, j, k, ndone_flag, ndone;
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
    
    
    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),
                                                                     sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* Let's determine which particles may be swallowed by whom, and the weights for feedback */
    i = FirstActiveParticle;
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
                if(blackhole_feed_evaluate(i, 0, &nexport, Send_count) < 0)
                    break;
        
        MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
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
        
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].GradRho[k];
#endif
            }
            BlackholeDataIn[j].mass_to_swallow_edd = BlackholeTempInfo[P[target].IndexMapToTempStruc].mass_to_swallow_edd;
            
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(BH_PHOTONMOMENTUM) 	|| defined(BH_BAL_WINDS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
#endif
            BlackholeDataIn[j].Density = BPP(place).DensAroundStar;
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            BlackholeDataIn[j].ID = P[place].ID;
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
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult",nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_feed_evaluate(j, 1, &dummy, &dummy);

        
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
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++) //
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
    #ifdef BH_REPOSITION_ON_POTMIN
            if(BPP(place).BH_MinPot > BlackholeDataOut[j].BH_MinPot)
            {
                BPP(place).BH_MinPot = BlackholeDataOut[j].BH_MinPot;
                for(k = 0; k < 3; k++)
                    BPP(place).BH_MinPotPos[k] = BlackholeDataOut[j].BH_MinPotPos[k];
            }
    #endif
    #if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BlackholeDataOut[j].BH_angle_weighted_kernel_sum;
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
    
}







/* do loop over neighbors to get quantities for accretion */
int blackhole_feed_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, listindex = 0;
    MyIDType id;
    MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass;
    double h_i2, r2, r, u, hinv, hinv3, wk, dwk, vrel, csnd;
    double dpos[3];
    
#if defined(UNIFIED_FEEDBACK) || defined(BH_ENFORCE_EDDINGTON_LIMIT)
#if (defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)) && (defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION))
    double meddington, medd_max_accretable, medd_markedswallow;
    double mass_to_swallow_edd, eddington_factor;
#endif
#endif
    
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    double norm,theta,BH_disk_hr,*Jgas_in_Kernel;
    double BH_angle_weighted_kernel_sum=0;
#endif
    
#ifdef BH_THERMALFEEDBACK
    double energy;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
    MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = BHPOTVALUEINIT;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk;
#endif
#if defined(BH_SWALLOWGAS) || defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
    int N_gas_toswallow=0;
    double w=0,p=0,mass_markedswallow=0,bh_mass_withdisk=0;
#endif
    
    
    //BlackholeTempInfo[target].mass_to_swallow_edd /
    /* these are the BH properties */
    if(mode == 0)
    {
        pos = P[target].Pos;
        rho = BPP(target).DensAroundStar;
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        h_i = PPP[target].Hsml;
        mass = P[target].Mass;
        bh_mass = BPP(target).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        velocity = P[target].Vel;
        id = P[target].ID;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = P[target].GradRho;
        BH_disk_hr = P[target].BH_disk_hr;
#endif
        mass_to_swallow_edd = BlackholeTempInfo[P[target].IndexMapToTempStruc].mass_to_swallow_edd;

    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        rho = BlackholeDataGet[target].Density;
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
        h_i = BlackholeDataGet[target].Hsml;
        mass = BlackholeDataGet[target].Mass;
        bh_mass = BlackholeDataGet[target].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;
#endif
        velocity = BlackholeDataGet[target].Vel;
        id = BlackholeDataGet[target].ID;
#if defined(BH_PHOTONMOMENTUM)  || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
#endif
        mass_to_swallow_edd = BlackholeDataGet[target].mass_to_swallow_edd;
    }
    
    if((mass<0)||(h_i<=0)) return -1;
    
    
    /* initialize variables before SPH loop is started */
    h_i2 = h_i * h_i;
    hinv = 1 / h_i;
    hinv3 = hinv * hinv * hinv;
#ifdef BH_ENFORCE_EDDINGTON_LIMIT
#if (defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)) && (defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION))
    meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON)) * (bh_mass/All.HubbleParam) * All.UnitTime_in_s;
    medd_max_accretable = All.BlackHoleEddingtonFactor * meddington * dt;
    
    eddington_factor = mass_to_swallow_edd / medd_max_accretable;   /* if <1 no problem, if >1, need to unset some swallowIDs */
#endif
#endif
#if defined(BH_SWALLOWGAS) || defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif
    
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    double m_to_swallow_thispart;
#endif
    
    /* Now start the actual SPH computation for this BH particle */
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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            BH_angle_weighted_kernel_sum = 0;
#endif
            
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                if(P[j].Mass > 0)
                {
                    for(k=0;k<3;k++) dpos[k] = pos[k] - P[j].Pos[k];
#ifdef PERIODIC		/*  find the closest image in the given box size  */
                    dpos[0]=NEAREST_X(dpos[0]); dpos[1]=NEAREST_Y(dpos[1]); dpos[2]=NEAREST_Z(dpos[2]);
#endif
                    r2=0; for(k=0;k<3;k++) r2+=dpos[k]*dpos[k];
                    
                    if(r2 < h_i2)
                    {
                        
                        vrel = 0;
                        for(k=0;k<3;k++) vrel += (P[j].Vel[k] - velocity[k])*(P[j].Vel[k] - velocity[k]);
                        vrel = sqrt(vrel) / All.cf_atime;       /* do this once and use below */
                        csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime) + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                        r = sqrt(r2);
                        
#ifdef BH_REPOSITION_ON_POTMIN
                        /* check if we've found a new potential minimum which is not moving too fast to 'jump' to */
                        if(P[j].Potential < minpot)
                        {
                            if(vrel <= csnd)
                            {
                                minpot = P[j].Potential;
                                for(k = 0; k < 3; k++) minpotpos[k] = P[j].Pos[k];
                            }
                        }
#endif
                        
                        
                        /* check_for_bh_merger.  Easy.  No Edd limit, just a pos and vel criteria. */
                        if((id != P[j].ID) && (P[j].Mass > 0) && (P[j].Type == 5))	/* we may have a black hole merger */
                        {
//                            r2=0;
//                            for(k=0;k<3;k++) r2+=(P[j].Pos[k] - pos[k])*(P[j].Pos[k] - pos[k]);
                            
//                            vrel = 0;
//                            for(k=0;k<3;k++) vrel += (P[j].Vel[k] - vel[k])*(P[j].Vel[k] - vel[k]);
//                            vrel = sqrt(vrel) / All.cf_atime;
                            csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime) + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                            
                            if(id != P[j].ID) /* check its not the same bh */
                            {
                                /* with feedback on, sound speeds get -very- low, would be silly to use;
                                 instead use the escape velocity and follow to pairing */
                                if(vrel > BH_CSND_FRAC_BH_MERGE * csnd)
                                {
                                    fprintf(FdBlackHolesDetails,
                                            "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n",
                                            ThisTask, All.Time, id, P[j].ID, vrel, csnd);
                                }
                                else
                                {
                                    if(P[j].SwallowID < id && P[j].ID < id) // makes it so only one swallows the other
                                        P[j].SwallowID = id;
                                }
                            }
                        } // if(P[j].Type == 5) //
                        
                        
                        
/* This is a similar loop to what we already did in blackhole_environment, but here we stochastially
 reduce GRAVCAPT events in order to (statistically) obey the eddington limit */
                        
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
#ifdef BH_GRAVCAPTURE_SWALLOWS
                        if(P[j].Type != 5)
#else
                            if((P[j].Type != 0)&&(P[j].Type != 5))
#endif
                            {
                                for(k = 0, vrel = 0; k < 3; k++)
                                    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);
                                vrel = sqrt(vrel) / All.cf_atime;
                                r = sqrt(r2);
                                csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(r*All.cf_atime)); /* escape velocity */
                                
                                if(vrel < csnd){ /* bound */
                                    if( All.ForceSoftening[5]*(1.0-vrel*vrel/(csnd*csnd))/r > 1.0 )
                                    { /* apocenter within 2.8*epsilon (softening length) */
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        /* only count gas and stars towards the Eddington limit */
                                        p=1/eddington_factor;   /*  >1 if below eddington limit (all particles swallowed), 
                                                                    <1 if above eddington limit (reduces accretion accordinginly) */

//                                        m_to_swallow_thispart=P[j].Mass;
//                                        if((P[j].Type != 1)||(All.ComovingIntegrationOn && (P[j].Type==0||P[j].Type==4)))
//                                        {
//                                            if(medd_max_accretable-mass_markedswallow <= 0)
//                                            {
//                                                p=0;
//                                            } else {
//#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
//                                                if(All.BAL_f_accretion>0) m_to_swallow_thispart *= All.BAL_f_accretion;
//#endif
//                                                if(m_to_swallow_thispart >  medd_max_accretable-mass_markedswallow)
//                                                    p = (medd_max_accretable-mass_markedswallow)/m_to_swallow_thispart;
//                                            }}
                                        
                                        
                                        w = get_random_number(P[j].ID);
                                        if(w < p) {
                                            printf("MARKING_BH_FOOD: j %d w %g p_acc %g facc %g TO_BE_SWALLOWED \n",j,w,p,m_to_swallow_thispart/P[j].Mass);
                                            if(P[j].SwallowID < id) {
                                                P[j].SwallowID = id;
//                                                if((P[j].Type != 1)||(All.ComovingIntegrationOn && (P[j].Type==0||P[j].Type==4)))
//                                                    mass_markedswallow += m_to_swallow_thispart;
                                            } /* P[j].SwallowID < id */
                                        } else { /* w < p */
                                            printf("REJECTED_BH_FOOD (based on eddington limit): j %d w %g p_acc %g facc %g TO_BE_SWALLOWED \n",j,w,p,m_to_swallow_thispart/P[j].Mass);
                                            
                                        }/* w < p */
#else
                                        /* simply mark all this to be accreted (can -greatly- exceed Eddington) */
                                        if(P[j].SwallowID < id)
                                            P[j].SwallowID = id;
#endif
                                    } /* if( All.ForceSoftening[5]*(1.0-vrel*vrel/(csnd*csnd))/sqrt(r2) > 1.0 ) */
                                } /* if(vrel < csnd) */
                            } /* type check */
#endif // BH_GRAVCAPTURE_SWALLOWS
                        
                        
                        

                        /* now is the more standard accretion only of gas, according to the mdot calculated before */
                        if(P[j].Type == 0)
                        {
                            /* here we have a gas particle */
                            u = r * hinv;
                            kernel_main(u,hinv3,hinv*hinv3,&wk,&dwk,-1);
                            
#ifdef BH_SWALLOWGAS
                            /* compute accretion probability */
                            if((bh_mass_withdisk - (mass+mass_markedswallow))>0)
                                p = (bh_mass_withdisk - (mass+mass_markedswallow)) * wk / rho;
                            else
                                p = 0;
#if defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                            p = 0;
#endif
#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                            if(All.BAL_f_accretion>0) p /= All.BAL_f_accretion;
#endif
                            w = get_random_number(P[j].ID);
                            if(w < p)
                            {
                                printf("MARKING_BH_FOOD: j %d w %g p %g TO_BE_SWALLOWED \n",j,w,p);
                                if(P[j].SwallowID < id) P[j].SwallowID = id;
                                N_gas_toswallow++;
#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                                mass_markedswallow += P[j].Mass*All.BAL_f_accretion;
#else
                                mass_markedswallow += P[j].Mass;
#endif
                            } // if(w < p)
#endif // BH_SWALLOWGAS
                            
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                            /* calculate the angle-weighting for the photon momentum */
                            if((mdot>0)&&(dt>0)&&(r>0))
                            {
                                /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                norm=0; for(k=0;k<3;k++) norm+=(dpos[k]/r)*Jgas_in_Kernel[k];
                                norm=fabs(norm); theta=acos(norm);
                                BH_angle_weighted_kernel_sum += bh_angleweight_localcoupling(j,BH_disk_hr,theta);
                            }
#endif
                            
#ifdef BH_THERMALFEEDBACK
#ifdef UNIFIED_FEEDBACK
                            meddington = (4*M_PI*GRAVITY*C*PROTONMASS/(All.BlackHoleRadiativeEfficiency*C*C*THOMPSON))*bh_mass*All.UnitTime_in_s/All.HubbleParam;
                            if(mdot > All.RadioThreshold * meddington)
#endif
                            {
                                energy = All.BlackHoleFeedbackFactor*All.BlackHoleRadiativeEfficiency * mdot*dt * pow(C/All.UnitVelocity_in_cm_per_s,2);
                                if(rho > 0)
                                    SphP[j].Injected_BH_Energy += (wk/rho) * energy * P[j].Mass;
                            }
#endif
                            
                        } // if(P[j].Type == 0)
                        
                        
                        
                        
                        
                        
                        

                        
                    } // if(r2 < h_i2)
                } // if(P[j].Mass > 0)
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
        } // mode==1
    } // while(startnode >= 0) (outer of the double-loop)
    
    
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BH_angle_weighted_kernel_sum;  /* need to correct target index */
#endif
#ifdef BH_REPOSITION_ON_POTMIN
        BPP(target).BH_MinPot = minpot;
        for(k = 0; k < 3; k++)
            BPP(target).BH_MinPotPos[k] = minpotpos[k];
#endif
    }
    else
    {
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        BlackholeDataResult[target].BH_angle_weighted_kernel_sum = BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
        BlackholeDataResult[target].BH_MinPot = minpot;
        for(k = 0; k < 3; k++)
            BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];
#endif
    }
    return 0;
} /* closes bh_evaluate routine */






//
//
//
//void blackhole_properties_loop(void)
//{
//    
//    int  i, n;
//    double fac, mdot, dt;
//    double medd;
//#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_GRAVACCRETION) || defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DYNFRICTION)
//    int k;
//#endif
//
//
//    
//#ifdef BH_GRAVACCRETION
//    double m_tmp_for_bhar;
//    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar,mdisk_for_bhar;
//    double f0_for_bhar;
//#endif
//#ifdef BH_SUBGRIDBHVARIABILITY
//    long nsubgridvar;
//    int jsub;
//    double varsg1,varsg2;
//    double omega_ri,n0_sgrid_elements,norm_subgrid,time_var_subgridvar;
//    gsl_rng *random_generator_forbh;
//#endif
//#ifdef BH_BONDI
//    double norm, soundspeed, bhvel, rho;
//#endif
//#ifdef KD_FRICTION
//    /* add a friction force for the black-holes, accounting for the environment */
//    double fac_friction, relvel, accgrv, accfrc;
//    double a_erf, lambda;
//#endif
//#ifdef BH_ALPHADISK_ACCRETION
//    double mdot_alphadisk;
//#endif
//#ifdef BH_ENFORCE_EDDINGTON_LIMIT
//    double meddington;
//#endif
//    
//    /* TODO: make this simply a loop over Nbh (active on loc proc) */
//    
//    for(i=0; i<N_active_loc_BHs; i++)
//    {
//        n = BlackholeTempInfo[i].index;
//        
//
////        printf("this task = %d    i = %d     n = %d  \n",ThisTask, i , n );
//
//            /* define the timestep */
//#ifndef WAKEUP
//        dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
//#else
//        dt = P[n].dt_step * All.Timebase_interval / All.cf_hubble_a;
//#endif
//            
//#ifdef BH_ENFORCE_EDDINGTON_LIMIT
//        /* define the Eddington limit for use later
//         * (Note: we take here a radiative efficiency
//         * as specified in the parameter file) */
//        meddington = (4 * M_PI * GRAVITY * C * PROTONMASS /
//                          (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON) ) *
//        (BPP(n).BH_Mass/All.HubbleParam) * All.UnitTime_in_s;
//#endif
//            
//        /* always initialize/default to zero accretion rate */
//        mdot=0;
//        BPP(n).BH_Mdot=0;
//        
//        normalize_temp_info_struct(i);
//
//
//            
//        
//#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
//        /* TODO:  move this to a separate function in blackhole_physics.c */
//        
//        /* pre-set quantities needed for long-range radiation pressure terms */
//        P[n].BH_disk_hr=1/3; P[n].GradRho[0]=P[n].GradRho[1]=0; P[n].GradRho[2]=1;
//        if(BlackholeTempInfo[i].Mgas_in_Kernel > 0)
//        {
//            /* estimate h/R surrounding the BH from the gas density gradients */
//            fac = 0; /* dummy variable */
//            for(k=0;k<3;k++)
//                fac += BlackholeTempInfo[i].GradRho_in_Kernel[k]*BlackholeTempInfo[i].GradRho_in_Kernel[k];
//            P[n].BH_disk_hr = P[n].DensAroundStar / (PPP[n].Hsml * sqrt(fac)) * 1.3;
//            /* 1.3 factor from integrating exponential disk
//             * with h/R=const over gaussian kernel, for width=1/3 (quintic kernel);
//             everything here is in code units, comes out dimensionless */
//            
//            /* use the gradrho vector as a surrogate to hold the orientation of the angular momentum */
//            fac=0;
//            for(k=0;k<3;k++)
//                fac += BlackholeTempInfo[i].Jgas_in_Kernel[k]*BlackholeTempInfo[i].Jgas_in_Kernel[k];
//            fac=sqrt(fac);
//            if(fac>0)
//                for(k=0;k<3;k++)
//                    P[n].GradRho[k] = BlackholeTempInfo[i].Jgas_in_Kernel[k]/fac;
//            /* now, the P[n].GradRho[k] field for the BH holds the orientation of the UNIT angular momentum vector
//             NOTE it is important that HARD-WIRED into the code, this blackhole calculation comes after the density calculation
//             but before the forcetree update and walk; otherwise, this won't be used correctly there */
//        }
//#endif // if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
//        
//        
//        // TODO:
//        // next come some HUGE blocks of code that are ignored most of the time
//        // want to move all of this to a single "evaluate_bh_acc" function that hides this
//        //
//#ifdef BH_GRAVACCRETION
//        /* calculate mdot: gravitational instability accretion rate from Hopkins & Quataert 2011 */
//        if(BlackholeTempInfo[i].Mgas_in_Kernel > 0)
//        {
//            m_tmp_for_bhar = BlackholeTempInfo[i].Mgas_in_Kernel + BlackholeTempInfo[i].Malt_in_Kernel;
//            r0_for_bhar = PPP[n].Hsml * All.cf_atime; /* convert to physical units */
//            j_tmp_for_bhar=0;
//            for(k=0;k<3;k++)
//                j_tmp_for_bhar += BlackholeTempInfo[i].Jalt_in_Kernel[k]*BlackholeTempInfo[i].Jalt_in_Kernel[k];
//            j_tmp_for_bhar=sqrt(j_tmp_for_bhar);
//            /* jx,y,z, is independent of 'a_scale' b/c ~ m*r*v, vphys=v/a, rphys=r*a */
//            
//            bh_mass = BPP(n).BH_Mass;
//#ifdef BH_ALPHADISK_ACCRETION
//            bh_mass += BPP(n).BH_Mass_AlphaDisk;
//#endif
//            fgas_for_bhar = BlackholeTempInfo[i].Mgas_in_Kernel / m_tmp_for_bhar;
//            fac = m_tmp_for_bhar * r0_for_bhar * sqrt(All.G*(m_tmp_for_bhar+bh_mass)/r0_for_bhar);
//            /* All.G is G in code (physical) units */
//            f_disk_for_bhar = fgas_for_bhar + (1.75*j_tmp_for_bhar/fac);
//            if(f_disk_for_bhar>1) f_disk_for_bhar=1;
//            
//            if((f_disk_for_bhar<=0)||(bh_mass <=0)||(fgas_for_bhar<=0)||(m_tmp_for_bhar<=0))
//            {
//                mdot = 0;
//            } else {
//                mdisk_for_bhar = m_tmp_for_bhar*f_disk_for_bhar * (All.UnitMass_in_g/(All.HubbleParam * 1.0e9*SOLAR_MASS)); /* mdisk/1e9msun */
//                bh_mass *= All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS); /* mbh/1e8msun */
//                r0_for_bhar *= All.UnitLength_in_cm/(All.HubbleParam * 3.086e20); /* r0/100pc */
//                f0_for_bhar = 0.31*f_disk_for_bhar*f_disk_for_bhar*pow(mdisk_for_bhar,-1./3.); /* dimensionless factor for equations */
//                fac = (10.0*(SOLAR_MASS/All.UnitMass_in_g)/(SEC_PER_YEAR/All.UnitTime_in_s)); /* basic units */
//                
//                mdot = All.BlackHoleAccretionFactor * fac * mdisk_for_bhar *
//                pow(f_disk_for_bhar,5./2.) * pow(bh_mass,1./6.) *
//                pow(r0_for_bhar,-3./2.) / (1 + f0_for_bhar/fgas_for_bhar);
//                
//                printf("BH GravAcc Eval :: mdot %g BHaccFac %g Norm %g fdisk %g bh_8 %g fgas %g f0 %g mdisk_9 %g r0_100 %g \n\n",
//                       mdot,All.BlackHoleAccretionFactor,fac,
//                       f_disk_for_bhar,bh_mass,fgas_for_bhar,f0_for_bhar,mdisk_for_bhar,r0_for_bhar);fflush(stdout);
//            } // if(f_disk_for_bhar<=0)
//        }
//#endif // ifdef BH_GRAVACCRETION
//        
//        
//#ifdef BH_BONDI
//        /* heres where we calculate the Bondi accretion rate, if that's going to be used */
//        bhvel = 0;
//#ifdef BH_USE_GASVEL_IN_BONDI
//        for(k=0;k<3;k++) bhvel += BlackholeTempInfo[i].BH_SurroundingGasVel[k]*BlackholeTempInfo[i].BH_SurroundingGasVel[k];
//#endif
//        rho = BPP(n).DensAroundStar * All.cf_a3inv; /* we want all quantities in physical units */
//        soundspeed = GAMMA*GAMMA_MINUS1 * BlackholeTempInfo[i].BH_InternalEnergy; // this is in physical units now
//        fac = pow(soundspeed+bhvel, 1.5);
//        if(fac > 0)
//        {
//            double AccretionFactor = All.BlackHoleAccretionFactor;
//#ifdef BH_VARIABLE_ACCRETION_FACTOR
//            /* variable-alpha model (Booth&Schaye 2009): now All.BlackHoleAccretionFactor is the slope of the density dependence */
//            AccretionFactor = 1.0;
//            if(rho > All.PhysDensThresh)
//                AccretionFactor = pow(rho/All.PhysDensThresh, All.BlackHoleAccretionFactor);
//#endif
//            mdot = 4. * M_PI * AccretionFactor * All.G * All.G * BPP(n).BH_Mass * BPP(n).BH_Mass * rho / fac;
//        }
//#endif // ifdef BH_BONDI
//        
//        
//#ifdef BH_GRAVCAPTURE_SWALLOWS
//        mdot = 0; /* force mdot=0 despite any earlier settings here */
//#endif
//        
//        
//#ifdef BH_ALPHADISK_ACCRETION
//        /* use the mass in the accretion disk from the previous timestep to determine the BH accretion rate */
//        mdot_alphadisk = mdot;
//        mdot = All.BlackHoleAccretionFactor *
//        (2.45 * (SOLAR_MASS/All.UnitMass_in_g)/(SEC_PER_YEAR/All.UnitTime_in_s)) * // normalization
//        pow( 0.1 , 8./7.) * // viscous disk 'alpha'
//        pow( BPP(n).BH_Mass*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , -5./14. ) * // mbh dependence
//        pow( BPP(n).BH_Mass_AlphaDisk*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , 10./7. ) * // m_disk dependence
//        pow( DMIN(0.2,DMIN(PPP[n].Hsml,All.ForceSoftening[5])*All.cf_atime*All.UnitLength_in_cm/(All.HubbleParam * 3.086e18)) , -25./14. ); // r_disk dependence
//        if(mdot<=0) mdot=0;
//        if(dt>0)
//        {
//#ifdef BH_BAL_WINDS
//            /* this is just here to prevent it from accidentally going to negative mass */
//            if(mdot > BPP(n).BH_Mass_AlphaDisk/(All.BAL_f_accretion*dt)) mdot = BPP(n).BH_Mass_AlphaDisk/(All.BAL_f_accretion*dt);
//#else
//            if(mdot > BPP(n).BH_Mass_AlphaDisk/dt) mdot = BPP(n).BH_Mass_AlphaDisk/dt;
//#endif
//        }
//#endif
//        
//        
//#ifdef BH_SUBGRIDBHVARIABILITY
//        /* account for sub-grid accretion rate variability */
//        if((mdot>0)&&(dt>0)&&(P[n].DensAroundStar>0))
//        {
//            omega_ri=sqrt(All.G*P[n].DensAroundStar*All.cf_a3inv); /* dynamical frequency in physical units */
//            n0_sgrid_elements=10.0; norm_subgrid=0.55*3.256/sqrt(n0_sgrid_elements);
//            nsubgridvar=(long)P[n].ID + (long)(All.Time/((All.TimeMax-All.TimeBegin)/1000.));
//            /* this line just allows 'resetting' the time constants every so often, while generally keeping them steady */
//            if(All.ComovingIntegrationOn)
//                fac=omega_ri * (evaluate_stellar_age_Gyr(0.001)/(0.001*All.UnitTime_in_Megayears/All.HubbleParam));
//            else
//                fac=omega_ri * All.Time; /* All.Time is physical time, this is good */
//            random_generator_forbh=gsl_rng_alloc(gsl_rng_ranlxd1);
//            gsl_rng_set(random_generator_forbh,nsubgridvar);
//            if(n0_sgrid_elements >= 1) {
//                for(jsub=1;jsub<=n0_sgrid_elements;jsub++) {
//                    varsg1=gsl_rng_uniform(random_generator_forbh);
//                    varsg2=gsl_ran_ugaussian(random_generator_forbh);
//                    time_var_subgridvar=fac*pow(omega_ri*dt,-((float)jsub)/n0_sgrid_elements) + 2.*M_PI*varsg1;
//                    mdot *= exp( norm_subgrid*cos(time_var_subgridvar)*varsg2 );
//                    /*
//                     printf("SUBGRIDVAR :: mdot %g x %g cosx %g om_ri %g All_t %g dt %g nsubgridvar %ld n0 %g norm %g jsub %d ru %g rg %g \n",
//                     mdot,x,cos(x),omega_ri,All.Time,dt,nsubgridvar,n0_sgrid_elements,norm_subgrid,jsub,varsg1,varsg2);fflush(stdout);
//                     */
//                }}
//            gsl_rng_free(random_generator_forbh);
//        } // if(mdot > 0)
//#endif
//        
//        
//#ifdef BH_ENFORCE_EDDINGTON_LIMIT
//        /* cap the maximum at the Eddington limit */
//        if(mdot > All.BlackHoleEddingtonFactor * meddington)
//            mdot = All.BlackHoleEddingtonFactor * meddington;
//#endif
//        /* alright, now we can FINALLY set the BH accretion rate */
//        BPP(n).BH_Mdot = mdot;
//        
//        // end of TODO, bhacc rate is finalized.
//        
//        
//        
//        /* dump the results so far to the 'blackhole_details' files */
//        fac=0; medd=0;
//#ifdef BH_ALPHADISK_ACCRETION
//        fac=BPP(n).BH_Mass_AlphaDisk;
//        medd=mdot_alphadisk;			// is this really what we want?!? I don't think so...
//#endif
//        fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g %g %g %g   %2.7f %2.7f %2.7f\n",
//                P[n].ID, All.Time, BPP(n).BH_Mass, fac, P[n].Mass, mdot, medd,
//                BPP(n).DensAroundStar*All.cf_a3inv, BlackholeTempInfo[i].BH_InternalEnergy,
//                P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]);
//        
//        
//#ifdef BH_DRAG
//        /* add a drag force for the black-holes, accounting for the accretion */
//        if((dt>0)&&(BPP(n).BH_Mass>0))
//        {
//            fac = BPP(n).BH_Mdot * dt / BPP(n).BH_Mass;
//#ifdef BH_STRONG_DRAG
//            /* make the force stronger to keep the BH from wandering */
//            fac = meddington * dt / BPP(n).BH_Mass;
//#endif
//            if(fac>1) fac=1;
//            for(k = 0; k < 3; k++)
//                P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac * BlackholeTempInfo[i].BH_SurroundingGasVel[k] / dt;
//        } // if((dt>0)&&(BPP(n).BH_Mass>0))
//#endif
//        
//        
//        
//#ifdef BH_DYNFRICTION
//        if(BlackholeTempInfo[i].DF_mmax_particles>0) /* found something in the kernel, we can proceed */
//        {
//            /* averaged value for colomb logarithm and integral over the distribution function */
//            /* fac_friction = log(lambda) * [erf(x) - 2*x*exp(-x^2)/sqrt(pi)]                  */
//            /*       lambda = b_max * v^2 / G / (M+m)                                          */
//            /*        b_max = Size of system (e.g. Rvir)                                       */
//            /*            v = Relative velocity of BH with respect to the environment          */
//            /*            M = Mass of BH                                                       */
//            /*            m = individual mass elements composing the large system (e.g. m<<M)  */
//            /*            x = v/sqrt(2)/sigma                                                  */
//            /*        sigma = width of the max. distr. of the host system                      */
//            /*                (e.g. sigma = v_disp / 3                                         */
//            bh_mass = BPP(n).BH_Mass;
//#ifdef BH_ALPHADISK_ACCRETION
//            bh_mass += BPP(n).BH_Mass_AlphaDisk;
//#endif
//            double bhvel_df=0; for(k=0;k<3;k++) bhvel_df += BlackholeTempInfo[i].DF_mean_vel[k]*BlackholeTempInfo[i].DF_mean_vel[k];
//            /* First term is approximation of the error function */
//            fac = 8 * (M_PI - 3) / (3 * M_PI * (4. - M_PI));
//            x = sqrt(bhvel_df) / (sqrt(2) * BlackholeTempInfo[i].DF_rms_vel);
//            double fac_friction =  x / fabs(x) * sqrt(1 - exp(-x * x * (4 / M_PI + fac * x * x) / (1 + fac * x * x))) - 2 * x / sqrt(M_PI) * exp(-x * x);
//            /* now the Coulomb logarithm */
//            fac = 50. * 3.086e21 / (All.UnitLength_in_cm/All.HubbleParam); /* impact parameter */
//            fac_friction *= log(1. + fac * bhvel_df / (All.G * bh_mass));
//            /* now we add a correction to only apply this force if M_BH is not >> <m_particles> */
//            fac_friction *= 1 / (1 + bh_mass / (5.*BlackholeTempInfo[i].DF_mmax_particles));
//            /* now the dimensional part of the force */
//            fac = (BlackholeTempInfo[i].Mgas_in_Kernel+BlackholeTempInfo[i].Malt_in_Kernel) /
//            ( (4*M_PI/3) * pow(PPP[n].Hsml*All.cf_atime,3) ); /* mean density of all mass inside kernel */
//            fac_friction *= 4*M_PI * All.G * All.G * fac * bh_mass / (bhvel_df*sqrt(bhvel_df));
//            /* now apply this to the actual acceleration */
//            if(fac_friction<0) fac_friction=0; if(isnan(fac_friction)) fac_friction=0;
//            for(k = 0; k < 3; k++)
//                P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac_friction * BlackholeTempInfo[i].DF_mean_vel[k];
//        }
//#endif
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        /* INCREMENT BH mass if mdot > 0 */
//        BPP(n).BH_Mass += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
//#ifdef BH_ALPHADISK_ACCRETION
//#ifdef BH_BAL_WINDS
//        BPP(n).BH_Mass_AlphaDisk += (mdot_alphadisk-BPP(n).BH_Mdot/All.BAL_f_accretion) * dt;
//        /* correct real particle mass for mass flux in accretion-disk wind */
//        P[n].Mass -= BPP(n).BH_Mdot*dt / All.BAL_f_accretion;		// again seems funny, must have been increased earlier
//#else
//        BPP(n).BH_Mass_AlphaDisk += (mdot_alphadisk-BPP(n).BH_Mdot) * dt;
//#endif
//        if(BPP(n).BH_Mass_AlphaDisk<0) BPP(n).BH_Mass_AlphaDisk=0;
//        if(P[n].Mass<0) P[n].Mass=0;
//#endif
//        
//#ifdef BH_BUBBLES
//        BPP(n).BH_Mass_bubbles += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
//#ifdef UNIFIED_FEEDBACK
//        if(BPP(n).BH_Mdot < All.RadioThreshold * meddington)
//            BPP(n).BH_Mass_radio += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
//#endif
//#endif
//    }// for(i=0; i<N_active_loc_BHs; i++)
//
//}
//
//

