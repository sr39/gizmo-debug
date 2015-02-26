/*! \file blackhole_feed.c
 *  \brief This is where particles are marked for gas accretion.
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clairity by parsing the existing code into
 *   smaller files and routines.  Some communication and black hole structures were modified
 *   to reduce memory usage.
 */

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
            BlackholeDataIn[j].mass_to_swallow_edd = BlackholeTempInfo[P[place].IndexMapToTempStruc].mass_to_swallow_edd;
            
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
    double h_i2, r2, r, u, hinv, hinv3, wk, dwk, vrel, vesc;
    double dpos[3];
    
#if defined(UNIFIED_FEEDBACK) || defined(BH_ENFORCE_EDDINGTON_LIMIT)
#if (defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)) && (defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION))
    double meddington, medd_max_accretable;
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
#if defined(BH_SWALLOWGAS)
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
#if (defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)) && (defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION))
        mass_to_swallow_edd = BlackholeTempInfo[P[target].IndexMapToTempStruc].mass_to_swallow_edd;
#endif
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
#if (defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)) && (defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION))
        mass_to_swallow_edd = BlackholeDataGet[target].mass_to_swallow_edd;
#endif
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
    
    eddington_factor = mass_to_swallow_edd / medd_max_accretable;   /* if <1 no problem, if >1, need to not set some swallowIDs */
#endif
#endif
    
#if defined(BH_SWALLOWGAS)
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
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
                        vesc = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime) + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                        r = sqrt(r2);
//                        if(P[j].Type==0)  printf("vrel=%g, vesc=%g, r=%g, cond2=%g, Type=%d\n", vrel, vesc, r, All.ForceSoftening[5]*(1.0-vrel*vrel/(vesc*vesc))/r, P[j].Type);

#ifdef BH_REPOSITION_ON_POTMIN
                        /* check if we've found a new potential minimum which is not moving too fast to 'jump' to */
                        if(P[j].Potential < minpot)
                        {
                            if(vrel <= vesc)
                            {
                                minpot = P[j].Potential;
                                for(k = 0; k < 3; k++) minpotpos[k] = P[j].Pos[k];
                            }
                        }
#endif
                        
                        
                        /* check_for_bh_merger.  Easy.  No Edd limit, just a pos and vel criteria. */
                        if((id != P[j].ID) && (P[j].Mass > 0) && (P[j].Type == 5))	/* we may have a black hole merger */
                        {
                            if(id != P[j].ID) /* check its not the same bh */
                            {
                                if(vrel > BH_CSND_FRAC_BH_MERGE * vesc)
                                {
                                    fprintf(FdBlackHolesDetails,
                                            "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n",
                                            ThisTask, All.Time, id, P[j].ID, vrel, vesc);
                                }
                                else
                                {
                                    printf("MARKING_BH_MERGER: P[j.]ID=%llu to be swallowed by id=%llu \n",
                                           (unsigned long long) P[j].ID, (unsigned long long) id);

                                    if(P[j].SwallowID < id && P[j].ID < id) // makes it so only one swallows the other
                                        P[j].SwallowID = id;
                                }
                            }
                        } // if(P[j].Type == 5) //
                        
                        
                        
/* This is a similar loop to what we already did in blackhole_environment, but here we stochastially
 reduce GRAVCAPT events in order to (statistically) obey the eddington limit */
                        
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
#ifdef BH_GRAVCAPTURE_SWALLOWS
//                        if(P[j].Type != 5)
                        if(P[j].Type == 0)		// enforce gas accretion only
#else
                            if((P[j].Type != 0)&&(P[j].Type != 5))
#endif
                            {                                
                                if(vrel < vesc){ /* bound */
                                    if( All.ForceSoftening[5]*(1.0-vrel*vrel/(vesc*vesc))/r > 1.0 )
                                    { /* apocenter within 2.8*epsilon (softening length) */
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        /* only count gas and stars towards the Eddington limit */
                                        if(P[j].Type == 1)
                                            p=10.0;
                                        else
                                            p=1/eddington_factor;   /*  >1 if below eddington limit (all particles swallowed),
                                                                        <1 if above eddington limit (reduces accretion accordinginly) */
                                        
                                        w = get_random_number(P[j].ID);
                                        if(w < p) {
                                            printf("MARKING_BH_FOOD: P[j.]ID=%llu to be swallowed by id=%llu \n",
                                                   (unsigned long long) P[j].ID, (unsigned long long) id);
                                            if(P[j].SwallowID < id) P[j].SwallowID = id;
                                        } else { /* w < p */
                                            printf("MARKING_BH_FOOD (should have been rejected): P[j.]ID=%llu to be swallowed by id=%llu \n",
                                                   (unsigned long long) P[j].ID, (unsigned long long) id);
                                            if(P[j].SwallowID < id)  P[j].SwallowID = id;
                                        }/* w < p */
#else
                                        /* simply mark all this to be accreted (can -greatly- exceed Eddington) */
                                        if(P[j].SwallowID < id)
                                            P[j].SwallowID = id;
#endif
                                    } /* if( All.ForceSoftening[5]*(1.0-vrel*vrel/(vesc*vesc))/sqrt(r2) > 1.0 ) */
                                } /* if(vrel < vesc) */
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
                            if((bh_mass_withdisk - (mass + mass_markedswallow))>0)
                                p = (bh_mass_withdisk - (mass + mass_markedswallow)) * wk / rho;
                            else
                                p = 0;
                            
#if defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                            p = 0;
#endif
                            
#if defined(BH_BAL_WINDS)
                            if(All.BAL_f_accretion>0) p /= All.BAL_f_accretion;
#endif
                            
                            w = get_random_number(P[j].ID);
                            if(w < p)
                            {
                                printf("MARKING_BH_FOOD: j %d w %g p %g TO_BE_SWALLOWED \n",j,w,p);
                                if(P[j].SwallowID < id) P[j].SwallowID = id;
                                N_gas_toswallow++;
#if defined(BH_BAL_WINDS)
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


