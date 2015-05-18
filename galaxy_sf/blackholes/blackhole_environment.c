#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"

/*! \file blackhole_environment.c
 *  \brief routines for evaluating black hole environment
 */
/*
 * This file was largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   by Paul Torrey (ptorrey@mit.edu) on 1/9/15 for clairity.  The main functional difference is that BlackholeTempInfo
 *   is now allocated only for N_active_loc_BHs, rather than NumPart (as was done before).  Some
 *   extra index gymnastics are required to follow this change through in the MPI comm routines.
 *   Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */



/* quantities that pass IN to the 'blackhole_environment_evaluate' routines */
static struct blackholedata_in
{
#if defined(BH_GRAVCAPTURE_GAS)
    MyDouble Mass;
#endif
    MyDouble Pos[3];
    MyFloat Vel[3];
    MyFloat Hsml;
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
}
*BlackholeDataIn, *BlackholeDataGet;



void blackhole_environment_loop(void)
{
    int i, j, k, nexport, nimport, place, ngrp, recvTask, dummy;
    int ndone_flag, ndone;
    MPI_Status status;
    
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackhole_temp_particle_data) +
                                                             sizemax(sizeof(struct blackholedata_in),
                                                                     sizeof(struct blackhole_temp_particle_data))));
    
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    /* Scan gas particles for angular momentum, boundedness, etc */
    i = FirstActiveParticle;	/* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local active BH particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(blackhole_environment_evaluate(i, 0, &nexport, Send_count) < 0)
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
            }
#if defined(BH_GRAVCAPTURE_GAS)
            BlackholeDataIn[j].Mass = P[place].Mass;
#endif
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].ID = P[place].ID;
            memcpy(BlackholeDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
        
        
        BlackholeDataPasserResult = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserResult", nexport * sizeof(struct blackhole_temp_particle_data));
        BlackholeDataPasserOut = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserOut", nimport * sizeof(struct blackhole_temp_particle_data));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_environment_evaluate(j, 1, &dummy, &dummy);
        
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
                    MPI_Sendrecv(&BlackholeDataPasserResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                 MPI_BYTE, recvTask, TAG_DENS_B,
                                 &BlackholeDataPasserOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                 MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
                }
            }
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_blackhole(&BlackholeDataPasserOut[j], P[place].IndexMapToTempStruc, 1);
        } // for(j = 0; j < nexport; j++)
        myfree(BlackholeDataPasserOut);
        myfree(BlackholeDataPasserResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
}




/* routine to return the values we need of the properties of the gas, stars, etc in the vicinity of the BH -- these all factor into the BHAR */
int blackhole_environment_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, k, n, listindex=0, mod_index;
    MyFloat *pos, h_i, *vel, hinv;
    MyIDType id;
    
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat hinv3, wk, dwk, u;
    u=wk=dwk=0;
#endif
    
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat mass, vrel, vbound, r2;
#endif
    
    double dP[3],dv[3],wt;
    struct blackhole_temp_particle_data out;
    memset(&out, 0, sizeof(struct blackhole_temp_particle_data));
    
    /* these are the BH properties */
    if(mode == 0)
    {
#if defined(BH_GRAVCAPTURE_GAS)
        mass = P[target].Mass;
#endif
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
    }
    else
    {
#if defined(BH_GRAVCAPTURE_GAS)
        mass = BlackholeDataGet[target].Mass;
#endif
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        mod_index = 0;                              /* this is not used for mode==1, but this avoids compiler error */
    }
    
    if(h_i < 0) return -1;
    hinv = 1./h_i;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    hinv3 = hinv*hinv*hinv;
#endif
    
    if(mode == 0)
    {
        startnode = All.MaxPart;  /* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;        /* open it */
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
                
                
                if( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
                {
                    wt = P[j].Mass;
                    dP[0] = P[j].Pos[0]-pos[0];
                    dP[1] = P[j].Pos[1]-pos[1];
                    dP[2] = P[j].Pos[2]-pos[2];
#ifdef PERIODIC
                    NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                    dv[0] = P[j].Vel[0]-vel[0];
                    dv[1] = P[j].Vel[1]-vel[1];
                    dv[2] = P[j].Vel[2]-vel[2];
#ifdef SHEARING_BOX
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[SHEARING_BOX_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[SHEARING_BOX_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
#ifdef BH_DYNFRICTION
                    for (k=0;k<3;k++)
                    {
                        out.DF_mean_vel[k] += wt*dv[k];
                        out.DF_rms_vel += wt*dv[k]*dv[k];
                        if(P[j].Mass>out.DF_mmax_particles) out.DF_mmax_particles=P[j].Mass;
                    }
#endif
                    if(P[j].Type==0) /* we found gas in BH's kernel */
                    {
                        out.Mgas_in_Kernel += wt;
                        out.BH_InternalEnergy += wt*SphP[j].InternalEnergy;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                        out.Jgas_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jgas_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jgas_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                        u=0;
                        for(k=0;k<3;k++) u+=dP[k]*dP[k];
                        u=sqrt(u)/h_i;
                        kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,1);
                        dwk /= u*h_i;
                        for(k=0;k<3;k++) out.GradRho_in_Kernel[k] += wt * dwk * fabs(dP[k]);
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
                        for(k=0;k<3;k++)
                        {
                            out.BH_SurroundingGasVel[k] += wt*dv[k];
                        }
#endif
                    } else { /* not gas, not BH */
                        out.Malt_in_Kernel += wt;
                        out.Jalt_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jalt_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jalt_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                    
                    
#if defined(BH_GRAVCAPTURE_GAS)
                    /* XM: I formally distinguish BH_GRAVCAPTURE_GAS and BH_GRAVCAPTURE_NONGAS. The former applies to
                     gas ONLY, as an accretion model. The later can be combined with any accreton model.
                     Currently, I only allow gas acretion to contribute to BH_Mdot (consistent with the energy radiating away).
                     For star particles, if there is an alpha-disk, they are captured to the disk. If not, they directly go
                     to the hole, without any contribution to BH_Mdot and feedback. This can be modified in the swallow loop
                     for other purposes. */
                    /* XM: The goal of the following part is to estimate BH_Mdot, which will be used to evaluate feedback strength.
                     Therefore, we only need it when we enable BH_GRAVCAPTURE_GAS as gas accretion model. */
                    if( (P[j].Mass > 0) && (P[j].Type == 0))
                    {
                        vrel = 0;
                        for(k=0;k<3;k++) vrel += (P[j].Vel[k] - vel[k])*(P[j].Vel[k] - vel[k]);
                        vrel = sqrt(vrel) / All.cf_atime;
                        
                        r2=0; for(k=0;k<3;k++) r2+=dP[k]*dP[k];
                        vbound = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime) + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                        
                        if(vrel < vbound) { /* bound */
                            if( All.ForceSoftening[5]*(1.0-vrel*vrel/(vbound*vbound))/sqrt(r2) > 1.0 ) { /* apocenter within 2.8*epsilon (softening length) */
                                
                                /* CAVEAT: when two BHs share some neighbours, this double counts the accretion */
                                if(P[j].SwallowID < id)
                                {
                                    out.mass_to_swallow_edd += P[j].Mass;
                                } /* P[j].SwallowID < id */
                            } /* if( All.ForceSoftening[5]*(1.0-vrel*vrel/(csnd*csnd))/sqrt(r2) > 1.0 ) */
                        } /* if(vrel < vbound) */
                    } /* type check */
#endif // BH_GRAVCAPTURE_GAS
                    
                    
                    
                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)
            
            if(mode == 0) /* local -> send directoy to local temp struct */
                out2particle_blackhole(&out, mod_index, 0);     /* target -> mod_index for reduced size struc */
            else
                BlackholeDataPasserResult[target] = out;        /* this is okay because target cycles over nimport only */
            
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];   /* non-local target o.k. */
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            } // if(listindex < NODELISTLENGTH)
        } // if(mode == 1)
    } // while(startnode >= 0)
    return 0;
}
