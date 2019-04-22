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
#if defined(SINGLE_STAR_STRICT_ACCRETION) || defined(NEWSINK)
    MyFloat SinkRadius;
#endif  
    MyDouble Pos[3];
    MyFloat Vel[3];
#ifdef BH_GRAVACCRETION
    MyFloat Jgas[3];
    MyFloat Jstar[3];
#endif
    MyFloat Hsml;
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
#ifdef BH_WAKEUP_GAS
    MyFloat TimeBin;
#endif
//#if !defined(SINGLE_STAR_STRICT_ACCRETION)
//    MyFloat SinkRadius;
//#endif
#if defined(NEWSINK_J_FEEDBACK)
    MyFloat Jsink[3];
#endif
}
*BlackholeDataIn, *BlackholeDataGet;



void blackhole_environment_loop(void)
{
    int i, j, k, nexport, nimport, place, ngrp, recvTask, dummy;
    int ndone_flag, ndone;
    MPI_Status status;
    
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackhole_temp_particle_data) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackhole_temp_particle_data))));
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
#ifdef SINGLE_STAR_STRICT_ACCRETION
            BlackholeDataIn[j].SinkRadius = P[place].SinkRadius;
#endif	    
#endif
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].ID = P[place].ID;
#if defined(NEWSINK)
	    BlackholeDataIn[j].SinkRadius = P[place].SinkRadius;
#endif
#ifdef BH_WAKEUP_GAS
	    BlackholeDataIn[j].TimeBin = P[place].TimeBin;
#endif
#if defined(NEWSINK_J_FEEDBACK)
                BlackholeDataIn[j].Jsink[0] = P[place].Jsink[0];BlackholeDataIn[j].Jsink[1] = P[place].Jsink[1];BlackholeDataIn[j].Jsink[2] = P[place].Jsink[2];
#endif
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
                                 recvTask, TAG_BH_A,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_A, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
   
        /*  DAA: PasserResult cycles over nimport!  PasserOut cycles over nexport! */
        BlackholeDataPasserResult = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserResult", nimport * sizeof(struct blackhole_temp_particle_data));
        BlackholeDataPasserOut = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserOut", nexport * sizeof(struct blackhole_temp_particle_data));
        
        memset( &BlackholeDataPasserResult[0], 0, nimport * sizeof(struct blackhole_temp_particle_data)  );
        memset( &BlackholeDataPasserOut[0],    0, nexport * sizeof(struct blackhole_temp_particle_data)  );

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
                                 MPI_BYTE, recvTask, TAG_BH_B,
                                 &BlackholeDataPasserOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                 MPI_BYTE, recvTask, TAG_BH_B, MPI_COMM_WORLD, &status);
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

    /* DAA: normalize/finalized some of the environment variables */
    for(i=0; i<N_active_loc_BHs; i++)
        normalize_temp_info_struct(i);
        
}




/* routine to return the values we need of the properties of the gas, stars, etc in the vicinity of the BH -- these all factor into the BHAR */
int blackhole_environment_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, k, n, listindex=0, mod_index;
    MyFloat *pos, h_i, *vel, hinv;
    MyIDType id;
    
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(NEWSINK)
    MyFloat hinv3, wk, dwk, u; u=wk=dwk=0;
#endif
    
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat mass, vrel, vbound, r2, sink_radius; sink_radius=0;
#endif
#if defined(NEWSINK)
    MyFloat csound_sq, hinv_gas1, hinv3_gas1, u_gas1, u_gas2,hinv_gas2, hinv3_gas2, int_zone_radius;
    integertime bh_timebin;
    int j2, n2;
    MyDouble Jpar[3], dt;
#if defined(NEWSINK_J_FEEDBACK)
    MyFloat Jsinktot, Jcrossdr[3], drcrossJcrossdr[3];
    MyFloat *Jsink;
#endif
#endif
    
    MyDouble dP[3],dv[3],wt;
    struct blackhole_temp_particle_data out;
    memset(&out, 0, sizeof(struct blackhole_temp_particle_data));
    
    /* these are the BH properties */
    if(mode == 0)
    {
#if defined(BH_GRAVCAPTURE_GAS)
        mass = P[target].Mass;
#if defined(SINGLE_STAR_STRICT_ACCRETION) || defined(NEWSINK)
        sink_radius = P[target].SinkRadius;
#endif	
#endif
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
#ifdef BH_WAKEUP_GAS
	bh_timebin = P[target].TimeBin;
#endif
#if defined(NEWSINK)
        int_zone_radius = P[target].Hsml * INT_ZONE_TO_HSML;
#if defined(NEWSINK_J_FEEDBACK)
        Jsink = P[target].Jsink;
#endif
#endif
    }
    else
    {
#if defined(BH_GRAVCAPTURE_GAS)
        mass = BlackholeDataGet[target].Mass;
#if defined(SINGLE_STAR_STRICT_ACCRETION) || defined(NEWSINK)
        sink_radius = BlackholeDataGet[target].SinkRadius;
#endif		
#endif
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        mod_index = 0;                              /* this is not used for mode==1, but this avoids compiler error */
#ifdef BH_WAKEUP_GAS
	bh_timebin = BlackholeDataGet[target].TimeBin;
#endif 
#if defined(NEWSINK)
        int_zone_radius = BlackholeDataGet[target].Hsml * INT_ZONE_TO_HSML;
#if defined(NEWSINK_J_FEEDBACK)
        Jsink = BlackholeDataGet[target].Jsink;
#endif
#endif
    }
    
    if(h_i < 0) return -1;
    hinv = 1./h_i;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(NEWSINK)
    hinv3 = hinv*hinv*hinv;
#endif
#if defined(NEWSINK)
        MyDouble dt_min = DMAX(DT_MIN_TOLERANCE_FACTOR * sqrt(pow(All.SofteningTable[5],3.0)/(All.G * mass)), 20.0*DMAX(All.MinSizeTimestep,All.Timebase_interval) );
#if defined(NEWSINK_J_FEEDBACK)
        Jsinktot = sqrt(Jsink[0]*Jsink[0] + Jsink[1]*Jsink[1] +Jsink[2]*Jsink[2]);
#endif
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
#if defined(NEWSINK)
            numngb = ngb_treefind_variable_targeted(pos, int_zone_radius, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG);	   
#else
            numngb = ngb_treefind_variable_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
#endif
        if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];

#ifdef BH_WAKEUP_GAS
		if (bh_timebin < P[j].LowestBHTimeBin) P[j].LowestBHTimeBin = bh_timebin;
//		SphP[j].wakeup = 1;
#endif		
                if( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
                {
                    wt = P[j].Mass;
                    dP[0] = P[j].Pos[0]-pos[0]; dP[1] = P[j].Pos[1]-pos[1]; dP[2] = P[j].Pos[2]-pos[2];
#ifdef BOX_PERIODIC
                    NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                    dv[0] = P[j].Vel[0]-vel[0]; dv[1] = P[j].Vel[1]-vel[1]; dv[2] = P[j].Vel[2]-vel[2];
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif

#ifdef BH_DYNFRICTION
#if (BH_DYNFRICTION == 1)    // DAA: dark matter + stars
                    if( !(P[j].Type==0) )
#if (BH_REPOSITION_ON_POTMIN == 2)
                    if( (P[j].Type != 5) )
#endif
#elif (BH_DYNFRICTION == 2)  // DAA: stars only
                    if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) )
#endif
                    {
                        double wtfac = wt;
#if (BH_REPOSITION_ON_POTMIN == 2)
                        double rfac = (dP[0]*dP[0] + dP[1]*dP[1] + dP[2]*dP[2]) * (10./(h_i*h_i) + 0.1/(All.ForceSoftening[5]*All.ForceSoftening[5]));
                        wtfac = wt / (1. + rfac); // simple function scaling ~ 1/r^2 for large r, to weight elements closer to the BH, so doesnt get 'pulled' by far-away elements //
#endif
                        if(P[j].Mass>out.DF_mmax_particles) out.DF_mmax_particles=P[j].Mass;
                        for (k=0;k<3;k++)
                        {
                            out.DF_mean_vel[k] += wt*dv[k];
#if (BH_REPOSITION_ON_POTMIN == 2)
                            out.DF_rms_vel += wt;
#else
                            out.DF_rms_vel += wt*dv[k]*dv[k];
#endif
                        }
                    }
#endif
                    
                    /* DAA: compute mass/angular momentum for GAS/STAR/DM components within BH kernel
                            this is done always now (regardless of the specific BH options used) */
                    if(P[j].Type==0)
                    {
                        /* we found gas in BH's kernel */
                        out.Mgas_in_Kernel += wt;
                        out.Sfr_in_Kernel += SphP[j].Sfr;
                        out.BH_InternalEnergy += wt*SphP[j].InternalEnergy;
                        out.Jgas_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jgas_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jgas_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
                        u=0;
                        for(k=0;k<3;k++) u+=dP[k]*dP[k];
                        u=sqrt(u)/h_i;
                        kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,1);
                        dwk /= u*h_i;
                        for(k=0;k<3;k++) out.GradRho_in_Kernel[k] += wt * dwk * fabs(dP[k]);
#endif
#if defined(BH_BONDI) || defined(BH_DRAG) || (BH_GRAVACCRETION == 5)
                        for(k=0;k<3;k++) {out.BH_SurroundingGasVel[k] += wt*dv[k];}
#endif
                    }
                    else if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) ) 
                    {
                        /* stars */
                        out.Mstar_in_Kernel += wt;
                        out.Jstar_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jstar_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jstar_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                    else 
                    { 
                        /* dark matter */
                        // DAA: Jalt_in_Kernel and Malt_in_Kernel are updated in normalize_temp_info_struct() to be TOTAL angular momentum and mass 
                        out.Malt_in_Kernel += wt;
                        out.Jalt_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jalt_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jalt_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                    
                    

#if defined(BH_GRAVCAPTURE_GAS)
                    /* XM: I formally distinguish BH_GRAVCAPTURE_GAS and BH_GRAVCAPTURE_NONGAS. The former applies to
                     gas ONLY, as an accretion model. The later can be combined with any accretion model.
                     Currently, I only allow gas accretion to contribute to BH_Mdot (consistent with the energy radiating away).
                     For star particles, if there is an alpha-disk, they are captured to the disk. If not, they directly go
                     to the hole, without any contribution to BH_Mdot and feedback. This can be modified in the swallow loop
                     for other purposes. The goal of the following part is to estimate BH_Mdot, which will be used to evaluate feedback strength.
                     Therefore, we only need it when we enable BH_GRAVCAPTURE_GAS as gas accretion model. */
                    if( (P[j].Mass > 0) && (P[j].Type == 0))
                    {
                        vrel=0; for(k=0;k<3;k++) {vrel += (P[j].Vel[k] - vel[k])*(P[j].Vel[k] - vel[k]);}
                        r2=0; for(k=0;k<3;k++) {r2+=dP[k]*dP[k];}
                        double dr_code = sqrt(r2); vrel = sqrt(vrel) / All.cf_atime; vbound = bh_vesc(j, mass, dr_code);
#ifdef SINGLE_STAR_STRICT_ACCRETION
			//			if(dr_code < DMAX(sink_radius, Get_Particle_Size(j)))
			if(dr_code < sink_radius)			
#endif			
                        if(vrel < vbound) { /* bound */
#ifdef SINGLE_STAR_STRICT_ACCRETION
                            double spec_mom=0; for(k=0;k<3;k++) {spec_mom += (P[j].Vel[k] - vel[k])*dP[k];} // delta_x.delta_v
                            spec_mom = (r2*vrel*vrel - spec_mom*spec_mom*All.cf_a2inv);  // specific angular momentum^2 = r^2(delta_v)^2 - (delta_v.delta_x)^2;
			    //                            if(spec_mom < All.G * (mass + P[j].Mass) * DMAX(Get_Particle_Size(j),sink_radius)) // check Bate 1995 angular momentum criterion (in addition to bounded-ness)
			    if(spec_mom < All.G * (mass + P[j].Mass) * sink_radius)
#endif
                            if( bh_check_boundedness(j,vrel,vbound,dr_code,sink_radius)==1 ) { /* apocenter within 2.8*epsilon (softening length) */
#ifdef SINGLE_STAR_FORMATION
			      double spec_energy = 0.5*(vrel*vrel - vbound*vbound); // specific energy of the 2-body system
			      if (spec_energy < P[j].SwallowEnergy) P[j].SwallowEnergy = spec_energy; // record the lowest energy we encounter, so we can accrete onto the sink that we are most tightly bound to if there are overlapping sinks
#endif
                                /* CAVEAT: when two BHs share some neighbours, this double counts the accretion. looks like this is true always since SwallowID=0 has just been initialized... only makes sense to check SwallowID if we update it... */
                                if(P[j].SwallowID < id) {out.mass_to_swallow_edd += P[j].Mass;} /* P[j].SwallowID < id */
                            } /* if( apocenter in tolerance range ) */
                        } /* if(vrel < vbound) */
                        
#if defined(NEWSINK)
                        if (dr_code <= int_zone_radius ) /*Check if gas in interaction radius*/
                        {
                            u=dr_code/h_i;
                            csound_sq = GAMMA*GAMMA_MINUS1 * SphP[j].InternalEnergyPred;
                            kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,-1);
                            wk = fabs(wk); /* Just to be safe */
                            out.intzone_gasmass += P[j].Mass; /* sum up all the mass in the sink radius*/ 
                            out.intzone_massweight_all += P[j].Mass / SphP[j].Density * wk; /* Volume weighted kernel sum of Eq 9 in Hubber 2013 */
#ifdef NEWSINK_BONDI			    
			    if(dr_code < All.ForceSoftening[5]) // Special switch to clean up gas that might want to hang around inside the softening radius
			     {
				 double bondi_mdot = All.G * All.G * mass * mass / pow(csound_sq + (dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]), 1.5) * P[j].Mass * wk;
				 out.t_rad_denom_sum += -DMAX(bondi_mdot,  -(dP[0]*dv[0] + dP[1]*dv[1] + dP[2]*dv[2]) * dr_code * P[j].Mass * wk);
				 out.min_bondi_mdot += 4. * M_PI * bondi_mdot;
				 out.gasmass_within_softening += P[j].Mass;
//#ifdef BH_OUTPUT_MOREINFO				 
//				 printf("Adding Bondi accretion rate %g for particle with r = %g, rho=%g and cs_eff=%g\n", 4*M_PI*bondi_mdot, dr_code, SphP[j].Density, sqrt(csound_sq + dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]));
//#endif				 
			     } else
#endif				
			    out.t_rad_denom_sum += (dP[0]*dv[0] + dP[1]*dv[1] + dP[2]*dv[2]) * dr_code * P[j].Mass * wk;// * DMAX(1, All.ForceSoftening[5]*All.ForceSoftening[5]/(DMAX(dr_code, P[j].Hsml)*DMAX(dr_code, P[j].Hsml))); /* sum in denominator of Eq 8 in Hbber 2013 */
			    
                            out.t_disc_num_sum += (sqrt(dr_code) * P[j].Mass * wk) / ( SphP[j].Density * csound_sq); /* sum in Eq 10 in Hubber 2013 without kernel weight*/
                            hinv_gas1 = 1.0/P[j].Hsml; hinv3_gas1 = hinv_gas1*hinv_gas1*hinv_gas1;
                            u_gas1=dr_code*hinv_gas1;
                            Jpar[0] = P[j].Mass*(dP[1]*dv[2] - dP[2]*dv[1]); /*angular momentum of particle relative to sink*/
                            Jpar[1] = P[j].Mass*(dP[2]*dv[0] - dP[0]*dv[2]);
                            Jpar[2] = P[j].Mass*(dP[0]*dv[1] - dP[1]*dv[0]);
                            out.gas_Erot_in_intzone += (Jpar[0]*Jpar[0] + Jpar[1]*Jpar[1] + Jpar[2]*Jpar[2])/(2.0*P[j].Mass*dr_code*dr_code); /* total rotational energy, L^2/(2 m r^2), not just bulk as in Eq 13 oh Hubber 2013 */			    
//                            out.gas_Egrav_in_intzone -= All.G * mass * P[j].Mass * 0.5 * (kernel_gravity(u, hinv, hinv3, -1) + kernel_gravity(u_gas1, hinv_gas1, hinv3_gas1, -1) ); /*Sink-gas interaction sum from Hubber 2013 Eq 14*/
			    out.gas_Egrav_in_intzone += grav_interaction_energy(dr_code, mass, P[j].Mass, All.ForceSoftening[5], P[j].Hsml);
                            //store properties of this neighbor particle
                            if (out.n_neighbor < NEWSINK_NEIGHBORMAX){
                                out.rgas[out.n_neighbor] = dr_code; //distance
                                out.xgas[out.n_neighbor] = dP[0]; //x coord
                                out.ygas[out.n_neighbor] = dP[1]; //y coord
                                out.zgas[out.n_neighbor] = dP[2]; //z coord
                                out.Hsmlgas[out.n_neighbor] = P[j].Hsml; //softening
                                out.mgas[out.n_neighbor] = P[j].Mass; //mass
                                out.gasID[out.n_neighbor] = P[j].ID; //unique ID, used for swallowing later
                                //get boundedness, this check is already done above but let's keep it just in case we change the code structure
                                if(vrel < vbound) { out.isbound[out.n_neighbor] = 1;}
#if defined(NEWSINK_EAT_SMALL_DT)
                                if ( out.isbound[out.n_neighbor]==1 ){ /*for bound gas get timestep of gas particle*/
#ifndef WAKEUP
                                    dt = (P[j].TimeBin ? (1 << P[j].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
                                    dt = P[j].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#ifdef NEWSINK_BONDI				    
                                    if (dt<dt_min){// || dr_code+P[j].Hsml < All.SofteningTable[5]*2.8) { /*Check if the timescale is too small or if the particle is inside the core of the gravitational softening radius. If yes, it gets eaten to avoid issues. */
#else
                                    if ( (dt<dt_min) || (sqrt(All.ErrTolIntAccuracy*r2*dr_code/(All.G * mass))<dt_min) || dr_code+P[j].Hsml < All.ForceSoftening[5]) { /*Check if the timescale is too small or if the particle is inside the core of the gravitational softening radius. If yes, it gets eaten to avoid issues. */
#endif					
                                        out.f_acc[out.n_neighbor] = 1.0;
#ifdef BH_OUTPUT_MOREINFO					
                                        printf("%d : Found gas with too small time step around BH with ID %d, gas id %d, gas mass %g, gas dt %g, BH_dt %g, dtmin %g, at distance %g while the interaction zone is %g \n", ThisTask, id, P[j].ID,P[j].Mass,dt,sqrt(All.ErrTolIntAccuracy*r2*dr_code/(All.G * mass)), dt_min, dr_code,int_zone_radius );
#endif					
                                    }
                                }
#endif
#if defined(NEWSINK_J_FEEDBACK)
                                /* We need a normalization factor for angular momentum feedback so we will go over all the neighbours*/
                                if (Jsinktot > 0 ){ /*Sink has angular momentum*/
                                    Jcrossdr[0] = -Jsink[2]*dP[1] + Jsink[1]*dP[2]; Jcrossdr[1] = Jsink[2]*dP[0] - Jsink[0]*dP[2];Jcrossdr[2] = -Jsink[1]*dP[0] + Jsink[0]*dP[1]; // L x dP cross product
                                    drcrossJcrossdr[0] = Jcrossdr[2]*dP[1] - Jcrossdr[1]*dP[2]; drcrossJcrossdr[1] = -Jcrossdr[2]*dP[0] + Jcrossdr[0]*dP[2];drcrossJcrossdr[2] = Jcrossdr[1]*dP[0] - Jcrossdr[0]*dP[1]; // dP x L x dP cross product
//                                    out.dv_ang_kick_norm[out.n_neighbor] = P[j].Mass * wk * sqrt(drcrossJcrossdr[0]*drcrossJcrossdr[0] + drcrossJcrossdr[1]*drcrossJcrossdr[1] + drcrossJcrossdr[2] * drcrossJcrossdr[2] );/*Normalization factor for angular momentum feedback kicks, see denominator of Eq 22 of Hubber 2013*/
				    out.dv_ang_kick_norm[out.n_neighbor] = P[j].Mass * sqrt(drcrossJcrossdr[0]*drcrossJcrossdr[0] + drcrossJcrossdr[1]*drcrossJcrossdr[1] + drcrossJcrossdr[2] * drcrossJcrossdr[2] ); 
                                }
#endif
                                out.n_neighbor +=1; //keep track of how many neighbors we have
                            }
                            else{
#ifdef BH_OUTPUT_MOREINFO				
                                printf("%d Gas neighbor number over limit for BH with ID %d Current neighbor number is %d\n", ThisTask, id, out.n_neighbor);
#endif				
                            }
                            // /* Start another cycle to get gravitational energy, this is very crude, should be replaced */
                            
                        } /* Sink radius check */
#endif // NEWSINK
                    } /* type check */
#endif // BH_GRAVCAPTURE_GAS
                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)
// #ifdef BH_OUTPUT_MOREINFO
// #if defined(NEWSINK)
    // printf("env n=%llu last wk is: %g \n", (unsigned long long) target, (MyFloat) wk);
    // printf("env n=%llu last cs sq is: %g \n", (unsigned long long) target, (MyFloat) csound_sq);
    // printf("env n=%llu BH h_i is: %g \n", (unsigned long long) target, (MyFloat) h_i);
    // printf("env n=%llu BH sink_radius is: %g \n", (unsigned long long) target, (MyFloat) sink_radius);
    // printf("env n=%llu BH int_zone_radius is: %g \n", (unsigned long long) target, (MyFloat) int_zone_radius);
    // printf("env n=%llu BH gasmass is: %g \n", (unsigned long long) target, (MyFloat) out.intzone_gasmass);
    // printf("env n=%llu BH t_disc_num_sum is: %g \n", (unsigned long long) target, (MyFloat) out.t_disc_num_sum);
    // printf("env n=%llu BH intzone_massweight_all is: %g \n", (unsigned long long) target, (MyFloat) out.intzone_massweight_all);
    // printf("env n=%llu BH t_rad_denom_sum is: %g \n", (unsigned long long) target, (MyFloat) out.t_rad_denom_sum);
    // printf("env n=%llu BH gas_Egrav_in_intzone is: %g \n", (unsigned long long) target, (MyFloat) out.gas_Egrav_in_intzone);
    // printf("env n=%llu BH gas_Erot_in_intzone is: %g \n", (unsigned long long) target, (MyFloat) out.gas_Erot_in_intzone);
    // printf("env n=%llu BH mass is: %g \n", (unsigned long long) target, (MyFloat) mass);
    // printf("env n=%llu last neighbor count is: %d while ngbnum is: %d\n", (unsigned long long) target, out.n_neighbor, numngb);
// #endif
// #endif
            if(mode == 0) /* local -> send directly to local temp struct */
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





/* -----------------------------------------------------------------------------------------------------
 * DAA: modified versions of blackhole_environment_loop and blackhole_environment_evaluate for a second
 * environment loop. Here we do a Bulge-Disk kinematic decomposition for gravitational torque accretion
 * ----------------------------------------------------------------------------------------------------- 
 */


#ifdef BH_GRAVACCRETION


void blackhole_environment_second_loop(void)
{
    int i, j, k, nexport, nimport, place, ngrp, recvTask, dummy;
    int ndone_flag, ndone;
    int mod_index;
    MPI_Status status;

    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackhole_temp_particle_data) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackhole_temp_particle_data))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* Scan gas particles for angular momentum, boundedness, etc */
    i = FirstActiveParticle;    /* first particle for this task */
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
                if(blackhole_environment_second_evaluate(i, 0, &nexport, Send_count) < 0)
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
            mod_index = P[place].IndexMapToTempStruc;  /* DAA: the index of the BlackholeTempInfo structure */
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
                BlackholeDataIn[j].Jgas[k] = BlackholeTempInfo[mod_index].Jgas_in_Kernel[k];    // DAA: Jgas is available after first environment loop
                BlackholeDataIn[j].Jstar[k] = BlackholeTempInfo[mod_index].Jstar_in_Kernel[k];
            }
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
                                 recvTask, TAG_BH_C,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_C, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);

        //DAA: PasserResult cycles over nimport, PasserOut cycles over nexport:
        BlackholeDataPasserResult = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserResult", nimport * sizeof(struct blackhole_temp_particle_data));
        BlackholeDataPasserOut = (struct blackhole_temp_particle_data *) mymalloc("BlackholeDataPasserOut", nexport * sizeof(struct blackhole_temp_particle_data));

 /* -------------- NOTE
 *         DAA: note that we actually don't need the full BlackholeDataPasserResult/BlackholeDataPasserOut structures here... 
 *         --> we could only pass Mbulge_in_Kernel !
 *  -------------- NOTE
 */

        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_environment_second_evaluate(j, 1, &dummy, &dummy);

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
                                             MPI_BYTE, recvTask, TAG_BH_D,
                                             &BlackholeDataPasserOut[Send_offset[recvTask]],
                                             Send_count[recvTask] * sizeof(struct blackhole_temp_particle_data),
                                             MPI_BYTE, recvTask, TAG_BH_D, MPI_COMM_WORLD, &status);
                            }
                        }
                    } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++)

        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            mod_index = P[place].IndexMapToTempStruc;
            /* DAA: we only need to update Mbulge_in_Kernel 
            out2particle_blackhole(&BlackholeDataPasserOut[j], P[place].IndexMapToTempStruc, 1); */
            BlackholeTempInfo[mod_index].MgasBulge_in_Kernel += BlackholeDataPasserOut[j].MgasBulge_in_Kernel;
            BlackholeTempInfo[mod_index].MstarBulge_in_Kernel += BlackholeDataPasserOut[j].MstarBulge_in_Kernel;
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




int blackhole_environment_second_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, n, listindex=0, mod_index;
    MyFloat *pos, h_i, *vel, *Jgas, *Jstar;
    MyIDType id;

    double dP[3],dv[3],J_tmp[3],wt;
    double MgasBulge_tmp, MstarBulge_tmp;
    MgasBulge_tmp=0; MstarBulge_tmp=0;

    /* these are the BH properties */
    if(mode == 0)
    {
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
        Jgas = BlackholeTempInfo[mod_index].Jgas_in_Kernel;      // DAA: Jgas/Jstar available after first environment loop
        Jstar = BlackholeTempInfo[mod_index].Jstar_in_Kernel;
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        mod_index = 0;                              /* this is not used for mode==1, but this avoids compiler error */
        Jgas = BlackholeDataGet[target].Jgas;
        Jstar = BlackholeDataGet[target].Jstar;
    }

    if(h_i < 0) return -1;

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
            numngb = ngb_treefind_variable_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
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
#ifdef BOX_PERIODIC           
                    NEAREST_XYZ(dP[0],dP[1],dP[2],-1); /*  find the closest image in the given box size  */
#endif
                    dv[0] = P[j].Vel[0]-vel[0];
                    dv[1] = P[j].Vel[1]-vel[1];
                    dv[2] = P[j].Vel[2]-vel[2];
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif

                    // DAA: jx,jy,jz, are independent of 'a' because ~ m*r*v, vphys=v/a, rphys=r*a //
                    J_tmp[0] = wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                    J_tmp[1] = wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                    J_tmp[2] = wt*(dP[0]*dv[1] - dP[1]*dv[0]);

                    if(P[j].Type==0)
                    {
                        /* GAS */
                        if( (J_tmp[0]*Jgas[0] + J_tmp[1]*Jgas[1] + J_tmp[2]*Jgas[2]) < 0 )  
                        {
                            /* DAA: assume the bulge component contains as many particles with positive azimuthal velocities  
                               as with negative azimuthal velocities relative to the angular momentum vector */
                            MgasBulge_tmp += 2.0 * wt;
                        }
                    }
                    else if( P[j].Type==4 || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn)) )
                    {
                        /* STARS */
                        if( (J_tmp[0]*Jstar[0] + J_tmp[1]*Jstar[1] + J_tmp[2]*Jstar[2]) < 0 )   
                        {
                            MstarBulge_tmp += 2.0 * wt;
                        }
                    }


                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)

            if(mode == 0) /* local -> send directly to local temp struct */
            {
                BlackholeTempInfo[mod_index].MgasBulge_in_Kernel = MgasBulge_tmp;
                BlackholeTempInfo[mod_index].MstarBulge_in_Kernel = MstarBulge_tmp;
            }
            else
            {
                BlackholeDataPasserResult[target].MgasBulge_in_Kernel = MgasBulge_tmp;
                BlackholeDataPasserResult[target].MstarBulge_in_Kernel = MstarBulge_tmp;
            }

        } // while(startnode >= 0)

        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];   /* non-local target o.k. */
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;  /* open it */
            } // if(listindex < NODELISTLENGTH)
        } // if(mode == 1)
    } // while(startnode >= 0)
    return 0;
}


#endif   //BH_GRAVACCRETION
