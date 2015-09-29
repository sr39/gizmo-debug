#include <stdio.h>
#include "../../proto.h"
#include "../../allvars.h"

/*! \file blackhole_util.c
 *  \brief util routines for memory (de)allocation and array setting for black holes
 */
/*
 * This file was largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   by Paul Torrey (ptorrey@mit.edu) for clairity.  It was rearranged and parsed into
 *   smaller files and routines. The main functional difference is that BlackholeTempInfo
 *   is now allocated only for N_active_loc_BHs, rather than NumPart (as was done before).  Some
 *   extra index gymnastics are required to follow this change through in the MPI comm routines.
 *   Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

/* function for allocating temp BH data struc needed for feedback routines*/
void blackhole_start(void)
{
    int i, j, Nbh;
    
    /* count the num BHs on this task */
    N_active_loc_BHs=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==5)
        {
            P[i].IndexMapToTempStruc = N_active_loc_BHs;         /* allows access via BlackholeTempInfo[P[i].IndexMapToTempStruc] */
            N_active_loc_BHs++;                     /* N_active_loc_BHs now set for BH routines */
        }
    }
    
    /* allocate the blackhole temp struct -- defined in blackhole.h */
    if(N_active_loc_BHs>0)
    {
        BlackholeTempInfo = mymalloc("BlackholeTempInfo", N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data));
    } else {
        BlackholeTempInfo = mymalloc("BlackholeTempInfo", 1 * sizeof(struct blackhole_temp_particle_data));    // allocate dummy, necessary?
    }
    
    Nbh=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==5)
        {
            BlackholeTempInfo[Nbh].index = i;               /* only meaningful field set here */
            BlackholeTempInfo[Nbh].BH_InternalEnergy = 0;
            BlackholeTempInfo[Nbh].accreted_Mass = 0;
            BlackholeTempInfo[Nbh].accreted_BH_Mass = 0;
            BlackholeTempInfo[Nbh].Mgas_in_Kernel=0;
            BlackholeTempInfo[Nbh].Malt_in_Kernel=0;
#ifdef BH_GRAVACCRETION_BTOD
            BlackholeTempInfo[Nbh].Mbulge_in_Kernel=0;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            BlackholeTempInfo[Nbh].BH_angle_weighted_kernel_sum=0;
#endif
#ifdef BH_DYNFRICTION
            BlackholeTempInfo[Nbh].DF_rms_vel=0;
            BlackholeTempInfo[Nbh].DF_mmax_particles=0;
#endif
            for(j=0;j<3;j++)
            {
                BlackholeTempInfo[Nbh].Jalt_in_Kernel[j]=0;
                BlackholeTempInfo[Nbh].accreted_momentum[j]=0;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_GRAVACCRETION)  
                BlackholeTempInfo[Nbh].Jgas_in_Kernel[j]=0;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeTempInfo[Nbh].GradRho_in_Kernel[j]=0;
#endif
#ifdef BH_DYNFRICTION
                BlackholeTempInfo[Nbh].DF_mean_vel[j]=0;
#endif
#if defined(BH_BONDI) || defined(BH_DRAG)
                BlackholeTempInfo[Nbh].BH_SurroundingGasVel[j]=0;
#endif
            }
            Nbh++;
        }
    }
    
    /* all future loops can now take the following form:
     for(i=0; i<N_active_loc_BHs; i++)
     {
     i_old = BlackholeTempInfo[i].index;
     ...
     }
     */
    
}


/* function for freeing temp BH data struc needed for feedback routines*/
void blackhole_end(void)
{
    int bin;
    double mdot, mdot_in_msun_per_year;
    double mass_real, total_mass_real, medd, total_mdoteddington;
    double mass_holes, total_mass_holes, total_mdot;
    
    /* sum up numbers to print for summary of the BH step (blackholes.txt) */
    mdot = mass_holes = mass_real = medd = 0;
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        if(TimeBinCount[bin])
        {
            mass_holes += TimeBin_BH_mass[bin];
            mass_real += TimeBin_BH_dynamicalmass[bin];
            mdot += TimeBin_BH_Mdot[bin];
            medd += TimeBin_BH_Medd[bin];
        }
    }
    MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        /* convert to solar masses per yr */
        mdot_in_msun_per_year = total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
        total_mdoteddington *= 1.0 / ((4 * M_PI * GRAVITY * C * PROTONMASS /
                                       (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON)) * All.UnitTime_in_s / All.HubbleParam);
        fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
                All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
                total_mass_real, total_mdoteddington);
        fflush(FdBlackHoles);
    }
    
    fflush(FdBlackHolesDetails);
#ifdef BH_OUTPUT_MOREINFO
    fflush(FdBhMergerDetails);
#ifdef BH_BAL_KICK
    fflush(FdBhWindDetails);
#endif
#endif
    myfree(BlackholeTempInfo);
}




/* simple routine to add quantities to BlackholeTempInfo */
void out2particle_blackhole(struct blackhole_temp_particle_data *out, int target, int mode)
{
    int k;
    ASSIGN_ADD(BlackholeTempInfo[target].BH_InternalEnergy,out->BH_InternalEnergy,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Mgas_in_Kernel,out->Mgas_in_Kernel,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Malt_in_Kernel,out->Malt_in_Kernel,mode);
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeTempInfo[target].Jalt_in_Kernel[k],out->Jalt_in_Kernel[k],mode);
#ifdef BH_DYNFRICTION
    ASSIGN_ADD(BlackholeTempInfo[target].DF_rms_vel,out->DF_rms_vel,mode);
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeTempInfo[target].DF_mean_vel[k],out->DF_mean_vel[k],mode);
    if(mode==0)
        BlackholeTempInfo[target].DF_mmax_particles = out->DF_mmax_particles;
    else
        if(out->DF_mmax_particles > BlackholeTempInfo[target].DF_mmax_particles)
            BlackholeTempInfo[target].DF_mmax_particles = out->DF_mmax_particles;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_GRAVACCRETION)  // DAA: need Jgas for GRAVACCRETION as well
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeTempInfo[target].Jgas_in_Kernel[k],out->Jgas_in_Kernel[k],mode);
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        ASSIGN_ADD(BlackholeTempInfo[target].GradRho_in_Kernel[k],out->GradRho_in_Kernel[k],mode);
#endif
    }
#endif
#if defined(BH_BONDI) || defined(BH_DRAG)
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeTempInfo[target].BH_SurroundingGasVel[k],out->BH_SurroundingGasVel[k],mode);
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    ASSIGN_ADD(BlackholeTempInfo[target].mass_to_swallow_edd, out->mass_to_swallow_edd, mode);
#endif
}



int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local)
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
            
            /* make sure we get all the particle types we need */
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_NONGAS) || defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_DYNFRICTION)
            if(P[p].Type < 0)
                continue;
#else
            if(P[p].Type != 0 && P[p].Type != 5)
                continue;
#endif
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
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
            
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist = hsml + 0.5 * current->len;;
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
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

