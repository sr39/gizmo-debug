#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../../allvars.h"
#include "../../proto.h"

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
    int i, Nbh;
    
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
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data));
    } else {
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", 1 * sizeof(struct blackhole_temp_particle_data));
    }
    
    memset( &BlackholeTempInfo[0], 0, N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data) );
    
    Nbh=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==5)
        {
            BlackholeTempInfo[Nbh].index = i;               /* only meaningful field set here */
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
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
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
        if((ThisTask == 0) && (total_mdot > 0) && (total_mass_real > 0))
        {
            /* convert to solar masses per yr */
            mdot_in_msun_per_year = total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
            total_mdoteddington *= 1.0 / bh_eddington_mdot(1);
            fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
                    All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
                    total_mass_real, total_mdoteddington);
        }
//#ifdef IO_REDUCED_MODE  DAA-IO: this is redundant
//        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
//#endif
//        {fflush(FdBlackHoles);} 
        fflush(FdBlackHoles);

//#ifndef IO_REDUCED_MODE   DAA-IO: BH_OUTPUT_MOREINFO overrides IO_REDUCED_MODE
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
        fflush(FdBlackHolesDetails);
#ifdef BH_OUTPUT_MOREINFO
        fflush(FdBhMergerDetails);
#ifdef BH_WIND_KICK
        fflush(FdBhWindDetails);
#endif
#endif
#endif
    }
    myfree(BlackholeTempInfo);
}




/* simple routine to add quantities to BlackholeTempInfo */
void out2particle_blackhole(struct blackhole_temp_particle_data *out, int target, int mode)
{
    int k;
    ASSIGN_ADD(BlackholeTempInfo[target].BH_InternalEnergy,out->BH_InternalEnergy,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Mgas_in_Kernel,out->Mgas_in_Kernel,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Mstar_in_Kernel,out->Mstar_in_Kernel,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Malt_in_Kernel,out->Malt_in_Kernel,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Sfr_in_Kernel,out->Sfr_in_Kernel,mode);
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeTempInfo[target].Jgas_in_Kernel[k],out->Jgas_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeTempInfo[target].Jstar_in_Kernel[k],out->Jstar_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeTempInfo[target].Jalt_in_Kernel[k],out->Jalt_in_Kernel[k],mode);
    }
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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeTempInfo[target].GradRho_in_Kernel[k],out->GradRho_in_Kernel[k],mode);
    }
#endif
#if defined(BH_BONDI) || defined(BH_DRAG) || (BH_GRAVACCRETION == 5)
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeTempInfo[target].BH_SurroundingGasVel[k],out->BH_SurroundingGasVel[k],mode);
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    ASSIGN_ADD(BlackholeTempInfo[target].mass_to_swallow_edd, out->mass_to_swallow_edd, mode);
#endif
#if defined(NEWSINK)
    ASSIGN_ADD(BlackholeTempInfo[target].gas_Erot_in_intzone, out->gas_Erot_in_intzone, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].gas_Egrav_in_intzone, out->gas_Egrav_in_intzone, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].t_rad_denom_sum, out->t_rad_denom_sum, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].t_disc_num_sum, out->t_disc_num_sum, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].intzone_massweight_all, out->intzone_massweight_all, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].intzone_gasmass, out->intzone_gasmass, mode);
    
    /*  Do an insertion sort for the gas particle properties  */
    int j, n=BlackholeTempInfo[target].n_neighbor;
    if ( (n + out->n_neighbor) > NEWSINK_NEIGHBORMAX){ // We are not supposed to have more neighbors than NEWSINK_NEIGHBORMAX
            printf("%d Gas neighbor number over limit of NEWSINK_NEIGHBORMAX for BH ID %d Current neighbor number is %d and we want to add %d more. We are keeping the closest ones.\n", ThisTask, BlackholeTempInfo[target].index ,n, out->n_neighbor);
    }
    //Regardless, we will collect the closest NEWSINK_NEIGHBORMAX particles
    for(j=0;j<(out->n_neighbor);j++){ //go over all the incoming particle data
            if (n==0){ //first one just gets stored
                BlackholeTempInfo[target].rgas[0] = out->rgas[j];
                BlackholeTempInfo[target].xgas[0] = out->xgas[j];
                BlackholeTempInfo[target].ygas[0] = out->ygas[j];
                BlackholeTempInfo[target].zgas[0] = out->zgas[j];
                BlackholeTempInfo[target].Hsmlgas[0] = out->Hsmlgas[j];
                BlackholeTempInfo[target].mgas[0] = out->mgas[j];
                BlackholeTempInfo[target].gasID[0] = out->gasID[j];
                BlackholeTempInfo[target].f_acc[0] = out->f_acc[j];
                BlackholeTempInfo[target].isbound[0] = out->isbound[j];
#if defined(NEWSINK_J_FEEDBACK)
                BlackholeTempInfo[target].dv_ang_kick_norm[0] = out->dv_ang_kick_norm[j];
#endif
                n++;
            }
            else{ //we already have some data, do an insertion sort
                k = n-1;
                while (k >= 0 && BlackholeTempInfo[target].rgas[k] > (out->rgas[j]) ) 
                {
                    if (k+1 < NEWSINK_NEIGHBORMAX){ //we can store this one, if we can't this will be just overwritten
                        BlackholeTempInfo[target].rgas[k+1] = BlackholeTempInfo[target].rgas[k];
                        BlackholeTempInfo[target].xgas[k+1] = BlackholeTempInfo[target].xgas[k];
                        BlackholeTempInfo[target].ygas[k+1] = BlackholeTempInfo[target].ygas[k];
                        BlackholeTempInfo[target].zgas[k+1] = BlackholeTempInfo[target].zgas[k];
                        BlackholeTempInfo[target].Hsmlgas[k+1] = BlackholeTempInfo[target].Hsmlgas[k];
                        BlackholeTempInfo[target].mgas[k+1] = BlackholeTempInfo[target].mgas[k];
                        BlackholeTempInfo[target].gasID[k+1] = BlackholeTempInfo[target].gasID[k];
                        BlackholeTempInfo[target].f_acc[k+1] = BlackholeTempInfo[target].f_acc[k];
                        BlackholeTempInfo[target].isbound[k+1] = BlackholeTempInfo[target].isbound[k];
#if defined(NEWSINK_J_FEEDBACK)
                        BlackholeTempInfo[target].dv_ang_kick_norm[k+1] = BlackholeTempInfo[target].dv_ang_kick_norm[k];
#endif
                    }
                    k--; 
                }
                if (k+1 < NEWSINK_NEIGHBORMAX){ //we can store this one, if we can't we will just throw it away
                    BlackholeTempInfo[target].rgas[k+1] = out->rgas[j];
                    BlackholeTempInfo[target].xgas[k+1] = out->xgas[j];
                    BlackholeTempInfo[target].ygas[k+1] = out->ygas[j];
                    BlackholeTempInfo[target].zgas[k+1] = out->zgas[j];
                    BlackholeTempInfo[target].Hsmlgas[k+1] = out->Hsmlgas[j];
                    BlackholeTempInfo[target].mgas[k+1] = out->mgas[j];
                    BlackholeTempInfo[target].gasID[k+1] = out->gasID[j];
                    BlackholeTempInfo[target].f_acc[k+1] = out->f_acc[j];
                    BlackholeTempInfo[target].isbound[k+1] = out->isbound[j];
#if defined(NEWSINK_J_FEEDBACK)
                    BlackholeTempInfo[target].dv_ang_kick_norm[k+1] = out->dv_ang_kick_norm[j];
#endif
                }
                if (n< NEWSINK_NEIGHBORMAX){n++;} //we have added another neighbor
            }
    }
    BlackholeTempInfo[target].n_neighbor = n; //update the number of neighbors stored
#endif

}



