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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeTempInfo[Nbh].GradRho_in_Kernel[j]=0;
                BlackholeTempInfo[Nbh].Jgas_in_Kernel[j]=0;
#endif
#ifdef BH_DYNFRICTION
                BlackholeTempInfo[Nbh].DF_mean_vel[j]=0;
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeTempInfo[target].Jgas_in_Kernel[k],out->Jgas_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeTempInfo[target].GradRho_in_Kernel[k],out->GradRho_in_Kernel[k],mode);
    }
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeTempInfo[target].BH_SurroundingGasVel[k],out->BH_SurroundingGasVel[k],mode);
#endif
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)

    ASSIGN_ADD(BlackholeTempInfo[target].mass_to_swallow_total, out->mass_to_swallow_total, mode);
    ASSIGN_ADD(BlackholeTempInfo[target].mass_to_swallow_edd, out->mass_to_swallow_edd, mode);

#endif
}
