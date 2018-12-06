/*! \file blackhole_local.h
 *  \brief comm structures to be used throughout blackhole routines
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


/* quantities that pass IN to the 'blackhole_evaluate' routines */
static struct blackholedata_in
{
    int orig_index;
    int trans_index;
    
    MyDouble Pos[3];
    MyFloat Density;
    MyFloat Mdot;
    MyFloat Dt;
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat BH_Mass;
#ifdef SINGLE_STAR_STRICT_ACCRETION
    MyFloat SinkRadius;
#endif  
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
    MyFloat Vel[3];
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
    MyFloat Jgas_in_Kernel[3];
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
    MyFloat BH_disk_hr;
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat mass_to_swallow_edd;
#endif
#if defined(NEWSINK)
#if !defined(SINGLE_STAR_STRICT_ACCRETION)
    MyFloat SinkRadius;
#endif
    /* properties of neighboring particles, used for preferential feeding */
    int n_neighbor; //number of neighbors currently stored in the arrays below
    MyFloat rgas[NEWSINK_NEIGHBORMAX]; /* Distance of gas from sink */
    MyFloat xgas[NEWSINK_NEIGHBORMAX]; /* x coordinate of gas from sink */
    MyFloat ygas[NEWSINK_NEIGHBORMAX]; /* y coordinate of gas from sink */
    MyFloat zgas[NEWSINK_NEIGHBORMAX]; /* z coordinate of gas from sink */
    MyFloat Hsmlgas[NEWSINK_NEIGHBORMAX]; /* gas smoothing length */
    MyFloat mgas[NEWSINK_NEIGHBORMAX]; /* Mass of gas particle */
    MyIDType gasID[NEWSINK_NEIGHBORMAX]; /* ID of gas particle */
    int isbound[NEWSINK_NEIGHBORMAX]; /* is it bound to the sink */
    MyFloat f_acc[NEWSINK_NEIGHBORMAX]; /* How much of the gas particle should be accreted */
#if defined(NEWSINK_J_FEEDBACK)
    MyDouble Jsink[3];
    MyFloat t_disc;
    MyDouble dv_ang_kick_norm[NEWSINK_NEIGHBORMAX]; /*Normalization term for angular momentum feedback kicks, see denominator of Eq 22 of Hubber 2013*/
#endif
#endif
}
*BlackholeDataIn, *BlackholeDataGet;


/* quantities that pass OUT of the 'blackhole_evaluate' routines */
static struct blackholedata_out
{
    
    int orig_index;
    int trans_index;
    
    MyLongDouble Mass;
    MyLongDouble BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    MyLongDouble BH_Mass_AlphaDisk;
#endif
    MyLongDouble accreted_Mass;
    MyLongDouble accreted_BH_Mass;
    MyLongDouble accreted_momentum[3];
#ifdef BH_REPOSITION_ON_POTMIN
    MyFloat BH_MinPotPos[3];
    MyFloat BH_MinPot;
#endif
#ifdef BH_COUNTPROGS
    int BH_CountProgs;
#endif
#ifdef GALSF
    MyFloat Accreted_Age;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_WIND_CONTINUOUS)
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#if defined(NEWSINK)
    MyFloat SinkRadius;
    MyFloat Mdot;
    MyFloat Dt;
    MyFloat f_acc[NEWSINK_NEIGHBORMAX]; /* How much of the gas particle should be accreted */
#if defined(NEWSINK_J_FEEDBACK)
    MyLongDouble accreted_J[3];
#endif
#endif
}
*BlackholeDataResult, *BlackholeDataOut;


