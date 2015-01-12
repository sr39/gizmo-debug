#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../../allvars.h"
#include "../../proto.h"

#include "../../kernel.h"


/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
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

#ifdef BLACK_HOLES

extern struct blackhole_temp_particle_data *BlackholeTempInfo;


/*  This is the master routine for the BH physics modules.
 *  It is called in calculate_non_standard_physics in run.c */
void blackhole_accretion(void)
{
    long i;
    int bin;
    double mdot, mdot_in_msun_per_year;
    double mass_real, total_mass_real, medd, total_mdoteddington;
    double mass_holes, total_mass_holes, total_mdot;
    
//#ifdef BH_GRAVACCRETION
//    double m_tmp_for_bhar;
//    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar;
//    double f0_for_bhar;
//#endif
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
    
    
    
    for(i = 0; i < NumPart; i++)
    {
        P[i].SwallowID = 0;
    }
        

        
    if(ThisTask == 0)  printf("Beginning black-hole accretion\n");
    blackhole_start();              /* allocates and cleans BlackholeTempInfo struct */

    
    
    /* this is the PRE-PASS loop.  BH does neighbor search, evaluates local gas properties */
    if(ThisTask == 0)  printf("Evaluating black-hole environment\n");
    blackhole_environment_loop();    /* populates BlackholeTempInfo based on surrounding gas (blackhole_environment.c).
                                        If using gravcap the desired mass accretion rate is calculated and set to BlackholeTempInfo.mass_to_swallow_edd
                                      */

    /*----------------------------------------------------------------------
     Now do a first set of local operations based on BH environment calculation:
     calculate mdot, dynamical friction, and other 'BH-centric' operations.
     No MPI comm necessary.
     ----------------------------------------------------------------------*/

    if(ThisTask == 0)  printf("Setting black-hole properties\n");
    blackhole_properties_loop();       /* do 'BH-centric' operations such as dyn-fric, mdot, etc. 
                                          This loop is at the end of this file.  */

    
    
    /*----------------------------------------------------------------------
     Now we perform a second pass over the black hole environment.
     Re-evaluate the decision to stochastically swallow gas if we exceed eddington.
     Use the above info to determine the weight functions for feedback
     ----------------------------------------------------------------------*/
    
    if(ThisTask == 0)  printf("Marking gas to swallow\n");
    blackhole_feed_loop();       /* BH mergers and gas/star/dm accretion events are evaluated
                                  - P[j].SwallowID's are set
                                  */
    
    
    
    /*----------------------------------------------------------------------
     Meow we do a THIRD pass over the particles, and
     this is where we can do the actual 'swallowing' operations
     (blackhole_evaluate_swallow), and 'kicking' operations
     ----------------------------------------------------------------------*/
    
    if(ThisTask == 0)  printf("Injecting feedback\n");
    blackhole_swallow_and_kick_loop();
    
    
    
    if(ThisTask == 0) printf("Doing whatever goes in the final loop\n");
    blackhole_final_loop();     /* this is causing problems with the alpha disk ?! */

    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
       Now do final operations on the results from the last pass
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/

    
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
    
    for(i = 0; i < NumPart; i++)
    {
        P[i].SwallowID = 0;
    }
    
    blackhole_end();            /* frees BlackholeTempInfo; cleans up */

    
}





int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
                           int *nexport, int *nsend_local)
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
#if defined(BH_REPOSITION_ON_POTMIN) || defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_NOGAS) || defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_DYNFRICTION)
            if(P[p].Type < 0)
                continue;
#else
            if(P[p].Type != 0 && P[p].Type != 5)
                continue;
#endif
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
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
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
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



#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
double bh_angleweight_localcoupling(int j, double hR, double theta)
{
#ifndef BH_PHOTONMOMENTUM
    // for now, if only BH_BAL_WINDS enabled, assume isotropic //
    return P[j].Hsml*P[j].Hsml;
#endif
    double b0,c0,f;
    // nathans 'B' and 'C' functions //
    b0=8.49403/(1.17286+hR);
    c0=64.4254/(2.5404+hR);
    f=1-(1+c0*exp(-b0*M_PI/2))/(1+c0*exp(-b0*(M_PI/2-theta)));
    return P[j].Hsml*P[j].Hsml * f;
    /* H^2 gives the fraction of the solid angle subtended by the particle (normalized by distance),
     the 'f' function gives the dForce/dOmega weighting from the radiative transfer calculations */
}


double bh_angleweight(double bh_lum_input, MyFloat bh_angle[3], double hR, double dx, double dy, double dz, int mode)
{
    double bh_lum = bh_lum_input;
    if(bh_lum <= 0) return 0;
    if(isnan(hR)) return 0;
    if(hR <= 0) return 0;
    double r2 = dx*dx+dy*dy+dz*dz;
    if(r2 <= 0) return 0;
    if(r2*All.UnitLength_in_cm*All.UnitLength_in_cm*All.cf_atime*All.cf_atime < 9.523e36) return 0; /* no force at < 1pc */
    
    double cos_theta = (dx*bh_angle[0] + dy*bh_angle[1] + dz*bh_angle[2])/sqrt(r2);
    if(cos_theta<0) cos_theta *= -1;
    if(isnan(cos_theta)) return 0;
    if(cos_theta <= 0) return 0;
    if(mode==1) bh_lum *= All.BlackHoleFeedbackFactor * All.BlackHoleRadiativeEfficiency * 4.597e20 * All.HubbleParam/All.UnitTime_in_s;
    if(cos_theta >= 1) return 1.441 * bh_lum;
    
    double hRe=hR; if(hRe<0.1) hRe=0.1; if(hRe>0.5) hRe=0.5;
    double y;
    y = -1.0 / (0.357-10.839*hRe+142.640*hRe*hRe-713.928*hRe*hRe*hRe+1315.132*hRe*hRe*hRe*hRe);
    y = 1.441 + (-6.42+9.92*hRe) * (exp(cos_theta*cos_theta*y)-exp(y)) / (1-exp(y));  // approximation to nathans fits
    //double A=5.57*exp(-hRe/0.52);double B=19.0*exp(-hRe/0.21);double C0=20.5+20.2/(1+exp((hRe-0.25)/0.035));
    //y = 1.441 + A*((1+C0*exp(-B*1.5708))/(1+C0*exp(-B*(1.5708-acos(cos_theta))))-1);
    // this is nathan's fitting function (fairly expensive with large number of exp calls and arc_cos
    //y=0.746559 - 9.10916*(-0.658128+exp(-0.418356*cos_theta*cos_theta));
    // this is normalized so the total flux is L/(4pi*r*r) and assumed monochromatic IR
    if(y>1.441) y=1.441; if(y<-5.0) y=-5.0;
    y*=2.3026; // so we can take exp, instead of pow //
    return exp(y) * bh_lum;
}
#endif /* end of #if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) */







void blackhole_properties_loop(void)
{
    int  i, n;
    double fac, mdot, dt;
    double medd;

//#ifdef BH_GRAVACCRETION
//    double m_tmp_for_bhar;
//    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar,mdisk_for_bhar;
//    double f0_for_bhar;
//#endif
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
#ifdef BH_ALPHADISK_ACCRETION
    double mdot_alphadisk;
#endif
    
    for(i=0; i<N_active_loc_BHs; i++)
    {
        n = BlackholeTempInfo[i].index;
        
        /* define the timestep */
#ifndef WAKEUP
        dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[n].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        
        /* always initialize/default to zero accretion rate */
        mdot=0;
        BPP(n).BH_Mdot=0;

        normalize_temp_info_struct(i);
        
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        set_blackhole_long_range_rp( i,  n);
#endif

        set_blackhole_mdot(i, n, dt);

#if defined(BH_DRAG) || defined(BH_DYNFRICTION)
        set_blackhole_drag(i, n, dt);
#endif
        
        set_blackhole_new_mass(i, n, dt);
        
        /* dump the results to the 'blackhole_details' files */
        fac=0; medd=0;
#ifdef BH_ALPHADISK_ACCRETION
        fac=BPP(n).BH_Mass_AlphaDisk;
        medd=mdot_alphadisk;			// is this really what we want?!? I don't think so...
#endif
        fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g %g %g %g   %2.7f %2.7f %2.7f\n",
                P[n].ID, All.Time, BPP(n).BH_Mass, fac, P[n].Mass, mdot, medd,
                BPP(n).DensAroundStar*All.cf_a3inv, BlackholeTempInfo[i].BH_InternalEnergy,
                P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]);

    }// for(i=0; i<N_active_loc_BHs; i++)
    
}










void normalize_temp_info_struct(int i)
{
#if defined(BH_DYNFRICTION) || ( defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG) )
    int k;
#endif
    /* for new quantities, divide out weights and convert to physical units */
    if(BlackholeTempInfo[i].Mgas_in_Kernel > 0)
    {
        BlackholeTempInfo[i].BH_InternalEnergy /= BlackholeTempInfo[i].Mgas_in_Kernel;
#ifdef BH_DYNFRICTION
        BlackholeTempInfo[i].DF_rms_vel /= BlackholeTempInfo[i].Mgas_in_Kernel;
        BlackholeTempInfo[i].DF_rms_vel = sqrt(BlackholeTempInfo[i].DF_rms_vel) / All.cf_atime;
        for(k=0;k<3;k++)
            BlackholeTempInfo[i].DF_mean_vel[k] /= BlackholeTempInfo[i].Mgas_in_Kernel * All.cf_atime;
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
        for(k=0;k<3;k++)
            BlackholeTempInfo[i].BH_SurroundingGasVel[k] /= BlackholeTempInfo[i].Mgas_in_Kernel * All.cf_atime;
#endif
    }
    else
    {
        BlackholeTempInfo[i].BH_InternalEnergy = 0;
    }
    
}


void set_blackhole_mdot(int i, int n, double dt)
{
    double mdot;
#ifdef BH_GRAVACCRETION
    int k;
    double m_tmp_for_bhar, mdisk_for_bhar, bh_mass;
    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar;
    double f0_for_bhar, fac;
#endif
//
//
//#ifdef BH_ALPHADISK_ACCRETION
//    double mdot_alphadisk;
//#endif
    
    
#ifdef BH_ENFORCE_EDDINGTON_LIMIT
    double meddington;
#endif
    
#ifdef BH_ENFORCE_EDDINGTON_LIMIT
    meddington = (4 * M_PI * GRAVITY * C * PROTONMASS /
                  (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON) ) *
    (BPP(n).BH_Mass/All.HubbleParam) * All.UnitTime_in_s;
#endif
    
#ifdef BH_GRAVACCRETION
    /* calculate mdot: gravitational instability accretion rate from Hopkins & Quataert 2011 */
    if(BlackholeTempInfo[i].Mgas_in_Kernel > 0)
    {
        m_tmp_for_bhar = BlackholeTempInfo[i].Mgas_in_Kernel + BlackholeTempInfo[i].Malt_in_Kernel;
        r0_for_bhar = PPP[n].Hsml * All.cf_atime; /* convert to physical units */
        j_tmp_for_bhar=0;
        for(k=0;k<3;k++)
            j_tmp_for_bhar += BlackholeTempInfo[i].Jalt_in_Kernel[k]*BlackholeTempInfo[i].Jalt_in_Kernel[k];
        j_tmp_for_bhar=sqrt(j_tmp_for_bhar);
        /* jx,y,z, is independent of 'a_scale' b/c ~ m*r*v, vphys=v/a, rphys=r*a */
        
        bh_mass = BPP(n).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass += BPP(n).BH_Mass_AlphaDisk;
#endif
        fgas_for_bhar = BlackholeTempInfo[i].Mgas_in_Kernel / m_tmp_for_bhar;
        fac = m_tmp_for_bhar * r0_for_bhar * sqrt(All.G*(m_tmp_for_bhar+bh_mass)/r0_for_bhar);
        /* All.G is G in code (physical) units */
        f_disk_for_bhar = fgas_for_bhar + (1.75*j_tmp_for_bhar/fac);
        if(f_disk_for_bhar>1) f_disk_for_bhar=1;
        
        if((f_disk_for_bhar<=0)||(bh_mass <=0)||(fgas_for_bhar<=0)||(m_tmp_for_bhar<=0))
        {
            mdot = 0;
        } else {
            mdisk_for_bhar = m_tmp_for_bhar*f_disk_for_bhar * (All.UnitMass_in_g/(All.HubbleParam * 1.0e9*SOLAR_MASS)); /* mdisk/1e9msun */
            bh_mass *= All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS); /* mbh/1e8msun */
            r0_for_bhar *= All.UnitLength_in_cm/(All.HubbleParam * 3.086e20); /* r0/100pc */
            f0_for_bhar = 0.31*f_disk_for_bhar*f_disk_for_bhar*pow(mdisk_for_bhar,-1./3.); /* dimensionless factor for equations */
            fac = (10.0*(SOLAR_MASS/All.UnitMass_in_g)/(SEC_PER_YEAR/All.UnitTime_in_s)); /* basic units */
            
            mdot = All.BlackHoleAccretionFactor * fac * mdisk_for_bhar *
            pow(f_disk_for_bhar,5./2.) * pow(bh_mass,1./6.) *
            pow(r0_for_bhar,-3./2.) / (1 + f0_for_bhar/fgas_for_bhar);
            
            printf("BH GravAcc Eval :: mdot %g BHaccFac %g Norm %g fdisk %g bh_8 %g fgas %g f0 %g mdisk_9 %g r0_100 %g \n\n",
                   mdot,All.BlackHoleAccretionFactor,fac,
                   f_disk_for_bhar,bh_mass,fgas_for_bhar,f0_for_bhar,mdisk_for_bhar,r0_for_bhar);fflush(stdout);
        } // if(f_disk_for_bhar<=0)
    }
#endif // ifdef BH_GRAVACCRETION
    
    
    
    
    
    
    
#ifdef BH_BONDI
    /* heres where we calculate the Bondi accretion rate, if that's going to be used */
    bhvel = 0;
#ifdef BH_USE_GASVEL_IN_BONDI
    for(k=0;k<3;k++) bhvel += BlackholeTempInfo[i].BH_SurroundingGasVel[k]*BlackholeTempInfo[i].BH_SurroundingGasVel[k];
#endif
    rho = BPP(n).DensAroundStar * All.cf_a3inv; /* we want all quantities in physical units */
    soundspeed = GAMMA*GAMMA_MINUS1 * BlackholeTempInfo[i].BH_InternalEnergy; // this is in physical units now
    fac = pow(soundspeed+bhvel, 1.5);
    if(fac > 0)
    {
        double AccretionFactor = All.BlackHoleAccretionFactor;
#ifdef BH_VARIABLE_ACCRETION_FACTOR
        /* variable-alpha model (Booth&Schaye 2009): now All.BlackHoleAccretionFactor is the slope of the density dependence */
        AccretionFactor = 1.0;
        if(rho > All.PhysDensThresh)
            AccretionFactor = pow(rho/All.PhysDensThresh, All.BlackHoleAccretionFactor);
#endif
        mdot = 4. * M_PI * AccretionFactor * All.G * All.G * BPP(n).BH_Mass * BPP(n).BH_Mass * rho / fac;
    }
#endif // ifdef BH_BONDI
    
    
    
    
    
    
    
#ifdef BH_GRAVCAPTURE_SWALLOWS
    mdot = 0; /* force mdot=0 despite any earlier settings here.  If this is set, we have to wait to swallow step to eval mdot. */
    mdot = BlackholeTempInfo[i].mass_to_swallow_edd / dt;       /* TODO: this can still greatly exceed eddington... */
#endif
    
    
    
    
    
    
    
#ifdef BH_ALPHADISK_ACCRETION
    /* use the mass in the accretion disk from the previous timestep to determine the BH accretion rate */
    BlackholeTempInfo[i].mdot_alphadisk = mdot;     /* if BH_GRAVCAPTURE_SWALLOWS is off, this gets the accretion rate */
    
    mdot = All.BlackHoleAccretionFactor *
    (2.45 * (SOLAR_MASS/All.UnitMass_in_g)/(SEC_PER_YEAR/All.UnitTime_in_s)) * // normalization
    pow( 0.1 , 8./7.) * // viscous disk 'alpha'
    pow( BPP(n).BH_Mass*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , -5./14. ) * // mbh dependence
    pow( BPP(n).BH_Mass_AlphaDisk*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , 10./7. ) * // m_disk dependence
    pow( DMIN(0.2,DMIN(PPP[n].Hsml,All.ForceSoftening[5])*All.cf_atime*All.UnitLength_in_cm/(All.HubbleParam * 3.086e18)) , -25./14. ); // r_disk dependence
    if(mdot<=0) mdot=0;
    if(dt>0)
    {
        if(mdot > BPP(n).BH_Mass_AlphaDisk/dt) mdot = BPP(n).BH_Mass_AlphaDisk/dt;
    }
#endif
    
    
    
    
    
    
    
#ifdef BH_SUBGRIDBHVARIABILITY
    /* account for sub-grid accretion rate variability */
    if((mdot>0)&&(dt>0)&&(P[n].DensAroundStar>0))
    {
        omega_ri=sqrt(All.G*P[n].DensAroundStar*All.cf_a3inv); /* dynamical frequency in physical units */
        n0_sgrid_elements=10.0; norm_subgrid=0.55*3.256/sqrt(n0_sgrid_elements);
        nsubgridvar=(long)P[n].ID + (long)(All.Time/((All.TimeMax-All.TimeBegin)/1000.));
        /* this line just allows 'resetting' the time constants every so often, while generally keeping them steady */
        if(All.ComovingIntegrationOn)
            fac=omega_ri * (evaluate_stellar_age_Gyr(0.001)/(0.001*All.UnitTime_in_Megayears/All.HubbleParam));
        else
            fac=omega_ri * All.Time; /* All.Time is physical time, this is good */
        random_generator_forbh=gsl_rng_alloc(gsl_rng_ranlxd1);
        gsl_rng_set(random_generator_forbh,nsubgridvar);
        if(n0_sgrid_elements >= 1) {
            for(jsub=1;jsub<=n0_sgrid_elements;jsub++) {
                varsg1=gsl_rng_uniform(random_generator_forbh);
                varsg2=gsl_ran_ugaussian(random_generator_forbh);
                time_var_subgridvar=fac*pow(omega_ri*dt,-((float)jsub)/n0_sgrid_elements) + 2.*M_PI*varsg1;
                mdot *= exp( norm_subgrid*cos(time_var_subgridvar)*varsg2 );
                /*
                 printf("SUBGRIDVAR :: mdot %g x %g cosx %g om_ri %g All_t %g dt %g nsubgridvar %ld n0 %g norm %g jsub %d ru %g rg %g \n",
                 mdot,x,cos(x),omega_ri,All.Time,dt,nsubgridvar,n0_sgrid_elements,norm_subgrid,jsub,varsg1,varsg2);fflush(stdout);
                 */
            }}
        gsl_rng_free(random_generator_forbh);
    } // if(mdot > 0)
#endif
    
    
    
    
    
    
#ifdef BH_ENFORCE_EDDINGTON_LIMIT
    /* cap the maximum at the Eddington limit */
    if(mdot > All.BlackHoleEddingtonFactor * meddington)
        mdot = All.BlackHoleEddingtonFactor * meddington;
#endif
    
    
    
    
    
    /* alright, now we can FINALLY set the BH accretion rate */
    BPP(n).BH_Mdot = mdot;
    

    
}



void set_blackhole_new_mass(int i, int n, double dt)
{
    
#ifdef BH_GRAVCAPTURE_SWALLOWS

#endif
    
    BPP(n).BH_Mass += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;        /* this is true BH mass, not P[i].Mass.  True for alpha and non-alpha disk */
    
    
#if defined(BH_ALPHADISK_ACCRETION) // && !defined(BH_GRAVCAPTURE_SWALLOWS)
    BPP(n).BH_Mass_AlphaDisk += (BlackholeTempInfo[i].mdot_alphadisk - BPP(n).BH_Mdot) * dt;
    if(BPP(n).BH_Mass_AlphaDisk<0) BPP(n).BH_Mass_AlphaDisk=0;
    if(P[n].Mass<0) P[n].Mass=0;
#endif
    
#ifdef BH_BUBBLES
    BPP(n).BH_Mass_bubbles += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
#ifdef UNIFIED_FEEDBACK
    if(BPP(n).BH_Mdot < All.RadioThreshold * meddington)
        BPP(n).BH_Mass_radio += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
#endif
#endif

    
}



#if defined(BH_DRAG) || defined(BH_DYNFRICTION)
void set_blackhole_drag(int i, int n, double dt)
{
    int k;
    double meddington, fac;

    meddington = (4 * M_PI * GRAVITY * C * PROTONMASS /
                  (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON) ) *
                  (BPP(n).BH_Mass/All.HubbleParam) * All.UnitTime_in_s;
    
#ifdef BH_DRAG
    /* add a drag force for the black-holes, accounting for the accretion */
    if((dt>0)&&(BPP(n).BH_Mass>0))
    {
        fac = BPP(n).BH_Mdot * dt / BPP(n).BH_Mass;
#ifdef BH_STRONG_DRAG
        /* make the force stronger to keep the BH from wandering */
        fac = meddington * dt / BPP(n).BH_Mass;
#endif
        if(fac>1) fac=1;
        for(k = 0; k < 3; k++)
            P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac * BlackholeTempInfo[i].BH_SurroundingGasVel[k] / dt;
    } // if((dt>0)&&(BPP(n).BH_Mass>0))
#endif
    
    
    
#ifdef BH_DYNFRICTION
    double bh_mass, x;
    
    if(BlackholeTempInfo[i].DF_mmax_particles>0) /* found something in the kernel, we can proceed */
    {
        /* averaged value for colomb logarithm and integral over the distribution function */
        /* fac_friction = log(lambda) * [erf(x) - 2*x*exp(-x^2)/sqrt(pi)]                  */
        /*       lambda = b_max * v^2 / G / (M+m)                                          */
        /*        b_max = Size of system (e.g. Rvir)                                       */
        /*            v = Relative velocity of BH with respect to the environment          */
        /*            M = Mass of BH                                                       */
        /*            m = individual mass elements composing the large system (e.g. m<<M)  */
        /*            x = v/sqrt(2)/sigma                                                  */
        /*        sigma = width of the max. distr. of the host system                      */
        /*                (e.g. sigma = v_disp / 3                                         */
        bh_mass = BPP(n).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass += BPP(n).BH_Mass_AlphaDisk;
#endif
        double bhvel_df=0; for(k=0;k<3;k++) bhvel_df += BlackholeTempInfo[i].DF_mean_vel[k]*BlackholeTempInfo[i].DF_mean_vel[k];
        /* First term is approximation of the error function */
        fac = 8 * (M_PI - 3) / (3 * M_PI * (4. - M_PI));
        x = sqrt(bhvel_df) / (sqrt(2) * BlackholeTempInfo[i].DF_rms_vel);
        double fac_friction =  x / fabs(x) * sqrt(1 - exp(-x * x * (4 / M_PI + fac * x * x) / (1 + fac * x * x))) - 2 * x / sqrt(M_PI) * exp(-x * x);
        /* now the Coulomb logarithm */
        fac = 50. * 3.086e21 / (All.UnitLength_in_cm/All.HubbleParam); /* impact parameter */
        fac_friction *= log(1. + fac * bhvel_df / (All.G * bh_mass));
        /* now we add a correction to only apply this force if M_BH is not >> <m_particles> */
        fac_friction *= 1 / (1 + bh_mass / (5.*BlackholeTempInfo[i].DF_mmax_particles));
        /* now the dimensional part of the force */
        fac = (BlackholeTempInfo[i].Mgas_in_Kernel+BlackholeTempInfo[i].Malt_in_Kernel) /
        ( (4*M_PI/3) * pow(PPP[n].Hsml*All.cf_atime,3) ); /* mean density of all mass inside kernel */
        fac_friction *= 4*M_PI * All.G * All.G * fac * bh_mass / (bhvel_df*sqrt(bhvel_df));
        /* now apply this to the actual acceleration */
        if(fac_friction<0) fac_friction=0; if(isnan(fac_friction)) fac_friction=0;
        for(k = 0; k < 3; k++)
            P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac_friction * BlackholeTempInfo[i].DF_mean_vel[k];
    }
#endif
    
    
    

    
}
#endif


#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
void set_blackhole_long_range_rp(int i, int n)
{
    int k;
    double fac;
    
    /* pre-set quantities needed for long-range radiation pressure terms */
    P[n].BH_disk_hr=1/3; P[n].GradRho[0]=P[n].GradRho[1]=0; P[n].GradRho[2]=1;
    if(BlackholeTempInfo[i].Mgas_in_Kernel > 0)
    {
        /* estimate h/R surrounding the BH from the gas density gradients */
        fac = 0; /* dummy variable */
        for(k=0;k<3;k++)
            fac += BlackholeTempInfo[i].GradRho_in_Kernel[k]*BlackholeTempInfo[i].GradRho_in_Kernel[k];
        P[n].BH_disk_hr = P[n].DensAroundStar / (PPP[n].Hsml * sqrt(fac)) * 1.3;
        /* 1.3 factor from integrating exponential disk
         * with h/R=const over gaussian kernel, for width=1/3 (quintic kernel);
         everything here is in code units, comes out dimensionless */
        
        /* use the gradrho vector as a surrogate to hold the orientation of the angular momentum */
        fac=0;
        for(k=0;k<3;k++)
            fac += BlackholeTempInfo[i].Jgas_in_Kernel[k]*BlackholeTempInfo[i].Jgas_in_Kernel[k];
        fac=sqrt(fac);
        if(fac>0)
            for(k=0;k<3;k++)
                P[n].GradRho[k] = BlackholeTempInfo[i].Jgas_in_Kernel[k]/fac;
        /* now, the P[n].GradRho[k] field for the BH holds the orientation of the UNIT angular momentum vector
         NOTE it is important that HARD-WIRED into the code, this blackhole calculation comes after the density calculation
         but before the forcetree update and walk; otherwise, this won't be used correctly there */
    }
}
#endif // if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)








void blackhole_final_loop(void)
{
    int i, k, n, bin;
    double  dt;
//#ifdef BH_GRAVACCRETION
//    double m_tmp_for_bhar;
//    double r0_for_bhar,j_tmp_for_bhar,fgas_for_bhar,f_disk_for_bhar,mdisk_for_bhar;
//    double f0_for_bhar;
//#endif
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



#endif // BLACK_HOLES
