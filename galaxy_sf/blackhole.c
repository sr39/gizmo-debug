#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES

#ifndef BH_CSND_FRAC_BH_MERGE
/* Relative velocity fraction (in units of soundspeed) for merging blackholes, default=1.0 */
#define BH_CSND_FRAC_BH_MERGE 1.0
#endif

#define BHPOTVALUEINIT 1.0e30

static int N_gas_swallowed, N_BH_swallowed;

/* quantities that pass IN to the 'blackhole_evaluate' routines */
static struct blackholedata_in
{
    MyDouble Pos[3];
    MyFloat Density;
    MyFloat Mdot;
    MyFloat Dt;
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
    MyFloat Vel[3];
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat Jgas_in_Kernel[3];
    MyFloat BH_disk_hr;
    MyFloat BH_angle_weighted_kernel_sum;
#endif
}
*BlackholeDataIn, *BlackholeDataGet;

/* quantities that pass OUT of the 'blackhole_evaluate' routines */
static struct blackholedata_out
{
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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat BH_angle_weighted_kernel_sum;
#endif
}
*BlackholeDataResult, *BlackholeDataOut;



/* this is a temporary structure for quantities used ONLY in the step below,
 which we calculate in a pre-pass before the 'final' loop on BH particles */
static struct blackholedata_topass
{
    MyFloat BH_InternalEnergy;
    MyLongDouble accreted_Mass;
    MyLongDouble accreted_BH_Mass;
    MyLongDouble accreted_momentum[3];
    MyLongDouble Mgas_in_Kernel;
    MyLongDouble Malt_in_Kernel;
    MyLongDouble Jalt_in_Kernel[3];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyLongDouble GradRho_in_Kernel[3];
    MyLongDouble Jgas_in_Kernel[3];
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_DYNFRICTION
    MyFloat DF_mean_vel[3];
    MyFloat DF_rms_vel;
    MyFloat DF_mmax_particles;
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
    MyFloat BH_SurroundingGasVel[3];
#endif
}
*BlackholeDataPasser,*BlackholeDataPasserOut,*BlackholeDataPasserResult;


void out2particle_blackhole(struct blackholedata_topass *out, int target, int mode);


/* simple routine to add quantities from blackhole_evaluate_PREPASS to BlackholeDataPasser */
void out2particle_blackhole(struct blackholedata_topass *out, int target, int mode)
{
    int k;
    ASSIGN_ADD(BlackholeDataPasser[target].BH_InternalEnergy,out->BH_InternalEnergy,mode);
    ASSIGN_ADD(BlackholeDataPasser[target].Mgas_in_Kernel,out->Mgas_in_Kernel,mode);
    ASSIGN_ADD(BlackholeDataPasser[target].Malt_in_Kernel,out->Malt_in_Kernel,mode);
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeDataPasser[target].Jalt_in_Kernel[k],out->Jalt_in_Kernel[k],mode);
#ifdef BH_DYNFRICTION
    ASSIGN_ADD(BlackholeDataPasser[target].DF_rms_vel,out->DF_rms_vel,mode);
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeDataPasser[target].DF_mean_vel[k],out->DF_mean_vel[k],mode);
    if(mode==0)
        BlackholeDataPasser[target].DF_mmax_particles = out->DF_mmax_particles;
    else
        if(out->DF_mmax_particles > BlackholeDataPasser[target].DF_mmax_particles)
            BlackholeDataPasser[target].DF_mmax_particles = out->DF_mmax_particles;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeDataPasser[target].Jgas_in_Kernel[k],out->Jgas_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeDataPasser[target].GradRho_in_Kernel[k],out->GradRho_in_Kernel[k],mode);
    }
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
    for(k=0;k<3;k++)
        ASSIGN_ADD(BlackholeDataPasser[target].BH_SurroundingGasVel[k],out->BH_SurroundingGasVel[k],mode);
#endif
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




/* Alright, now this is the 'master' routine for the BH physics modules */
void blackhole_accretion(void)
{
    int i, j, k, n, bin;
    int ndone_flag, ndone;
    int ngrp, recvTask, place, nexport, nimport, dummy;
    int Ntot_gas_swallowed, Ntot_BH_swallowed;
    double fac, mdot, meddington, dt, mdot_in_msun_per_year;
    double mass_real, total_mass_real, medd, total_mdoteddington;
    double mass_holes, total_mass_holes, total_mdot;
    MPI_Status status;
    double x,bh_mass;
    x=bh_mass=0;
    
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
#ifdef BH_ALPHADISK_ACCRETION
    double mdot_alphadisk;
#endif
    
    if(ThisTask == 0)
    {
        printf("Beginning black-hole accretion\n");
        fflush(stdout);
    }
    
    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
     First, we enter a loop to calculate properties of the gas surrounding
     the BH -- we now ALWAYS do this, independent of the modules set
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/
    
    /* allocate buffers to arrange communication */
    BlackholeDataPasser = (struct blackholedata_topass *) mymalloc("BlackholeDataPasser",NumPart * sizeof(struct blackholedata_topass));
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),
                                                                     sizeof(struct blackholedata_out))));
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
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(blackhole_evaluate_PREPASS(i, 0, &nexport, Send_count) < 0)
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
        BlackholeDataPasserResult = (struct blackholedata_topass *) mymalloc("BlackholeDataPasserResult", nexport * sizeof(struct blackholedata_topass));
        BlackholeDataPasserOut = (struct blackholedata_topass *) mymalloc("BlackholeDataPasserOut", nimport * sizeof(struct blackholedata_topass));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_evaluate_PREPASS(j, 1, &dummy, &dummy);
        
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
                                 Recv_count[recvTask] * sizeof(struct blackholedata_topass),
                                 MPI_BYTE, recvTask, TAG_DENS_B,
                                 &BlackholeDataPasserOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_topass),
                                 MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
                }
            }
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_blackhole(&BlackholeDataPasserOut[j], place, 1);
        } // for(j = 0; j < nexport; j++)
        myfree(BlackholeDataPasserOut);
        myfree(BlackholeDataPasserResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    
    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
     now that the PRE-PASS loop is done, do a first set of
     operations on the relevant quantities, calculating mdot,
     dynamical friction forces, and other 'BH-centric' operations
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/
    
    for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    {
        if(P[n].Type == 5)
        {
            /* define the timestep */
#ifndef WAKEUP
            dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            dt = P[n].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            /* define the Eddington limit for use later (Note: we take here a radiative efficiency as specified in the parameter file) */
            meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON)) * (BPP(n).BH_Mass/All.HubbleParam) * All.UnitTime_in_s;
            
            /* always initialize/default to zero accretion rate */
            mdot=0;
            BPP(n).BH_Mdot=0;
            
            /* for new quantities, divide out weights and convert to physical units */
            if(BlackholeDataPasser[n].Mgas_in_Kernel > 0)
            {
                BlackholeDataPasser[n].BH_InternalEnergy /= BlackholeDataPasser[n].Mgas_in_Kernel;
#ifdef BH_DYNFRICTION
                BlackholeDataPasser[n].DF_rms_vel /= BlackholeDataPasser[n].Mgas_in_Kernel;
                BlackholeDataPasser[n].DF_rms_vel = sqrt(BlackholeDataPasser[n].DF_rms_vel) / All.cf_atime;
                for(k=0;k<3;k++)
                    BlackholeDataPasser[n].DF_mean_vel[k] /= BlackholeDataPasser[n].Mgas_in_Kernel * All.cf_atime;
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
                for(k=0;k<3;k++)
                    BlackholeDataPasser[n].BH_SurroundingGasVel[k] /= BlackholeDataPasser[n].Mgas_in_Kernel * All.cf_atime;
#endif
            }
            else
            {
                BlackholeDataPasser[n].BH_InternalEnergy = 0;
            }
            
            
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            /* pre-set quantities needed for long-range radiation pressure terms */
            P[n].BH_disk_hr=1/3; P[n].GradRho[0]=P[n].GradRho[1]=0; P[n].GradRho[2]=1;
            if(BlackholeDataPasser[n].Mgas_in_Kernel > 0)
            {
                /* estimate h/R surrounding the BH from the gas density gradients */
                fac = 0; /* dummy variable */
                for(k=0;k<3;k++)
                    fac += BlackholeDataPasser[n].GradRho_in_Kernel[k]*BlackholeDataPasser[n].GradRho_in_Kernel[k];
                P[n].BH_disk_hr = P[n].DensAroundStar / (PPP[n].Hsml * sqrt(fac)) * 1.3;
                /* 1.3 factor from integrating exponential disk with h/R=const over gaussian kernel, for width=1/3 (quintic kernel);
                 everything here is in code units, comes out dimensionless */
                
                /* use the gradrho vector as a surrogate to hold the orientation of the angular momentum */
                fac=0;
                for(k=0;k<3;k++)
                    fac += BlackholeDataPasser[n].Jgas_in_Kernel[k]*BlackholeDataPasser[n].Jgas_in_Kernel[k];
                fac=sqrt(fac);
                if(fac>0)
                    for(k=0;k<3;k++)
                        P[n].GradRho[k] = BlackholeDataPasser[n].Jgas_in_Kernel[k]/fac;
                /* now, the P[n].GradRho[k] field for the BH holds the orientation of the UNIT angular momentum vector
                 NOTE it is important that HARD-WIRED into the code, this blackhole calculation comes after the density calculation
                 but before the forcetree update and walk; otherwise, this won't be used correctly there */
            }
#endif // if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)

            
            
#ifdef BH_GRAVACCRETION
            /* calculate mdot: gravitational instability accretion rate from Hopkins & Quataert 2011 */
            if(BlackholeDataPasser[n].Mgas_in_Kernel > 0)
            {
                m_tmp_for_bhar = BlackholeDataPasser[n].Mgas_in_Kernel + BlackholeDataPasser[n].Malt_in_Kernel;
                r0_for_bhar = PPP[n].Hsml * All.cf_atime; /* convert to physical units */
                j_tmp_for_bhar=0; for(k=0;k<3;k++) j_tmp_for_bhar += BlackholeDataPasser[n].Jalt_in_Kernel[k]*BlackholeDataPasser[n].Jalt_in_Kernel[k];
                j_tmp_for_bhar=sqrt(j_tmp_for_bhar);
                /* jx,y,z, is independent of 'a_scale' b/c ~ m*r*v, vphys=v/a, rphys=r*a */
                
                bh_mass = BPP(n).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
                bh_mass += BPP(n).BH_Mass_AlphaDisk;
#endif
                fgas_for_bhar = BlackholeDataPasser[n].Mgas_in_Kernel / m_tmp_for_bhar;
                fac = m_tmp_for_bhar * r0_for_bhar * sqrt(All.G*(m_tmp_for_bhar+bh_mass)/r0_for_bhar); /* All.G is G in code (physical) units */
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
            for(k=0;k<3;k++) bhvel += BlackholeDataPasser[n].BH_SurroundingGasVel[k]*BlackholeDataPasser[n].BH_SurroundingGasVel[k];
#endif
            rho = BPP(n).DensAroundStar * All.cf_a3inv; /* we want all quantities in physical units */
            soundspeed = GAMMA*GAMMA_MINUS1 * BlackholeDataPasser[n].BH_InternalEnergy; // this is in physical units now
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
            mdot = 0; /* force mdot=0 despite any earlier settings here */
#endif
            
            
#ifdef BH_ALPHADISK_ACCRETION
            /* use the mass in the accretion disk from the previous timestep to determine the BH accretion rate */
            mdot_alphadisk = mdot;
            mdot = All.BlackHoleAccretionFactor *
                (2.45 * (SOLAR_MASS/All.UnitMass_in_g)/(SEC_PER_YEAR/All.UnitTime_in_s)) * // normalization
                pow( 0.1 , 8./7.) * // viscous disk 'alpha'
                pow( BPP(n).BH_Mass*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , -5./14. ) * // mbh dependence
                pow( BPP(n).BH_Mass_AlphaDisk*All.UnitMass_in_g / (All.HubbleParam * 1.0e8*SOLAR_MASS) , 10./7. ) * // m_disk dependence
                pow( DMIN(0.2,DMIN(PPP[n].Hsml,All.ForceSoftening[5])*All.cf_atime*All.UnitLength_in_cm/(All.HubbleParam * 3.086e18)) , -25./14. ); // r_disk dependence
            if(mdot<=0) mdot=0;
            if(dt>0)
            {
#ifdef BH_BAL_WINDS
                /* this is just here to prevent it from accidentally going to negative mass */
                if(mdot > BPP(n).BH_Mass_AlphaDisk/(All.BAL_f_accretion*dt)) mdot = BPP(n).BH_Mass_AlphaDisk/(All.BAL_f_accretion*dt);
#else
                if(mdot > BPP(n).BH_Mass_AlphaDisk/dt) mdot = BPP(n).BH_Mass_AlphaDisk/dt;
#endif
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
            
            
            /* dump the results so far to the 'blackhole_details' files */
            fac=0; medd=0;
#ifdef BH_ALPHADISK_ACCRETION
            fac=BPP(n).BH_Mass_AlphaDisk;
            medd=mdot_alphadisk;
#endif
            fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g %g %g %g   %2.7f %2.7f %2.7f\n",
                    P[n].ID, All.Time, BPP(n).BH_Mass, fac, P[n].Mass, mdot, medd,
                    BPP(n).DensAroundStar*All.cf_a3inv, BlackholeDataPasser[n].BH_InternalEnergy,
                    P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]);
            
            
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
                    P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac * BlackholeDataPasser[n].BH_SurroundingGasVel[k] / dt;
            } // if((dt>0)&&(BPP(n).BH_Mass>0))
#endif
            
            
            
#ifdef BH_DYNFRICTION
            if(BlackholeDataPasser[n].DF_mmax_particles>0) /* found something in the kernel, we can proceed */
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
                double bhvel_df=0; for(k=0;k<3;k++) bhvel_df += BlackholeDataPasser[n].DF_mean_vel[k]*BlackholeDataPasser[n].DF_mean_vel[k];
                /* First term is approximation of the error function */
                fac = 8 * (M_PI - 3) / (3 * M_PI * (4. - M_PI));
                x = sqrt(bhvel_df) / (sqrt(2) * BlackholeDataPasser[n].DF_rms_vel);
                double fac_friction =  x / fabs(x) * sqrt(1 - exp(-x * x * (4 / M_PI + fac * x * x) / (1 + fac * x * x))) - 2 * x / sqrt(M_PI) * exp(-x * x);
                /* now the Coulomb logarithm */
                fac = 50. * 3.086e21 / (All.UnitLength_in_cm/All.HubbleParam); /* impact parameter */
                fac_friction *= log(1. + fac * bhvel_df / (All.G * bh_mass));
                /* now we add a correction to only apply this force if M_BH is not >> <m_particles> */
                fac_friction *= 1 / (1 + bh_mass / (5.*BlackholeDataPasser[n].DF_mmax_particles));
                /* now the dimensional part of the force */
                fac = (BlackholeDataPasser[n].Mgas_in_Kernel+BlackholeDataPasser[n].Malt_in_Kernel) /
                        ( (4*M_PI/3) * pow(PPP[n].Hsml*All.cf_atime,3) ); /* mean density of all mass inside kernel */
                fac_friction *= 4*M_PI * All.G * All.G * fac * bh_mass / (bhvel_df*sqrt(bhvel_df));
                /* now apply this to the actual acceleration */
                if(fac_friction<0) fac_friction=0; if(isnan(fac_friction)) fac_friction=0;
                for(k = 0; k < 3; k++)
                    P[n].GravAccel[k] += All.cf_atime*All.cf_atime * fac_friction * BlackholeDataPasser[n].DF_mean_vel[k];
            }
#endif
            
            
            /* INCREMENT BH mass if mdot > 0 */
            BPP(n).BH_Mass += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot * dt;
#ifdef BH_ALPHADISK_ACCRETION
#ifdef BH_BAL_WINDS
            BPP(n).BH_Mass_AlphaDisk += (mdot_alphadisk-BPP(n).BH_Mdot/All.BAL_f_accretion) * dt;
            /* correct real particle mass for mass flux in accretion-disk wind */
            P[n].Mass -= BPP(n).BH_Mdot*dt / All.BAL_f_accretion;
#else
            BPP(n).BH_Mass_AlphaDisk += (mdot_alphadisk-BPP(n).BH_Mdot) * dt;
#endif
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
        } // if(P[n].Type == 5)
    }// for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n]) if(P[n].Type == 5)
    
    
    
    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
     Now let's do ANOTHER pass over the particles, and
     invoke the functions that calculate when to stochastically swallow gas
     and deal with black hole mergers (blackhole_evaluate), as well 
     as using the above info to determine the weight functions for feedback
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/
    
    if(ThisTask == 0)
    {
        printf("Start swallowing of gas particles and black holes\n");
        fflush(stdout);
    }
    N_gas_swallowed = N_BH_swallowed = 0;
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
                if(blackhole_evaluate(i, 0, &nexport, Send_count) < 0)
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
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
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
            blackhole_evaluate(j, 1, &dummy, &dummy);
        
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
            BlackholeDataPasser[place].BH_angle_weighted_kernel_sum += BlackholeDataOut[j].BH_angle_weighted_kernel_sum;
#endif
        }

        myfree(BlackholeDataOut);
        myfree(BlackholeDataResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    
    
    
    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
     Ok, now we do a THIRD pass over the particles, and
     this is where we can do the actual 'swallowing' operations
     (blackhole_evaluate_swallow), and 'kicking' operations
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/
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
                if(P[i].SwallowID == 0)
                    if(blackhole_evaluate_swallow(i, 0, &nexport, Send_count) < 0)
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
        
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            for(k = 0; k < 3; k++)
            {
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].GradRho[k];
#endif
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
            BlackholeDataIn[j].BH_angle_weighted_kernel_sum = BlackholeDataPasser[place].BH_angle_weighted_kernel_sum;
#endif
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            BlackholeDataIn[j].ID = P[place].ID;
            BlackholeDataIn[j].Mass = P[place].Mass;
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
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult", nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_evaluate_swallow(j, 1, &dummy, &dummy);
        
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
        }
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            BlackholeDataPasser[place].accreted_Mass += BlackholeDataOut[j].Mass;
            BlackholeDataPasser[place].accreted_BH_Mass += BlackholeDataOut[j].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BPP(place).BH_Mass_AlphaDisk += BlackholeDataOut[j].BH_Mass_AlphaDisk;
#endif
#ifdef BH_BUBBLES
            BPP(place).b7.dBH_accreted_BHMass_bubbles += BlackholeDataOut[j].BH_Mass_bubbles;
#ifdef UNIFIED_FEEDBACK
            BPP(place).b8.dBH_accreted_BHMass_radio += BlackholeDataOut[j].BH_Mass_radio;
#endif
#endif
            for(k = 0; k < 3; k++)
                BlackholeDataPasser[place].accreted_momentum[k] += BlackholeDataOut[j].accreted_momentum[k];
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

    /* active particle loop is done, free the pre-pass quantities we used to calculate all this! */
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("Accretion done: %d gas particles swallowed, %d BH particles swallowed\n",
               Ntot_gas_swallowed, Ntot_BH_swallowed);
        fflush(stdout);
    }
    
    
    /*----------------------------------------------------------------------
     ------------------------------------------------------------------------
       Now do final operations on the results from the last pass
     ------------------------------------------------------------------------
     ----------------------------------------------------------------------*/
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
    
    for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    {
        if(P[n].Type == 5)
        {
            if(((BlackholeDataPasser[n].accreted_Mass>0)||(BlackholeDataPasser[n].accreted_BH_Mass>0)) && P[n].Mass > 0)
                {
                    for(k = 0; k < 3; k++)
                    {
                        /* correct momentum for accretion (so we conserve total momentum) */
#ifndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        P[n].Vel[k] = (P[n].Vel[k]*P[n].Mass + BlackholeDataPasser[n].accreted_momentum[k]) / (BlackholeDataPasser[n].accreted_Mass + P[n].Mass);
#else
                        P[n].Vel[k] = (P[n].Vel[k]*P[n].Mass + BlackholeDataPasser[n].accreted_momentum[k]) / (BlackholeDataPasser[n].accreted_BH_Mass + P[n].Mass);
#endif
                    } //for(k = 0; k < 3; k++)
                    P[n].Mass += BlackholeDataPasser[n].accreted_Mass;
                    BPP(n).BH_Mass += BlackholeDataPasser[n].accreted_BH_Mass;
#ifdef BH_BUBBLES
                    BPP(n).BH_Mass_bubbles += BPP(n).b7.BH_accreted_BHMass_bubbles;
#ifdef UNIFIED_FEEDBACK
                    BPP(n).BH_Mass_radio += BPP(n).b8.BH_accreted_BHMass_radio;
#endif
#endif
                } // if(((BlackholeDataPasser[n].accreted_Mass>0)||(BlackholeDataPasser[n].accreted_BH_Mass>0)) && P[n].Mass > 0)
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
            /* now that its been added to the BH, use the actual swallowed mass to define the BHAR */
#ifdef BH_ALPHADISK_ACCRETION
            if(dt>0)
                BPP(n).BH_Mdot += BlackholeDataPasser[n].accreted_BH_Mass / dt;
#else
            if(dt>0)
                BPP(n).BH_Mdot = BlackholeDataPasser[n].accreted_BH_Mass / dt;
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
        } // if(P[n].Type == 5) //
    } // for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n]) //
    
#ifdef BH_BUBBLES
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    MPI_Allreduce(&num_activebh, &total_num_activebh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("The total number of active BHs is: %d\n", total_num_activebh);
        fflush(stdout);
    }
    if(total_num_activebh > 0)
    {
        bh_dmass = mymalloc("bh_dmass", num_activebh * sizeof(double));
        tot_bh_dmass = mymalloc("tot_bh_dmass", total_num_activebh * sizeof(double));
        bh_posx = mymalloc("bh_posx", num_activebh * sizeof(float));
        bh_posy = mymalloc("bh_posy", num_activebh * sizeof(float));
        bh_posz = mymalloc("bh_posz", num_activebh * sizeof(float));
        tot_bh_posx = mymalloc("tot_bh_posx", total_num_activebh * sizeof(float));
        tot_bh_posy = mymalloc("tot_bh_posy", total_num_activebh * sizeof(float));
        tot_bh_posz = mymalloc("tot_bh_posz", total_num_activebh * sizeof(float));
        bh_id = mymalloc("bh_id", num_activebh * sizeof(MyIDType));
        tot_bh_id = mymalloc("tot_bh_id", total_num_activebh * sizeof(MyIDType));
        
        for(n = 0; n < num_activebh; n++)
        {
            bh_dmass[n] = 0.0;
            bh_posx[n] = 0.0;
            bh_posy[n] = 0.0;
            bh_posz[n] = 0.0;
            bh_id[n] = 0;
        }
        for(n = 0; n < total_num_activebh; n++)
        {
            tot_bh_dmass[n] = 0.0;
            tot_bh_posx[n] = 0.0;
            tot_bh_posy[n] = 0.0;
            tot_bh_posz[n] = 0.0;
            tot_bh_id[n] = 0;
        }
        for(n = FirstActiveParticle, l = 0; n >= 0; n = NextActiveParticle[n])
            if(P[n].Type == 5)
            {
                if(BPP(n).BH_Mass_bubbles > 0
                   && BPP(n).BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * BPP(n).BH_Mass_ini)
                {
#ifndef UNIFIED_FEEDBACK
                    bh_dmass[l] = BPP(n).BH_Mass_bubbles - BPP(n).BH_Mass_ini;
#else
                    bh_dmass[l] = BPP(n).BH_Mass_radio - BPP(n).BH_Mass_ini;
                    BPP(n).BH_Mass_radio = BPP(n).BH_Mass;
#endif
                    BPP(n).BH_Mass_ini = BPP(n).BH_Mass;
                    BPP(n).BH_Mass_bubbles = BPP(n).BH_Mass;
                    bh_posx[l] = P[n].Pos[0];
                    bh_posy[l] = P[n].Pos[1];
                    bh_posz[l] = P[n].Pos[2];
                    bh_id[l] = P[n].ID;
                    l++;
                }
            }
        common_num_activebh = mymalloc("common_num_activebh", NTask * sizeof(int));
        disp = mymalloc("disp", NTask * sizeof(int));
        MPI_Allgather(&num_activebh, 1, MPI_INT, common_num_activebh, 1, MPI_INT, MPI_COMM_WORLD);
        for(k = 1, disp[0] = 0; k < NTask; k++)
            disp[k] = disp[k - 1] + common_num_activebh[k - 1];
        
        MPI_Allgatherv(bh_dmass, num_activebh, MPI_DOUBLE, tot_bh_dmass, common_num_activebh, disp, MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allgatherv(bh_posx, num_activebh, MPI_FLOAT, tot_bh_posx, common_num_activebh, disp, MPI_FLOAT,MPI_COMM_WORLD);
        MPI_Allgatherv(bh_posy, num_activebh, MPI_FLOAT, tot_bh_posy, common_num_activebh, disp, MPI_FLOAT,MPI_COMM_WORLD);
        MPI_Allgatherv(bh_posz, num_activebh, MPI_FLOAT, tot_bh_posz, common_num_activebh, disp, MPI_FLOAT,MPI_COMM_WORLD);
#ifndef LONGIDS
        MPI_Allgatherv(bh_id, num_activebh, MPI_UNSIGNED, tot_bh_id, common_num_activebh, disp, MPI_UNSIGNED,MPI_COMM_WORLD);
#else
        MPI_Allgatherv(bh_id, num_activebh, MPI_UNSIGNED_LONG_LONG, tot_bh_id, common_num_activebh, disp,MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
#endif
        for(l = 0; l < total_num_activebh; l++)
        {
            bh_center[0] = tot_bh_posx[l];
            bh_center[1] = tot_bh_posy[l];
            bh_center[2] = tot_bh_posz[l];
            if(tot_bh_dmass[l] > 0)
                bh_bubble(tot_bh_dmass[l], bh_center, tot_bh_id[l]);
        }
        myfree(disp);
        myfree(common_num_activebh);
        myfree(tot_bh_id);
        myfree(bh_id);
        myfree(tot_bh_posz);
        myfree(tot_bh_posy);
        myfree(tot_bh_posx);
        myfree(bh_posz);
        myfree(bh_posy);
        myfree(bh_posx);
        myfree(tot_bh_dmass);
        myfree(bh_dmass);
    }
    myfree(Ngblist);
#endif //  BH_BUBBLES
    
    
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
    
    myfree(BlackholeDataPasser);
    fflush(FdBlackHolesDetails);
}





/* do loop over neighbors to get quantities for accretion */
int blackhole_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, listindex = 0;
    MyIDType id;
    MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass;
    double h_i2, r2, r, u, hinv, hinv3, wk, dwk, vrel, csnd;
    double dpos[3];
    
#if defined(UNIFIED_FEEDBACK) || defined(BH_ENFORCE_EDDINGTON_LIMIT)
    double meddington;
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
    double medd_max_accretable,m_to_swallow_thispart;
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
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
#endif
    }
    if((mass<0)||(h_i<=0)) return -1;
    
    
    /* initialize variables before SPH loop is started */
    h_i2 = h_i * h_i;
    hinv = 1 / h_i;
    hinv3 = hinv * hinv * hinv;
#ifdef BH_ENFORCE_EDDINGTON_LIMIT
    meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (All.BlackHoleRadiativeEfficiency * C * C * THOMPSON)) * (bh_mass/All.HubbleParam) * All.UnitTime_in_s;
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
	medd_max_accretable = All.BlackHoleEddingtonFactor * meddington * dt;
#endif
#endif
#if defined(BH_SWALLOWGAS) || defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif
    
    /* Now start the actual SPH computation for this particle */
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
                        
#ifdef BH_REPOSITION_ON_POTMIN
                        /* check if we've found a new potential minimum which is not moving too fast to 'jump' to */
                        if(P[j].Potential < minpot)
                        {
                            vrel = 0;
                            for(k=0;k<3;k++) vrel += (P[j].Vel[k] - velocity[k])*(P[j].Vel[k] - velocity[k]);
                            vrel = sqrt(vrel) / All.cf_atime;
                            /* escape velocity from BH at distance r, plus 10 km/s ionized gas sound speed 'floor' */
                            csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime) + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                            if(vrel <= csnd)
                            {
                                minpot = P[j].Potential;
                                for(k = 0; k < 3; k++) minpotpos[k] = P[j].Pos[k];
                            }
                        }
#endif
                        
                        if(P[j].Type == 5)	/* we (may) have a black hole merger */
                        {
                            if(id != P[j].ID) /* check its not the same bh, duh! */
                            {
                                /* compute relative velocity of BHs */
                                for(k = 0, vrel = 0; k < 3; k++)
                                    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);
                                vrel = sqrt(vrel) / All.cf_atime;
                                
                                /* with feedback on, sound speeds get -very- low, would be silly to use; 
                                    instead use the escape velocity and follow to pairing */
                                csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime));
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
                        
                        
#if defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS)
                        /* this will determine accretion from all particle types, based only on the
                         ability of the particle to be gravitationally captured by the BH */
#ifdef BH_GRAVCAPTURE_SWALLOWS
                        if(P[j].Type != 5)
#else
                        if((P[j].Type != 0)&&(P[j].Type != 5))
#endif
                        {
                            /* compute relative velocity */
                            for(k = 0, vrel = 0; k < 3; k++)
                                vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);
                            vrel = sqrt(vrel) / All.cf_atime;
                            r = sqrt(r2);
                            csnd = sqrt(2.0*All.G*(mass+P[j].Mass)/(r*All.cf_atime)); /* escape velocity */
                            
                            if(vrel < csnd)
                            { /* bound */
                                if( All.ForceSoftening[5]*(1.0-vrel*vrel/(csnd*csnd))/r > 1.0 )
                                { /* apocenter within 2.8*epsilon (softening length) */
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                    /* only count gas and stars towards the Eddington limit */
                                    p=1; m_to_swallow_thispart=P[j].Mass;
                                    if((P[j].Type != 1)||(All.ComovingIntegrationOn && (P[j].Type==0||P[j].Type==4)))
                                    {
                                        if(medd_max_accretable-mass_markedswallow <= 0)
                                        {
                                            p=0;
                                        } else {
#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                                            if(All.BAL_f_accretion>0) m_to_swallow_thispart *= All.BAL_f_accretion;
#endif
                                            if(m_to_swallow_thispart >  medd_max_accretable-mass_markedswallow)
                                                p = (medd_max_accretable-mass_markedswallow)/m_to_swallow_thispart;
                                        }}
                                    w = get_random_number(P[j].ID);
                                    if(w < p) {
                                        printf("MARKING_BH_FOOD: j %d w %g p_acc %g facc %g TO_BE_SWALLOWED \n",j,w,p,m_to_swallow_thispart/P[j].Mass);
                                        if(P[j].SwallowID < id) {
                                            P[j].SwallowID = id;
                                            if((P[j].Type != 1)||(All.ComovingIntegrationOn && (P[j].Type==0||P[j].Type==4)))
                                                mass_markedswallow += m_to_swallow_thispart;
                                        } /* P[j].SwallowID < id */
                                    } /* w < p */
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
                            r = sqrt(r2);
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
                            /* compute random number, uniform in [0,1] */
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
        BlackholeDataPasser[target].BH_angle_weighted_kernel_sum = BH_angle_weighted_kernel_sum;
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




int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, bin, listindex = 0;
    MyIDType id;
    MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3];
    MyFloat *pos, *velocity, h_i, bh_mass, hinv, hinv3;
    MyFloat mdot,dt;
    MyFloat dir[3],norm,mom=0;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat f_accreted,v_kick;
    f_accreted=0;
    double BH_angle_weighted_kernel_sum, mom_wt, e_wind, m_wind;
    MyFloat theta,*Jgas_in_Kernel,BH_disk_hr,kernel_zero,dwk;
    kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1);
#endif
#ifdef BH_BUBBLES
    MyLongDouble accreted_BH_mass_bubbles = 0;
    MyLongDouble accreted_BH_mass_radio = 0;
#endif
#ifdef GALSF
    double accreted_age = 1;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk, accreted_BH_mass_alphadisk;
#endif
    
    if(mode == 0)
    {
        pos = P[target].Pos;
        velocity = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        bh_mass = BPP(target).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = P[target].GradRho;
        BH_disk_hr = P[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeDataPasser[target].BH_angle_weighted_kernel_sum;
#endif
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        velocity = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        bh_mass = BlackholeDataGet[target].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;
#endif
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
        BH_angle_weighted_kernel_sum = BlackholeDataGet[target].BH_angle_weighted_kernel_sum;
#endif
    }
    
    accreted_mass = 0;
    accreted_BH_mass = 0;
#ifdef BH_ALPHADISK_ACCRETION
    accreted_BH_mass_alphadisk = 0;
#endif
    accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
#ifdef BH_COUNTPROGS
    int accreted_BH_progs = 0;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    mom = All.BlackHoleFeedbackFactor *
        All.BlackHoleRadiativeEfficiency * mdot * dt * (C / All.UnitVelocity_in_cm_per_s);
    mom_wt = 0;
#endif
    
    hinv=h_i;
    hinv3=hinv*hinv*hinv;
    
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
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                
                
                /* check for BH-BH mergers */
                if(P[j].SwallowID == id)
                {
                    if(P[j].Type == 5)
                    {
                        fprintf(FdBlackHolesDetails,
                                "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
                                ThisTask, All.Time, id, P[j].ID, bh_mass, BPP(j).BH_Mass);
                        
                        accreted_mass += FLT(P[j].Mass);
                        accreted_BH_mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
#ifdef BH_BUBBLES
                        accreted_BH_mass_bubbles += FLT(BPP(j).BH_Mass_bubbles - BPP(j).BH_Mass_ini);
#ifdef UNIFIED_FEEDBACK
                        accreted_BH_mass_radio += FLT(BPP(j).BH_Mass_radio - BPP(j).BH_Mass_ini);
#endif
#endif
#ifndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
#else
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(BPP(j).BH_Mass * P[j].Vel[k]);
#endif
#ifdef BH_COUNTPROGS
                        accreted_BH_progs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin;
                        TimeBin_BH_mass[bin] -= BPP(j).BH_Mass;
                        TimeBin_BH_dynamicalmass[bin] -= P[j].Mass;
                        TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0)
                            TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;
                        P[j].Mass = 0;
                        BPP(j).BH_Mass = 0;
                        BPP(j).BH_Mdot = 0;
#ifdef BH_BUBBLES
                        BPP(j).BH_Mass_bubbles = 0;
                        BPP(j).BH_Mass_ini = 0;
#ifdef UNIFIED_FEEDBACK
                        BPP(j).BH_Mass_radio = 0;
#endif
#endif
#ifdef GALSF
                        accreted_age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5)
                } // if(P[j].SwallowID == id)

                
                /* now, do the accretion/swallowing of other particle types */
#if !defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                if(P[j].Type == 0)
#endif
                {
                    if((P[j].SwallowID == id)&&(P[j].Type != 5))
                    {
                        
#if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                        /* here we do the BAL wind model if and only if there is a swallowing event (since mdot is not continuous) */
                        printf("BAL kick: j %d Type(j) %d f_acc %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",j,P[j].Type,
                               All.BAL_f_accretion,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2],
                               All.BAL_v_outflow);fflush(stdout);
                        
                        dir[0]=dir[1]=dir[2]=0;
                        if(P[j].Type==0) f_accreted=All.BAL_f_accretion; else f_accreted=1;
                        if(P[j].Type==0) v_kick=All.BAL_v_outflow; else v_kick=0.1*All.BAL_v_outflow;
                        accreted_mass += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass);
#ifdef BH_GRAVCAPTURE_SWALLOWS
                        accreted_BH_mass += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass);
#endif // #ifdef BH_GRAVCAPTURE_SWALLOWS
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(f_accreted*(1. - All.BlackHoleRadiativeEfficiency)*P[j].Mass * P[j].Vel[k]);
                        P[j].Mass *= (1-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue *= (1-f_accreted);
#endif
                        for(k = 0; k < 3; k++) dir[k]=P[j].Pos[k]-pos[k];
                        for(k = 0, norm = 0; k < 3; k++) norm += dir[k]*dir[k];
                        if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm);}
                        for(k = 0; k < 3; k++)
                        {
                            P[j].Vel[k] += v_kick*All.cf_atime*dir[k]/norm;
                            if(P[j].Type==0) SphP[j].VelPred[k] += v_kick*All.cf_atime*dir[k]/norm;
                        }
#else // #if defined(BH_BAL_WINDS) && defined(BH_GRAVCAPTURE_SWALLOWS) && !defined(BH_GRAVCAPTURE_NOGAS)
                        /* in this case, no kick, so just zero out the mass and 'get rid of' the
                         particle (preferably by putting it somewhere irrelevant) */
                        accreted_mass += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass);
#ifdef BH_GRAVCAPTURE_SWALLOWS
#ifdef BH_ALPHADISK_ACCRETION
                        /* mass goes into the alpha disk, before going into the BH */
                        accreted_BH_mass_alphadisk += FLT(P[j].Mass);
#else
                        /* mass goes directly to the BH, not just the parent particle */
                        accreted_BH_mass += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass);
#endif
#endif // #ifdef BH_GRAVCAPTURE_SWALLOWS
#ifndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT((1. - All.BlackHoleRadiativeEfficiency) * P[j].Mass * P[j].Vel[k]);
#endif // ndef BH_IGNORE_ACCRETED_GAS_MOMENTUM
                        P[j].Mass = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue = 0;
#endif
#endif // ifdef BH_BAL_WINDS else
                        
                        N_gas_swallowed++;
                        
                    } // if(P[j].SwallowID == id)
                } // if(P[j].Type == 0)
                
                
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0)
                {
                    if(P[j].Type==0)
                    {
                        if((P[j].Mass>0)&&(P[j].SwallowID==0)) // not swallowed!
                        {
                            for(norm=0,k=0;k<3;k++)
                            {
                                dir[k] = (pos[k]-P[j].Pos[k]);
                                norm += dir[k]*dir[k];
                            }
                            if(norm>0)
                            {
                                norm=sqrt(norm); for(k=0;k<3;k++) dir[k]/=norm;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                                /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                for(norm=0,k=0;k<3;k++) norm += dir[k]*Jgas_in_Kernel[k];
                                theta = acos(fabs(norm));
#ifdef BH_PHOTONMOMENTUM
                                /* now we get the weight function based on what we calculated earlier */
                                mom_wt = All.BH_FluxMomentumFactor * bh_angleweight_localcoupling(j,BH_disk_hr,theta) / BH_angle_weighted_kernel_sum;
                                if(BH_angle_weighted_kernel_sum<=0) mom_wt=0;
                                /* add initial L/c optical/UV coupling to the gas at the dust sublimation radius */
                                v_kick = mom_wt * mom / P[j].Mass;
                                for(k = 0; k < 3; k++)
                                {
                                    P[j].Vel[k] += v_kick*All.cf_atime*dir[k];
                                    SphP[j].VelPred[k] += v_kick*All.cf_atime*dir[k];
                                }
#endif // BH_PHOTONMOMENTUM
#if defined(BH_BAL_WINDS) && (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))
                                mom_wt = bh_angleweight_localcoupling(j,BH_disk_hr,theta) / BH_angle_weighted_kernel_sum;
                                m_wind = mom_wt * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * mdot*dt; /* mass to couple */
                                if(BH_angle_weighted_kernel_sum<=0) m_wind=0;
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
                                for(e_wind=0,k=0;k<3;k++)
                                {
                                    // relative wind-particle velocity (in code units) including BH-particle motion;
                                    norm = All.cf_atime*All.BAL_v_outflow*dir[k] + velocity[k]-P[j].Vel[k];
                                    // momentum conservation gives the following change in velocities
                                    P[j].Vel[k] += norm * m_wind/P[j].Mass;
                                    SphP[j].VelPred[k] += norm * m_wind/P[j].Mass;
                                    // and the shocked wind energy is given by
                                    e_wind += (norm/All.cf_atime)*(norm/All.cf_atime);
                                }
                                e_wind *= 0.5*m_wind;
                                /* now add wind shock energy to particle */
                                e_wind *= 1 / P[j].Mass;
                                SphP[j].InternalEnergy += e_wind;
                                SphP[j].InternalEnergyPred += e_wind;
#endif //defined(BH_BAL_WINDS) && (!defined(BH_GRAVCAPTURE_SWALLOWS) || defined(BH_GRAVCAPTURE_NOGAS))
#endif // defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                            } // norm > 0
                        } // (P[j].Mass>0)&&(P[j].SwallowID==0)
                    } // P[j].Type==0
                } // (mom>0)&&(BH_angle_weighted_kernel_sum>0)
                
                
                
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
        }
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        BlackholeDataPasser[target].accreted_Mass = accreted_mass;
        BlackholeDataPasser[target].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BPP(target).BH_Mass_AlphaDisk += accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++)
            BlackholeDataPasser[target].accreted_momentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
        BPP(target).b7.dBH_accreted_BHMass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
        BPP(target).b8.dBH_accreted_BHMass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
        BPP(target).BH_CountProgs += accreted_BH_progs;
#endif
#ifdef KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER
        for(k = 0; k < 3; k++)
            BlackholeDataPasser[target].BH_SwallowPos[k] = bh_swallowpos[k];
#endif
#ifdef GALSF
        if(P[target].StellarAge > accreted_age)
            P[target].StellarAge = accreted_age;
#endif
    }
    else
    {
        BlackholeDataResult[target].Mass = accreted_mass;
        BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeDataResult[target].BH_Mass_AlphaDisk = accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++)
            BlackholeDataResult[target].accreted_momentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
        BlackholeDataResult[target].BH_Mass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
        BlackholeDataResult[target].BH_Mass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
        BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
#ifdef GALSF
        BlackholeDataResult[target].Accreted_Age = accreted_age;
#endif
    }
    
    return 0;
} /* closes bh_evaluate_swallow */




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



#ifdef BH_BUBBLES
void bh_bubble(double bh_dmass, MyFloat center[3], MyIDType BH_id)
{
    double phi, theta;
    double dx, dy, dz, rr, r2, dE;
    double E_bubble, totE_bubble;
    double BubbleDistance = 0.0, BubbleRadius = 0.0, BubbleEnergy = 0.0;
    double ICMDensity;
    double Mass_bubble, totMass_bubble;
    double u_to_temp_fac;
    MyDouble pos[3];
    int numngb, tot_numngb, startnode, numngb_inbox;
    int n, i, j, dummy;
    
#ifdef CR_BUBBLES
    double tinj = 0.0, instant_reheat = 0.0;
    double sum_instant_reheat = 0.0, tot_instant_reheat = 0.0;
#endif
    
    u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS /
    BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
    
    if(All.ComovingIntegrationOn)
    {
        
        BubbleDistance = All.BubbleDistance;
        BubbleRadius = All.BubbleRadius;
        
        /*switch to comoving if it is assumed that Rbub should be constant with redshift */
        
        /* BubbleDistance = All.BubbleDistance / All.Time;
         BubbleRadius = All.BubbleRadius / All.Time; */
    }
    else
    {
        BubbleDistance = All.BubbleDistance;
        BubbleRadius = All.BubbleRadius;
    }
    
    BubbleEnergy = All.RadioFeedbackFactor * All.BlackHoleRadiativeEfficiency * bh_dmass * All.UnitMass_in_g / All.HubbleParam * pow(C, 2);	/*in cgs units */
    
    phi = 2 * M_PI * get_random_number(BH_id);
    theta = acos(2 * get_random_number(BH_id + 1) - 1);
    rr = pow(get_random_number(BH_id + 2), 1. / 3.) * BubbleDistance;
    
    pos[0] = sin(theta) * cos(phi);
    pos[1] = sin(theta) * sin(phi);
    pos[2] = cos(theta);
    
    for(i = 0; i < 3; i++)
        pos[i] *= rr;
    
    for(i = 0; i < 3; i++)
        pos[i] += center[i];
    
    
    /* First, let's see how many particles are in the bubble of the default radius */
    
    numngb = 0;
    E_bubble = 0.;
    Mass_bubble = 0.;
    
    startnode = All.MaxPart;
    do
    {
        numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);
        
        for(n = 0; n < numngb_inbox; n++)
        {
            j = Ngblist[n];
            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];
            
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
            dx = NEAREST_X(dx);
            dy = NEAREST_Y(dy);
            dz = NEAREST_Z(dz);
#endif
            r2 = dx * dx + dy * dy + dz * dz;
            
            if(r2 < BubbleRadius * BubbleRadius)
            {
                if(P[j].Type == 0)
                {
                    numngb++;
                    
                    E_bubble += P[j].Mass * Particle_Internal_energy_i(j);
                    Mass_bubble += P[j].Mass;
                }
            }
        }
    }
    while(startnode >= 0);
    
    
    MPI_Allreduce(&numngb, &tot_numngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&E_bubble, &totE_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Mass_bubble, &totMass_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    totE_bubble *= All.UnitEnergy_in_cgs;
    
    if(totMass_bubble > 0)
    {
        if(ThisTask == 0)
        {
            printf("found %d particles in bubble with energy %g and total mass %g \n",
                   tot_numngb, totE_bubble, totMass_bubble);
            fflush(stdout);
        }
        
        
        /*calculate comoving density of ICM inside the bubble */
        
        ICMDensity = totMass_bubble / (4.0 * M_PI / 3.0 * pow(BubbleRadius, 3));
        
        if(All.ComovingIntegrationOn)
            ICMDensity = ICMDensity * All.cf_a3inv;	/*now physical */
        
        /*Rbub=R0*[(Ejet/Ejet,0)/(rho_ICM/rho_ICM,0)]^(1./5.) - physical */
        
        rr = rr / BubbleDistance;
        
        BubbleRadius =
        All.BubbleRadius * pow((BubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)),
                               1. / 5.);
        
        BubbleDistance =
        All.BubbleDistance * pow((BubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)),
                                 1. / 5.);
        
        if(All.ComovingIntegrationOn)
        {
            /*switch to comoving if it is assumed that Rbub should be constant with redshift */
            /* BubbleRadius = BubbleRadius / All.Time;
             BubbleDistance = BubbleDistance / All.Time; */
        }
        
        /*recalculate pos */
        rr = rr * BubbleDistance;
        
        pos[0] = sin(theta) * cos(phi);
        pos[1] = sin(theta) * sin(phi);
        pos[2] = cos(theta);
        
        for(i = 0; i < 3; i++)
            pos[i] *= rr;
        
        for(i = 0; i < 3; i++)
            pos[i] += center[i];
        
        /* now find particles in Bubble again,
         and recalculate number, mass and energy */
        
        numngb = 0;
        E_bubble = 0.;
        Mass_bubble = 0.;
        tot_numngb = 0;
        totE_bubble = 0.;
        totMass_bubble = 0.;
        
        startnode = All.MaxPart;
        
        do
        {
            numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = Ngblist[n];
                dx = pos[0] - P[j].Pos[0];
                dy = pos[1] - P[j].Pos[1];
                dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
                dx = NEAREST_X(dx);
                dy = NEAREST_Y(dy);
                dz = NEAREST_Z(dz);
#endif
                r2 = dx * dx + dy * dy + dz * dz;
                
                if(r2 < BubbleRadius * BubbleRadius)
                {
                    if(P[j].Type == 0 && P[j].Mass > 0)
                    {
                        numngb++;
                        
                        E_bubble += P[j].Mass * Particle_Internal_energy_i(j);
                        Mass_bubble += P[j].Mass;
                    }
                }
            }
        }
        while(startnode >= 0);
        
        
        MPI_Allreduce(&numngb, &tot_numngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&E_bubble, &totE_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Mass_bubble, &totMass_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        totE_bubble *= All.UnitEnergy_in_cgs;
        
        if(totMass_bubble > 0)
        {
            if(ThisTask == 0)
            {
                printf("found %d particles in bubble of rescaled radius with energy %g and total mass %g \n",
                       tot_numngb, totE_bubble, totMass_bubble);
                printf("energy shall be increased by: (Eini+Einj)/Eini = %g \n",
                       (BubbleEnergy + totE_bubble) / totE_bubble);
                fflush(stdout);
            }
        }
        
        /* now find particles in Bubble again, and inject energy */
        
#ifdef CR_BUBBLES
        sum_instant_reheat = 0.0;
        tot_instant_reheat = 0.0;
#endif
        
        startnode = All.MaxPart;
        
        do
        {
            numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = Ngblist[n];
                dx = pos[0] - P[j].Pos[0];
                dy = pos[1] - P[j].Pos[1];
                dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
                dx = NEAREST_X(dx);
                dy = NEAREST_Y(dy);
                dz = NEAREST_Z(dz);
#endif
                r2 = dx * dx + dy * dy + dz * dz;
                
                if(r2 < BubbleRadius * BubbleRadius)
                {
                    if(P[j].Type == 0 && P[j].Mass > 0)
                    {
                        /* energy we want to inject in this particle */
                        
                        if(All.StarformationOn)
                            dE = ((BubbleEnergy / All.UnitEnergy_in_cgs) / totMass_bubble) * P[j].Mass;
                        else
                            dE = (BubbleEnergy / All.UnitEnergy_in_cgs) / tot_numngb;
                        
                        if(u_to_temp_fac * dE / P[j].Mass > 5.0e9)
                            dE = 5.0e9 * P[j].Mass / u_to_temp_fac;
                        
#ifndef CR_BUBBLES
                        SphP[j].InternalEnergy += dE / P[j].Mass;
#else
                        
                        tinj = 10.0 * All.HubbleParam * All.cf_hubble_a / All.UnitTime_in_Megayears;
                        
                        instant_reheat =
                        CR_Particle_SupernovaFeedback(&SphP[j], dE / P[j].Mass * All.CR_AGNEff, tinj);
                        
                        if(instant_reheat > 0)
                        {
                            SphP[j].InternalEnergy += instant_reheat;
                        }
                        
                        if(All.CR_AGNEff < 1)
                        {
                            SphP[j].InternalEnergy += (1 - All.CR_AGNEff) * dE / P[j].Mass;
                        }
                        
                        
                        sum_instant_reheat += instant_reheat * P[j].Mass;
#endif
                        
                    }
                }
            }
        }
        while(startnode >= 0);
        
#ifdef CR_BUBBLES
        MPI_Allreduce(&sum_instant_reheat, &tot_instant_reheat, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if(ThisTask == 0)
        {
            printf("Total BubbleEnergy %g Thermalized Energy %g \n", BubbleEnergy,
                   tot_instant_reheat * All.UnitEnergy_in_cgs);
            fflush(stdout);
            
        }
#endif
    }
    else
    {
        if(ThisTask == 0)
        {
            printf("No particles in bubble found! \n");
            fflush(stdout);
        }
        
    }
    
}
#endif /* end of BH_BUBBLE */




/* routine to return the values we need of the properties of the gas, stars, etc in the vicinity of the BH -- these all factor into the BHAR */
int blackhole_evaluate_PREPASS(int target, int mode, int *nexport, int *nSend_local)
{
    /* initialize variables before SPH loop is started */
    int startnode, numngb, j, k, n, id, listindex=0;
    MyFloat *pos, h_i, *vel, u, wk, dwk, hinv, hinv3;
    wk=dwk=u=0;
    double dP[3],dv[3],wt;
    struct blackholedata_topass out;
    memset(&out, 0, sizeof(struct blackholedata_topass));
    
    if(mode == 0)
    {
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        vel = BlackholeDataGet[target].Vel;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
    }
    
    if(h_i < 0) return -1;
    hinv = 1./h_i;
    hinv3 = hinv*hinv*hinv;
    
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
                    for (k=0;k<3;k++)
                    {
                        dP[k]=P[j].Pos[k]-pos[k];
#ifdef PERIODIC
                        dP[k]=NEAREST(dP[k]);
#endif
                        dv[k]=P[j].Vel[k]-vel[k];
                    }
#ifdef BH_DYNFRICTION
                    for (k=0;k<3;k++)
                    {
                        out.DF_mean_vel[k] += wt*dv[k];
                        out.DF_rms_vel += wt*dv[k]*dv[k];
                        if(P[j].Mass>out.DF_mmax_particles) out.DF_mmax_particles=P[j].Mass;
                    }
#endif
                    if(P[j].Type==0) /* gas */
                    {
                        out.Mgas_in_Kernel += wt;
                        out.BH_InternalEnergy += wt*SphP[j].InternalEnergy;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                        out.Jgas_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jgas_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jgas_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                        u=0; for(k=0;k<3;k++) u+=dP[k]*dP[k]; u=sqrt(u)/h_i;
                        kernel_main(u,hinv3,hinv3*hinv,&wk,&dwk,1);
                        dwk /= u*h_i;
                        for(k=0;k<3;k++)
                        {
                            out.GradRho_in_Kernel[k] += wt * dwk * fabs(dP[k]);
                        }
#endif
#if defined(BH_USE_GASVEL_IN_BONDI) || defined(BH_DRAG)
                        for(k=0;k<3;k++)
                        {
                            out.BH_SurroundingGasVel[k] += wt*dv[k];
                        }
#endif
                    } else { /* not gas */
                        out.Malt_in_Kernel += wt;
                        out.Jalt_in_Kernel[0] += wt*(dP[1]*dv[2] - dP[2]*dv[1]);
                        out.Jalt_in_Kernel[1] += wt*(dP[2]*dv[0] - dP[0]*dv[2]);
                        out.Jalt_in_Kernel[2] += wt*(dP[0]*dv[1] - dP[1]*dv[0]);
                    }
                } // ( (P[j].Mass > 0) && (P[j].Type != 5) && (P[j].ID != id) )
            } // for(n = 0; n < numngb; n++)
            
            if(mode == 0)
                out2particle_blackhole(&out, target, 0);
            else
                BlackholeDataPasserResult[target] = out;
            
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            } // if(listindex < NODELISTLENGTH)
        } // if(mode == 1)
    } // while(startnode >= 0)
    return 0;
}





#endif // BLACK_HOLES
