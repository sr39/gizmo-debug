#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"


void reconstruct_timebins(void)
{
    int i, bin;
    long long glob_sum;
    
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        TimeBinCount[bin] = 0;
        TimeBinCountSph[bin] = 0;
        FirstInTimeBin[bin] = -1;
        LastInTimeBin[bin] = -1;
#ifdef GALSF
        TimeBinSfr[bin] = 0;
#endif
#ifdef BLACK_HOLES
        TimeBin_BH_mass[bin] = 0;
        TimeBin_BH_dynamicalmass[bin] = 0;
        TimeBin_BH_Mdot[bin] = 0;
        TimeBin_BH_Medd[bin] = 0;
#endif
    }
    
    for(i = 0; i < NumPart; i++)
    {
        bin = P[i].TimeBin;
        
        if(TimeBinCount[bin] > 0)
        {
            PrevInTimeBin[i] = LastInTimeBin[bin];
            NextInTimeBin[i] = -1;
            NextInTimeBin[LastInTimeBin[bin]] = i;
            LastInTimeBin[bin] = i;
        }
        else
        {
            FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
            PrevInTimeBin[i] = NextInTimeBin[i] = -1;
        }
        TimeBinCount[bin]++;
        if(P[i].Type == 0)
            TimeBinCountSph[bin]++;
        
#ifdef GALSF
        if(P[i].Type == 0)
            TimeBinSfr[bin] += SphP[i].Sfr;
#endif
#ifdef BLACK_HOLES
        if(P[i].Type == 5)
        {
            TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
            TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
            TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
            TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
        }
#endif
    }
    
    make_list_of_active_particles();
    
    for(i = FirstActiveParticle, NumForceUpdate = 0; i >= 0; i = NextActiveParticle[i])
    {
        NumForceUpdate++;
        if(i >= NumPart)
        {
            printf("Bummer i=%d\n", i);
            terminate("inconsistent list");
        }
    }
    
    sumup_large_ints(1, &NumForceUpdate, &glob_sum);
    
    GlobNumForceUpdate = glob_sum;
}





void drift_particle(int i, integertime time1)
{
    int j;
    double dt_drift;
    
    integertime time0 = P[i].Ti_current;
    
    if(time1 < time0)
    {
        printf("i=%d time0=%d time1=%d\n", i, (int)time0, (int)time1);
        terminate("no prediction into past allowed");
    }
    
    if(time1 == time0)
        return;
    
    if(All.ComovingIntegrationOn)
        dt_drift = get_drift_factor(time0, time1);
    else
        dt_drift = (time1 - time0) * All.Timebase_interval;
    
    for(j = 0; j < 3; j++)
        P[i].Pos[j] += P[i].Vel[j] * dt_drift;
    
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(P[i].Type>0)
    {
        if(dt_drift>0)
        {
            PPP[i].Hsml *= exp(P[i].Particle_DivVel * dt_drift / NUMDIMS);
            if(P[i].Hsml < All.ForceSoftening[P[i].Type])
                P[i].Hsml = All.ForceSoftening[P[i].Type];
        }
    }
#endif
    
#ifdef DISTORTIONTENSORPS
    do_phase_space_drift(i, dt_drift);
#endif
    
    if((P[i].Type == 0) && (P[i].Mass > 0))
        {
            double dt_gravkick, dt_hydrokick, dt_entr;
            
            if(All.ComovingIntegrationOn)
            {
                dt_entr = dt_hydrokick = (time1 - time0) * All.Timebase_interval / All.cf_hubble_a;
                dt_gravkick = get_gravkick_factor(time0, time1);
            }
            else
            {
                dt_entr = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
            }
            
            
            double divVel = P[i].Particle_DivVel;
#ifdef PMGRID
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick +
                    SphP[i].HydroAccel[j]*All.cf_atime * dt_hydrokick; /* make sure v is in code units */
#else
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick +
                    SphP[i].HydroAccel[j]*All.cf_atime * dt_hydrokick; /* make sure v is in code units */
#endif
            
#if defined (VS_TURB) || defined (AB_TURB)
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += SphP[i].TurbAccel[j] * dt_gravkick;
#endif
            
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
            bp_cr_diff_update(i, dt_entr);
#endif
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P[i].Mass = DMAX(P[i].Mass + SphP[i].DtMass * dt_entr, 0.5 * SphP[i].MassTrue);
#endif
            
            SphP[i].Density *= exp(-divVel * dt_drift);
            double etmp = SphP[i].InternalEnergyPred + SphP[i].DtInternalEnergy * dt_entr;
            if(etmp<0.5*SphP[i].InternalEnergyPred) {SphP[i].InternalEnergyPred *= 0.5;} else {SphP[i].InternalEnergyPred=etmp;}
            if(SphP[i].InternalEnergyPred<All.MinEgySpec) SphP[i].InternalEnergyPred=All.MinEgySpec;
            
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
            SphP[i].EgyWtDensity *= exp(-divVel * dt_drift);
#endif
            
            
#ifdef EOS_DEGENERATE
            {
                double maxfac = dt_entr;
                double xnuc, tmpfac;
                for(j = 0; j < EOS_NSPECIES; j++)
                {
                    xnuc = SphP[i].xnucPred[j] + SphP[i].dxnuc[j] * dt_entr;
                    if(xnuc > 1.0)
                    {
                        tmpfac = (1.0 - xnuc) / SphP[i].dxnuc[j];
                        if(tmpfac < maxfac)
                            maxfac = tmpfac;
                    }
                    if(xnuc < 0.0)
                    {
                        tmpfac = (0.0 - xnuc) / SphP[i].dxnuc[j];
                        if(tmpfac < maxfac)
                            maxfac = tmpfac;
                    }
                }
                
                if(maxfac > 0)
                {
                    for(j = 0; j < EOS_NSPECIES; j++)
                    {
                        SphP[i].xnucPred[j] += SphP[i].dxnuc[j] * maxfac;
                    }
                }
            }
#endif
            
            SphP[i].Pressure = get_pressure(i);
            
            PPP[i].Hsml *= exp(divVel * dt_drift / NUMDIMS);
            
            if(PPP[i].Hsml < All.MinGasHsml)
                PPP[i].Hsml = All.MinGasHsml;
            
            if(PPP[i].Hsml > All.MaxHsml)
                PPP[i].Hsml = All.MaxHsml;
            
            drift_sph_extra_physics(i, time0, time1, dt_entr);
        }
    
    P[i].Ti_current = time1;
}


void check_particle_for_temperature_minimum(int i)
{
    if(All.MinEgySpec)
    {
        if(SphP[i].InternalEnergy < All.MinEgySpec)
        {
            SphP[i].InternalEnergy = All.MinEgySpec;
            SphP[i].DtInternalEnergy = 0;
            //SphP[i].dInternalEnergy = 0;
        }
    }
}



void move_particles(integertime time1)
{
    int i;
    for(i = 0; i < NumPart; i++)
        drift_particle(i, time1);
}


/* return the pressure of particle i */
double get_pressure(int i)
{
    MyFloat press = 0;
    
#if !defined(EOS_DEGENERATE)
    /* this is the 'normal' pressure */
    press = Get_Particle_Pressure(i);
#endif
    
    
#ifdef GALSF_EFFECTIVE_EQS
    /* modify pressure to 'interpolate' between effective EOS and isothermal */
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh)
        press = All.FactorForSofterEQS * press +
        (1 - All.FactorForSofterEQS) * All.cf_afac1 * GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;
#endif
    
    
#ifdef EOS_DEGENERATE
    /* call tabulated eos with physical units */
    struct eos_result res;
    eos_calc_egiven(SphP[i].Density * All.UnitDensity_in_cgs, SphP[i].xnucPred, SphP[i].InternalEnergyPred, &SphP[i].temp, &res);
    press = res.p.v / All.UnitPressure_in_cgs;
    SphP[i].dpdr = (res.p.drho + res.temp * gsl_pow_2(res.p.dtemp / (SphP[i].Density * All.UnitDensity_in_cgs)) / res.e.dtemp);
#endif
    
    
#ifdef COSMIC_RAYS
    int CRpop;
#if defined( CR_UPDATE_PARANOIA )
    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
        CR_Particle_Update(SphP + i, CRpop);
#endif
#ifndef CR_NOPRESSURE
    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
        press += CR_Comoving_Pressure(SphP + i, CRpop);
#endif
#endif
    
#ifdef BP_REAL_CRs
    press += SphP[i].CRpPressure;
#endif
    
    
#ifdef TRUELOVE_CRITERION_PRESSURE
    /* add an extra pressure term to prevent artificial fragmentation */
    MyFloat xJeans;
    xJeans=(1.83*2.0/GAMMA)*All.G*PPP[i].Hsml*PPP[i].Hsml*SphP[i].Density*SphP[i].Density;
    //xJeans*=1.5; /* above is NJeans=5, this is NJeans=9 */
    if(All.ComovingIntegrationOn) xJeans *= All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
#endif
    
    return press;
}



void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
    double dt_mag = dt_entr;
#ifdef DIVBCLEANING_DEDNER
    SphP[i].PhiPred += SphP[i].DtPhi * dt_mag;
#endif
    int j;
    for(j = 0; j < 3; j++)
        SphP[i].BPred[j] += SphP[i].DtB[j] * dt_mag;
#endif
}




/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
    int i, j;
    double boxsize[3];
    
    for(j = 0; j < 3; j++)
        boxsize[j] = All.BoxSize;
    
#ifdef LONG_X
    boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
    boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
    boxsize[2] *= LONG_Z;
#endif
    
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
        {
            while(P[i].Pos[j] < 0)
                P[i].Pos[j] += boxsize[j];
            
            while(P[i].Pos[j] >= boxsize[j])
                P[i].Pos[j] -= boxsize[j];
        }
}
#endif




/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */



double INLINE_FUNC Particle_density_for_energy_i(int i)
{
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    return SphP[i].EgyWtDensity;
#else
    return SphP[i].Density;
#endif
}

double INLINE_FUNC Get_Particle_Pressure(int i)
{
#ifndef EOS_DEGENERATE
    return GAMMA_MINUS1 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i);
#endif
}

double INLINE_FUNC Particle_Internal_energy_i(int i)
{
    return SphP[i].InternalEnergy;
}

double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{
#ifndef EOS_DEGENERATE
    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
#else
    return sqrt(SphP[i].dpdr);
#endif
}

