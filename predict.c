#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"

/* Routines for the drift/predict step */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially in detail (although the actual algorithm 
 * structure remains essentially the same) 
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

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
    {
        P[i].Pos[j] += P[i].Vel[j] * dt_drift;
    }
#ifdef ONEDIM
    P[i].Pos[1]=P[i].Pos[2]=0;
#endif
#ifdef TWODIMS
    P[i].Pos[2]=0;
#endif
    
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
            
#if defined(TURB_DRIVING)
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += SphP[i].TurbAccel[j] * dt_gravkick;
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
            
            /* check for reflecting boundaries: if so, do the reflection! */
#if defined(REFLECT_BND_X) || defined(REFLECT_BND_Y) || defined(REFLECT_BND_Z)
            double box_upper[3]; box_upper[0]=box_upper[1]=box_upper[2]=1;
#ifdef PERIODIC
            box_upper[0]=boxSize_X; box_upper[1]=boxSize_Y; box_upper[2]=boxSize_Z;
#endif
            for(j = 0; j < 3; j++)
            {
                /* skip the non-reflecting boundaries */
#ifndef REFLECT_BND_X
                if(j==0) continue;
#endif
#ifndef REFLECT_BND_Y
                if(j==1) continue;
#endif
#ifndef REFLECT_BND_Z
                if(j==2) continue;
#endif
                if(P[i].Pos[j] <= 0)
                {
                    if(P[i].Vel[j]<0) {P[i].Vel[j]=-P[i].Vel[j]; SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;}
                    P[i].Pos[j]=(0+((double)P[i].ID)*1.e-6)*box_upper[j];
                }
                if(P[i].Pos[j] >= box_upper[j])
                {
                    if(P[i].Vel[j]>0) {P[i].Vel[j]=-P[i].Vel[j]; SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;}
                    P[i].Pos[j]=box_upper[j]*(1-((double)P[i].ID)*1.e-6);
                }
            }
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
            
            PPP[i].Hsml *= exp(divVel * dt_drift / NUMDIMS);
            
            if(PPP[i].Hsml < All.MinHsml)
                PPP[i].Hsml = All.MinHsml;
            
            if(PPP[i].Hsml > All.MaxHsml)
                PPP[i].Hsml = All.MaxHsml;
            
            drift_sph_extra_physics(i, time0, time1, dt_entr);

        
            SphP[i].Pressure = get_pressure(i);
#ifdef GAMMA_ENFORCE_ADIABAT
            SphP[i].InternalEnergyPred = SphP[i].Pressure / (SphP[i].Density * GAMMA_MINUS1);
#endif
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
    
#ifdef GAMMA_ENFORCE_ADIABAT
    press = GAMMA_ENFORCE_ADIABAT * pow(SphP[i].Density, GAMMA);
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
    SphP[i].dp_drho = (res.p.drho + res.temp * gsl_pow_2(res.p.dtemp / (SphP[i].Density * All.UnitDensity_in_cgs)) / res.e.dtemp);
#endif
    
    
#ifdef COSMIC_RAYS
    press += Get_Particle_CosmicRayPressure(i);
#endif
    
    
#ifdef TRUELOVE_CRITERION_PRESSURE
    /* add an extra pressure term to suppress fragmentation at/below the explicit resolution scale */
    MyFloat xJeans;
#ifdef HYDRO_SPH
    /* robertson & kravtsov modification of Bate & Burkert formulation */
    xJeans=(1.83*2.0/GAMMA)*All.G*PPP[i].Hsml*PPP[i].Hsml*SphP[i].Density*SphP[i].Density;
    //xJeans*=1.5; /* above is NJeans=5, this is NJeans=9 */
#else
    /* standard finite-volume formulation of this */
    double NJeans = 2; // set so that resolution = lambda_Jeans/NJeans //
    double h_eff = Get_Particle_Size(i);
    xJeans = NJeans*NJeans/(M_PI*GAMMA) * All.G * h_eff*h_eff * SphP[i].Density*SphP[i].Density;
#endif
    if(All.ComovingIntegrationOn) xJeans *= All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
#endif
    
    return press;
}



void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
    int k;
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    for(k=0;k<3;k++) {SphP[i].BPred[k] += SphP[i].DtB[k] * dt_entr * BphysVolphys_to_BcodeVolCode;} // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1;
    double dtphi_code = (PhiphysVolphys_to_PhicodeVolCode) * SphP[i].DtPhi;
    SphP[i].PhiPred += dtphi_code  * dt_entr;
    double t_damp = Get_Particle_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        SphP[i].PhiPred *= exp( -dt_entr * t_damp );
    }
#endif
#endif
#ifdef COSMIC_RAYS
    double etmp = SphP[i].CosmicRayEnergyPred + SphP[i].DtCosmicRayEnergy * dt_entr;
    if(etmp<0.5*SphP[i].CosmicRayEnergyPred) {SphP[i].CosmicRayEnergyPred *= 0.5;} else {SphP[i].CosmicRayEnergyPred=etmp;}
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
    {
        boxsize[j] = All.BoxSize;
    }
    
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
    {
        for(j = 0; j < 3; j++)
        {
            while(P[i].Pos[j] < 0)
            {
                P[i].Pos[j] += boxsize[j];
#ifdef SHEARING_BOX
                if(j==0) {
                    P[i].Vel[SHEARING_BOX_PHI_COORDINATE]-=Shearing_Box_Vel_Offset;
                    if(P[i].Type==0) SphP[i].VelPred[SHEARING_BOX_PHI_COORDINATE]-=Shearing_Box_Vel_Offset;}
#if (SHEARING_BOX != 1)
                /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                // PFH: needs update for shearing coordinates, if desired //
                //if(j==0) {P[i].Pos[SHEARING_BOX_PHI_COORDINATE]+= Shearing_Box_Vel_Offset * All.Time;
#endif
#endif
            }
            
            while(P[i].Pos[j] >= boxsize[j])
            {
                P[i].Pos[j] -= boxsize[j];
#ifdef SHEARING_BOX
                if(j==0) {
                    P[i].Vel[SHEARING_BOX_PHI_COORDINATE]+=Shearing_Box_Vel_Offset;
                    if(P[i].Type==0) SphP[i].VelPred[SHEARING_BOX_PHI_COORDINATE]+=Shearing_Box_Vel_Offset;}
#endif
            }
        }
    }
}
#endif




/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */


/* this function returns the effective (grid-equivalent) particle 'size'; useful for things like 
    time-stepping and limiter functions */
double INLINE_FUNC Get_Particle_Size(int i)
{
#ifdef ONEDIM
    return 2 / PPP[i].NumNgb * PPP[i].Hsml;
#else
#ifdef TWODIMS
    return sqrt( (M_PI/2.0) / PPP[i].NumNgb ) * PPP[i].Hsml;
#else
    return pow( 4.0 * M_PI / (3.0 * PPP[i].NumNgb) , 1./3. ) * PPP[i].Hsml;
#endif
#endif
}



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
    return GAMMA_MINUS1 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i);
}

#ifdef COSMIC_RAYS
double INLINE_FUNC Get_Particle_CosmicRayPressure(int i)
{
    return GAMMA_COSMICRAY_MINUS1 * (SphP[i].CosmicRayEnergyPred * SphP[i].Density) / P[i].Mass; // cosmic ray pressure = (4/3-1) * e_cr = 1/3 * (E_cr/Vol) //

}
#endif


double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{
#ifdef EOS_DEGENERATE
    return sqrt(SphP[i].dp_drho);
#endif
#ifdef COSMIC_RAYS
    return sqrt(GAMMA*GAMMA_MINUS1 * SphP[i].InternalEnergyPred + GAMMA_COSMICRAY*GAMMA_COSMICRAY_MINUS1 * SphP[i].CosmicRayEnergyPred);
#endif
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
}

#ifdef MAGNETIC
double INLINE_FUNC Get_Particle_BField(int i_particle_id, int k_vector_component)
{
    return SphP[i_particle_id].BPred[k_vector_component] * SphP[i_particle_id].Density / P[i_particle_id].Mass;
}

/* this function is needed to control volume fluxes of the normal components of B and phi in the 
    -bad- situation where the meshless method 'faces' do not properly close (usually means you are 
    using boundary conditions that you should not) */
double Get_DtB_FaceArea_Limiter(int i)
{
#ifdef HYDRO_SPH
    return 1;
#else
    /* define some variables */
    double dt_entr = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
    int j; double dB[3],dBmag=0,Bmag=0;
    /* check the magnitude of the predicted change in B-fields, vs. B-magnitude */
    for(j=0;j<3;j++)
    {
        dB[j] = SphP[i].DtB[j] * dt_entr / All.cf_atime; /* converts to code units of Vol_code*B_code = Vol_phys*B_phys/a */
        dBmag += dB[j]*dB[j];
        Bmag += SphP[i].BPred[j]*SphP[i].BPred[j];
    }
    dBmag = sqrt(dBmag); Bmag = sqrt(Bmag);
    /* also make sure to check the actual pressure, since if P>>B, we will need to allow larger changes in B per timestep */
    double P_BV_units = sqrt(2.*SphP[i].Pressure*All.cf_a3inv)*P[i].Mass/SphP[i].Density * All.cf_afac3 / All.cf_a2inv;
    /* the above should be in CODE Bcode*Vol_code units! */
    double Bmag_max = DMAX(Bmag, DMIN( P_BV_units, 10.*Bmag ));
    /* now check how accurately the cell is 'closed': the face areas are ideally zero */
    double area_sum = fabs(SphP[i].Face_Area[0])+fabs(SphP[i].Face_Area[1])+fabs(SphP[i].Face_Area[2]);
    /* but this needs to be normalized to the 'expected' area given Hsml */
#ifdef ONEDIM
    double area_norm = 2.;
#else
#ifdef TWODIMS
    double area_norm = 2.*M_PI * PPP[i].Hsml * All.cf_atime;
#else
    double area_norm = 4.*M_PI * PPP[i].Hsml*PPP[i].Hsml * All.cf_atime*All.cf_atime;
#endif
#endif
    /* ok, with that in hand, define an error tolerance based on this */
    if(area_norm>0)
    {
        double area_norm_min_threshold = 0.001;
        double area_norm_weight = 200.0;
#ifdef PM_HIRES_REGION_CLIPPING
        area_norm_min_threshold *= 0.01;
        area_norm_weight *= 2.5; // can be as low as 1.0 (PFH) //
#endif
        if(area_sum/area_norm > area_norm_min_threshold)
        {
            double tol = (All.CourantFac/0.2) * DMAX( 0.01, area_norm/(area_norm_weight * area_sum) );
            tol *= Bmag_max; /* give the limiter dimensions */
            if(dBmag > tol) {return tol/dBmag;} /* now actually check if we exceed this */
        }
    }
    return 1;
#endif
}


#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Particle_PhiField(int i_particle_id)
{
    return SphP[i_particle_id].PhiPred * SphP[i_particle_id].Density / P[i_particle_id].Mass;
}

double INLINE_FUNC Get_Particle_PhiField_DampingTimeInv(int i_particle_id)
{
    /* this timescale should always be returned as a -physical- time */
#ifdef HYDRO_SPH
    /* PFH: add simple damping (-phi/tau) term */
    double damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (SphP[i_particle_id].MaxSignalVel*All.cf_afac3 / (All.cf_atime*Get_Particle_Size(i_particle_id)));
    /* PFH: add div_v term from Tricco & Price to DtPhi */
    // damping_tinv += 0.5 * P[i_particle_id].Particle_DivVel*All.cf_a2inv; // this is not needed for the form of Vphi we evolve //
#else
    double damping_tinv;
    if((All.StarformationOn)||(All.ComovingIntegrationOn))
    {
        // only see a small performance drop from fastestwavespeed above to maxsignalvel below, despite the fact that below is purely local (so allows more flexible adapting to high dynamic range)
        damping_tinv = 0.0;
        
        if(PPP[i_particle_id].Hsml > 0)
        {
            double vsig2 = 0.5 * All.cf_afac3 * fabs(SphP[i_particle_id].MaxSignalVel);
            double phi_B_eff = 0.0;
            if(vsig2 > 0) {phi_B_eff = Get_Particle_PhiField(i_particle_id) / (All.cf_atime * vsig2);}
            double vsig1 = 0.0;
            if(SphP[i_particle_id].Density > 0)
            {
                vsig1 = All.cf_afac3 *
                sqrt( Particle_effective_soundspeed_i(i_particle_id)*Particle_effective_soundspeed_i(i_particle_id) +
                     (All.cf_afac1 / All.cf_atime) *
                     (Get_Particle_BField(i_particle_id,0)*Get_Particle_BField(i_particle_id,0) +
                      Get_Particle_BField(i_particle_id,1)*Get_Particle_BField(i_particle_id,1) +
                      Get_Particle_BField(i_particle_id,2)*Get_Particle_BField(i_particle_id,2) +
                      phi_B_eff*phi_B_eff) / SphP[i_particle_id].Density );
            }
            double vsig_max = DMAX( DMAX(vsig1,vsig2) , 0.1 * All.FastestWaveSpeed );            
            damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (vsig_max / (All.cf_atime*Get_Particle_Size(i_particle_id)));
        }
    } else {
        damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveSpeed / Get_Particle_Size(i_particle_id); // fastest wavespeed has units of [vphys]
        //double damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveDecay * All.cf_a2inv; // no improvement over fastestwavespeed; decay has units [vphys/rphys]
    }
#endif
    return damping_tinv;
}

#endif // dedner
#endif // magnetic
