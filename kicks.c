#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

void apply_long_range_kick(integertime, integertime);

void do_first_halfstep_kick(void)
{
    int i;
    integertime ti_step, tstart=0, tend=0;
    
#if defined (VS_TURB) || defined (AB_TURB) || defined (TURB_DRIVING)
    do_turb_driving_step_first_half();
#endif
    
#ifdef PMGRID
    if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif
    
    /* collisionless particles only need an update if they are active; however, to 
        maintain manifest conservation in the hydro, need to check -ALL- sph particles every timestep */
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Mass > 0)
        {
            /* 'full' kick for active particles */
            if(TimeBinActive[P[i].TimeBin])
            {
                ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
                tstart = P[i].Ti_begstep;	/* beginning of step */
                tend = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
            }
            do_the_kick(i, tstart, tend, P[i].Ti_current, 0);
        }
    } // for(i = 0; i < NumPart; i++) // 
}


void do_second_halfstep_kick(void)
{
    int i;
    integertime ti_step, tstart=0, tend=0;
    
#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep + ti_step / 2;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif

    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Mass > 0)
        {
            /* 'full' kick for active particles */
            if(TimeBinActive[P[i].TimeBin])
            {
                ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
                tstart = P[i].Ti_begstep + ti_step / 2;	/* midpoint of step */
                tend = P[i].Ti_begstep + ti_step;	/* end of step */
            }
            do_the_kick(i, tstart, tend, P[i].Ti_current, 1);
            if(TimeBinActive[P[i].TimeBin])
            {
                if(P[i].Type == 0 && P[i].Mass > 0)
                {
#ifdef EOS_DEGENERATE
                    for(j = 0; j < 3; j++)
                        SphP[i].xnucPred[j] = SphP[i].xnuc[j];
#endif
                    SphP[i].Pressure = get_pressure(i);
                    set_predicted_sph_quantities_for_extra_physics(i);
                }
            }
        }
    } // for(i = 0; i < NumPart; i++) //
    
#if defined (VS_TURB) || defined (AB_TURB) || defined (TURB_DRIVING)
    do_turb_driving_step_second_half();
#endif
}


#ifdef PMGRID
void apply_long_range_kick(integertime tstart, integertime tend)
{
    int i, j;
    double dt_gravkick, dvel[3];
    
    if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
    else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;
    
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Mass > 0)
            for(j = 0; j < 3; j++)	/* do the kick, only collisionless particles */
            {
                dvel[j] = P[i].GravPM[j] * dt_gravkick;
                P[i].Vel[j] += dvel[j];
                P[i].dp[j] += P[i].Mass * dvel[j];
            }
#ifdef DISTORTIONTENSORPS
        do_long_range_phase_space_kick(i, dt_gravkick);
#endif
    }
}
#endif


void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode)
{
    int j;
    double dp[3], dt_entr, dt_gravkick, dt_hydrokick;
    double mass_old, mass_pred, mass_new;
    mass_old = mass_pred = mass_new = P[i].Mass;    
    
    /* First, we do the pure hydro update for gas (because we use a total energy equation, its much easier to 
        do this than to deal with the gravitational force first and then have to subtract it back out) */

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* need to do the slightly more complicated update scheme to maintain exact mass conservation */
    double dMass, d_inc = 0.5; // fraction of delta_conserved to couple per kick step (each 'kick' is 1/2-timestep) //
    //double dv[3], v_old[3], dMass, ent_old=0;
    if(P[i].Type==0)
    {
        //ent_old = SphP[i].InternalEnergy;
        //for(j=0;j<3;j++) v_old[j] = P[i].Vel[j];
        if(SphP[i].dMass != 0)
        {
            // update the --conserved-- variables of each particle //
            if(mode != 0)
            {
                dMass = (tend - tstart) * All.Timebase_interval / All.cf_hubble_a * SphP[i].DtMass;
                if(dMass >= SphP[i].dMass) dMass = SphP[i].dMass; // try to get close to what the time-integration scheme would give //
                SphP[i].dMass -= dMass;
            } else {
                dMass = SphP[i].dMass;
            }
            if(fabs(dMass) > 0.9*SphP[i].MassTrue) dMass *= 0.9*SphP[i].MassTrue/fabs(dMass); // limiter to prevent madness //

            
            // load and update the particle masses //
            // particle mass update here, from hydro fluxes //
            mass_old = SphP[i].MassTrue;
            mass_pred = P[i].Mass;
            mass_new = mass_old + dMass;
            SphP[i].MassTrue = mass_new;
            // UNITS: remember all time derivatives (DtX, dX) are in -physical- units; as are mass, entropy/internal energy, but -not- velocity //
            /*
            double e_old = mass_old * SphP[i].InternalEnergy;
            for(j = 0; j< 3; j++) e_old += 0.5*mass_old * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // physical //
            
            // do the momentum-space kick //
            for(j = 0; j < 3; j++)
            {
                dp[j] = d_inc * SphP[i].dMomentum[j]; // dv[j] = SphP[i].HydroAccel[j] * dt_hydrokick; //non-conservative// //
                // now update the velocity based on the total momentum change //
                P[i].Vel[j] = (mass_old*P[i].Vel[j] + dp[j]*All.cf_atime) / mass_new; // call after tabulating dP[j] //
            }
            
            // kick for gas internal energy/entropy //
            e_old += d_inc * SphP[i].dInternalEnergy; // increment of total (thermal+kinetic) energy //
            for(j = 0; j< 3; j++) e_old -= 0.5*mass_new * (P[i].Vel[j]/All.cf_atime)*(P[i].Vel[j]/All.cf_atime); // subtract off the new kinetic energy //
            SphP[i].InternalEnergy = e_old / mass_new; // obtain the new internal energy per unit mass //
            check_particle_for_temperature_minimum(i); // if we've fallen below the minimum temperature, force the 'floor' //
            */
             
            // at the end of this kick, need to re-zero the dInternalEnergy, and other
            // conserved-variable SPH quantities set in the hydro loop, to avoid double-counting them
            if(mode==0)
            {
                //SphP[i].dInternalEnergy = 0;
                //SphP[i].dMomentum[0] = SphP[i].dMomentum[1] = SphP[i].dMomentum[2] = 0;
                SphP[i].dMass = 0;
            }
        }
    } // if(P[i].Type==0) //
#endif
    
    /* only enter the 'normal' kick loop below for genuinely active particles */
    if(TimeBinActive[P[i].TimeBin])
    {
        /* get the timestep (physical units for dt_entr and dt_hydrokick) */
        if(All.ComovingIntegrationOn)
        {
            dt_entr = dt_hydrokick = (tend - tstart) * All.Timebase_interval / All.cf_hubble_a;
            dt_gravkick = get_gravkick_factor(tstart, tend);
        }
        else
        {
            dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
        }
        
        if(P[i].Type==0)
        {
            double grav_acc[3], dEnt_Gravity = 0;
            for(j = 0; j < 3; j++)
            {
                grav_acc[j] = All.cf_a2inv * P[i].GravAccel[j];
#ifdef PMGRID
                grav_acc[j] += All.cf_a2inv * P[i].GravPM[j];
#endif
            }
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
            for(j=0;j<3;j++) {dEnt_Gravity += -(SphP[i].GravWorkTerm[j] * All.cf_atime * dt_hydrokick) * grav_acc[j];}
#endif
            double dEnt = SphP[i].InternalEnergy + SphP[i].DtInternalEnergy * dt_hydrokick + dEnt_Gravity;
            
#ifndef HYDRO_SPH
            /* if we're using a Riemann solver, we include an energy/entropy-type switch to ensure
                that we don't corrupt the temperature evolution of extremely cold, adiabatic flows */
            double e_thermal,e_kinetic,e_potential;
            e_potential=0; for(j=0;j<3;j++) {e_potential += grav_acc[j]*grav_acc[j];}
            e_potential = P[i].Mass * sqrt(e_potential) * (KERNEL_CORE_SIZE*PPP[i].Hsml*All.cf_atime); // = M*|a_grav|*h (physical)
            e_kinetic = 0.5 * P[i].Mass * All.cf_a2inv * SphP[i].MaxKineticEnergyNgb;
            e_thermal = DMAX(0.5*SphP[i].InternalEnergy, dEnt) * P[i].Mass;
            int do_entropy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            if(0.01*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#else
            if(0.005*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#endif
            if(0.01*e_potential > e_thermal) {do_entropy=1;}
            // also check for flows which are totally dominated by the adiabatic component of their temperature evolution //
            // double mach = fabs(SphP[i].MaxSignalVel/Particle_effective_soundspeed_i(i) - 2.0); //
            // if(mach < 1.1) {do_entropy=1;} // (actually, this switch tends to do more harm than good!) //
            do_entropy = 0; // seems unstable in tests like interacting blastwaves... //
            if(do_entropy)
            {
                /* use the pure-SPH entropy equation, which is exact up to the mass flux, for adiabatic flows */
                SphP[i].DtInternalEnergy = (SphP[i].Pressure/SphP[i].Density) * P[i].Particle_DivVel*All.cf_a2inv;
                if(All.ComovingIntegrationOn) SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a;
                dEnt = SphP[i].InternalEnergy + SphP[i].DtInternalEnergy * dt_hydrokick; /* gravity term not included here, as it makes this unstable */
            }
#endif
            
            if(dEnt < 0.5*SphP[i].InternalEnergy) {SphP[i].InternalEnergy *= 0.5;} else {SphP[i].InternalEnergy = dEnt;}
            check_particle_for_temperature_minimum(i); /* if we've fallen below the minimum temperature, force the 'floor' */
        }
        
        /* now, kick for non-SPH quantities (accounting for momentum conservation if masses are changing) */
        for(j = 0; j < 3; j++)
        {
            dp[j] = 0;
            if(P[i].Type==0)
                dp[j] += mass_pred * SphP[i].HydroAccel[j] * All.cf_atime * dt_hydrokick; // convert to code units
#ifdef RT_RAD_PRESSURE
            if(P[i].Type==0)
                dp[j] += mass_pred * SphP[i].RadAccel[j] * dt_hydrokick;
#endif
            dp[j] += mass_pred * P[i].GravAccel[j] * dt_gravkick;
#ifdef RELAXOBJECT
            dp[j] -= mass_pred * P[i].Vel[j] * All.RelaxFac * dt_gravkick;
#endif
            P[i].Vel[j] += dp[j] / mass_new; /* correctly accounts for mass change if its allowed */
        }

 
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
                if(P[i].Vel[j]<0) {P[i].Vel[j]=-P[i].Vel[j]; SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0; dp[j]+=2*P[i].Vel[j]*mass_new;}
                P[i].Pos[j]=(0+((double)P[i].ID)*1.e-6)*box_upper[j];
            }
            if(P[i].Pos[j] >= box_upper[j])
            {
                if(P[i].Vel[j]>0) {P[i].Vel[j]=-P[i].Vel[j]; SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0; dp[j]+=2*P[i].Vel[j]*mass_new;}
                P[i].Pos[j]=box_upper[j]*(1-((double)P[i].ID)*1.e-6);
            }
        }
#endif
        /* any other gas-specific kicks (e.g. B-fields, radiation) go here */
        if(P[i].Type==0)
        {
            do_sph_kick_for_extra_physics(i, tstart, tend, dt_entr);

            /* after completion of a full step, set the predicted values of SPH quantities
             * to the current values. They will then predicted further along in drift operations */
            if(mode==1)
            {
                for(j = 0; j < 3; j++)
                    SphP[i].VelPred[j] = P[i].Vel[j];//(mass_old*v_old[j] + dp[j]) / mass_new;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                P[i].Mass = SphP[i].MassTrue; //mass_old + SphP[i].DtMass * dt_hydrokick;
#endif
                SphP[i].InternalEnergyPred = SphP[i].InternalEnergy; //ent_old + SphP[i].DtInternalEnergy * dt_entr;
            }
        }
        
        /* set the momentum shift so we know how to move the tree! */
        for(j=0;j<3;j++) {P[i].dp[j] += dp[j];}
#ifdef DISTORTIONTENSORPS
        /* momentum-space correction for following phase-space distribution (call after momentum-space kicks) */
        do_the_phase_space_kick(i, dt_gravkick);
#endif
        
    } // if(TimeBinActive[P[i].TimeBin]) //
}


void set_predicted_sph_quantities_for_extra_physics(int i)
{
#if defined(MAGNETIC)
    int j1;
    for(j1 = 0; j1 < 3; j1++)
        SphP[i].BPred[j1] = SphP[i].B[j1];
#if defined(DIVBCLEANING_DEDNER)
    SphP[i].PhiPred = SphP[i].Phi;
#endif
#endif
}


void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
    int j; j=0;
#ifdef MAGNETIC
    for(j = 0; j < 3; j++)
        SphP[i].B[j] += SphP[i].DtB[j] * dt_entr;
#ifdef DIVBCLEANING_DEDNER
    SphP[i].Phi += SphP[i].DtPhi * dt_entr;
#endif
#endif
#ifdef NUCLEAR_NETWORK
    for(j = 0; j < EOS_NSPECIES; j++)
        SphP[i].xnuc[j] += SphP[i].dxnuc[j] * dt_entr * All.UnitTime_in_s;
    
    network_normalize(SphP[i].xnuc, &SphP[i].InternalEnergy, &All.nd, &All.nw);
#endif
}

