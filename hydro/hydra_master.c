#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef COSMIC_RAYS
#include "../cosmic_rays/cosmic_rays.h"
#endif
#ifdef MACHNUM
#include "../cosmic_rays/machfinder.h"
#endif
#ifdef JD_DPP
#include "../cosmic_rays/cr_electrons.h"
#endif
#define NDEBUG
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif

#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/*! \file hydra_master.c
 *  \brief Computation of SPH forces and rate of entropy/internal energy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */


/* some very useful notes on the hydro variables in comoving integrations:
 
 v_code = a * v_peculiar/physical (canonical momentum)
 r_code = r_physical / a (comoving coordinates)
 m_code = m_physical
 rho_code = rho_physical * a^3 (from length/mass scaling)
 InternalEnergy_code = InternalEnergy_physical
 Pressure_code =
    InternalEnergy_code * rho_code * (gamma-1) = Pressure_physical * a^3 (energy SPH)
    -- the distinction between these cases should be taken care of in the factors
        All.cf_afac1/2/3, which will correctly assign between the two --
 B_code = a*a * B_physical (comoving magnetic fields)
 Phi_code = B_code*v_code (damping field for Dedner divergence cleaning)
    (note: spec egy of phi field is: phi*phi/(2*mu0*rho*ch*ch); compare Bfield is B*B/(mu0*rho);
    so [phi]~[B]*[ch], where ch is the signal velocity used in the damping equation;
 
 -- Time derivatives (rate of change from hydro forces) here are all 
        assumed to end up in *physical* units ---
 HydroAccel, dMomentum are assumed to end up in *physical* units
    (note, this is different from GADGET's convention, where 
     HydroAccel is in units of (Pcode/rhocode)/rcode)
 DtInternalEnergy and dInternalEnergy are assumed to end up in *physical* units
 DtMass and dMass are assumed to end up in *physical* units
*/


#ifdef MACHNUM
double fac_mu, fac_vsic_fix, fac_egy;
#else
static double fac_mu, fac_vsic_fix;
#endif
#ifdef MAGNETIC
static double fac_magnetic_pressure;
#endif


/* --------------------------------------------------------------------------------- */
/* define the kernel structure -- purely for handy purposes to clean up notation */
/* --------------------------------------------------------------------------------- */
/* structure to hold fluxes being passed from the hydro sub-routine */
struct Conserved_var_Riemann
{
    MyDouble rho;
    MyDouble p;
    MyDouble v[3];
};
struct kernel_hydra
{
    double dx, dy, dz;
    double r, vsig, sound_i, sound_j;
    double dvx, dvy, dvz, vdotr2;
    double wk_i, wk_j, dwk_i, dwk_j;
    double h_i, h_j, dwk_ij, rho_ij_inv;
#ifdef HYDRO_SPH
    double p_over_rho2_i;
#endif
#if defined(HYDRO_SPH) || defined(CONDUCTION_EXPLICIT) || defined(TURB_DIFF_ENERGY)
    double spec_egy_u_i;
#endif
#ifdef MAGNETIC
    double mj_r, mf_Ind, b2_i, b2_j;
#ifdef MAGFORCE
    double mf_i, mf_j;
#endif
#if defined(MAGNETIC_SIGNALVEL)
    double alfven2_i, alfven2_j;
#endif
#if defined(MAGNETIC_DISSIPATION)
    double mf_dissInd, mf_dissEnt;
#endif
#endif // MAGNETIC //
};
#ifndef HYDRO_SPH
#include "reimann.h"
#endif


/* --------------------------------------------------------------------------------- */
/* inputs to the routine: put here what's needed to do the calculation! */
/* --------------------------------------------------------------------------------- */
struct hydrodata_in
{
    /* basic hydro variables */
    MyDouble Pos[3];
    MyFloat Vel[3];
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat Density;
    MyFloat Pressure;
    MyFloat ConditionNumber;
    MyIDType ID;
    int Timestep;
#ifdef HYDRO_SPH
    MyFloat DhsmlHydroSumFactor;
    MyFloat alpha;
#endif
    
    /* matrix of the conserved variable gradients: rho, u, vx, vy, vz */
    struct
    {
        MyDouble Density[3];
        MyDouble Pressure[3];
        MyDouble Velocity[3][3];
#ifdef MAGNETIC
        MyDouble B[3][3];
#ifdef DIVBCLEANING_DEDNER
        MyDouble Phi[3];
#endif
#endif
    } Gradients;
    MyFloat NV_T[3][3];
    
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    MyFloat EgyWtRho;
    MyFloat InternalEnergyPred;
#endif

#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    MyFloat Metallicity[NUM_METAL_SPECIES];
#endif
    
    
#ifdef TURB_DIFFUSION
    MyFloat TD_DiffCoeff;
#endif
    
#ifdef CONDUCTION_EXPLICIT
    MyFloat Kappa_Conduction;
#endif
    
#ifdef BP_REAL_CRs
    MyFloat CRpPressure;
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
    MyFloat CRpE[BP_REAL_CRs];
    MyFloat CRpN[BP_REAL_CRs];
#endif
#endif
    
#ifdef MAGNETIC
#ifdef MAGFORCE
    MyFloat BPred[3];
#endif
#if defined(TRICCO_RESISTIVITY_SWITCH)
    MyFloat Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
    MyFloat PhiPred;
#endif
#endif // MAGNETIC //
    
#ifdef EOS_DEGENERATE
    MyFloat dp_drho;
#endif
    
#ifndef DONOTUSENODELIST
    int NodeList[NODELISTLENGTH];
#endif
}
*HydroDataIn, *HydroDataGet;



/* --------------------------------------------------------------------------------- */
/* outputs: this is what the routine needs to return to the particles to set their final values */
/* --------------------------------------------------------------------------------- */
struct hydrodata_out
{
    MyLongDouble Acc[3];
    //MyLongDouble dMomentum[3]; //???
    MyLongDouble DtInternalEnergy;
    //MyLongDouble dInternalEnergy; //???
    MyFloat MaxSignalVel;
    MyFloat MaxKineticEnergyNgb;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    MyLongDouble DtMass;
    //MyLongDouble dMass; //???
#endif
    
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    MyFloat Dyield[NUM_METAL_SPECIES];
#endif
    
#if defined(MAGNETIC)
    MyFloat DtB[3];
#ifdef DIVBFORCE3
    MyFloat magacc[3];
    MyFloat magcorr[3];
#endif
#if defined(HYDRO_SPH) && defined(DIVBCLEANING_DEDNER)
    MyFloat GradPhi[3];
#endif
#endif // MAGNETIC //
    
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
    MyFloat DtCRpE[BP_REAL_CRs];
    MyFloat DtCRpN[BP_REAL_CRs];
#endif
#if  defined(CR_SHOCK)
    MyFloat CR_EnergyChange[NUMCRPOP];
    MyFloat CR_BaryonFractionChange[NUMCRPOP];
#endif
}
*HydroDataResult, *HydroDataOut;




/* --------------------------------------------------------------------------------- */
/* this subroutine actually loads the particle data into the structure to share between nodes */
/* --------------------------------------------------------------------------------- */
static inline void particle2in_hydra(struct hydrodata_in *in, int i);
static inline void out2particle_hydra(struct hydrodata_out *out, int i, int mode);
static inline void particle2in_hydra(struct hydrodata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = SphP[i].VelPred[k];
    }
    in->Hsml = PPP[i].Hsml;
    in->Mass = P[i].Mass;
    in->Density = SphP[i].Density;
    in->Pressure = SphP[i].Pressure;
    in->ID = P[i].ID;
    in->Timestep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
    in->ConditionNumber = SphP[i].ConditionNumber;
#ifdef HYDRO_SPH
    in->DhsmlHydroSumFactor = SphP[i].DhsmlHydroSumFactor;
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    in->alpha = SphP[i].alpha_limiter * SphP[i].alpha;
#else
    in->alpha = SphP[i].alpha_limiter;
#endif
#endif
    
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    in->EgyWtRho = SphP[i].EgyWtDensity;
    in->InternalEnergyPred = SphP[i].InternalEnergyPred;
#endif
    
    int j;
    for(j=0;j<3;j++)
        for(k=0;k<3;k++)
            in->NV_T[j][k] = SphP[i].NV_T[j][k];
    
    /* matrix of the conserved variable gradients: rho, u, vx, vy, vz */
    for(k=0;k<3;k++)
    {
        in->Gradients.Density[k] = SphP[i].Gradients.Density[k];
        in->Gradients.Pressure[k] = SphP[i].Gradients.Pressure[k];
        in->Gradients.Velocity[0][k] = SphP[i].Gradients.Velocity[0][k];
        in->Gradients.Velocity[1][k] = SphP[i].Gradients.Velocity[1][k];
        in->Gradients.Velocity[2][k] = SphP[i].Gradients.Velocity[2][k];
#ifdef MAGNETIC
        in->Gradients.B[0][k] = SphP[i].Gradients.B[0][k];
        in->Gradients.B[1][k] = SphP[i].Gradients.B[1][k];
        in->Gradients.B[2][k] = SphP[i].Gradients.B[2][k];
#ifdef DIVBCLEANING_DEDNER
        in->Gradients.Phi[k] = SphP[i].Gradients.Phi[k];
#endif
#endif
    }
    

#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    for(k=0;k<NUM_METAL_SPECIES;k++)
        in->Metallicity[k] = P[i].Metallicity[k];
#endif
    
#ifdef TURB_DIFFUSION
    in->TD_DiffCoeff = SphP[i].TD_DiffCoeff;
#endif
    
#ifdef CONDUCTION_EXPLICIT
    in->Kappa_Conduction = SphP[i].Kappa_Conduction;
#endif
    
#ifdef EOS_DEGENERATE
    in->dp_drho = SphP[i].dp_drho;
#endif
    
#ifdef BP_REAL_CRs
    in->CRpPressure = SphP[i].CRpPressure;
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
    int Nbin;
    for( Nbin = 0;Nbin < BP_REAL_CRs; Nbin++ )
    {
        in->CRpE[Nbin] = SphP[i].CRpE[Nbin];
        in->CRpN[Nbin] = SphP[i].CRpN[Nbin];
    }
#endif
#endif
    
#ifdef MAGNETIC
#ifdef MAGFORCE
    for(k = 0; k < 3; k++)
        in->BPred[k] = SphP[i].BPred[k];
#endif
#if defined(TRICCO_RESISTIVITY_SWITCH)
    in->Balpha = SphP[i].Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
    in->PhiPred = SphP[i].PhiPred;
#endif
#endif // MAGNETIC // 
    
}



/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_hydra(struct hydrodata_out *out, int i, int mode)
{
    int k;
    /* these are zero-d out at beginning of hydro loop so should always be added */
    for(k = 0; k < 3; k++)
    {
        SphP[i].HydroAccel[k] += out->Acc[k];
        //SphP[i].dMomentum[k] += out->dMomentum[k]; //???
    }
    SphP[i].DtInternalEnergy += out->DtInternalEnergy;
    //SphP[i].dInternalEnergy += out->dInternalEnergy; //???
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    SphP[i].DtMass += out->DtMass;
    //SphP[i].dMass += out->dMass; //???
#endif
    if(SphP[i].MaxSignalVel < out->MaxSignalVel)
        SphP[i].MaxSignalVel = out->MaxSignalVel;

    if(SphP[i].MaxKineticEnergyNgb < out->MaxKineticEnergyNgb)
        SphP[i].MaxKineticEnergyNgb = out->MaxKineticEnergyNgb;
    
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    for(k=0;k<NUM_METAL_SPECIES;k++)
        P[i].Metallicity[k] += out->Dyield[k] / P[i].Mass;
#endif
    
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
    int Nbin;
    for( Nbin = 0; Nbin < BP_REAL_CRs; Nbin++ )
    {
        ASSIGN_ADD(SphP[i].DtCRpE[Nbin], out->DtCRpE[Nbin], mode);
        ASSIGN_ADD(SphP[i].DtCRpN[Nbin], out->DtCRpN[Nbin], mode);
    }
#endif
    
#if defined(MAGNETIC)
    for(k = 0; k < 3; k++)
    {
        ASSIGN_ADD(SphP[i].DtB[k], out->DtB[k], mode);
#ifdef DIVBFORCE3
        ASSIGN_ADD(SphP[i].magacc[k], out->magacc[k], mode);
        ASSIGN_ADD(SphP[i].magcorr[k], out->magcorr[k], mode);
#endif
#if defined(HYDRO_SPH) && defined(DIVBCLEANING_DEDNER)
        ASSIGN_ADD(SphP[i].Gradients.Phi[k], out->GradPhi[k], mode);
#endif
    }
#endif // MAGNETIC //
}


/* --------------------------------------------------------------------------------- */
/* need to link to the file "hydra_evaluate" which actually contains the computation part of the loop! */
/* --------------------------------------------------------------------------------- */
#include "hydra_evaluate.h"

/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/* This will perform final operations and corrections on the output from the 
    hydro routines, AFTER the neighbors have all been checked and summed */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
void hydro_final_operations_and_cleanup(void)
{
    int i,k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            double dt;
            dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#ifdef CR_SHOCK
            rShockEnergy = 0.0;
            if(SphP[i].DtInternalEnergy > 0.0)
                rShockEnergy = SphP[i].DtInternalEnergy * dt / (All.cf_atime*All.cf_atime * fac_egy);
#endif
            
            /* we calculated the flux of conserved variables: these are used in the kick operation. But for
             intermediate drift operations, we need the primive variables, so reduce to those here 
             (remembering that v_phys = v_code/All.cf_atime, for the sake of doing the unit conversions to physical) */
            for(k=0;k<3;k++)
            {
                SphP[i].DtInternalEnergy -= (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].HydroAccel[k];
                /* we solved for total energy flux (and remember, HydroAccel is still momentum -- keep units straight here!) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[i].DtInternalEnergy += 0.5 * (SphP[i].VelPred[k]/All.cf_atime) * (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].DtMass;
                SphP[i].HydroAccel[k] -= (SphP[i].VelPred[k]/All.cf_atime) * SphP[i].DtMass; /* we solved for momentum flux */
#endif
                SphP[i].HydroAccel[k] /= P[i].Mass; /* we solved for momentum flux */
            }
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            SphP[i].DtInternalEnergy -= SphP[i].InternalEnergyPred * SphP[i].DtMass;
#endif
            SphP[i].DtInternalEnergy /= P[i].Mass;
            /* ok, now: HydroAccel = dv/dt, DtInternalEnergy = du/dt (energy per unit mass) */
            
            // need to explicitly include adiabatic correction from the hubble-flow (for drifting) here //
            if(All.ComovingIntegrationOn) SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a; //???
            // = du/dlna -3*(gamma-1)*u ; then dlna/dt = H(z) =  All.cf_hubble_a //

#ifdef EOS_DEGENERATE
            /* DtInternalEnergy stores the energy change rate in internal units */
            SphP[i].DtInternalEnergy *= All.UnitEnergy_in_cgs / All.UnitTime_in_s;
#endif
            
            
            
#ifdef MACHNUM
            /* Estimates the Mach number of particle i for non-radiative runs, or the Mach number,
             density jump and specific energy jump in case of cosmic rays! */
#if (CR_SHOCK == 2)
            GetMachNumberCR(SphP + i);
#else
            GetMachNumber(SphP + i);
#endif // COSMIC_RAYS //
#ifdef MACHSTATISTIC
            GetShock_DtEnergy(SphP + i);
#endif
#endif // MACHNUM //
#ifdef CR_SHOCK
            if(rShockEnergy > 0.0)
            {
                /* Feed fraction "All.CR_ShockEfficiency" into CR and see what amount of energy instantly gets rethermalized
                 * for this, we need the physical time step, which is Delta t_p = Delta t_c / All.cf_hubble_a */
                rNonRethermalizedEnergy = CR_Particle_ShockInject(SphP + i, rShockEnergy, dt);
                /* Fraction of total energy that went and remained in CR is rNonRethermalizedEnergy / rShockEnergy, hence, we conserve energy if we do */
                SphP[i].DtInternalEnergy *= (1.0 - rNonRethermalizedEnergy / rShockEnergy);
                assert(rNonRethermalizedEnergy >= 0.0);
                assert(rNonRethermalizedEnergy <= (rShockEnergy * All.CR_ShockEfficiency));
            }
#endif // CR_SHOCK //
#if (!defined(COOLING) && !defined(CR_SHOCK) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION))) || (defined(CR_SHOCK) &&  (!defined(COOLING) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION))))
            if((P[i].TimeBin)&&(dt>0))	/* upon start-up, we need to protect against dt==0 */
            {
                for(int CRpop = 0; CRpop < NUMCRPOP; CRpop++)
                {
                    utherm = 0.0;
                    utherm = CR_Particle_ThermalizeAndDissipate(SphP + i, dt, CRpop);
                    SphP[i].DtInternalEnergy += utherm * fac_egy / (dt * All.cf_hubble_a);
                }
            }
#endif
            
                        
            
#if defined(HYDRO_SPH) && defined(MAGNETIC)
            double tmpb,phiphi;
#ifdef DIVBFORCE3
            // PFH: check if magnitude of correction > force itself; in which case limit it to equal value //
            phiphi = SphP[i].magcorr[0]*SphP[i].magcorr[0] + SphP[i].magcorr[1]*SphP[i].magcorr[1] + SphP[i].magcorr[2]*SphP[i].magcorr[2];
            tmpb = SphP[i].magacc[0]*SphP[i].magacc[0] + SphP[i].magacc[1]*SphP[i].magacc[1] + SphP[i].magacc[2]*SphP[i].magacc[2];
            if(phiphi > DIVBFORCE3 * tmpb)
            {
                tmpb = sqrt(DIVBFORCE3 * tmpb / phiphi); // save sqrt for here (speeds this up)
                for(k = 0; k < 3; k++)
                    SphP[i].magcorr[k] *= tmpb;
            }
            // add the corrected mhd acceleration to the hydro acceleration //
            for(k = 0; k < 3; k++)
                SphP[i].HydroAccel[k] += (SphP[i].magacc[k] - SphP[i].magcorr[k]);
#endif
#ifdef DIVBCLEANING_DEDNER
            /* full correct form of D(phi)/Dt = -ch*ch*div.dot.B - phi/tau - (1/2)*phi*div.dot.v */
            /* PFH: here's the div.dot.B term: make sure div.dot.B def'n matches appropriate grad_phi conjugate pair: recommend direct diff div.dot.B */
            tmpb = 0.5 * SphP[i].MaxSignalVel / fac_mu; // has units of v_code now
            phiphi = tmpb * tmpb * All.DivBcleanHyperbolicSigma * SphP[i].divB;
            // phiphi above now has units of [Bcode]*[vcode]^2/[rcode]=(Bcode*vcode)*vcode/rcode; needs to have units of [Phicode]*[vcode]/[rcode]
            // [GradPhi]=[Phicode]/[rcode] = [DtB] = [Bcode]*[vcode]/[rcode] IFF [Phicode]=[Bcode]*[vcode]; this also makes the above self-consistent //
            // (implicitly, this gives the correct evolution in comoving, adiabatic coordinates where the sound speed is the relevant speed at which
            //   the 'damping wave' propagates. another choice (provided everything else is self-consistent) is fine, it just makes different assumptions
            //   about the relevant 'desired' timescale for damping wave propagation in the expanding box) //
            
            /* PFH: add simple damping (-phi/tau) term */
            phiphi += SphP[i].PhiPred * 0.5 * SphP[i].MaxSignalVel / (KERNEL_CORE_SIZE*PPP[i].Hsml*fac_mu) *
                All.DivBcleanParabolicSigma; // need to be sure this translates into vcode/rcode units
            
            /* PFH: add div_v term from Tricco & Price to DtPhi */
            phiphi += SphP[i].PhiPred * 0.5 * P[i].Particle_DivVel*All.cf_a2inv;
            
            /* multiply by negative sign and cosmological correction term */
            //SphP[i].DtPhi = - phiphi * All.cf_atime * All.cf_atime;	/* Compensate for the 1/Ha^2 in dt_mag */
            SphP[i].DtPhi = -phiphi; // should be in units of [Phi]*[v_code]/[r_code]; from that point dt_mag will correctly integrate it //
            phiphi = sqrt(SphP[i].Gradients.Phi[0]*SphP[i].Gradients.Phi[0] + SphP[i].Gradients.Phi[1]*SphP[i].Gradients.Phi[1] + SphP[i].Gradients.Phi[2]*SphP[i].Gradients.Phi[2]);
            tmpb = sqrt(SphP[i].DtB[0]*SphP[i].DtB[0] + SphP[i].DtB[1]*SphP[i].DtB[1] + SphP[i].DtB[2]*SphP[i].DtB[2]);
            if(phiphi > All.DivBcleanQ * tmpb && tmpb != 0)
                for(k = 0; k < 3; k++)
                    SphP[i].Gradients.Phi[k] *=  tmpb * All.DivBcleanQ / phiphi;
#ifdef MAGNETIC_DISSIPATION
            SphP[i].DtInternalEnergy -= (SphP[i].BPred[0] * SphP[i].Gradients.Phi[0] + SphP[i].BPred[1] * SphP[i].Gradients.Phi[1] +
                                  SphP[i].BPred[2] * SphP[i].Gradients.Phi[2]) * fac_magnetic_pressure / SphP[i].Density;
#endif
            SphP[i].DtB[0] += SphP[i].Gradients.Phi[0];
            SphP[i].DtB[1] += SphP[i].Gradients.Phi[1];
            SphP[i].DtB[2] += SphP[i].Gradients.Phi[2];
#endif /* End DEDNER */
#endif /* End Magnetic */
            
            
#ifdef RT_RAD_PRESSURE
            /* radiative accelerations */
            if(All.Time != All.TimeBegin)
                for(k = 0; k < 3; k++)
                {
                    SphP[i].RadAccel[k] = 0.0;
                    for(j = 0; j < N_BINS; j++)
                        SphP[i].RadAccel[k] += SphP[i].n_gamma[j] * nu[j];
                    SphP[i].RadAccel[k] *= SphP[i].n[k] / P[i].Mass * ELECTRONVOLT_IN_ERGS /
                    All.UnitEnergy_in_cgs * All.HubbleParam / (C / All.UnitVelocity_in_cm_per_s) / dt / SphP[i].Density;
                }
#endif
            
            
#ifdef GALSF_SUBGRID_WINDS
            /* if we have winds, we decouple particles briefly if delaytime>0 */
            if(SphP[i].DelayTime > 0)
            {
                for(k = 0; k < 3; k++)
                    SphP[i].HydroAccel[k] = 0;//SphP[i].dMomentum[k] = 0;
                SphP[i].DtInternalEnergy = 0; //SphP[i].dInternalEnergy = 0;
                double windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
                windspeed *= fac_mu;
                double hsml_c = pow(All.WindFreeTravelDensFac * All.PhysDensThresh / (SphP[i].Density * All.cf_a3inv), (1. / 3.));
                SphP[i].MaxSignalVel = hsml_c * DMAX((2 * windspeed), SphP[i].MaxSignalVel);
            }
#endif
            
            
#ifdef BND_PARTICLES
            /* this flag signals all particles with id=0 are frozen (boundary particles) */
            if(P[i].ID == 0)
            {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[i].DtMass = 0;//SphP[i].dMass = 0;//???
#endif
                SphP[i].DtInternalEnergy = 0;//SphP[i].dInternalEnergy = 0;//???
                for(k = 0; k < 3; k++) SphP[i].HydroAccel[k] = 0;//SphP[i].dMomentum[k] = 0;//???
#ifdef SPH_BND_DTB
                for(k = 0; k < 3; k++) SphP[i].DtB[k] = 0;
#endif
#ifdef SPH_BND_BFLD
                for(k = 0; k < 3; k++) SphP[i].B[k] = 0;
#endif
            }
#endif
        
        } // closes P[i].Type==0 check and so closes loop over particles i
    } // for (loop over active particles) //
    
    
    
#ifdef NUCLEAR_NETWORK
    if(ThisTask == 0)
    {
        printf("Doing nuclear network.\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    tstart = my_second();
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type == 0)
        {
            /* evaluate network here, but do it only for high enough temperatures */
            if(SphP[i].temp > All.NetworkTempThreshold)
            {
                nuc_particles++;
                network_integrate(SphP[i].temp, SphP[i].Density * All.UnitDensity_in_cgs, SphP[i].xnuc,
                                  SphP[i].dxnuc, dt*All.UnitTime_in_s, &dedt_nuc, NULL, &All.nd, &All.nw);
                SphP[i].DtInternalEnergy += dedt_nuc * All.UnitEnergy_in_cgs / All.UnitTime_in_s;
            }
            else
            {
                for(k = 0; k < EOS_NSPECIES; k++)
                {
                    SphP[i].dxnuc[k] = 0;
                }
            }
        }
    tend = my_second();
    timenetwork += timediff(tstart, tend);
    MPI_Allreduce(&nuc_particles, &nuc_particles_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("Nuclear network done for %d particles.\n", nuc_particles_sum);
    }
    timewait1 += timediff(tend, my_second());
#endif
}




/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the driver routine for the calculation of hydrodynamical
 *  force, fluxes, etc. */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
void hydro_force(void)
{
    int i, j, k, ngrp, ndone, ndone_flag;
    int recvTask, place;
    double timeall=0, timecomp1=0, timecomp2=0, timecommsumm1=0, timecommsumm2=0, timewait1=0, timewait2=0, timenetwork=0;
    double timecomp, timecomm, timewait, tstart, tend, t0, t1;
    int save_NextParticle;
    long long n_exported = 0;
    /* need to zero out all numbers that can be set -EITHER- by an active particle in the domain, 
     or by one of the neighbors we will get sent */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type==0)
        {
            SphP[i].MaxSignalVel = -1.e10;
            SphP[i].MaxKineticEnergyNgb = -1.e10;
            SphP[i].DtInternalEnergy = 0;//SphP[i].dInternalEnergy = 0;//???
            for(k=0;k<3;k++)
                SphP[i].HydroAccel[k] = 0;//SphP[i].dMomentum[k] = 0;//???
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            SphP[i].DtMass = 0;//SphP[i].dMass = 0;//???
#endif
#ifdef WAKEUP
            SphP[i].wakeup = 0;
#endif
        }
    
    /* --------------------------------------------------------------------------------- */
    // Global factors for comoving integration of hydro //
    fac_mu = 1 / (All.cf_afac3 * All.cf_atime);
    // code_vel * fac_mu = sqrt[code_pressure/code_density] = code_soundspeed //
    // note also that signal_vel in forms below should be in units of code_soundspeed //
    fac_vsic_fix = All.cf_hubble_a * All.cf_afac1;
#ifdef MAGNETIC
    fac_magnetic_pressure = MU0_1 * All.cf_afac1 / All.cf_atime;
    // code_Bfield*code_Bfield * fac_magnetic_pressure = code_pressure //
    // -- use this to get alfven velocities, etc, as well as comoving units for magnetic integration //
#endif
#ifdef MACHNUM
    fac_egy = All.cf_afac1;
#endif
    
    /* --------------------------------------------------------------------------------- */
    /* allocate buffers to arrange communication */
    int NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct hydrodata_in) +
                                                             sizeof(struct hydrodata_out) +
                                                             sizemax(sizeof(struct hydrodata_in),
                                                                     sizeof(struct hydrodata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    CPU_Step[CPU_HYDMISC] += measure_time();
    t0 = my_second();
    NextParticle = FirstActiveParticle;	/* beginn with this index */
    
    do
    {
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
        tstart = my_second();
        
#ifdef OMP_NUM_THREADS
        pthread_t mythreads[OMP_NUM_THREADS - 1];
        int threadid[OMP_NUM_THREADS - 1];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        pthread_mutex_init(&mutex_nexport, NULL);
        pthread_mutex_init(&mutex_partnodedrift, NULL);
        TimerFlag = 0;
        for(j = 0; j < OMP_NUM_THREADS - 1; j++)
        {
            threadid[j] = j + 1;
            pthread_create(&mythreads[j], &attr, hydro_evaluate_primary, &threadid[j]);
        }
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            hydro_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
        
#ifdef OMP_NUM_THREADS
        for(j = 0; j < OMP_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
#endif
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            NextParticle = save_NextParticle;
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle)
                    break;
                
                if(ProcessedFlag[NextParticle] != 1)
                    break;
                
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle)
            {
                /* in this case, the buffer is too small to process even a single particle */
                endrun(115508);
            }
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1)
                        k = j + 1;
                    
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            
                            DataIndexTable[j] = DataIndexTable[k];
                            DataNodeList[j] = DataNodeList[k];
                            DataIndexTable[j].IndexGet = j;
                            new_export++;
                            
                            DataIndexTable[k].Index = old_index;
                            k++;
                            break;
                        }
                }
                else
                    new_export++;
            
            Nexport = new_export;
        }
        
        n_exported += Nexport;
        for(j = 0; j < NTask; j++)
            Send_count[j] = 0;
        for(j = 0; j < Nexport; j++)
            Send_count[DataIndexTable[j].Task]++;
        
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        tstart = my_second();
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        tend = my_second();
        timewait1 += timediff(tstart, tend);
        
        for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            Nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        HydroDataGet = (struct hydrodata_in *) mymalloc("HydroDataGet", Nimport * sizeof(struct hydrodata_in));
        HydroDataIn = (struct hydrodata_in *) mymalloc("HydroDataIn", Nexport * sizeof(struct hydrodata_in));
        
        /* prepare particle data for export */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_hydra(&HydroDataIn[j], place);		// MADE D_IND CHANGE IN HERE
#ifndef DONOTUSENODELIST
            memcpy(HydroDataIn[j].NodeList,
                   DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif
            
        }
        
        /* exchange particle data */
        tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
                                 recvTask, TAG_HYDRO_A,
                                 &HydroDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
                                 recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm1 += timediff(tstart, tend);
        
        myfree(HydroDataIn);
        HydroDataResult = (struct hydrodata_out *) mymalloc("HydroDataResult", Nimport * sizeof(struct hydrodata_out));
        HydroDataOut = (struct hydrodata_out *) mymalloc("HydroDataOut", Nexport * sizeof(struct hydrodata_out));
        report_memory_usage(&HighMark_sphhydro, "SPH_HYDRO");
        
        /* now do the particles that were sent to us */
        tstart = my_second();
        NextJ = 0;
#ifdef OMP_NUM_THREADS
        for(j = 0; j < OMP_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, hydro_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            hydro_evaluate_secondary(&mainthreadid);
        }
        
#ifdef OMP_NUM_THREADS
        for(j = 0; j < OMP_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
        
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        tend = my_second();
        timecomp2 += timediff(tstart, tend);
        
        if(NextParticle < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        tstart = my_second();
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        tend = my_second();
        timewait2 += timediff(tstart, tend);
        
        /* get the result */
        tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct hydrodata_out),
                                 MPI_BYTE, recvTask, TAG_HYDRO_B,
                                 &HydroDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct hydrodata_out),
                                 MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm2 += timediff(tstart, tend);
        
        /* add the result to the local particles */
        tstart = my_second();
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_hydra(&HydroDataOut[j], place, 1);
        }
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        
        myfree(HydroDataOut);
        myfree(HydroDataResult);
        myfree(HydroDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    
    /* --------------------------------------------------------------------------------- */
    /* do final operations on results */
    /* --------------------------------------------------------------------------------- */
    hydro_final_operations_and_cleanup();
    
    
    /* --------------------------------------------------------------------------------- */
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_HYDCOMPUTE] += timecomp;
    CPU_Step[CPU_HYDWAIT] += timewait;
    CPU_Step[CPU_HYDCOMM] += timecomm;
    CPU_Step[CPU_HYDNETWORK] += timenetwork;
    CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}



/* --------------------------------------------------------------------------------- */
/* one of the core sub-routines used to do the MPI version of the hydro evaluation
 (don't put actual operations here!!!) */
/* --------------------------------------------------------------------------------- */
void *hydro_evaluate_primary(void *p)
{
    int thread_id = *(int *) p;
    int i, j;
    int *exportflag, *exportnodecount, *exportindex, *ngblist;
    
    ngblist = Ngblist + thread_id * NumPart;
    exportflag = Exportflag + thread_id * NTask;
    exportnodecount = Exportnodecount + thread_id * NTask;
    exportindex = Exportindex + thread_id * NTask;
    /* Note: exportflag is local to each thread */
    for(j = 0; j < NTask; j++)
        exportflag[j] = -1;
    
    while(1)
    {
        int exitFlag = 0;
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
            if(BufferFullFlag != 0 || NextParticle < 0)
            {
                exitFlag = 1;
            }
            else
            {
                i = NextParticle;
                ProcessedFlag[i] = 0;
                NextParticle = NextActiveParticle[NextParticle];
            }
        }
        UNLOCK_NEXPORT;
        if(exitFlag)
            break;
        
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            if(hydro_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
                break;		/* export buffer has filled up */
        }
        ProcessedFlag[i] = 1;	/* particle successfully finished */
    }
    return NULL;
}



/* --------------------------------------------------------------------------------- */
/* one of the core sub-routines used to do the MPI version of the hydro evaluation
 (don't put actual operations here!!!) */
/* --------------------------------------------------------------------------------- */
void *hydro_evaluate_secondary(void *p)
{
    int thread_id = *(int *) p;
    int j, dummy, *ngblist;
    ngblist = Ngblist + thread_id * NumPart;
    while(1)
    {
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
            j = NextJ;
            NextJ++;
        }
        UNLOCK_NEXPORT;
        
        if(j >= Nimport)
            break;
        
        hydro_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }
    return NULL;
}

