#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#define NDEBUG

/*! \file hydra_master.c
 *  \brief This contains the "second hydro loop", where the hydro fluxes are computed.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* some very useful notes on the hydro variables in comoving integrations:
 
 v_code = a * v_peculiar/physical (canonical momentum)
 r_code = r_physical / a (comoving coordinates)
 m_code = m_physical
 rho_code = rho_physical * a^3 (from length/mass scaling)
 InternalEnergy_code = InternalEnergy_physical
 Pressure_code =
    InternalEnergy_code * rho_code * (gamma-1) = Pressure_physical * a^3 (energy SPH)
    -- the distinction between these cases and e.g. entropy sph (now depricated) 
        should be taken care of in the factors
        All.cf_afac1/2/3, which will correctly assign between the two --
 B_code = a*a * B_physical (comoving magnetic fields)
 Phi_code = B_code*v_code = a^3 * Phi_physical (damping field for Dedner divergence cleaning)
    (note: spec egy of phi field is: phi*phi/(2*mu0*rho*ch*ch); compare Bfield is B*B/(mu0*rho);
    so [phi]~[B]*[ch], where ch is the signal velocity used in the damping equation);

 -- Time derivatives (rate of change from hydro forces) here are all
        assumed to end up in *physical* units ---
 HydroAccel, dMomentum are assumed to end up in *physical* units
    (note, this is different from GADGET's convention, where 
     HydroAccel is in units of (Pcode/rhocode)/rcode)
 DtInternalEnergy and dInternalEnergy are assumed to end up in *physical* units
 DtMass and dMass are assumed to end up in *physical* units

 -----------------------------------------
 
 // All.cf_atime = a = 1/(1+z), the cosmological scale factor //
 All.cf_atime = All.Time;
 // All.cf_a2inv is just handy //
 All.cf_a2inv = 1 / (All.Time * All.Time);
 // All.cf_a3inv * Density_code = Density_physical //
 All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
 // Pressure_code/Density_code = All.cf_afac1 * Pressure_physical/Density_physical //
 All.cf_afac1 = 1;
 // All.cf_afac2 * Pressure_code/Density_code * 1/r_code = Pressure_physical/Density_physical * 1/r_physical //
 All.cf_afac2 = 1 / (All.Time * All.cf_afac1);
 // All.cf_afac3 * cs_code = All.cf_afac3 * sqrt(Pressure_code/Density_code) = sqrt(Pressure_phys/Density_phys) = cs_physical //
 All.cf_afac3 = 1 / sqrt(All.cf_afac1);
 // time units: proper time dt_phys = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
 code time unit in comoving is dlna, so dt_phys = dt_code / All.cf_hubble_a   //
 All.cf_hubble_a = hubble_function(All.Time); // hubble_function(a) = H(a) = H(z) //
 // dt_code * v_code/r_code = All.cf_hubble_a2 * dt_phys * v_phys/r_phys //
 All.cf_hubble_a2 = All.Time * All.Time * hubble_function(All.Time);
 
 
 -----------------------------------------
 A REMINDER ABOUT GIZMO/GADGET VELOCITY UNITS:: (direct quote from Volker)
 
 The IC file should contain the *peculiar* velocity divided by sqrt(a),
 not the *physical* velocity. Let "x" denote comoving
 coordinates and "r=a*x" physical coordinates. Then I call
 
 comoving velocity: dx/dt
 physical velocity: dr/dt = H(a)*r + a*dx/dt
 peculiar velocity: v = a * dx/dt
 
 The physical velocity is hence the peculiar velocity plus the Hubble flow.
 
 The internal velocity variable is not given by dx/d(ln a). Rather, it is given by
 the canonical momentum p = a^2 * dx/dt.
 The IC-file and snapshot files of gadget/GIZMO don't
 contain the variable "p" directly because of historical reasons.
 Instead, they contain the velocity variable
 u = v/sqrt(a) = sqrt(a) * dx/dt = p / a^(3/2), which is just what the
 manual says. (The conversion between u and p is done on the fly when
 reading or writing snapshot files.)
 
 Also note that d(ln a)/dt is equal to the
 Hubble rate, i.e.: d(ln a)/dt = H(a) = H_0 * sqrt(omega_m/a^3 + omega_v
 + (1 - omega_m - omega_v)/a^2).
 
 Best wishes,
 Volker
 
 -----------------------------------------
*/


static double fac_mu, fac_vsic_fix;
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
    MyDouble u;
    MyDouble cs;
#ifdef MAGNETIC
    MyDouble B[3];
    MyDouble B_normal_corrected;
#ifdef DIVBCLEANING_DEDNER
    MyDouble phi;
#endif
#endif
#ifdef COSMIC_RAYS
    MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS];
#ifdef COSMIC_RAYS_M1
    MyDouble CosmicRayFlux[N_CR_PARTICLE_BINS][3];
#endif
#ifdef COSMIC_RAYS_ALFVEN
    MyDouble CosmicRayAlfvenEnergy[N_CR_PARTICLE_BINS][2];
#endif
#endif
};


struct kernel_hydra
{
    double dp[3];
    double r, vsig, sound_i, sound_j;
    double dv[3], vdotr2;
    double wk_i, wk_j, dwk_i, dwk_j;
    double h_i, h_j, dwk_ij, rho_ij_inv;
    double spec_egy_u_i;
#ifdef HYDRO_SPH
    double p_over_rho2_i;
#endif
#ifdef MAGNETIC
    double b2_i, b2_j;
    double alfven2_i, alfven2_j;
#ifdef HYDRO_SPH
    double mf_i, mf_j;
#endif
#endif // MAGNETIC //
};
#ifndef HYDRO_SPH
#include "reimann.h"
#endif


/* ok here we define some important variables for our generic communication
    and flux-exchange structures. these can be changed, and vary across the code, but need to be set! */

#define MASTER_FUNCTION_NAME hydro_force_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int MASTER_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_hydra    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_hydra  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if((P[i].Type==0)&&(P[i].Mass>0)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */




/* --------------------------------------------------------------------------------- */
/* inputs to the routine: put here what's needed to do the calculation! */
/* --------------------------------------------------------------------------------- */
struct INPUT_STRUCT_NAME
{
    /* basic hydro variables */
    MyDouble Pos[3];
    MyFloat Vel[3];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    MyFloat ParticleVel[3];
#endif
    MyFloat Hsml;
    MyFloat Mass;
    MyFloat Density;
    MyFloat Pressure;
    MyFloat ConditionNumber;
    MyFloat InternalEnergyPred;
    MyFloat SoundSpeed;
    integertime Timestep;
    MyFloat DhsmlNgbFactor;
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
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_M1)
        MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS][3];
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        MyDouble Metallicity[NUM_METAL_SPECIES][3];
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        MyDouble InternalEnergy[3];
#endif
#ifdef DOGRAD_SOUNDSPEED
        MyDouble SoundSpeed[3];
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_COMPGRAD_EDDINGTON_TENSOR)
        MyDouble Rad_E_gamma_ET[N_RT_FREQ_BINS][3];
#endif
    } Gradients;
    MyFloat NV_T[3][3];
    
#if defined(KERNEL_CRK_FACES)
    MyFloat Tensor_CRK_Face_Corrections[16];
#endif
#ifdef HYDRO_PRESSURE_SPH
    MyFloat EgyWtRho;
#endif

#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    MyFloat Metallicity[NUM_METAL_SPECIES];
#endif

#ifdef CHIMES_TURB_DIFF_IONS 
    MyDouble ChimesNIons[CHIMES_TOTSIZE]; 
#endif 
    
#ifdef RT_SOLVER_EXPLICIT
    MyDouble Rad_E_gamma[N_RT_FREQ_BINS];
    MyDouble Rad_Kappa[N_RT_FREQ_BINS];
    MyDouble RT_DiffusionCoeff[N_RT_FREQ_BINS];
#if defined(RT_EVOLVE_FLUX) || defined(HYDRO_SPH)
    MyDouble ET[N_RT_FREQ_BINS][6];
#endif
#ifdef RT_EVOLVE_FLUX
    MyDouble Rad_Flux[N_RT_FREQ_BINS][3];
#endif
#ifdef RT_INFRARED
    MyDouble Radiation_Temperature;
#endif
#if defined(RT_EVOLVE_INTENSITIES)
    MyDouble Rad_Intensity_Pred[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS];
#endif
#endif
    
#ifdef TURB_DIFFUSION
    MyFloat TD_DiffCoeff;
#endif
    
#ifdef CONDUCTION
    MyFloat Kappa_Conduction;
#endif

#ifdef MHD_NON_IDEAL
    MyFloat Eta_MHD_OhmicResistivity_Coeff;
    MyFloat Eta_MHD_HallEffect_Coeff;
    MyFloat Eta_MHD_AmbiPolarDiffusion_Coeff;
#endif
    
#ifdef VISCOSITY
    MyFloat Eta_ShearViscosity;
    MyFloat Zeta_BulkViscosity;
#endif
    
#ifdef MAGNETIC
    MyFloat BPred[3];
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
    MyFloat Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
    MyFloat PhiPred;
#endif
#endif // MAGNETIC //

#ifdef COSMIC_RAYS
    MyDouble CosmicRayPressure[N_CR_PARTICLE_BINS];
    MyDouble CosmicRayDiffusionCoeff[N_CR_PARTICLE_BINS];
#ifdef COSMIC_RAYS_M1
    MyDouble CosmicRayFlux[N_CR_PARTICLE_BINS][3];
#endif
#ifdef COSMIC_RAYS_ALFVEN
    MyDouble CosmicRayAlfvenEnergy[N_CR_PARTICLE_BINS][2];
#endif
#endif
    
#ifdef GALSF_SUBGRID_WINDS
    MyDouble DelayTime;
#endif
    
#ifdef EOS_ELASTIC
    int CompositionType;
    MyFloat Elastic_Stress_Tensor[3][3];
#endif
    
#ifndef DONOTUSENODELIST
    int NodeList[NODELISTLENGTH];
#endif
}
*DATAIN_NAME, *DATAGET_NAME;



/* --------------------------------------------------------------------------------- */
/* outputs: this is what the routine needs to return to the particles to set their final values */
/* --------------------------------------------------------------------------------- */
struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Acc[3];
    //MyLongDouble dMomentum[3]; //manifest-indiv-timestep-debug//
    MyLongDouble DtInternalEnergy;
    //MyLongDouble dInternalEnergy; //manifest-indiv-timestep-debug//
    MyFloat MaxSignalVel;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    MyFloat MaxKineticEnergyNgb;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    MyLongDouble DtMass;
    MyLongDouble dMass;
    MyLongDouble GravWorkTerm[3];
#endif
    
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    MyFloat Dyield[NUM_METAL_SPECIES];
#endif

#ifdef CHIMES_TURB_DIFF_IONS 
    MyDouble ChimesIonsYield[CHIMES_TOTSIZE]; 
#endif 
    
#if defined(RT_SOLVER_EXPLICIT)
#if defined(RT_EVOLVE_ENERGY)
    MyFloat Dt_Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#if defined(RT_EVOLVE_FLUX)
    MyFloat Dt_Rad_Flux[N_RT_FREQ_BINS][3];
#endif
#if defined(RT_INFRARED)
    MyFloat Dt_Rad_E_gamma_T_weighted_IR;
#endif
#if defined(RT_EVOLVE_INTENSITIES)
    MyFloat Dt_Rad_Intensity[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS];
#endif
#endif
    
#if defined(MAGNETIC)
    MyDouble Face_Area[3];
    MyFloat DtB[3];
    MyFloat divB;
#if defined(DIVBCLEANING_DEDNER)
#ifdef HYDRO_MESHLESS_FINITE_VOLUME // mass-based phi-flux
    MyFloat DtPhi;
#endif
    MyFloat DtB_PhiCorr[3];
#endif
#endif // MAGNETIC //
    
#ifdef COSMIC_RAYS
    MyDouble DtCosmicRayEnergy[N_CR_PARTICLE_BINS];
#ifdef COSMIC_RAYS_ALFVEN
    MyDouble DtCosmicRayAlfvenEnergy[N_CR_PARTICLE_BINS][2];
#endif
#endif

}
*DATARESULT_NAME, *DATAOUT_NAME;




/* --------------------------------------------------------------------------------- */
/* this subroutine actually loads the particle data into the structure to share between nodes */
/* --------------------------------------------------------------------------------- */
static inline void particle2in_hydra(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration);
static inline void particle2in_hydra(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = SphP[i].VelPred[k];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        in->ParticleVel[k] = SphP[i].ParticleVel[k];
#endif
    }
    in->Hsml = PPP[i].Hsml;
    in->Mass = P[i].Mass;
    in->Density = SphP[i].Density;
    in->Pressure = SphP[i].Pressure;
    in->InternalEnergyPred = SphP[i].InternalEnergyPred;
    in->SoundSpeed = Particle_effective_soundspeed_i(i);
    in->Timestep = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0);
    in->ConditionNumber = SphP[i].ConditionNumber;
#ifdef MHD_CONSTRAINED_GRADIENT
    /* since it is not used elsewhere, we can use the sign of the condition number as a bit 
        to conveniently indicate the status of the parent particle flag, for the constrained gradients */
    if(SphP[i].FlagForConstrainedGradients == 0) {in->ConditionNumber *= -1;}
#endif
#ifdef BH_WIND_SPAWN
    if(P[i].ID == All.AGNWindID) {in->ConditionNumber *= -1;} /* as above, use sign of condition number as a bitflag to indicate if this is, or is not, a wind particle */
#endif
    in->DhsmlNgbFactor = PPP[i].DhsmlNgbFactor;
#ifdef HYDRO_SPH
    in->DhsmlHydroSumFactor = SphP[i].DhsmlHydroSumFactor;
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    in->alpha = SphP[i].alpha_limiter * SphP[i].alpha;
#else
    in->alpha = SphP[i].alpha_limiter;
#endif
#endif
    
#ifdef HYDRO_PRESSURE_SPH
    in->EgyWtRho = SphP[i].EgyWtDensity;
#endif
#if defined(KERNEL_CRK_FACES)
    for(k=0;k<16;k++) {in->Tensor_CRK_Face_Corrections[k] = SphP[i].Tensor_CRK_Face_Corrections[k];}
#endif

    int j;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {in->NV_T[j][k] = SphP[i].NV_T[j][k];}}

    
    /* matrix of the conserved variable gradients: rho, u, vx, vy, vz */
    for(k=0;k<3;k++)
    {
        in->Gradients.Density[k] = SphP[i].Gradients.Density[k];
        in->Gradients.Pressure[k] = SphP[i].Gradients.Pressure[k];
        for(j=0;j<3;j++) {in->Gradients.Velocity[j][k] = SphP[i].Gradients.Velocity[j][k];}
#ifdef MAGNETIC
        for(j=0;j<3;j++) {in->Gradients.B[j][k] = SphP[i].Gradients.B[j][k];}
#ifdef DIVBCLEANING_DEDNER
        in->Gradients.Phi[k] = SphP[i].Gradients.Phi[k];
#endif
#endif
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_M1)
        for(j=0;j<N_CR_PARTICLE_BINS;j++) {in->Gradients.CosmicRayPressure[j][k] = SphP[i].Gradients.CosmicRayPressure[j][k];}
#endif
#if defined(TURB_DIFF_METALS) && !defined(TURB_DIFF_METALS_LOWORDER)
        for(j=0;j<NUM_METAL_SPECIES;j++) {in->Gradients.Metallicity[j][k] = SphP[i].Gradients.Metallicity[j][k];}
#endif
#ifdef DOGRAD_INTERNAL_ENERGY
        in->Gradients.InternalEnergy[k] = SphP[i].Gradients.InternalEnergy[k];
#endif
#ifdef DOGRAD_SOUNDSPEED
        in->Gradients.SoundSpeed[k] = SphP[i].Gradients.SoundSpeed[k];
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_COMPGRAD_EDDINGTON_TENSOR)
        for(j=0;j<N_RT_FREQ_BINS;j++) {in->Gradients.Rad_E_gamma_ET[j][k] = SphP[i].Gradients.Rad_E_gamma_ET[j][k];}
#endif
    }

#ifdef RT_SOLVER_EXPLICIT
    for(k=0;k<N_RT_FREQ_BINS;k++)
    {
#ifdef RT_EVOLVE_ENERGY
        in->Rad_E_gamma[k] = SphP[i].Rad_E_gamma_Pred[k];
#else
        in->Rad_E_gamma[k] = SphP[i].Rad_E_gamma[k];
#endif
        in->Rad_Kappa[k] = SphP[i].Rad_Kappa[k];
        in->RT_DiffusionCoeff[k] = rt_diffusion_coefficient(i,k);
#if defined(RT_EVOLVE_FLUX) || defined(HYDRO_SPH)
        {int k_dir; for(k_dir=0;k_dir<6;k_dir++) in->ET[k][k_dir] = SphP[i].ET[k][k_dir];}
#endif
#ifdef RT_EVOLVE_FLUX
        {int k_dir; for(k_dir=0;k_dir<3;k_dir++) in->Rad_Flux[k][k_dir] = SphP[i].Rad_Flux_Pred[k][k_dir];}
#endif
#if defined(RT_EVOLVE_INTENSITIES)
        {int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {in->Rad_Intensity_Pred[k][k_dir] = SphP[i].Rad_Intensity_Pred[k][k_dir];}}
#endif
    }
#ifdef RT_INFRARED
        in->Radiation_Temperature = SphP[i].Radiation_Temperature;
#endif
#endif

#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->Metallicity[k] = P[i].Metallicity[k];}
#endif

#ifdef CHIMES_TURB_DIFF_IONS  
    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) {in->ChimesNIons[k] = SphP[i].ChimesNIons[k]; }
#endif 
    
#ifdef TURB_DIFFUSION
    in->TD_DiffCoeff = SphP[i].TD_DiffCoeff;
#endif
    
#ifdef CONDUCTION
    in->Kappa_Conduction = SphP[i].Kappa_Conduction;
#endif
    
#ifdef MHD_NON_IDEAL
    in->Eta_MHD_OhmicResistivity_Coeff = SphP[i].Eta_MHD_OhmicResistivity_Coeff;
    in->Eta_MHD_HallEffect_Coeff = SphP[i].Eta_MHD_HallEffect_Coeff;
    in->Eta_MHD_AmbiPolarDiffusion_Coeff = SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff;
#endif
    

#ifdef VISCOSITY
    in->Eta_ShearViscosity = SphP[i].Eta_ShearViscosity;
    in->Zeta_BulkViscosity = SphP[i].Zeta_BulkViscosity;
#endif
    
#ifdef MAGNETIC
    for(k = 0; k < 3; k++) {in->BPred[k] = Get_Particle_BField(i,k);}
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
    in->Balpha = SphP[i].Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
    in->PhiPred = Get_Particle_PhiField(i);
#endif
#endif // MAGNETIC //
    
#ifdef COSMIC_RAYS
    for(j=0;j<N_CR_PARTICLE_BINS;j++)
    {
        in->CosmicRayPressure[j] = Get_Particle_CosmicRayPressure(i,j);
        in->CosmicRayDiffusionCoeff[j] = SphP[i].CosmicRayDiffusionCoeff[j];
#ifdef COSMIC_RAYS_M1
        for(k=0;k<3;k++) {in->CosmicRayFlux[j][k] = SphP[i].CosmicRayFluxPred[j][k];}
#endif
#ifdef COSMIC_RAYS_ALFVEN
        for(k=0;k<2;k++) {in->CosmicRayAlfvenEnergy[j][k] = SphP[i].CosmicRayAlfvenEnergyPred[j][k];}
#endif
    }
#endif

#ifdef EOS_ELASTIC
    in->CompositionType = SphP[i].CompositionType;
    {int k_v; for(k=0;k<3;k++) {for(k_v=0;k_v<3;k_v++) {in->Elastic_Stress_Tensor[k][k_v] = SphP[i].Elastic_Stress_Tensor_Pred[k][k_v];}}}
#endif
    
#ifdef GALSF_SUBGRID_WINDS
    in->DelayTime = SphP[i].DelayTime;
#endif

}



/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_hydra(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration);
static inline void out2particle_hydra(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int k;
    /* these are zero-d out at beginning of hydro loop so should always be added */
    for(k = 0; k < 3; k++)
    {
        SphP[i].HydroAccel[k] += out->Acc[k];
        //SphP[i].dMomentum[k] += out->dMomentum[k]; //manifest-indiv-timestep-debug//
    }
    SphP[i].DtInternalEnergy += out->DtInternalEnergy;
    //SphP[i].dInternalEnergy += out->dInternalEnergy; //manifest-indiv-timestep-debug//

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    SphP[i].DtMass += out->DtMass;
    SphP[i].dMass += out->dMass;
    for(k=0;k<3;k++) {SphP[i].GravWorkTerm[k] += out->GravWorkTerm[k];}
#endif
    if(SphP[i].MaxSignalVel < out->MaxSignalVel) {SphP[i].MaxSignalVel = out->MaxSignalVel;}
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    if(SphP[i].MaxKineticEnergyNgb < out->MaxKineticEnergyNgb) {SphP[i].MaxKineticEnergyNgb = out->MaxKineticEnergyNgb;}
#endif
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
    for(k=0;k<NUM_METAL_SPECIES;k++) {SphP[i].Dyield[k] += out->Dyield[k];}
#endif

#ifdef CHIMES_TURB_DIFF_IONS  
    for (k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++) 
      SphP[i].ChimesNIons[k] = DMAX(SphP[i].ChimesNIons[k] + out->ChimesIonsYield[k], 0.5 * SphP[i].ChimesNIons[k]); 
#endif 
    
#if defined(RT_SOLVER_EXPLICIT)
#if defined(RT_EVOLVE_ENERGY)
    for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[i].Dt_Rad_E_gamma[k] += out->Dt_Rad_E_gamma[k];}
#endif
#if defined(RT_EVOLVE_FLUX)
    for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Dt_Rad_Flux[k][k_dir] += out->Dt_Rad_Flux[k][k_dir];}}
#endif
#if defined(RT_INFRARED)
    SphP[i].Dt_Rad_E_gamma_T_weighted_IR += out->Dt_Rad_E_gamma_T_weighted_IR;
#endif
#if defined(RT_EVOLVE_INTENSITIES)
    for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {SphP[i].Dt_Rad_Intensity[k][k_dir] += out->Dt_Rad_Intensity[k][k_dir];}}
#endif
#endif
    
#if defined(MAGNETIC)
    /* can't just do DtB += out-> DtB, because for SPH methods, the induction equation is solved in the density loop; need to simply add it here */
    for(k=0;k<3;k++) {SphP[i].DtB[k] += out->DtB[k]; SphP[i].Face_Area[k] += out->Face_Area[k];}
    SphP[i].divB += out->divB;
#if defined(DIVBCLEANING_DEDNER)
#ifdef HYDRO_MESHLESS_FINITE_VOLUME // mass-based phi-flux
    SphP[i].DtPhi += out->DtPhi;
#endif
    for(k=0;k<3;k++) {SphP[i].DtB_PhiCorr[k] += out->DtB_PhiCorr[k];}
#endif // Dedner //
#endif // MAGNETIC //

#ifdef COSMIC_RAYS
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        SphP[i].DtCosmicRayEnergy[k] += out->DtCosmicRayEnergy[k];
#ifdef COSMIC_RAYS_ALFVEN
        int kAlf; for(kAlf=0;kAlf<2;kAlf++) {SphP[i].DtCosmicRayAlfvenEnergy[k][kAlf] += out->DtCosmicRayAlfvenEnergy[k][kAlf];}
#endif
    }
#endif
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
            dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* signal velocity needs to include rate of gas flow -over- the resolution element, which can be non-zero here */
            double v2_p = SphP[i].MaxSignalVel*SphP[i].MaxSignalVel;
            for(k=0;k<3;k++) {v2_p += (SphP[i].VelPred[k]-SphP[i].ParticleVel[k])*(SphP[i].VelPred[k]-SphP[i].ParticleVel[k]);}
            SphP[i].MaxSignalVel = sqrt(v2_p);
#endif

#if defined(MAGNETIC)
            /* need to subtract out the source terms proportional to the (non-zero) B-field divergence; to stabilize the scheme */
            for(k = 0; k < 3; k++)
            {
#ifndef HYDRO_SPH
                /* this part of the induction equation has to do with advection of div-B, it is not present in SPH */
                SphP[i].DtB[k] -= SphP[i].divB * SphP[i].VelPred[k]/All.cf_atime;
#endif
                SphP[i].HydroAccel[k] -= SphP[i].divB * Get_Particle_BField(i,k)*All.cf_a2inv;
                SphP[i].DtInternalEnergy -= SphP[i].divB * (SphP[i].VelPred[k]/All.cf_atime) * Get_Particle_BField(i,k)*All.cf_a2inv;
            }

            double magnorm_closure = Get_DtB_FaceArea_Limiter(i);

#if defined(DIVBCLEANING_DEDNER) && !defined(HYDRO_SPH)
            // ok now deal with the divB correction forces and damping fields //
            double tolerance_for_correction,db_vsig_h_norm;
            tolerance_for_correction = 10.0;
            db_vsig_h_norm = 0.1; // can be as low as 0.03 //
#ifdef PM_HIRES_REGION_CLIPPING
            tolerance_for_correction = 0.5; // could be as high as 0.75 //
            //db_vsig_h_norm = 0.03;
#endif

            double DtB_PhiCorr=0,DtB_UnCorr=0,db_vsig_h=0,PhiCorr_Norm=1.0;
            for(k=0; k<3; k++)
            {
                DtB_UnCorr += SphP[i].DtB[k] * SphP[i].DtB[k]; // physical units //
                db_vsig_h = db_vsig_h_norm * (SphP[i].BPred[k]*All.cf_atime) * (0.5*SphP[i].MaxSignalVel*All.cf_afac3) / (Get_Particle_Size(i)*All.cf_atime);
                DtB_UnCorr += db_vsig_h * db_vsig_h;
                DtB_PhiCorr += SphP[i].DtB_PhiCorr[k] * SphP[i].DtB_PhiCorr[k];
            }
            
            /* take a high power of these: here we'll use 4, so it works like a threshold */
            DtB_UnCorr*=DtB_UnCorr; DtB_PhiCorr*=DtB_PhiCorr; tolerance_for_correction *= tolerance_for_correction;
            /* now re-normalize the correction term if its unacceptably large */
            if((DtB_PhiCorr > 0)&&(!isnan(DtB_PhiCorr))&&(DtB_UnCorr>0)&&(!isnan(DtB_UnCorr))&&(tolerance_for_correction>0)&&(!isnan(tolerance_for_correction)))
            {
                
                if(DtB_PhiCorr > tolerance_for_correction * DtB_UnCorr) {PhiCorr_Norm *= tolerance_for_correction * DtB_UnCorr / DtB_PhiCorr;}
                for(k=0; k<3; k++)
                {
                    SphP[i].DtB[k] += PhiCorr_Norm * SphP[i].DtB_PhiCorr[k];
                    SphP[i].DtInternalEnergy += PhiCorr_Norm * SphP[i].DtB_PhiCorr[k] * Get_Particle_BField(i,k)*All.cf_a2inv;
                }
            }
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME // mass-based phi-flux
            SphP[i].DtPhi *= magnorm_closure;
#else
            SphP[i].DtPhi = 0;
#endif
            if((!isnan(SphP[i].divB))&&(PPP[i].Hsml>0)&&(SphP[i].divB!=0)&&(SphP[i].Density>0))
            {
                double tmp_ded = 0.5 * SphP[i].MaxSignalVel / (fac_mu*All.cf_atime); // has units of v_physical now
                /* do a check to make sure divB isn't something wildly divergent (owing to particles being too close) */
                double b2_max = 0.0;
                for(k=0;k<3;k++) {b2_max += Get_Particle_BField(i,k)*Get_Particle_BField(i,k);}
                b2_max = 100.0 * fabs( sqrt(b2_max) * All.cf_a2inv * P[i].Mass / (SphP[i].Density*All.cf_a3inv) * 1.0 / (PPP[i].Hsml*All.cf_atime) );
                if(fabs(SphP[i].divB) > b2_max) {SphP[i].divB *= b2_max / fabs(SphP[i].divB);}
                /* ok now can apply this to get the growth rate of phi */
                // SphP[i].DtPhi -= tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * SphP[i].divB;
                SphP[i].DtPhi -= tmp_ded * tmp_ded * All.DivBcleanHyperbolicSigma * SphP[i].divB * SphP[i].Density*All.cf_a3inv; // mass-based phi-flux
            }
#endif
#endif // MAGNETIC
            
            
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
            
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_DISABLE_STREAMING) && !defined(COSMIC_RAYS_ALFVEN)
            /* energy transfer from CRs to gas due to the streaming instability (mediated by high-frequency Alfven waves, but they thermalize quickly
                (note this is important; otherwise build up CR 'traps' where the gas piles up and cools but is entirely supported by CRs in outer disks) */
            double vstream_0 = Get_CosmicRayStreamingVelocity(i), vA=10.*vstream_0;
#ifdef MAGNETIC /* account for the fact that the loss term is always [or below] the Alfven speed, regardless of the bulk streaming speed */
            vA=0; for(k=0;k<3;k++) {vA += Get_Particle_BField(i,k)*Get_Particle_BField(i,k);}
            vA = All.cf_afac3 * sqrt(All.cf_afac1 * vA/ (All.cf_atime * SphP[i].Density));
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
            vA /= Get_Gas_Ionized_Fraction(i); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
#endif
            vstream_0 = DMIN(vA, vstream_0);
            for(k=0;k<N_CR_PARTICLE_BINS;k++)
            {
                double L_cr = Get_CosmicRayGradientLength(i,k), v_st_eff = SphP[i].CosmicRayDiffusionCoeff[k] / (GAMMA_COSMICRAY * L_cr + MIN_REAL_NUMBER); // maximum possible streaming speed from combined diffusivity
                double cr_stream_cool = fabs(GAMMA_COSMICRAY_MINUS1 * DMIN(v_st_eff, vstream_0) / L_cr); // if upper-limit to streaming is less than nominal 'default' v_stream/loss term, this should be lower too
                SphP[i].DtCosmicRayEnergy[k] -= SphP[i].CosmicRayEnergyPred[k] * cr_stream_cool; SphP[i].DtInternalEnergy += SphP[i].CosmicRayEnergyPred[k] * cr_stream_cool;
            }
#endif // CRs
            
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            SphP[i].DtInternalEnergy -= SphP[i].InternalEnergyPred * SphP[i].DtMass;
#endif
#ifdef MAGNETIC
#ifndef HYDRO_SPH
            for(k=0;k<3;k++)
            {
                SphP[i].DtInternalEnergy += -Get_Particle_BField(i,k)*All.cf_a2inv * SphP[i].DtB[k];
            }
#endif
            for(k=0;k<3;k++) {SphP[i].DtB[k] *= magnorm_closure;}
#endif
            SphP[i].DtInternalEnergy /= P[i].Mass;
            /* ok, now: HydroAccel = dv/dt, DtInternalEnergy = du/dt (energy per unit mass) */
            
            /* zero out hydrodynamic PdV work terms if the particle is at the maximum smoothing, these will be incorrect */
            if(PPP[i].Hsml >= 0.99*All.MaxHsml) {SphP[i].DtInternalEnergy = 0;}
            
            // need to explicitly include adiabatic correction from the hubble-flow (for drifting) here //
            if(All.ComovingIntegrationOn) SphP[i].DtInternalEnergy -= 3*(GAMMA(i)-1) * SphP[i].InternalEnergyPred * All.cf_hubble_a;
            // = du/dlna -3*(gamma-1)*u ; then dlna/dt = H(z) =  All.cf_hubble_a //
            
            
#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_EVOLVE_FLUX) && !defined(RT_RADPRESSURE_IN_HYDRO) //#elif defined(RT_COMPGRAD_EDDINGTON_TENSOR) /* // -- moved for OTVET+FLD to drift-kick operation to deal with limiters more accurately -- // */
            /* calculate the radiation pressure force */
            double radacc[3],fluxcorr; radacc[0]=radacc[1]=radacc[2]=0;  int kfreq;
            for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++)
            {
                double vol_inv = SphP[i].Density*All.cf_a3inv/P[i].Mass, f_kappa_abs = rt_absorb_frac_albedo(i,kfreq), vel_i[3]={0}, vdot_h[3]={0}, flux_i[3]={0}, flux_mag=0, erad_i=0, flux_corr=1, work_band=0;
                erad_i = SphP[i].Rad_E_gamma_Pred[kfreq]*vol_inv;
                for(k=0;k<3;k++) {flux_i[k]=SphP[i].Rad_Flux_Pred[kfreq][k]*vol_inv; vel_i[k]=SphP[i].VelPred[k]/All.cf_atime; flux_mag+=flux_i[k]*flux_i[k];}
                eddington_tensor_dot_vector(SphP[i].ET[kfreq],vel_i,vdot_h); for(k=0;k<3;k++) {vdot_h[k] = erad_i * (vel_i[k] + vdot_h[k]);} // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
                double flux_thin = erad_i * C_LIGHT_CODE_REDUCED; if(flux_mag>0) {flux_mag=sqrt(flux_mag);} else {flux_mag=1.e-20*flux_thin;}
                flux_corr = DMIN(1., flux_thin/flux_mag);
#if defined(RT_ENABLE_R15_GRADIENTFIX)
                flux_corr = flux_thin/flux_mag; // set to maximum (optically thin limit)
#endif
                for(k=0;k<3;k++)
                {
                    radacc[k] += (SphP[i].Rad_Kappa[kfreq]/C_LIGHT_CODE_REDUCED) * (flux_corr*flux_i[k] - vdot_h[k]); // note these 'vdoth' terms shouldn't be included in FLD, since its really assuming the entire right-hand-side of the flux equation reaches equilibrium with the pressure tensor, which gives the expression in rt_utilities
                    work_band += radacc[k] * vel_i[k] * P[i].Mass; // PdV work done by photons [absorbed ones are fully-destroyed, so their loss of energy and momentum is already accounted for by their deletion in this limit -- note that we have to be careful about the RSOL factors here! //
                }
                SphP[i].Dt_Rad_E_gamma[kfreq] += (2.*f_kappa_abs-1.)*work_band;
                SphP[i].DtInternalEnergy -= 2.*f_kappa_abs*work_band / P[i].Mass;
            }
            for(k=0;k<3;k++)
            {
#ifdef RT_RAD_PRESSURE_OUTPUT
                SphP[i].Rad_Accel[k] = radacc[k]; // physical units, as desired
#else
                SphP[i].HydroAccel[k] += radacc[k]; // physical units, as desired
#endif
            } 
#endif
#ifdef RT_RADPRESSURE_IN_HYDRO
            int kfreq; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {
                double fac = (1./3.) * return_flux_limiter(i,kfreq) * SphP[i].Rad_E_gamma_Pred[kfreq] * P[i].Particle_DivVel*All.cf_a2inv * (1.-2.*rt_absorb_frac_albedo(i,kfreq));
                SphP[i].Dt_Rad_E_gamma[kfreq] -= fac; SphP[i].DtInternalEnergy += fac / P[i].Mass; /* exact energy conservation */
            }
#endif
            
            
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME)) /* update the metal masses from exchange */
            for(k=0;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] = DMAX(P[i].Metallicity[k] + SphP[i].Dyield[k] / P[i].Mass , 0.01*P[i].Metallicity[k]);}
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
            
            
#ifdef BOX_BND_PARTICLES
            /* this flag signals all particles with id=0 are frozen (boundary particles) */
            if(P[i].ID == 0)
            {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[i].DtMass = 0;
                SphP[i].dMass = 0;
                for(k = 0; k < 3; k++) SphP[i].GravWorkTerm[k] = 0;
#endif
                SphP[i].DtInternalEnergy = 0;//SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
                for(k = 0; k < 3; k++) SphP[i].HydroAccel[k] = 0;//SphP[i].dMomentum[k] = 0;//manifest-indiv-timestep-debug//
#ifdef MAGNETIC
                for(k = 0; k < 3; k++) SphP[i].DtB[k] = 0;
#ifdef DIVBCLEANING_DEDNER
                for(k = 0; k < 3; k++) SphP[i].DtB_PhiCorr[k] = 0;
                SphP[i].DtPhi = 0;
#endif
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
            if(SphP[i].Temperature > All.NetworkTempThreshold)
            {
                nuc_particles++;
                network_integrate(SphP[i].Temperature, SphP[i].Density * All.UnitDensity_in_cgs, SphP[i].xnuc,
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
    MPI_Allreduce(&nuc_particles, &nuc_particles_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("Nuclear network done for %d particles.\n", nuc_particles_sum);
    }
    timewait1 += timediff(tend, my_second());
#endif
}




/* this function exists to loop over the hydro variables and do any needed 'pre-processing' before they enter the primary hydro force loop */
void hydro_force_initial_operations_preloop(void)
{
    // Set global factors for comoving integration of hydro //
    fac_mu = 1 / (All.cf_afac3 * All.cf_atime); // code_vel * fac_mu = sqrt[code_pressure/code_density] = code_soundspeed //
    fac_vsic_fix = All.cf_hubble_a * All.cf_afac1; // note also that signal_vel in forms below should be in units of code_soundspeed //
#ifdef MAGNETIC
    fac_magnetic_pressure = All.cf_afac1 / All.cf_atime; // code_Bfield*code_Bfield * fac_magnetic_pressure = code_pressure -- use this to get alfven velocities, etc, as well as comoving units for magnetic integration //
#endif
    
    /* need to zero out all numbers that can be set -EITHER- by an active particle in the domain, or by one of the neighbors we will get sent */
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type==0)
        {
            SphP[i].MaxSignalVel = -1.e10;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            SphP[i].MaxKineticEnergyNgb = -1.e10;
#endif
            SphP[i].DtInternalEnergy = 0; //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
            for(k=0;k<3;k++) {SphP[i].HydroAccel[k] = 0;} //SphP[i].dMomentum[k] = 0;//manifest-indiv-timestep-debug//
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            SphP[i].DtMass = 0; SphP[i].dMass = 0; for(k=0;k<3;k++) SphP[i].GravWorkTerm[k] = 0;
#endif
#if defined(TURB_DIFF_METALS) || (defined(METALS) && defined(HYDRO_MESHLESS_FINITE_VOLUME))
            for(k=0;k<NUM_METAL_SPECIES;k++) {SphP[i].Dyield[k] = 0;}
#endif
#if defined(RT_SOLVER_EXPLICIT)
#if defined(RT_EVOLVE_ENERGY)
            for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[i].Dt_Rad_E_gamma[k] = 0;}
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Dt_Rad_Flux[k][k_dir] = 0;}}
#endif
#if defined(RT_INFRARED)
            SphP[i].Dt_Rad_E_gamma_T_weighted_IR = 0;
#endif
#if defined(RT_EVOLVE_FLUX)
            for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Dt_Rad_Flux[k][k_dir] = 0;}}
#endif
#if defined(RT_EVOLVE_INTENSITIES)
            for(k=0;k<N_RT_FREQ_BINS;k++) {int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {SphP[i].Dt_Rad_Intensity[k][k_dir] = 0;}}
#endif
#endif
#ifdef MAGNETIC
            SphP[i].divB = 0; for(k=0;k<3;k++) {SphP[i].Face_Area[k] = 0;}
#ifdef DIVBCLEANING_DEDNER
            for(k=0;k<3;k++) {SphP[i].DtB_PhiCorr[k] = 0;}
#endif
#ifndef HYDRO_SPH
            for(k=0;k<3;k++) {SphP[i].DtB[k] = 0;}
#ifdef DIVBCLEANING_DEDNER
            SphP[i].DtPhi = 0;
#endif
#endif
#endif // magnetic //
#ifdef COSMIC_RAYS
            for(k=0;k<N_CR_PARTICLE_BINS;k++)
            {
                SphP[i].DtCosmicRayEnergy[k] = 0;
#ifdef COSMIC_RAYS_ALFVEN
                int kAlf; for(kAlf=0;kAlf<2;kAlf++) {SphP[i].DtCosmicRayAlfvenEnergy[k][kAlf] = 0;}
#endif
            }
#endif
#ifdef WAKEUP
            PPPZ[i].wakeup = 0;
#endif
        }
}





/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the driver routine for the calculation of hydrodynamical
 *  force, fluxes, etc. */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
void hydro_force(void)
{
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    hydro_force_initial_operations_preloop(); /* do initial pre-processing operations as needed before main hydro force loop */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    hydro_final_operations_and_cleanup(); /* do final operations on results */
    /* collect timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_HYDCOMPUTE] += timecomp; CPU_Step[CPU_HYDWAIT] += timewait; CPU_Step[CPU_HYDCOMM] += timecomm;
    CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */



