/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the 'core' of the hydro force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*   In this routine, we find the gas particle neighbors, and do the loop over 
*  neighbors to calculate the hydro fluxes. The actual flux calculation, 
*  and the returned values, should be in PHYSICAL (not comoving) units */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, k, n, startnode, numngb, kernel_mode, listindex = 0;
    double hinv_i,hinv3_i,hinv4_i,hinv_j,hinv3_j,hinv4_j,V_i,V_j,dt_hydrostep,r2,rinv,rinv_soft,u;
    double v_hll,k_hll,b_hll; v_hll=k_hll=0,b_hll=1;
    struct kernel_hydra kernel;
    struct hydrodata_in local;
    struct hydrodata_out out;
    struct Conserved_var_Riemann Fluxes;
    memset(&out, 0, sizeof(struct hydrodata_out));
    memset(&kernel, 0, sizeof(struct kernel_hydra));
    memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifndef HYDRO_SPH
    struct Input_vec_Riemann Riemann_vec;
    struct Riemann_outputs Riemann_out;
    double face_area_dot_vel, face_vel_i=0, face_vel_j=0;
    double Face_Area_Vec[3], Face_Area_Norm = 0;
    face_area_dot_vel = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_MASS
    double epsilon_entropic_eos_big = 0.5; // can be anything from (small number=more diffusive, less accurate entropy conservation) to ~1.1-1.3 (least diffusive, most noisy)
    double epsilon_entropic_eos_small = 1.e-3; // should be << epsilon_entropic_eos_big
    if(All.ComovingIntegrationOn) {epsilon_entropic_eos_big = 0.6; epsilon_entropic_eos_small=1.e-2;}
#endif
#if defined(RT_EVOLVE_NGAMMA_IN_HYDRO)
    double Fluxes_E_gamma[N_RT_FREQ_BINS];
#endif
    
    if(mode == 0)
    {
        particle2in_hydra(&local, target); // this setup allows for all the fields we need to define (don't hard-code here)
    }
    else
    {
        local = HydroDataGet[target]; // this setup allows for all the fields we need to define (don't hard-code here)
    }
    
    /* --------------------------------------------------------------------------------- */
    /* pre-define Particle-i based variables (so we save time in the loop below) */
    /* --------------------------------------------------------------------------------- */
    kernel.sound_i = local.SoundSpeed;
    kernel.spec_egy_u_i = local.InternalEnergyPred;
    kernel.h_i = local.Hsml;
    kernel_hinv(kernel.h_i, &hinv_i, &hinv3_i, &hinv4_i);
    hinv_j=hinv3_j=hinv4_j=0;
    V_i = local.Mass / local.Density;
    dt_hydrostep = local.Timestep * All.Timebase_interval / All.cf_hubble_a; /* (physical) timestep */
    out.MaxSignalVel = kernel.sound_i;
    kernel_mode = 0; /* need dwk and wk */
    double cnumcrit2 = ((double)CONDITION_NUMBER_DANGER)*((double)CONDITION_NUMBER_DANGER) - local.ConditionNumber*local.ConditionNumber;
    //define units used for upwind instead of time-centered formulation//
    //double cs_t_to_comoving_x = All.cf_afac3 / All.cf_atime; /* convert to code (comoving) length units */
    //double delta_halfstep_i=0,delta_halfstep_j=0;
    
#if defined(HYDRO_SPH)
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    kernel.p_over_rho2_i = local.Pressure / (local.EgyWtRho*local.EgyWtRho);
#else 
    kernel.p_over_rho2_i = local.Pressure / (local.Density*local.Density);
#endif
#endif
    
#ifdef MAGNETIC
    kernel.b2_i = local.BPred[0]*local.BPred[0] + local.BPred[1]*local.BPred[1] + local.BPred[2]*local.BPred[2];
#if defined(HYDRO_SPH)
    double magfluxv[3]; magfluxv[0]=magfluxv[1]=magfluxv[2]=0;
    kernel.mf_i = local.Mass * fac_magnetic_pressure / (local.Density * local.Density);
    kernel.mf_j = local.Mass * fac_magnetic_pressure;
    // PFH: comoving factors here to convert from B*B/rho to P/rho for accelerations //
    double mm_i[3][3], mm_j[3][3];
    for(k = 0; k < 3; k++)
    {
        for(j = 0; j < 3; j++)
            mm_i[k][j] = local.BPred[k] * local.BPred[j];
    }
    for(k = 0; k < 3; k++)
        mm_i[k][k] -= 0.5 * kernel.b2_i;
#endif
#ifdef MAGNETIC_SIGNALVEL
    kernel.alfven2_i = kernel.b2_i * fac_magnetic_pressure / local.Density;
#ifdef ALFVEN_VEL_LIMITER
    kernel.alfven2_i = DMIN(kernel.alfven2_i, ALFVEN_VEL_LIMITER * kernel.sound_i*kernel.sound_i);
#endif
    double vcsa2_i = kernel.sound_i*kernel.sound_i + kernel.alfven2_i;
#endif
#endif // MAGNETIC //
    
    
    /* --------------------------------------------------------------------------------- */
    /* Now start the actual SPH computation for this particle */
    /* --------------------------------------------------------------------------------- */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
#ifndef DONOTUSENODELIST
        startnode = HydroDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
#else
        startnode = All.MaxPart;	/* root node */
#endif
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            /* --------------------------------------------------------------------------------- */
            /* get the neighbor list */
            /* --------------------------------------------------------------------------------- */
            numngb = ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
                                       exportnodecount, exportindex, ngblist);
            if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n];
                
                /* check if I need to compute this pair-wise interaction from "i" to "j", or skip it and 
                    let it be computed from "j" to "i" */
                int TimeStep_J = (P[j].TimeBin ? (1 << P[j].TimeBin) : 0);
#ifndef SHEARING_BOX // (shearing box means the fluxes at the boundaries are not actually symmetric, so can't do this) //
                if(local.Timestep > TimeStep_J) continue; /* compute from particle with smaller timestep */
                /* use relative positions to break degeneracy */
                if(local.Timestep == TimeStep_J)
                {
                    int n0=0; if(local.Pos[n0] == P[j].Pos[n0]) {n0++; if(local.Pos[n0] == P[j].Pos[n0]) n0++;}
                    if(local.Pos[n0] < P[j].Pos[n0]) continue;
                }
#endif
                if(P[j].Mass <= 0) continue;
                if(SphP[j].Density <= 0) continue;
#ifdef GALSF_SUBGRID_WINDS
                if(SphP[j].DelayTime > 0) continue; /* no hydro forces for decoupled wind particles */
#endif
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC  /* find the closest image in the given box size  */
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                kernel.h_j = PPP[j].Hsml;
                
                /* force applied for all particles inside each-others kernels! */
                if((r2 >= kernel.h_i * kernel.h_i) && (r2 >= kernel.h_j * kernel.h_j)) continue;
                if(r2 <= 0) continue;
                
                /* --------------------------------------------------------------------------------- */
                /* ok, now we definitely have two interacting particles */
                /* --------------------------------------------------------------------------------- */
                
                /* --------------------------------------------------------------------------------- */
                /* calculate a couple basic properties needed: separation, velocity difference (needed for timestepping) */
                kernel.r = sqrt(r2);
                rinv = 1 / kernel.r;
                /* we require a 'softener' to prevent numerical madness in interpolating functions */
                rinv_soft = 1.0 / sqrt(r2 + 0.0001*kernel.h_i*kernel.h_i);
#ifdef SHEARING_BOX
                /* in a shearing box, need to set dv appropriately for the shearing boundary conditions */
                MyDouble VelPred_j[3];
                for(k=0;k<3;k++) {VelPred_j[k]=SphP[j].VelPred[k];}
                if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[SHEARING_BOX_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#else
                /* faster to just set a pointer directly */
                MyDouble *VelPred_j = SphP[j].VelPred;
#endif
                kernel.dv[0] = local.Vel[0] - VelPred_j[0];
                kernel.dv[1] = local.Vel[1] - VelPred_j[1];
                kernel.dv[2] = local.Vel[2] - VelPred_j[2];
                kernel.rho_ij_inv = 2.0 / (local.Density + SphP[j].Density);
                
                /* --------------------------------------------------------------------------------- */
                /* sound speed, relative velocity, and signal velocity computation */
                kernel.sound_j = Particle_effective_soundspeed_i(j);
                kernel.vsig = kernel.sound_i + kernel.sound_j;
#ifdef COSMIC_RAYS
                double CosmicRayPressure_j = Get_Particle_CosmicRayPressure(j); /* compute this for use below */
#endif
#ifdef MAGNETIC
                double BPred_j[3];
                for(k=0;k<3;k++) {BPred_j[k]=Get_Particle_BField(j,k);} /* defined j b-field in appropriate units for everything */
#ifdef DIVBCLEANING_DEDNER
                double PhiPred_j = Get_Particle_PhiField(j); /* define j phi-field in appropriate units */
#endif
#ifdef MAGNETIC_SIGNALVEL
                kernel.b2_j = BPred_j[0]*BPred_j[0] + BPred_j[1]*BPred_j[1] + BPred_j[2]*BPred_j[2];
                kernel.alfven2_j = kernel.b2_j * fac_magnetic_pressure / SphP[j].Density;
#ifdef ALFVEN_VEL_LIMITER
                kernel.alfven2_j = DMIN(kernel.alfven2_j, ALFVEN_VEL_LIMITER * kernel.sound_j*kernel.sound_j);
#endif
                double vcsa2_j = kernel.sound_j*kernel.sound_j + kernel.alfven2_j;
                double Bpro2_j = (BPred_j[0]*kernel.dp[0] + BPred_j[1]*kernel.dp[1] + BPred_j[2]*kernel.dp[2]) / kernel.r;
                Bpro2_j *= Bpro2_j;
                double magneticspeed_j = sqrt(0.5 * (vcsa2_j + sqrt(DMAX((vcsa2_j*vcsa2_j -
                        4 * kernel.sound_j*kernel.sound_j * Bpro2_j*fac_magnetic_pressure/SphP[j].Density), 0))));
                double Bpro2_i = (local.BPred[0]*kernel.dp[0] + local.BPred[1]*kernel.dp[1] + local.BPred[2]*kernel.dp[2]) / kernel.r;
                Bpro2_i *= Bpro2_i;
                double magneticspeed_i = sqrt(0.5 * (vcsa2_i + sqrt(DMAX((vcsa2_i*vcsa2_i -
                        4 * kernel.sound_i*kernel.sound_i * Bpro2_i*fac_magnetic_pressure/local.Density), 0))));
                kernel.vsig = magneticspeed_i + magneticspeed_j;
#endif
#endif
                kernel.vdotr2 = kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2];
                // hubble-flow correction: need in -code- units, hence extra a2 appearing here //
                if(All.ComovingIntegrationOn) kernel.vdotr2 += All.cf_hubble_a2 * r2;
                if(kernel.vdotr2 < 0)
                {
#ifdef HYDRO_SPH
                    kernel.vsig -= 3 * fac_mu * kernel.vdotr2 * rinv;
#else
                    kernel.vsig -= fac_mu * kernel.vdotr2 * rinv;
#endif
                }
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
                double KE = kernel.dv[0]*kernel.dv[0] + kernel.dv[1]*kernel.dv[1] + kernel.dv[2]*kernel.dv[2];
                if(KE > out.MaxKineticEnergyNgb) out.MaxKineticEnergyNgb = KE;
                if(TimeBinActive[P[j].TimeBin])
                {
                    if(KE > SphP[j].MaxKineticEnergyNgb) SphP[j].MaxKineticEnergyNgb = KE;
                }
#endif
                
                /* --------------------------------------------------------------------------------- */
                /* calculate the kernel functions (centered on both 'i' and 'j') */
                if(kernel.r < kernel.h_i)
                {
                    u = kernel.r * hinv_i;
                    kernel_main(u, hinv3_i, hinv4_i, &kernel.wk_i, &kernel.dwk_i, kernel_mode);
                }
                else
                {
                    kernel.dwk_i = 0;
                    kernel.wk_i = 0;
                }
                if(kernel.r < kernel.h_j)
                {
                    kernel_hinv(kernel.h_j, &hinv_j, &hinv3_j, &hinv4_j);
                    u = kernel.r * hinv_j;
                    kernel_main(u, hinv3_j, hinv4_j, &kernel.wk_j, &kernel.dwk_j, kernel_mode);
                }
                else
                {
                    kernel.dwk_j = 0;
                    kernel.wk_j = 0;
                }
                
                /* --------------------------------------------------------------------------------- */
                /* with the overhead numbers above calculated, we now 'feed into' the "core" 
                    hydro computation (SPH, meshless godunov, etc -- doesn't matter, should all take the same inputs) 
                    the core code is -inserted- here from the appropriate .h file, depending on the mode 
                    the code has been compiled in */
                /* --------------------------------------------------------------------------------- */
#ifdef HYDRO_SPH
#include "hydra_core_sph.h"
#else
#include "hydra_core_meshless.h"
#endif

#ifndef HYDRO_SPH
/* the following macros are useful for all the diffusion operations below: this is the diffusion term associated
    with the HLL reimann problem solution. This adds numerical diffusion (albeit limited to the magnitude of the 
    physical diffusion coefficients), but stabilizes the relevant equations */
#ifdef MAGNETIC
                double bhat[3]={Riemann_out.Face_B[0],Riemann_out.Face_B[1],Riemann_out.Face_B[2]};
                double bhat_mag=bhat[0]*bhat[0]+bhat[1]*bhat[1]+bhat[2]*bhat[2];
                if(bhat_mag>0) {bhat_mag=1./sqrt(bhat_mag); bhat[0]*=bhat_mag; bhat[1]*=bhat_mag; bhat[2]*=bhat_mag;}
                v_hll = 0.5*fabs(face_vel_i-face_vel_j) + DMAX(magneticspeed_i,magneticspeed_j);
#define B_dot_grad_weights(grad_i,grad_j) {if(bhat_mag<=0) {b_hll=1;} else {double q_tmp_sum=0,b_tmp_sum=0; for(k=0;k<3;k++) {\
                                           double q_tmp=0.5*(grad_i[k]+grad_j[k]); q_tmp_sum+=q_tmp*q_tmp; b_tmp_sum+=bhat[k]*q_tmp;}\
                                           if((b_tmp_sum!=0)&&(q_tmp_sum>0)) {b_hll=fabs(b_tmp_sum)/sqrt(q_tmp_sum);} else {b_hll=0;}}}
#else
                v_hll = 0.5*fabs(face_vel_i-face_vel_j) + DMAX(kernel.sound_i,kernel.sound_j);
#define B_dot_grad_weights(grad_i,grad_j) {b_hll=1;}
#endif
#define HLL_correction(ui,uj,wt,kappa) (k_hll = v_hll * (wt) * kernel.r * All.cf_atime / fabs(kappa),\
                                        k_hll = (0.2 + k_hll) / (0.2 + k_hll + k_hll*k_hll),\
                                        -0.5*k_hll*Face_Area_Norm*v_hll*((ui)-(uj)))
#endif
                
                
#ifdef CONDUCTION
#include "conduction.h"
#endif

#ifdef VISCOSITY
#include "viscosity.h"
#endif
                
#ifdef TURB_DIFFUSION
#include "turbulent_diffusion.h"
#endif
                
#ifdef COSMIC_RAYS
#include "../galaxy_sf/cosmic_ray_diffusion.h"
#endif
                
#ifdef RT_DIFFUSION_EXPLICIT
#include "../radiation/rt_diffusion_explicit.h"
#endif
                
                
                /* --------------------------------------------------------------------------------- */
                /* now we will actually assign the hydro variables for the evolution step */
                /* --------------------------------------------------------------------------------- */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                double dmass_holder = Fluxes.rho * dt_hydrostep;
                double dmass_limiter = 0.01 * DMAX(0,DMIN(DMIN(local.Mass,SphP[j].MassTrue),P[j].Mass));
                if(fabs(dmass_holder) > dmass_limiter) {dmass_holder *= dmass_limiter / fabs(dmass_holder);}
                out.dMass += dmass_holder;
                out.DtMass += Fluxes.rho;
                SphP[j].dMass -= dmass_holder; 
                double gravwork[3]; gravwork[0]=Fluxes.rho*kernel.dp[0]; gravwork[1]=Fluxes.rho*kernel.dp[1]; gravwork[2]=Fluxes.rho*kernel.dp[2];
                for(k=0;k<3;k++) {out.GravWorkTerm[k] += gravwork[k];}
#endif
                for(k=0;k<3;k++)
                {
                    out.Acc[k] += Fluxes.v[k];
                    //out.dMomentum[k] += Fluxes.v[k] * dt_hydrostep; //manifest-indiv-timestep-debug//
                    //SphP[j].dMomentum[k] -= Fluxes.v[k] * dt_hydrostep; //manifest-indiv-timestep-debug//
                }
                out.DtInternalEnergy += Fluxes.p;
#if defined(RT_EVOLVE_NGAMMA_IN_HYDRO)
                for(k=0;k<N_RT_FREQ_BINS;k++) {out.Dt_E_gamma[k] += Fluxes_E_gamma[k];}
#endif
#ifdef MAGNETIC
                for(k=0;k<3;k++) {out.DtB[k]+=Fluxes.B[k];}
                out.divB += Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                out.DtPhi += Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                for(k=0;k<3;k++) {out.DtInternalEnergy+=magfluxv[k]*local.Vel[k]/All.cf_atime;}
#else
                for(k=0;k<3;k++) {out.Face_Area[k] += Face_Area_Vec[k];}
                double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_i);
                out.DtInternalEnergy += 0.5 * kernel.b2_i*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                //out.DtPhi += (Riemann_out.phi_normal_mean - local.PhiPred*All.cf_a3inv) * wt_face_sum; // now use mass-based phi-flux
                /*
                double phi_normal_full = Riemann_out.phi_normal_mean + Riemann_out.phi_normal_db;
                for(k=0;k<3;k++) {out.DtB_PhiCorr[k] += phi_normal_full * Face_Area_Vec[k];}
                */ 
                for(k=0; k<3; k++)
                {
                    out.DtB_PhiCorr[k] += Riemann_out.phi_normal_db * Face_Area_Vec[k];
                    out.DtB[k] += Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                    out.DtInternalEnergy += Riemann_out.phi_normal_mean * Face_Area_Vec[k] * local.BPred[k]*All.cf_a2inv;
                }
#endif
#endif
#endif // magnetic //
                
#ifdef COSMIC_RAYS
                out.DtCosmicRayEnergy += Fluxes.CosmicRayPressure;
#endif
                
                //out.dInternalEnergy += Fluxes.p * dt_hydrostep; //manifest-indiv-timestep-debug//
                //SphP[j].dInternalEnergy -= Fluxes.p * dt_hydrostep; //manifest-indiv-timestep-debug//
                
                /* if this is particle j's active timestep, you should sent them the time-derivative
                 information as well, for their subsequent drift operations */
#ifndef SHEARING_BOX
                if(TimeBinActive[P[j].TimeBin])
                {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                    SphP[j].DtMass -= Fluxes.rho;
                    for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] -= gravwork[k];}
#endif
                    for(k=0;k<3;k++) {SphP[j].HydroAccel[k] -= Fluxes.v[k];}
                    SphP[j].DtInternalEnergy -= Fluxes.p;
#if defined(RT_EVOLVE_NGAMMA_IN_HYDRO)
                    for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[j].Dt_E_gamma[k] -= Fluxes_E_gamma[k];}
#endif
#ifdef MAGNETIC
                    for(k=0;k<3;k++) {SphP[j].DtB[k]-=Fluxes.B[k];}
                    SphP[j].divB -= Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                    SphP[j].DtPhi -= Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                    for(k=0;k<3;k++) {SphP[j].DtInternalEnergy-=magfluxv[k]*VelPred_j[k]/All.cf_atime;}
#else
                    for(k=0;k<3;k++) {SphP[j].Face_Area[k] -= Face_Area_Vec[k];}
                    double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_j);
                    SphP[j].DtInternalEnergy -= 0.5 * kernel.b2_j*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                    //SphP[j].DtPhi -= (Riemann_out.phi_normal_mean - PhiPred_j*All.cf_a3inv) * wt_face_sum; // mass-based phi-flux
                    /*
                    for(k=0;k<3;k++) {SphP[j].DtB_PhiCorr[k] -= phi_normal_full * Face_Area_Vec[k];;}
                    */
                    for(k=0; k<3; k++)
                    {
                        SphP[j].DtB_PhiCorr[k] -= Riemann_out.phi_normal_db * Face_Area_Vec[k];
                        SphP[j].DtB[k] -= Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                        SphP[j].DtInternalEnergy -= Riemann_out.phi_normal_mean * Face_Area_Vec[k] * BPred_j[k]*All.cf_a2inv;
                    }
#endif
#endif
#endif // magnetic //

#ifdef COSMIC_RAYS
                    SphP[j].DtCosmicRayEnergy -= Fluxes.CosmicRayPressure;
#endif
                }
#endif

                /* if we have mass fluxes, we need to have metal fluxes if we're using them (or any other passive scalars) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                if(dmass_holder != 0)
                {
#ifdef METALS
                    if(Fluxes.rho > 0)
                    {
                        /* particle i gains mass from particle j */
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            out.Dyield[k] += (P[j].Metallicity[k] - local.Metallicity[k]) * dmass_holder;
                    } else {
                        /* particle j gains mass from particle i */
                        dmass_holder /= -P[j].Mass;
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            P[j].Metallicity[k] += (local.Metallicity[k] - P[j].Metallicity[k]) * dmass_holder;
                    }
#endif
                }
#endif
                
                /* --------------------------------------------------------------------------------- */
                /* don't forget to save the signal velocity for time-stepping! */
                /* --------------------------------------------------------------------------------- */
                if(kernel.vsig > out.MaxSignalVel) out.MaxSignalVel = kernel.vsig;
                if(TimeBinActive[P[j].TimeBin])
                    if(kernel.vsig > SphP[j].MaxSignalVel) SphP[j].MaxSignalVel = kernel.vsig;
#ifdef WAKEUP
                if(kernel.vsig > WAKEUP * SphP[j].MaxSignalVel) SphP[j].wakeup = 1;
#endif
                
                
            } // for(n = 0; n < numngb; n++) //
        } // while(startnode >= 0) //
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = HydroDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1) //
#endif
    } // while(startnode >= 0) //
    
    /* Now collect the result at the right place */
    if(mode == 0)
        out2particle_hydra(&out, target, 0);
    else
        HydroDataResult[target] = out;
    
    return 0;
}

