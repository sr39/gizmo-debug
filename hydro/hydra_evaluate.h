/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*   In this routine, we find the SPH neighbors, and do the loop over 
*  neighbors to calculate the hydro fluxes. The actual flux calculation, 
*  and the returned values, should be in PHYSICAL (not comoving) units */
/* --------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------- */
int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, k, n, startnode, numngb, kernel_mode, listindex = 0;
    double hinv_i,hinv3_i,hinv4_i,hinv_j,hinv3_j,hinv4_j,V_i,V_j,dt_hydrostep,r2,rinv,rinv_soft,u;
    struct kernel_hydra kernel;
    struct hydrodata_in local;
    struct hydrodata_out out;
    struct Conserved_var_Riemann Fluxes;
    memset(&out, 0, sizeof(struct hydrodata_out));
    memset(&kernel, 0, sizeof(struct kernel_hydra));
#ifndef HYDRO_SPH
    struct Input_vec_Riemann Riemann_vec;
    struct Riemann_outputs Riemann_out;
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
#ifdef EOS_DEGENERATE
    kernel.sound_i = sqrt(local.dp_drho);
#else
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    kernel.sound_i = sqrt(GAMMA*local.Pressure/local.EgyWtRho);
#else
    kernel.sound_i = sqrt(GAMMA*local.Pressure/local.Density);
#endif
#endif
    kernel.h_i = local.Hsml;
    kernel_hinv(kernel.h_i, &hinv_i, &hinv3_i, &hinv4_i);
    hinv_j=hinv3_j=hinv4_j=0;
    V_i = local.Mass / local.Density;
    kernel_mode = -1; /* only need wk */
    dt_hydrostep = local.Timestep * All.Timebase_interval / All.cf_hubble_a; /* (physical) timestep */
    out.MaxSignalVel = kernel.sound_i;
#if defined(HYDRO_SPH) || defined(CONDUCTION_EXPLICIT) || defined(TURB_DIFFUSION)
    kernel_mode = 0; /* need dwk and wk */
#endif
    double cnumcrit2 = ((double)CONDITION_NUMBER_DANGER)*((double)CONDITION_NUMBER_DANGER) - local.ConditionNumber*local.ConditionNumber;
    double cs_t_to_comoving_x = All.cf_afac3 / All.cf_atime; /* convert to code (comoving) length units */
    double delta_halfstep_i=0,delta_halfstep_j=0;
    
    double rho_for_egy;
    rho_for_egy = local.Density;
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    rho_for_egy = local.EgyWtRho;
#endif
    
#if defined(HYDRO_SPH)
    kernel.p_over_rho2_i = local.Pressure / (rho_for_egy*rho_for_egy);
#endif
#if defined(HYDRO_SPH) || defined(CONDUCTION_EXPLICIT) || defined(TURB_DIFFUSION)
    kernel.spec_egy_u_i = local.Pressure / (GAMMA_MINUS1 * rho_for_egy);
#endif
    
#ifdef MAGNETIC
    kernel.b2_i = local.BPred[0]*local.BPred[0] + local.BPred[1]*local.BPred[1] + local.BPred[2]*local.BPred[2];
#ifdef MAGFORCE
    // PFH: comoving factors here to convert from B*B/rho to P/rho for accelerations //
    kernel.mf_i = fac_magnetic_pressure / (local.Density * local.Density);
    kernel.mf_j = fac_magnetic_pressure;
#endif
#ifdef MAGFORCE
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
                if(local.Timestep > TimeStep_J) continue; /* compute from particle with smaller timestep */
                if((local.Timestep == TimeStep_J) && (P[j].ID < local.ID)) continue; /* use ID to break degeneracy */
                if(P[j].Mass <= 0) continue;
#ifdef GALSF_SUBGRID_WINDS
                if(SphP[j].DelayTime > 0) continue; /* no hydro forces for decoupled wind particles */
#endif
                kernel.dx = local.Pos[0] - P[j].Pos[0];
                kernel.dy = local.Pos[1] - P[j].Pos[1];
                kernel.dz = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC  /* find the closest image in the given box size  */
                kernel.dx = NEAREST_X(kernel.dx);
                kernel.dy = NEAREST_Y(kernel.dy);
                kernel.dz = NEAREST_Z(kernel.dz);
#endif
                r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
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
                kernel.dvx = local.Vel[0] - SphP[j].VelPred[0];
                kernel.dvy = local.Vel[1] - SphP[j].VelPred[1];
                kernel.dvz = local.Vel[2] - SphP[j].VelPred[2];
                kernel.rho_ij_inv = 2.0 / (local.Density + SphP[j].Density);
                
                /* --------------------------------------------------------------------------------- */
                /* sound speed, relative velocity, and signal velocity computation */
                kernel.sound_j = Particle_effective_soundspeed_i(j);
#ifdef MAGNETIC_SIGNALVEL
                kernel.b2_j = SphP[j].BPred[0]*SphP[j].BPred[0] + SphP[j].BPred[1]*SphP[j].BPred[1] + SphP[j].BPred[2]*SphP[j].BPred[2];
                kernel.alfven2_j = kernel.b2_j * fac_magnetic_pressure / SphP[j].Density;
#ifdef ALFVEN_VEL_LIMITER
                kernel.alfven2_j = DMIN(kernel.alfven2_j, ALFVEN_VEL_LIMITER * kernel.sound_j*kernel.sound_j);
#endif
                double vcsa2_j = kernel.sound_j*kernel.sound_j + kernel.alfven2_j;
                double Bpro2_j = (SphP[j].BPred[0]*kernel.dx + SphP[j].BPred[1]*kernel.dy + SphP[j].BPred[2]*kernel.dz) / kernel.r;
                Bpro2_j *= Bpro2_j;
                double magneticspeed_j = sqrt(0.5 * (vcsa2_j + sqrt(DMAX((vcsa2_j*vcsa2_j -
                        4 * kernel.sound_j*kernel.sound_j * Bpro2_j*fac_magnetic_pressure/SphP[j].Density), 0))));
                double Bpro2_i = (local.BPred[0]*kernel.dx + local.BPred[1]*kernel.dy + local.BPred[2]*kernel.dz) / kernel.r;
                Bpro2_i *= Bpro2_i;
                double magneticspeed_i = sqrt(0.5 * (vcsa2_i + sqrt(DMAX((vcsa2_i*vcsa2_i -
                        4 * kernel.sound_i*kernel.sound_i * Bpro2_i*fac_magnetic_pressure/local.Density), 0))));
                kernel.vsig = magneticspeed_i + magneticspeed_j;
#endif
                kernel.vdotr2 = kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz;
                // hubble-flow correction: need in -code- units, hence extra a2 appearing here //
                if(All.ComovingIntegrationOn) kernel.vdotr2 += All.cf_hubble_a2 * r2;
                kernel.vsig = kernel.sound_i + kernel.sound_j;
                if(kernel.vdotr2 < 0)
                {
#ifdef HYDRO_SPH
                    kernel.vsig -= 3 * fac_mu * kernel.vdotr2 * rinv;
#else
                    kernel.vsig -= fac_mu * kernel.vdotr2 * rinv;
#endif
                }
                
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
                
#ifdef CONDUCTION_EXPLICIT
#include "conduction.h"
#endif
                
#ifdef TURB_DIFFUSION
#include "turbulent_diffusion.h"
#endif
                
                /* --------------------------------------------------------------------------------- */
                /* now we will actually assign the hydro variables for the evolution step */
                /* --------------------------------------------------------------------------------- */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                //out.dMass += Fluxes.rho * dt_hydrostep; //???
                out.DtMass += Fluxes.rho;
                //SphP[j].dMass -= Fluxes.rho * dt_hydrostep; //???
#endif
                for(k=0;k<3;k++)
                {
                    //out.dMomentum[k] += Fluxes.v[k] * dt_hydrostep; //???
                    out.Acc[k] += Fluxes.v[k];
                    //SphP[j].dMomentum[k] -= Fluxes.v[k] * dt_hydrostep; //???
                }
                //out.dInternalEnergy += Fluxes.p * dt_hydrostep; //???
                out.DtInternalEnergy += Fluxes.p;
                //SphP[j].dInternalEnergy -= Fluxes.p * dt_hydrostep; //???
                
                /* if this is particle j's active timestep, you should sent them the time-derivative
                 information as well, for their subsequent drift operations */
                if(TimeBinActive[P[j].TimeBin])
                {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                    SphP[j].DtMass -= Fluxes.rho;
#endif
                    for(k=0;k<3;k++) SphP[j].HydroAccel[k] -= Fluxes.v[k];
                    SphP[j].DtInternalEnergy -= Fluxes.p;
                }

                /* if we have mass fluxes, we need to have metal fluxes if we're using them (or any other passive scalars) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                if(Fluxes.rho != 0)
                {
                    double dmass;
                    dmass = Fluxes.rho * dt_hydrostep;
#ifdef METALS
                    if(Fluxes.rho > 0)
                    {
                        /* particle i gains mass from particle j */
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            out.Dyield[k] += (P[j].Metallicity[k] - local.Metallicity[k]) * dmass;
                    } else {
                        /* particle j gains mass from particle i */
                        dmass /= -P[j].Mass;
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            P[j].Metallicity[k] += (local.Metallicity[k] - P[j].Metallicity[k]) * dmass;
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

