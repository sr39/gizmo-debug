/* --------------------------------------------------------------------------------- */
/* ... turbulent diffusion (sub-grid) models ...
 *
 *  The basic equations here follow the Smagorinky eddy diffusion model. For SPH, the
 *    discretization comes from Wadsley 2008 & Shen 2010. However, some caution is needed,
 *    for SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *    operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used, which
 *    greatly minimizes artificial (numerical) diffusion.
 *  In either case, since we solve the diffusion equations explicitly, a stronger timestep
 *    restriction is necessary (since the equations are parabolic); this is in timestep.c.
 *    This is very important (many implementations of these equations in the literature
 *    do not include the appropriate timestep and flux-limiters; that makes the equations
 *    numerically unstable (you can get an answer, but it might be wrong, independent of resolution)
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
#ifdef HYDRO_SPH
    double diffusion_wt,diffusion_wt_dt,dv2_ij,diffusion_wt_dt_m,diffusion_m_min,diff_vi_dot_r,diffusion_du_ij;
    /* out.dA_dt +=  diffusion_wt*(local.A-P[j].A) // this is the template for all turbulent diffusion terms */
    diffusion_wt = (local.TD_DiffCoeff + SphP[j].TD_DiffCoeff) * P[j].Mass * kernel.rho_ij_inv * kernel.dwk_ij / kernel.r;
    /* note, units of TD_DiffCoeff have already been corrected so that this combination has physical units */
    diffusion_wt_dt = diffusion_wt * dt_hydrostep;
    if(fabs(diffusion_wt_dt)>0.01)
    {
        diffusion_wt *= 0.01/fabs(diffusion_wt_dt);
        diffusion_wt_dt = diffusion_wt * dt_hydrostep;
    }
#ifdef TURB_DIFF_ENERGY  // turbulent thermal energy diffusion //
    diffusion_du_ij = kernel.spec_egy_u_i - SphP[j].InternalEnergyPred;
    Fluxes.p += local.Mass * diffusion_wt * diffusion_du_ij;
#endif
#ifdef TURB_DIFF_VELOCITY // turbulent 'velocity diffusion': this is a turbulent effective viscosity
    double diff_vdotr2_phys = kernel.vdotr2;
    if(All.ComovingIntegrationOn) diff_vdotr2_phys -= All.cf_hubble_a2 * r2;
        if(diff_vdotr2_phys < 0)
        {
            diffusion_wt_dt_m = KERNEL_CORE_SIZE * 0.5 * (kernel.h_i + kernel.h_j);
            dv2_ij = local.Mass * fac_mu*fac_mu * diffusion_wt * kernel.vdotr2 / ((r2 + 0.0001*diffusion_wt_dt_m*diffusion_wt_dt_m) * All.cf_atime);
            Fluxes.v[0] += dv2_ij * kernel.dp[0]; /* momentum (physical units) */
            Fluxes.v[1] += dv2_ij * kernel.dp[1];
            Fluxes.v[2] += dv2_ij * kernel.dp[2];
            diff_vi_dot_r = local.Vel[0]*kernel.dp[0] + local.Vel[1]*kernel.dp[1] + local.Vel[2]*kernel.dp[2];
            Fluxes.p += dv2_ij * (diff_vi_dot_r - 0.5*diff_vdotr2_phys) / All.cf_atime; /* remember, this is -total- energy now */
        }
#endif
#ifdef TURB_DIFF_METALS // turbulent diffusion of metals (passive scalar mixing) //
    if((local.Mass>0)&&(P[j].Mass>0))
    {
        diffusion_wt_dt_m = 0.25 * (P[j].Mass + local.Mass) * diffusion_wt_dt;
        diffusion_m_min = 0.01 * DMIN(local.Mass,P[j].Mass);
        if(diffusion_wt_dt_m > 0) if(diffusion_wt_dt_m >  diffusion_m_min) diffusion_wt_dt_m =  diffusion_m_min;
            if(diffusion_wt_dt_m < 0) if(diffusion_wt_dt_m < -diffusion_m_min) diffusion_wt_dt_m = -diffusion_m_min;
                for(k=0;k<NUM_METAL_SPECIES;k++)
                {
                    out.Dyield[k] += diffusion_wt_dt_m * (local.Metallicity[k] - P[j].Metallicity[k]);
                    P[j].Metallicity[k] -= diffusion_wt_dt_m * (local.Metallicity[k] - P[j].Metallicity[k]) / P[j].Mass;
                }
    }
#endif
    
#else // ends SPH portion of these routines
    
    if((local.Mass>0)&&(P[j].Mass>0)&&(local.TD_DiffCoeff>0)&&(SphP[j].TD_DiffCoeff))
    {
        double diffusion_wt,diff_tmp;
        double wt_i,wt_j;
        // wt_i = wt_j = 0.5;
        wt_i = PPP[j].Hsml / (PPP[j].Hsml + local.Hsml); wt_j = 1.-wt_i; // this is consistent with our second-order face location //
        diffusion_wt = wt_i*local.TD_DiffCoeff + wt_j*SphP[j].TD_DiffCoeff; // arithmetic mean
        //diffusion_wt = 2.0 * (local.TD_DiffCoeff * SphP[j].TD_DiffCoeff) / (local.TD_DiffCoeff + SphP[j].TD_DiffCoeff); // geometric mean
        
        diffusion_wt *= All.cf_atime; // based on units TD_DiffCoeff is defined with, this makes it physical for a dimensionless quantity gradient below
        //diffusion_wt *= (wt_i*local.Density + wt_j*SphP[j].Density) * All.cf_a3inv; // mean density; should really use solution to the Riemann problem here;
        diffusion_wt *= Riemann_out.Face_Density;
        
        /* calculate the implied mass flux 'across the boundary' from the sub-grid model,
         so that we can apply a slope-limiter; this is needed for stability */
        double massflux = 0.0;
        for(k=0;k<3;k++) {massflux+=Face_Area_Vec[k]*Face_Area_Vec[k];}
        massflux = fabs( sqrt(massflux) * diffusion_wt /
                        (0.25*DMIN(local.Hsml,PPP[j].Hsml) * All.cf_atime)
                        * dt_hydrostep / (DMIN(local.Mass,P[j].Mass)) );
        if(massflux > 0.1) {diffusion_wt *= 0.1/massflux;}
        
#ifdef TURB_DIFF_ENERGY  // turbulent thermal energy diffusion //
        double diff_u = 0.0;
        for(k=0;k<3;k++)
        {
            diff_u += Face_Area_Vec[k] * (wt_i*local.Gradients.InternalEnergy[k]
                                        + wt_j*SphP[j].Gradients.InternalEnergy[k]);
        }
        double c_max = 2.0 * Face_Area_Norm * (local.InternalEnergyPred-SphP[j].InternalEnergyPred) * rinv; // inter-particle gradient times tolerance //
        diff_u = -diffusion_wt * MINMOD(c_max,diff_u);

        double conduction_wt = dt_hydrostep * diff_u; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double du_ij_cond = 0.5*DMIN(DMIN(0.5*fabs(local.Mass*local.InternalEnergyPred-P[j].Mass*SphP[j].InternalEnergyPred),
                                              local.Mass*local.InternalEnergyPred),
                                         P[j].Mass*SphP[j].InternalEnergyPred);
            if(fabs(conduction_wt)>du_ij_cond) {conduction_wt *= du_ij_cond/fabs(conduction_wt);}
            Fluxes.p += conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
#endif
#ifdef TURB_DIFF_VELOCITY // turbulent 'velocity diffusion': this is a turbulent effective viscosity
        double diff_vdotr2_phys = kernel.vdotr2;
        if(All.ComovingIntegrationOn) diff_vdotr2_phys -= All.cf_hubble_a2 * r2;
            if(diff_vdotr2_phys < 0)
            {
                int k_v;
                double diffusion_wt_v = diffusion_wt / All.cf_atime; // needed because Vcode = a * Vphys in gradient below
                for(k_v=0;k_v<3;k_v++)
                {
                    double diff_tmp = 0.0;
                    for(k=0;k<3;k++)
                    {
                        diff_tmp += Face_Area_Vec[k] * (wt_i*local.Gradients.Velocity[k_v][k]
                                                      + wt_j*SphP[j].Gradients.Velocity[k_v][k]);
                    }
                    double c_max = 2.0 * Face_Area_Norm * (local.Vel[k_v]-VelPred_j[k_v]) * rinv; // inter-particle gradient times tolerance //
                    diff_tmp = -diffusion_wt_v * MINMOD(c_max,diff_tmp);
                    Fluxes.v[k_v] += diff_tmp;
                    Fluxes.p += diff_tmp * 0.5*(local.Vel[k_v] + VelPred_j[k_v])/All.cf_atime;
                }
            }
#endif
#ifdef TURB_DIFF_METALS // turbulent diffusion of metals (passive scalar mixing) //
        int k_species;
        double diffusion_wt_z = diffusion_wt * dt_hydrostep;
        for(k_species=0;k_species<NUM_METAL_SPECIES;k_species++)
        {
            double diff_tmp = 0.0;
            for(k=0;k<3;k++)
            {
                diff_tmp += Face_Area_Vec[k] * (wt_i*local.Gradients.Metallicity[k_species][k]
                                              + wt_j*SphP[j].Gradients.Metallicity[k_species][k]);
            }
            double c_max = 2.0 * Face_Area_Norm * (local.Metallicity[k_species]-P[j].Metallicity[k_species]) * rinv; // inter-particle gradient times tolerance //
            diff_tmp = -diffusion_wt_z * MINMOD(c_max,diff_tmp);

            double zlim = 0.25*DMIN(DMIN(0.5*fabs(local.Mass*local.Metallicity[k_species]-P[j].Mass*P[j].Metallicity[k_species]),
                                              local.Mass*local.Metallicity[k_species]),
                                                P[j].Mass*P[j].Metallicity[k_species]);
            if(fabs(diff_tmp)>zlim) {diff_tmp*=zlim/fabs(diff_tmp);}
            out.Dyield[k_species] += diff_tmp;
            P[j].Metallicity[k_species] -= diff_tmp / P[j].Mass;
        }
#endif
    }
    
#endif
}
