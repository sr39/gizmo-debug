/* --------------------------------------------------------------------------------- */
/* ... turbulent diffusion (sub-grid) models ... */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
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
    diffusion_du_ij = kernel.spec_egy_u_i - Particle_Internal_energy_i(j);
    Fluxes.p += local.Mass * diffusion_wt * diffusion_du_ij;
#endif
#ifdef TURB_DIFF_VELOCITY // turbulent 'velocity diffusion': this is a turbulent effective viscosity
    double diff_vdotr2_phys = kernel.vdotr2;
    if(All.ComovingIntegrationOn) diff_vdotr2_phys -= All.cf_hubble_a2 * r2;
    if(diff_vdotr2_phys < 0)
    {
        diffusion_wt_dt_m = KERNEL_CORE_SIZE * 0.5 * (kernel.h_i + kernel.h_j);
        dv2_ij = local.Mass * fac_mu*fac_mu * diffusion_wt * kernel.vdotr2 / ((r2 + 0.0001*diffusion_wt_dt_m*diffusion_wt_dt_m) * All.cf_atime);
        Fluxes.v[0] += dv2_ij * kernel.dx; /* momentum (physical units) */
        Fluxes.v[1] += dv2_ij * kernel.dy;
        Fluxes.v[2] += dv2_ij * kernel.dz;
        diff_vi_dot_r = local.Vel[0]*kernel.dx + local.Vel[1]*kernel.dy + local.Vel[2]*kernel.dz;
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
}
