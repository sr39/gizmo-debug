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
#ifdef TURB_DIFF_METALS // turbulent diffusion of metals (passive scalar mixing) //
    
#ifdef HYDRO_SPH
    /* out.dA_dt +=  diffusion_wt*(local.A-P[j].A) // this is the template for all turbulent diffusion terms */
    double diffusion_wt = (local.TD_DiffCoeff + SphP[j].TD_DiffCoeff) * P[j].Mass * kernel.rho_ij_inv * kernel.dwk_ij / kernel.r * All.cf_a2inv;
    /* note, units of TD_DiffCoeff have already been corrected so that this combination has physical units */
    double diffusion_wt_dt = diffusion_wt * dt_hydrostep;
    if(fabs(diffusion_wt_dt)>0.01)
    {
        diffusion_wt *= 0.01/fabs(diffusion_wt_dt);
        diffusion_wt_dt = diffusion_wt * dt_hydrostep;
    }
    if((local.Mass>0)&&(P[j].Mass>0))
    {
        double diffusion_wt_dt_m = 2. * local.Mass * diffusion_wt_dt;
        double diffusion_m_min = 0.05 * DMIN(local.Mass,P[j].Mass);
        if(diffusion_wt_dt_m > 0) if(diffusion_wt_dt_m >  diffusion_m_min) diffusion_wt_dt_m =  diffusion_m_min;
            if(diffusion_wt_dt_m < 0) if(diffusion_wt_dt_m < -diffusion_m_min) diffusion_wt_dt_m = -diffusion_m_min;
                for(k=0;k<NUM_METAL_SPECIES;k++)
                {
                    out.Dyield[k] += diffusion_wt_dt_m * (local.Metallicity[k] - P[j].Metallicity[k]);
                    P[j].Metallicity[k] -= diffusion_wt_dt_m * (local.Metallicity[k] - P[j].Metallicity[k]) / P[j].Mass;
                }
    }
    
#else // ends SPH portion of these routines
    
    if((local.Mass>0)&&(P[j].Mass>0)&&((local.TD_DiffCoeff>0)||(SphP[j].TD_DiffCoeff>0)))
    {
        double wt_i=0.5, wt_j=0.5, cmag, d_scalar;
        double diffusion_wt = wt_i*local.TD_DiffCoeff + wt_j*SphP[j].TD_DiffCoeff; // arithmetic mean
        diffusion_wt *= Riemann_out.Face_Density;
        double diffusion_wt_physical = diffusion_wt;
        diffusion_wt /= All.cf_atime; // based on units TD_DiffCoeff is defined with, this makes it physical for a dimensionless quantity gradient below
        /* calculate implied mass flux 'across the boundary' to prevent excessively large coefficients */
        double massflux = fabs( Face_Area_Norm * diffusion_wt_physical / (DMIN(kernel.h_i,kernel.h_j)*All.cf_atime) * dt_hydrostep / (DMIN(local.Mass,P[j].Mass)) );
        if(massflux > 0.25) {diffusion_wt *= 0.25/massflux;}
        
        int k_species;
        double diffusion_wt_z = diffusion_wt * dt_hydrostep;
        double rho_i = local.Density*All.cf_a3inv, rho_j = SphP[j].Density*All.cf_a3inv, rho_ij = 0.5*(rho_i+rho_j);
        for(k_species=0;k_species<NUM_METAL_SPECIES;k_species++)
        {
            cmag = 0.0;
            d_scalar = local.Metallicity[k_species]-P[j].Metallicity[k_species];
            for(k=0;k<3;k++)
            {
                double grad_ij = wt_i*local.Gradients.Metallicity[k_species][k] + wt_j*SphP[j].Gradients.Metallicity[k_species][k];
                double grad_direct = d_scalar * kernel.dp[k] * rinv*rinv;
                grad_ij = MINMOD(grad_ij , grad_direct);
                cmag += Face_Area_Vec[k] * grad_ij;
            }
            double hll_corr = rho_ij * HLL_correction(local.Metallicity[k_species], P[j].Metallicity[k_species], rho_ij, diffusion_wt_physical) / (-diffusion_wt);
            double cmag_corr = cmag + hll_corr;
            cmag = MINMOD(1.5*cmag, cmag_corr);
            double f_direct = Face_Area_Norm*d_scalar*rinv;
            if((d_scalar*cmag < 0) && (fabs(f_direct) > fabs(cmag))) {cmag = 0;}
            
            cmag *= -diffusion_wt_z;
            if(fabs(cmag) > 0)
            {
                double zlim = 0.5*DMIN(DMIN(0.5*fabs(DMIN(local.Mass,P[j].Mass)*d_scalar),local.Mass*local.Metallicity[k_species]),P[j].Mass*P[j].Metallicity[k_species]);
                if(fabs(cmag)>zlim) {cmag*=zlim/fabs(cmag);}
                out.Dyield[k_species] += cmag;
                P[j].Metallicity[k_species] -= cmag / P[j].Mass;
            }
        }
    }
    
#endif
#endif
}
