/* --------------------------------------------------------------------------------- */
/* ... real cosmic ray diffusion/streaming evaluation ...
 *
 * This does the cosmic ray diffusion/streaming operations.
 *    Coefficients are calculated in gradients.c. Here they are used to stream/diffuse.
 * The SPH formulation of this is possible (copy from conduction.h, if desired), but
 *    known to be very noisy and systematically problematic, because of the inconsistent
 *    second-derivative operator. So MFM/MFV is strongly recommended, and much more
 *    accurate.
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    double cosmo_unit = All.cf_a3inv;
    double scalar_i = local.CosmicRayPressure * cosmo_unit;
    double scalar_j = CosmicRayPressure_j * cosmo_unit;
    double kappa_i = local.CosmicRayDiffusionCoeff; // physical units
    double kappa_j = SphP[j].CosmicRayDiffusionCoeff;
    
    if((kappa_i>0)&&(kappa_j>0)&&(local.Mass>0)&&(P[j].Mass>0))
    {
        double d_scalar = scalar_i - scalar_j;
#ifdef HYDRO_SPH
        // SPH: use the sph 'effective areas' oriented along the lines between particles and direct-difference gradients
        double Face_Area_Norm = local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density);
        double diffusion_wt = -Face_Area_Norm * (d_scalar*rinv) * All.cf_atime; // multiplies implied gradient by face area: should give physical units //
#ifdef MAGNETIC
        kappa_i *= Bpro2_i; kappa_j *= Bpro2_j;
#endif
        double cmag = diffusion_wt * (2.*kappa_i*kappa_j/(kappa_i+kappa_j)); // geometric-weighted kappa (see Cleary & Monaghan '99)
        double check_for_stability_sign = 1;
#else
        
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        double *grad_i = local.Gradients.CosmicRayPressure;
        double *grad_j = SphP[j].Gradients.CosmicRayPressure;
        double flux_wt = 1;
        
        double diffusion_wt = 0.5*(kappa_i+kappa_j);
        int do_isotropic = 1;
        double b_hll=1, cmag=0, wt_i=0.5, wt_j=0.5;
        double grad_ij[3];
        for(k=0;k<3;k++)
        {
            double q_grad = wt_i*grad_i[k] + wt_j*grad_j[k];
            double q_direct = d_scalar * kernel.dp[k] * rinv*rinv;
#ifdef MAGNETIC
            grad_ij[k] = MINMOD_G(q_grad , q_direct);
            if(q_grad*q_direct < 0) {if(fabs(q_direct) > 2.*fabs(q_grad)) {grad_ij[k] = 0.0;}}
#else
            grad_ij[k] = MINMOD(q_grad , q_direct);
#endif
        }
        
#ifdef MAGNETIC
        if(bhat_mag > 0)
        {
            do_isotropic = 0;
            double B_interface_dot_grad_T = 0.0, grad_mag = 0.0;
            for(k=0;k<3;k++)
            {
                B_interface_dot_grad_T += bhat[k] * grad_ij[k];
                grad_mag += grad_ij[k]*grad_ij[k];
            }
            for(k=0;k<3;k++) {cmag += bhat[k] * Face_Area_Vec[k];}
            cmag *= B_interface_dot_grad_T;
            if(grad_mag > 0) {grad_mag = sqrt(grad_mag);} else {grad_mag=1;}
            b_hll = B_interface_dot_grad_T / grad_mag;
            b_hll *= b_hll;
        }
#endif
        if(do_isotropic) {for(k=0;k<3;k++) {cmag += Face_Area_Vec[k] * grad_ij[k];}}
        cmag /= All.cf_atime; // cmag has units of u/r -- convert to physical
        
        /* obtain HLL correction terms for Reimann problem solution */
        double hll_corr = b_hll * flux_wt * HLL_correction(scalar_i,scalar_j, flux_wt, diffusion_wt) / (-diffusion_wt);
        double cmag_corr = cmag + hll_corr;
        cmag = MINMOD(HLL_DIFFUSION_COMPROMISE_FACTOR*cmag, cmag_corr);
        /* slope-limiter to ensure heat always flows from hot to cold */
        double d_scalar_b = b_hll * d_scalar;
        double f_direct = Face_Area_Norm*d_scalar_b*rinv/All.cf_atime;
        double check_for_stability_sign = d_scalar*cmag;
        if((check_for_stability_sign < 0) && (fabs(f_direct) > 0.005*fabs(cmag))) {cmag = 0;}
        cmag *= -diffusion_wt; /* multiply through coefficient to get flux */
        
#endif // end of SPH/NOT SPH check
        
        /* follow that with a flux limiter as well */
        diffusion_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(diffusion_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double CR_egy_i = local.CosmicRayPressure*V_i / GAMMA_COSMICRAY_MINUS1; // (E_cr = Volume * (Pressure/(GAMMA_CR-1))) - this is physical units //
            double CR_egy_j = CosmicRayPressure_j*V_j / GAMMA_COSMICRAY_MINUS1;
            double du_ij_cond = 1.0*DMIN(CR_egy_i, CR_egy_j);
            if(check_for_stability_sign<0) {du_ij_cond *= 1.e-2;}
            if(fabs(diffusion_wt)>du_ij_cond) {diffusion_wt *= du_ij_cond/fabs(diffusion_wt);}
            Fluxes.CosmicRayPressure += diffusion_wt / dt_hydrostep;
        } // if(diffusion_wt > 0)
        
    } // close check that kappa and particle masses are positive
}

