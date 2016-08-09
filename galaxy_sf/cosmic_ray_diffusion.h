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
    double kappa_i = fabs(local.CosmicRayDiffusionCoeff); // physical units
    double kappa_j = fabs(SphP[j].CosmicRayDiffusionCoeff);
    
    if(((kappa_i>0)||(kappa_j>0))&&(local.Mass>0)&&(P[j].Mass>0))
    {
        double d_scalar = scalar_i - scalar_j;
        
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        double *grad_i = local.Gradients.CosmicRayPressure;
        double *grad_j = SphP[j].Gradients.CosmicRayPressure;
        double flux_wt = 1;
        
        double diffusion_wt = 0.5*(kappa_i+kappa_j);
        int do_isotropic = 1;
        double b_hll=1, cmag=0, wt_i=0.5, wt_j=0.5, grad_dot_x_ij = 0;
        double grad_ij[3];
        for(k=0;k<3;k++)
        {
            double q_grad = wt_i*grad_i[k] + wt_j*grad_j[k];
            double q_direct = d_scalar * kernel.dp[k] * rinv*rinv;
            grad_dot_x_ij += q_grad * kernel.dp[k];
#ifdef MAGNETIC
            grad_ij[k] = MINMOD_G(q_grad , q_direct);
            if(q_grad*q_direct < 0) {if(fabs(q_direct) > 5.*fabs(q_grad)) {grad_ij[k] = 0.0;}}
#else
            grad_ij[k] = MINMOD(q_grad , q_direct);
#endif
            // negative coefficient is used here as shorthand for a particle being a local extremum in CR density.
            //  in this case we use a zeroth-order estimate for the flux: more diffusive, but needed to get the gradients resolved
            if((local.CosmicRayDiffusionCoeff<0)||(SphP[j].CosmicRayDiffusionCoeff<0)) {grad_ij[k] = q_direct;}
        }

        double grad_mag = 0.0;
        for(k=0;k<3;k++) {grad_mag += grad_ij[k]*grad_ij[k];}
        if(grad_mag > 0) {grad_mag = sqrt(grad_mag);} else {grad_mag=MIN_REAL_NUMBER;}
        //double Flux_squared_for_streamloss = grad_mag*grad_mag;

#ifdef MAGNETIC
        if(bhat_mag > 0)
        {
            do_isotropic = 0;
            double B_interface_dot_grad_T = 0.0;
            for(k=0;k<3;k++)
            {
                B_interface_dot_grad_T += bhat[k] * grad_ij[k];
                cmag += bhat[k] * Face_Area_Vec[k];
            }
            cmag *= B_interface_dot_grad_T;
            b_hll = B_interface_dot_grad_T / grad_mag;
            b_hll *= b_hll;
            //Flux_squared_for_streamloss *= b_hll;
        }
#endif
        if(do_isotropic) {for(k=0;k<3;k++) {cmag += Face_Area_Vec[k] * grad_ij[k];}}
        //double cmag_0 = cmag;
        cmag /= All.cf_atime; // cmag has units of u/r -- convert to physical
        
        double check_for_stability_sign = 1; /* if we are using the zeroth-order method, no HLL flux, etc is needed */
        if((local.CosmicRayDiffusionCoeff>=0)&&(SphP[j].CosmicRayDiffusionCoeff>=0))
        {
            /* obtain HLL correction terms for Reimann problem solution */
            double d_scalar_tmp = d_scalar - grad_dot_x_ij;
            double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp);
            double hll_corr = b_hll * flux_wt * HLL_correction(d_scalar_hll, 0, flux_wt, diffusion_wt) / (-diffusion_wt);
            double cmag_corr = cmag + hll_corr;
            cmag = MINMOD(HLL_DIFFUSION_COMPROMISE_FACTOR*cmag, cmag_corr);
            /* slope-limiter to ensure heat always flows from hot to cold */
            double d_scalar_b = b_hll * d_scalar;
            double f_direct = Face_Area_Norm*d_scalar_b*rinv/All.cf_atime;
            check_for_stability_sign = d_scalar*cmag;
            if((check_for_stability_sign < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag))) {cmag = 0;}
        }
        cmag *= -diffusion_wt; /* multiply through coefficient to get flux */
        //if(fabs(cmag_0) != 0) {Flux_squared_for_streamloss *= cmag*cmag / (cmag_0*cmag_0);}
        
        /* follow that with a flux limiter as well */
        diffusion_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(diffusion_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double CR_egy_i = local.CosmicRayPressure*V_i / GAMMA_COSMICRAY_MINUS1; // (E_cr = Volume * (Pressure/(GAMMA_CR-1))) - this is physical units //
            double CR_egy_j = CosmicRayPressure_j*V_j / GAMMA_COSMICRAY_MINUS1;
            double prefac_duij = 0.25, flux_multiplier = 1;
            if((local.CosmicRayDiffusionCoeff<0)||(SphP[j].CosmicRayDiffusionCoeff<0)) {prefac_duij = 0.05;}
            double du_ij_cond = prefac_duij * DMAX(DMIN( fabs(CR_egy_i-CR_egy_j) , DMAX(CR_egy_i , CR_egy_j)) , DMIN(CR_egy_i , CR_egy_j));
            if(diffusion_wt > 0) {du_ij_cond=DMIN(du_ij_cond,0.5*CR_egy_j);} else {du_ij_cond=DMIN(du_ij_cond,0.5*CR_egy_i);} // prevent flux from creating negative values //
            if(fabs(diffusion_wt)>du_ij_cond) {flux_multiplier = du_ij_cond/fabs(diffusion_wt);}
            diffusion_wt *= flux_multiplier;            
            Fluxes.CosmicRayPressure += diffusion_wt / dt_hydrostep;
            /*
            double x_dot_A=0; for(k=0;k<3;k++) {x_dot_A+=kernel.dp[k]*Face_Area_Vec[k];}
            Flux_squared_for_streamloss *= (GAMMA_COSMICRAY_MINUS1/GAMMA_COSMICRAY) * (1./(1.*NUMDIMS)) *
                (1./(0.5*(kappa_i+kappa_j) * 0.5*(scalar_i+scalar_j))) * (flux_multiplier*flux_multiplier) * (fabs(x_dot_A)*All.cf_atime); // rescale, convert to physical units
            if((Flux_squared_for_streamloss>0)&&(!isnan(Flux_squared_for_streamloss))) {Streaming_Loss_Term = Flux_squared_for_streamloss;}
            */
        } // if(diffusion_wt > 0)
        
    } // close check that kappa and particle masses are positive
}

