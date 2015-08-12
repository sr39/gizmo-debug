/* --------------------------------------------------------------------------------- */
/* ... real conduction evaluation ...
 *
 * For SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *  operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used.
 *  In either case, since we solve the conduction equations explicitly, a stronger timestep
 *  restriction is necessary (since the equations are parabolic); this is in timestep.c
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    double scalar_i = local.InternalEnergyPred;
    double scalar_j = SphP[j].InternalEnergyPred;
    double kappa_i = local.Kappa_Conduction; // physical units
    double kappa_j = SphP[j].Kappa_Conduction;
    
    if((kappa_i>0)&&(kappa_j>0)&&(local.Mass>0)&&(P[j].Mass>0))
    {
        double d_scalar = scalar_i - scalar_j;
        double rho_i, rho_j, rho_ij;
        rho_i = local.Density*All.cf_a3inv; rho_j = SphP[j].Density*All.cf_a3inv; rho_ij = 0.5*(rho_i+rho_j); // physical units
        
#ifdef HYDRO_SPH
        // SPH: use the sph 'effective areas' oriented along the lines between particles and direct-difference gradients
        double Face_Area_Norm = local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density);
        double diffusion_wt = -Face_Area_Norm * (d_scalar*rinv) * All.cf_atime; // multiplies implied gradient by face area: should give physical units //
#ifdef MAGNETIC
        kappa_i *= Bpro2_i; kappa_j *= Bpro2_j;
#endif
        double cmag = diffusion_wt * (2.*kappa_i*kappa_j/(kappa_i+kappa_j)); // geometric-weighted kappa (see Cleary & Monaghan '99)
        
#else
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        double *grad_i = local.Gradients.InternalEnergy;
        double *grad_j = SphP[j].Gradients.InternalEnergy;
        double flux_wt = rho_ij;
        
        double diffusion_wt = 0.5*(kappa_i+kappa_j);
        if(diffusion_wt > 0) {diffusion_wt = kappa_i*kappa_j / diffusion_wt;} else {diffusion_wt = 0;}
        
        int do_isotropic = 1;
        double b_hll=1, cmag=0, c_max=0, wt_i=0.5, wt_j=0.5;
        double grad_ij[3];
        for(k=0;k<3;k++)
        {
            double q_grad = wt_i*grad_i[k] + wt_j*grad_j[k];
            double q_direct = d_scalar * kernel.dp[k] * rinv*rinv;
            grad_ij[k] = MINMOD(q_grad , q_direct);
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
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k];
                cmag += bhat[k] * Face_Area_Vec[k];
            }
            cmag *= B_interface_dot_grad_T;
            if(grad_mag > 0) {grad_mag = sqrt(grad_mag);} else {grad_mag=1;}
            b_hll = B_interface_dot_grad_T / grad_mag;
            b_hll *= b_hll;
        }
#endif
        if(do_isotropic)
        {
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k];
                cmag += Face_Area_Vec[k] * grad_ij[k];
            }
        }
        c_max *= rinv*rinv;
        cmag /= All.cf_atime; c_max /= All.cf_atime; // c_max and cmag have units of u/r -- convert to physical
        
        /* obtain HLL correction terms for Reimann problem solution */
        cmag += flux_wt * HLL_correction(scalar_i,scalar_j, flux_wt, diffusion_wt) / (-diffusion_wt);
        
        /* slope-limiter to ensure heat always flows from hot to cold */
        double d_scalar_b = b_hll * d_scalar;
        cmag = MINMOD(MINMOD(MINMOD(cmag , c_max*d_scalar_b), fabs(c_max)*d_scalar_b) , Face_Area_Norm*d_scalar_b*rinv/All.cf_atime);
        
        /* now multiply through the coefficient to get the actual flux */
        cmag *= -diffusion_wt;
        
#endif // end of SPH/NOT SPH check
        
        /* follow that with a flux limiter as well */
        diffusion_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(diffusion_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double du_ij_cond = 0.5*DMIN(DMIN(0.5*fabs(DMIN(local.Mass,P[j].Mass)*d_scalar),
                                              local.Mass*scalar_i), P[j].Mass*scalar_j);
            if(fabs(diffusion_wt)>du_ij_cond) {diffusion_wt *= du_ij_cond/fabs(diffusion_wt);}
            Fluxes.p += diffusion_wt / dt_hydrostep;
        } // if(diffusion_wt > 0)
        
    } // close check that kappa and particle masses are positive
}
