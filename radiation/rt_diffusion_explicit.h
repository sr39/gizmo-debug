/* --------------------------------------------------------------------------------- */
/* ... explicit radiation diffusion/flux transport evaluation ...
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
double c_light = RT_SPEEDOFLIGHT_REDUCTION * (C/All.UnitVelocity_in_cm_per_s);

#if !defined(RT_EVOLVE_FLUX) /* this means we just solve the diffusion equation for the eddington tensor, done in the loop below */
{
    int k_freq;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_E_gamma[k_freq] = 0;
        double kappa_ij = 0.5 * (local.DiffusionCoeff[k_freq] + rt_diffusion_coefficient(j,k_freq)); // physical
        if((kappa_ij>0)&&(local.Mass>0)&&(P[j].Mass>0))
        {
            double scalar_i = local.E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin (1/code volume) //
            double scalar_j = SphP[j].E_gamma_Pred[k_freq] / V_j;
            
            double d_scalar = (scalar_i - scalar_j); // units (1/code volume)
            double conduction_wt = kappa_ij * All.cf_a3inv/All.cf_atime;  // weight factor and conversion to physical units
            double cmag=0., grad_norm=0, grad_dot_x_ij=0.0;
            for(k=0;k<3;k++)
            {
                /* the flux is determined by the energy density gradient */
                double grad = 0.5 * (local.Gradients.E_gamma_ET[k_freq][k] + SphP[j].Gradients.E_gamma_ET[k_freq][k]); // (1/(code volume*code length))
                double grad_direct = d_scalar * kernel.dp[k] * rinv*rinv; // (1/(code volume*code length))
                grad_dot_x_ij += grad * kernel.dp[k];
                grad = MINMOD_G( grad , grad_direct );
#if defined(GALSF) || defined(COOLING) || defined(BLACKHOLES)
                double grad_direct_vs_abs_fac = 2.0;
#else
                double grad_direct_vs_abs_fac = 5.0;
#endif
                if(grad*grad_direct < 0) {if(fabs(grad_direct) > grad_direct_vs_abs_fac*fabs(grad)) {grad = 0.0;}}
                cmag += Face_Area_Vec[k] * grad;
                grad_norm += grad*grad;
            }
            double A_dot_grad_alignment = cmag*cmag / (Face_Area_Norm*Face_Area_Norm * grad_norm); // dimensionless
            cmag *= -conduction_wt; // multiplies through the coefficient to get actual flux (physical) //

            /* here we add the HLL-like correction term. this greatly reduces noise and improves the stability of the diffusion.
            	however it comes at the cost of (significant) additional numerical diffusion */
            double v_eff_light = DMIN(c_light , kappa_ij / (Get_Particle_Size(j)*All.cf_atime)); // physical
            double c_hll = 0.5*fabs(face_vel_i-face_vel_j) + v_eff_light;
            double q = 0.5 * c_hll * kernel.r * All.cf_atime / fabs(1.e-37 + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q); // physical
            double d_scalar_tmp = d_scalar - grad_dot_x_ij; // (1/code volume)
            double d_scalar_hll = MINMOD(d_scalar , d_scalar_tmp) * All.cf_a3inv; // physical
            double hll_tmp = -A_dot_grad_alignment * q * Face_Area_Norm * c_hll * d_scalar_hll; // physical
            
            /* add asymptotic-preserving correction so that numerical flux doesn't dominate in optically thick limit */
            double tau_c_j = Get_Particle_Size(j)*All.cf_atime * SphP[j].Kappa_RT[k_freq]*SphP[j].Density*All.cf_a3inv; // = L_particle / (lambda_mean_free_path) = L*kappa*rho (physical) //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            hll_tmp *= hll_corr;
            
            double thold_hll = 2.0*fabs(cmag);
            if(fabs(hll_tmp)>thold_hll) {hll_tmp*=thold_hll/fabs(hll_tmp);}
            double cmag_corr = cmag + hll_tmp;
            cmag = MINMOD(HLL_DIFFUSION_COMPROMISE_FACTOR*cmag, cmag_corr);
            /* flux-limiter to ensure flow is always down the local gradient */
            double f_direct = -conduction_wt * (1./9.) * Face_Area_Norm*d_scalar*rinv; // physical
            double check_for_stability_sign = f_direct*cmag;
            if((check_for_stability_sign < 0) && (fabs(f_direct) > HLL_DIFFUSION_OVERSHOOT_FACTOR*fabs(cmag))) {cmag = 0;}

            // prevent super-luminal local fluxes //
            double R_flux = fabs(cmag) / (3. * Face_Area_Norm * c_light * (fabs(d_scalar)*All.cf_a3inv) + 1.e-37); // physical
            R_flux = (1. + 12.*R_flux) / (1. + 12.*R_flux*(1.+R_flux)); // 12 arbitrary but >>1 gives good behavior here //
#ifndef FREEZE_HYDRO
            cmag *= R_flux;
#endif
            cmag *= dt_hydrostep; // all in physical units //
            if(fabs(cmag) > 0)
            {
                // enforce a flux limiter for stability (to prevent overshoot) //
                double thold_hll = 0.25 * DMIN(fabs(scalar_i*V_i-scalar_j*V_j),DMAX(fabs(scalar_i*V_i),fabs(scalar_j*V_j))); // physical
                if(check_for_stability_sign<0) {thold_hll *= 1.e-2;}
                if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                Fluxes_E_gamma[k_freq] += cmag / dt_hydrostep;
            } // if(conduction_wt > 0)
            
        } // close check that kappa and particle masses are positive
    }
}


#else /* RT_EVOLVE_FLUX is ON, so we don't solve a diffusion equation, but a system of two standard advection-like equations */


{
    int k_freq;
    double c_hll = 0.5*fabs(face_vel_i-face_vel_j) + c_light;
    double thold_hll_0 = 0.1, thold_hll_0RE=1.1;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_E_gamma[k_freq] = 0;
        Fluxes_Flux[k_freq][0]=Fluxes_Flux[k_freq][1]=Fluxes_Flux[k_freq][2]=0;
        double scalar_i = local.E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin (1/code volume)//
        double scalar_j = SphP[j].E_gamma_Pred[k_freq] / V_j;
        if((scalar_i>0)&&(scalar_j>0)&&(local.Mass>0)&&(P[j].Mass>0)&&(dt_hydrostep>0))
        {
            double d_scalar = scalar_i - scalar_j;
            double cmag=0., cmag_flux[3]={0}, grad_norm=0;
            for(k=0;k<3;k++)
            {
                /* the flux is already known (its explicitly evolved, rather than determined by the gradient of the energy density */
                double grad = 0.5*(local.Flux[k_freq][k]/V_i + SphP[j].Flux_Pred[k_freq][k]/V_j);
                grad_norm += grad*grad;
                cmag += Face_Area_Vec[k] * grad; /* remember, our 'flux' variable is a volume-integral */
                
                int k_xyz=k, k_et_al, k_et_loop[3];
                if(k_xyz==0) {k_et_loop[0]=0; k_et_loop[1]=3; k_et_loop[2]=5;}
                if(k_xyz==1) {k_et_loop[0]=3; k_et_loop[1]=1; k_et_loop[2]=4;}
                if(k_xyz==2) {k_et_loop[0]=5; k_et_loop[1]=4; k_et_loop[2]=2;}
                for(k_et_al=0;k_et_al<3;k_et_al++)
                {
                    cmag_flux[k_xyz] += c_light*c_light * Face_Area_Vec[k_et_al] * 0.5*(scalar_i*local.ET[k_freq][k_et_loop[k_et_al]] + scalar_j*SphP[j].ET[k_freq][k_et_loop[k_et_al]]);
                }
            }
            double A_dot_grad_alignment = cmag*cmag / (Face_Area_Norm*Face_Area_Norm * grad_norm); // dimensionless

            /* add asymptotic-preserving correction so that numerical flux doesn't unphysically dominate in optically thick limit */
            double kappa_ij = 0.5 * (local.DiffusionCoeff[k_freq] + rt_diffusion_coefficient(j,k_freq)); // physical
            double v_eff_light = DMIN(c_light , kappa_ij / (Get_Particle_Size(j)*All.cf_atime)); // physical
            c_hll = 0.5*fabs(face_vel_i-face_vel_j) + v_eff_light; // physical
            double tau_c_j = Get_Particle_Size(j)*All.cf_atime * SphP[j].Kappa_RT[k_freq]*SphP[j].Density*All.cf_a3inv; // = L_particle / (lambda_mean_free_path) = L*kappa*rho //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            /* q below is a limiter to try and make sure the diffusion speed given by the hll flux doesn't exceed the diffusion speed in the diffusion limit */
            double q = 0.5 * c_hll * kernel.r * All.cf_atime / fabs(1.e-37 + kappa_ij); q = (0.2 + q) / (0.2 + q + q*q); // physical
            double renormerFAC = DMIN(1.,fabs(A_dot_grad_alignment * q * hll_corr));
            
            double hll_tmp = -Face_Area_Norm * c_hll * d_scalar; /* simple HLL term for frame moving at 1/2 inter-particle velocity: here not limited */
            hll_tmp *= renormerFAC;
            
            double thold_hll = thold_hll_0*fabs(cmag);
            if(fabs(hll_tmp)>thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
            double cmag_corr = cmag + hll_tmp;
            cmag = MINMOD(thold_hll_0RE*cmag, cmag_corr);
            
            /* flux-limiter to ensure flow is always down the local gradient [no 'uphill' flow] */
            double f_direct = -Face_Area_Norm * c_hll * d_scalar * renormerFAC;
            if((f_direct*cmag < 0) && (fabs(f_direct) > fabs(cmag))) {cmag = 0;}
            
            cmag *= dt_hydrostep; // all in physical units //
            if(fabs(cmag) > 0)
            {
                // enforce a flux limiter for stability (to prevent overshoot) //
                thold_hll = 0.25 * DMIN(fabs(scalar_i*V_i-scalar_j*V_j),DMAX(fabs(scalar_i*V_i),fabs(scalar_j*V_j)));
                if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                Fluxes_E_gamma[k_freq] += cmag / dt_hydrostep;
            } // if(conduction_wt > 0)
            
            for(k=0;k<3;k++)
            {
                double flux_i=local.Flux[k_freq][k], flux_j=SphP[j].Flux_Pred[k_freq][k];
                double d_flux =  flux_i/V_i - flux_j/V_j;
                hll_tmp = -Face_Area_Norm * c_hll * d_flux * renormerFAC;
                if(fabs(hll_tmp) > 0)
                {
                    thold_hll = thold_hll_0 * fabs(cmag_flux[k]);
                    if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
                    double cmag_tmp = cmag_flux[k] + hll_tmp;
                    cmag_flux[k] = MINMOD(thold_hll_0RE*cmag_flux[k], cmag_tmp);
                }
                cmag_flux[k] *= dt_hydrostep;
                if(fabs(cmag_flux[k]) > 0)
                {
                    thold_hll = 1.e-37 + 2.0 * DMIN(fabs(flux_i-flux_j),DMAX(fabs(flux_i),fabs(flux_j)));
                    if(fabs(cmag_flux[k])>thold_hll) {cmag_flux[k] *= thold_hll/fabs(cmag_flux[k]);}
                    Fluxes_Flux[k_freq][k] += cmag_flux[k] / dt_hydrostep;
                }
            }
        } // close check that energy and masses are positive
    }
}

#endif
