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
        double kappa_i = local.DiffusionCoeff[k_freq];
        double kappa_j = rt_diffusion_coefficient(j,k_freq);
        if((kappa_i>0)&&(kappa_j>0)&&(local.Mass>0)&&(P[j].Mass>0))
        {
            double scalar_i = local.E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin //
            double scalar_j = SphP[j].E_gamma_Pred[k_freq] / V_j;
            
            double d_scalar = scalar_i - scalar_j;
	    	double conduction_wt = 0.5*(kappa_i+kappa_j) * All.cf_a3inv/All.cf_atime;  // weight factor and conversion to physical units
#ifdef HYDRO_SPH
            double Face_Area_Norm = local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density);
            double d_ET[6];
            for(k=0;k<6;k++) {d_ET[k] = scalar_i*local.ET[k_freq][k] - scalar_j*SphP[j].ET[k_freq][k];}
            double flux_tmp = conduction_wt * (d_ET[0]*kernel.dp[0]*kernel.dp[0] + d_ET[1]*kernel.dp[1]*kernel.dp[1] + d_ET[2]*kernel.dp[2]*kernel.dp[2] +
                                               2.*(d_ET[3]*kernel.dp[0]*kernel.dp[1] + d_ET[4]*kernel.dp[1]*kernel.dp[2] +
                                                   d_ET[5]*kernel.dp[0]*kernel.dp[2])) / (1.e-37 + kernel.r * kernel.r * kernel.r);
            double cmag = -Face_Area_Norm * flux_tmp;
            cmag = MINMOD(cmag , -conduction_wt*Face_Area_Norm*d_scalar*rinv);
#else
            double cmag=0., c_max=0.;
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k];
                /* the flux is determined by the energy density gradient */
                double grad = 0.5 * (local.Gradients.E_gamma_ET[k_freq][k] + SphP[j].Gradients.E_gamma_ET[k_freq][k]);
                double grad_direct = d_scalar * kernel.dp[k] * rinv*rinv;
                grad = MINMOD( grad , grad_direct );
                cmag += Face_Area_Vec[k] * grad;
            }
            /* here we add the HLL-like correction term. this greatly reduces noise and improves the stability of the diffusion.
            	however it comes at the cost of (significant) additional numerical diffusion */
            double c_hll = 0.5*fabs(face_vel_i-face_vel_j) + c_light;
			double q = 0.5 * c_hll * kernel.r * All.cf_atime / fabs(1.e-37 + 0.5*(kappa_i+kappa_j));
			q = (0.2 + q) / (0.2 + q + q*q);
            double hll_tmp = -q * Face_Area_Norm * c_hll * d_scalar/(-conduction_wt); 

            /* add asymptotic-preserving correction so that numerical flux doesn't dominate in optically thick limit */
            double tau_c_j = Get_Particle_Size(j)*All.cf_atime * SphP[j].Kappa_RT[k_freq]*SphP[j].Density*All.cf_a3inv; // = L_particle / (lambda_mean_free_path) = L*kappa*rho //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            hll_tmp *= hll_corr;
            
            double thold_hll = 2.0*fabs(cmag);
            if(fabs(hll_tmp)>thold_hll) hll_tmp*=thold_hll/fabs(hll_tmp);
            cmag += hll_tmp;
            /* slope-limiter to ensure flow is always down the local gradient */
            c_max *= rinv*rinv;
            cmag = MINMOD(MINMOD(MINMOD(cmag , c_max*d_scalar), fabs(c_max)*d_scalar) , Face_Area_Norm*d_scalar*rinv);
            cmag *= -conduction_wt; // multiplies through the coefficient to get actual flux //
#endif
			// prevent super-luminal local fluxes //
			double R_flux = fabs(cmag) / (3. * Face_Area_Norm * c_light * fabs(d_scalar) + 1.e-37);
			R_flux = (1. + 12.*R_flux) / (1. + 12.*R_flux*(1.+R_flux)); // 12 arbitrary but >>1 gives good behavior here //
			cmag *= R_flux;
			cmag *= dt_hydrostep; // all in physical units //
			if(fabs(cmag) > 0)
			{
				// enforce a flux limiter for stability (to prevent overshoot) //
				double thold_hll = 0.25 * fabs(scalar_i*V_i-scalar_j*V_j); 
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
    double thold_hll_0 = 2.0;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_E_gamma[k_freq] = 0;
        Fluxes_Flux[k_freq][0]=Fluxes_Flux[k_freq][1]=Fluxes_Flux[k_freq][2]=0;
        double scalar_i = local.E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin //
        double scalar_j = SphP[j].E_gamma_Pred[k_freq] / V_j;
        if((scalar_i>0)&&(scalar_j>0)&&(local.Mass>0)&&(P[j].Mass>0)&&(dt_hydrostep>0))
        {
            
            double d_scalar = scalar_i - scalar_j;
            double cmag=0., c_max=0., cmag_flux[3]={0};
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k]; /* for flux-limiter function */
                /* the flux is already known (its explicitly evolved, rather than determined by the gradient of the energy density */
                cmag += Face_Area_Vec[k] * 0.5*(local.Flux[k_freq][k]/V_i + SphP[j].Flux_Pred[k_freq][k]/V_j); /* remember, our 'flux' variable is a volume-integral */
                int k_xyz=k, k_et_al, k_et_loop[3];
                if(k_xyz==0) {k_et_loop[0]=0; k_et_loop[1]=3; k_et_loop[2]=5;}
                if(k_xyz==1) {k_et_loop[0]=3; k_et_loop[1]=1; k_et_loop[2]=4;}
                if(k_xyz==2) {k_et_loop[0]=5; k_et_loop[1]=4; k_et_loop[2]=2;}
                for(k_et_al=0;k_et_al<3;k_et_al++)
                {
                    cmag_flux[k_xyz] += c_light*c_light * Face_Area_Vec[k_et_al] * 0.5*(scalar_i*local.ET[k_freq][k_et_loop[k_et_al]] + scalar_j*SphP[j].ET[k_freq][k_et_loop[k_et_al]]);
                }
            }
            c_max *= rinv*rinv;
            
            /* add asymptotic-preserving correction so that numerical flux doesn't unphysically dominate in optically thick limit */
            double tau_c_j = Get_Particle_Size(j)*All.cf_atime * SphP[j].Kappa_RT[k_freq]*SphP[j].Density*All.cf_a3inv; // = L_particle / (lambda_mean_free_path) = L*kappa*rho //
            double hll_corr = 1./(1. + 1.5*DMAX(tau_c_i[k_freq],tau_c_j));
            /* q below is a limiter to try and make sure the diffusion speed given by the hll flux doesn't exceed the diffusion speed in the diffusion limit */
            double q = 0.5 * c_hll * kernel.r * All.cf_atime / fabs(1.e-37 + 0.5*(local.DiffusionCoeff[k_freq]+rt_diffusion_coefficient(j,k_freq))); q = (0.2 + q) / (0.2 + q + q*q);
            double hll_tmp = -Face_Area_Norm * c_hll * d_scalar; /* simple HLL term for frame moving at 1/2 inter-particle velocity: here not limited */
            hll_tmp *= q * hll_corr;
            double thold_hll = thold_hll_0*fabs(cmag);
            if(fabs(hll_tmp)>thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
            cmag += hll_tmp;
            /* slope-limiter to ensure flow is always down the local gradient [no 'uphill' flow] */
            cmag = MINMOD(cmag, -Face_Area_Norm * c_hll * d_scalar);
            cmag *= dt_hydrostep; // all in physical units //
            if(fabs(cmag) > 0)
            {
                // enforce a flux limiter for stability (to prevent overshoot) //
                thold_hll = 0.25 * fabs(scalar_i*V_i-scalar_j*V_j);
                if(fabs(cmag)>thold_hll) {cmag *= thold_hll/fabs(cmag);}
                Fluxes_E_gamma[k_freq] += cmag / dt_hydrostep;
            } // if(conduction_wt > 0)
            
            for(k=0;k<3;k++)
            {
                double d_flux = local.Flux[k_freq][k]/V_i - SphP[j].Flux_Pred[k_freq][k]/V_j;
                hll_tmp = -Face_Area_Norm * c_hll * d_flux;
                thold_hll = thold_hll_0*fabs(cmag_flux[k]);
                if(fabs(hll_tmp)>thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
                cmag_flux[k] += hll_tmp;
                Fluxes_Flux[k_freq][k] += cmag_flux[k];
            }
        } // close check that energy and masses are positive
    }
}

#endif
