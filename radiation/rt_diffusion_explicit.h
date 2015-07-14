/* --------------------------------------------------------------------------------- */
/* ... explicit radiation diffusion evaluation ...
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
    int k_freq;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        Fluxes_E_gamma[k_freq] = 0;
        double kappa_i = local.Kappa_RT[k_freq];
        double kappa_j = SphP[j].Kappa_RT[k_freq];
        if((kappa_i>0)&&(kappa_j>0)&&(local.Mass>0)&&(P[j].Mass>0))
        {
            double scalar_i = local.E_gamma[k_freq] / V_i; // volumetric photon number density in this frequency bin //
            double scalar_j = SphP[j].E_gamma[k_freq] / V_j;
            
            double d_scalar = scalar_i - scalar_j;
	    	double conduction_wt = 0.5*(kappa_i+kappa_j) * All.cf_a3inv/All.cf_atime;  // weight factor and conversion to physical units
#ifdef HYDRO_SPH
            conduction_wt *= d_scalar * P[j].Mass * (0.5*(kernel.dwk_i+kernel.dwk_j)) / (kernel.r * local.Density * SphP[j].Density);
#else
            double cmag=0., c_max=0.;
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k];
                cmag += Face_Area_Vec[k] * 0.5*(local.Gradients.E_gamma_ET[k_freq][k] + SphP[j].Gradients.E_gamma_ET[k_freq][k]);
            }
            /* slope-limiter to ensure flow is always down the local gradient */
            c_max *= rinv*rinv;
            cmag = MINMOD(MINMOD(MINMOD(cmag , c_max*d_scalar), fabs(c_max)*d_scalar) , Face_Area_Norm*d_scalar*rinv);
            conduction_wt *= -cmag; // multiplies through the coefficient to get actual flux //
#endif
			Fluxes_E_gamma[k_freq] += conduction_wt;
        } // close check that kappa and particle masses are positive
    }
}
