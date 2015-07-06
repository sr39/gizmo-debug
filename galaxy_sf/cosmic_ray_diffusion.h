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
#ifdef COSMIC_RAYS
    if((local.CosmicRayDiffusionCoeff>0)&&(SphP[j].CosmicRayDiffusionCoeff>0)&&(local.Mass>0)&&(P[j].Mass>0))
    {
#ifdef HYDRO_SPH
        
        /* First, we have the usual SPH anisotropic conduction definition */
        kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);
        double conduction_wt = P[j].Mass * kernel.dwk_ij / (kernel.r * local.Density * SphP[j].Density) * All.cf_atime; // physical units //
        double CR_pressure_i = local.CosmicRayPressure * All.cf_a3inv; // physical units
        double CR_pressure_j = CosmicRayPressure_j * All.cf_a3inv; // physical units
        double du_ij_cond = CR_pressure_i - CR_pressure_j;
        conduction_wt *= du_ij_cond; // multiply by specific energy difference //
#ifdef MAGNETIC
        // account for suppression of conduction along field lines //
#ifndef MAGNETIC_SIGNALVEL
        double Bpro2_i=0;
        double Bpro2_j=0;
        if(kernel.b2_i>0)
        {
            Bpro2_i = local.BPred[0]*kernel.dp[0] + local.BPred[1]*kernel.dp[1] + local.BPred[2]*kernel.dp[2];
            Bpro2_i *= Bpro2_i / (kernel.b2_i * kernel.r*kernel.r);
        }
        if(kernel.b2_j>0)
        {
            Bpro2_j = BPred_j[0]*kernel.dp[0] + BPred_j[1]*kernel.dp[1] + BPred_j[2]*kernel.dp[2];
            Bpro2_j *= Bpro2_j / (kernel.b2_j * kernel.r*kernel.r);
        }
#else
        // these are already calculated in MAGNETIC_SIGNALVEL above
        Bpro2_i /= kernel.b2_i;
        Bpro2_j /= kernel.b2_j;
#endif
        conduction_wt *= 2.0 * local.CosmicRayDiffusionCoeff*Bpro2_i*SphP[j].CosmicRayDiffusionCoeff*Bpro2_j /
        (local.CosmicRayDiffusionCoeff*Bpro2_i + SphP[j].CosmicRayDiffusionCoeff*Bpro2_j);
#else
        conduction_wt *= 2.0 * local.CosmicRayDiffusionCoeff*SphP[j].CosmicRayDiffusionCoeff/(local.CosmicRayDiffusionCoeff + SphP[j].CosmicRayDiffusionCoeff);
#endif
        // compute actual change that will occur this timestep //
        conduction_wt *= dt_hydrostep; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a limiter for stability (to prevent artificial oscillations) //
            double CR_egy_i = local.CosmicRayPressure*(V_i/All.cf_a3inv) / GAMMA_COSMICRAY_MINUS1; // (E_cr = Volume * (Pressure/(GAMMA_CR-1)))
            double CR_egy_j = CR_pressure_j*(V_j/All.cf_a3inv) / GAMMA_COSMICRAY_MINUS1;
            double du_ij_cond = 0.5*DMIN(DMIN(0.5*fabs(CR_egy_i-CR_egy_j),CR_egy_i),CR_egy_j);
            if(fabs(conduction_wt)>du_ij_cond) {conduction_wt *= du_ij_cond/fabs(conduction_wt);}
            // now apply time rate of change to particle 'i'
            Fluxes.CosmicRayPressure += conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
        
#else
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        
        double conduction_wt;
        double wt_i,wt_j;
        wt_i = wt_j = 0.5;
        //wt_i = PPP[j].Hsml / (PPP[j].Hsml + local.Hsml); wt_j = 1.-wt_i; // this is consistent with our second-order face location //
        conduction_wt = wt_i*local.CosmicRayDiffusionCoeff + wt_j*SphP[j].CosmicRayDiffusionCoeff; // arithmetic mean
        //conduction_wt = 2.0 * (CosmicRayDiffusionCoeff * SphP[j].CosmicRayDiffusionCoeff) / (CosmicRayDiffusionCoeff + SphP[j].CosmicRayDiffusionCoeff); // geometric mean
        conduction_wt *= All.cf_a3inv / All.cf_atime; // based on units CosmicRayDiffusionCoeff is defined with, this makes it physical for a dimensionless quantity gradient below
        
        double cmag = 0.0;
#ifdef MAGNETIC
        double B_interface[3];
        /* should use the solution in the appropriate face of the Riemann problem for interface values */
        for(k=0;k<3;k++) {B_interface[k] = Riemann_out.Face_B[k];}
        double B_interface_mag = 0.0;
        double B_interface_dot_grad_T = 0.0;
        for(k=0;k<3;k++)
        {
            B_interface_dot_grad_T += B_interface[k] * (wt_i*local.Gradients.CosmicRayPressure[k]
                                                        + wt_j*SphP[j].Gradients.CosmicRayPressure[k]);
            B_interface_mag += B_interface[k] * B_interface[k];
        }
        if(B_interface_mag > 0)
        {
            for(k=0;k<3;k++)
            {
                cmag += B_interface[k] * Face_Area_Vec[k];
            }
            cmag *= B_interface_dot_grad_T / B_interface_mag;
        } else {
            /* no magnetic field; use isotropic conduction equation */
            for(k=0;k<3;k++)
            {
                cmag += Face_Area_Vec[k] * (wt_i*local.Gradients.CosmicRayPressure[k]
                                            + wt_j*SphP[j].Gradients.CosmicRayPressure[k]);
            }
        }
#else
        for(k=0;k<3;k++)
        {
            cmag += Face_Area_Vec[k] * (wt_i*local.Gradients.CosmicRayPressure[k]
                                        + wt_j*SphP[j].Gradients.CosmicRayPressure[k]);
        }
#endif
        /* slope-limiter to ensure cosmic rays always flows from high pressure to low */
        double c_max = 0.0;
        for(k=0;k<3;k++) {c_max += Face_Area_Vec[k] * kernel.dp[k];}
        c_max *= rinv*rinv;
        double d_cr_p = local.CosmicRayPressure-CosmicRayPressure_j;
        cmag = MINMOD(MINMOD(MINMOD(cmag , c_max*d_cr_p), fabs(c_max)*d_cr_p) , 2.0*Face_Area_Norm*d_cr_p*rinv);
        /* now multiply through the coefficient to get the actual flux */
        cmag *= -conduction_wt;
        
        /* follow that with a flux limiter as well */
        conduction_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double CR_egy_i = local.CosmicRayPressure*V_i / GAMMA_COSMICRAY_MINUS1; // (E_cr = Volume * (Pressure/(GAMMA_CR-1)))
            double CR_egy_j = CosmicRayPressure_j*V_j / GAMMA_COSMICRAY_MINUS1;
            double du_ij_cond = 0.5*DMIN(DMIN(0.5*fabs(CR_egy_i-CR_egy_j),CR_egy_i),CR_egy_j);
            if(fabs(conduction_wt)>du_ij_cond) {conduction_wt *= du_ij_cond/fabs(conduction_wt);}
            Fluxes.CosmicRayPressure += conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
        
#endif // end of SPH/NOT SPH check
        
    } // close check that CosmicRayDiffusionCoeff and particle masses are positive
#endif
}
