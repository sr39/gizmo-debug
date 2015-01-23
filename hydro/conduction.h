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
    if((local.Kappa_Conduction>0)&&(SphP[j].Kappa_Conduction>0)&&(local.Mass>0)&&(P[j].Mass>0))
    {
#ifdef HYDRO_SPH
        
        /* First, we have the usual SPH anisotropic conduction definition */
        kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);
        double conduction_wt = P[j].Mass * kernel.dwk_ij / (kernel.r * local.Density * SphP[j].Density) * All.cf_atime; // physical units //
        double du_ij_cond = kernel.spec_egy_u_i - SphP[j].InternalEnergyPred;
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
        conduction_wt *= 2.0 * local.Kappa_Conduction*Bpro2_i*SphP[j].Kappa_Conduction*Bpro2_j /
        (local.Kappa_Conduction*Bpro2_i + SphP[j].Kappa_Conduction*Bpro2_j);
#else
        conduction_wt *= 2.0 * local.Kappa_Conduction*SphP[j].Kappa_Conduction/(local.Kappa_Conduction + SphP[j].Kappa_Conduction);
        // this uses geometric-weighted kappa (as advocated by Cleary & Monaghan '99 for stability):
        //    equally valid (slightly more accurate, but less stable) is to use arithmetic mean: = (local.Kappa_Conduction + SphP[j].Kappa_Conduction)
#endif
        // compute actual change that will occur this timestep //
        conduction_wt *= dt_hydrostep; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a limiter for stability (to prevent artificial oscillations) //
            double du_ij_cond = 0.25*DMIN(DMIN(0.5*fabs(local.Mass*kernel.spec_egy_u_i-P[j].Mass*SphP[j].InternalEnergyPred),local.Mass*kernel.spec_egy_u_i),P[j].Mass*SphP[j].InternalEnergyPred);
            if(fabs(conduction_wt)>du_ij_cond) {conduction_wt *= du_ij_cond/fabs(conduction_wt);}
            // now apply time rate of change to particle 'i'
            Fluxes.p += conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
        
#else
        
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //
        
        double conduction_wt;
        double wt_i,wt_j;
        // wt_i = wt_j = 0.5;
        wt_i = PPP[j].Hsml / (PPP[j].Hsml + local.Hsml); wt_j = 1.-wt_i; // this is consistent with our second-order face location //
        conduction_wt = wt_i*local.Kappa_Conduction + wt_j*SphP[j].Kappa_Conduction; // arithmetic mean
        //conduction_wt = 2.0 * (local.Kappa_Conduction * SphP[j].Kappa_Conduction) / (local.Kappa_Conduction + SphP[j].Kappa_Conduction); // geometric mean
        conduction_wt *= All.cf_atime; // based on units TD_DiffCoeff is defined with, this makes it physical for a dimensionless quantity gradient below
        /* if we use -DIFFUSIVITIES-, we need a density here; if we use -CONDUCTIVITITIES-, no density */
        // conduction_wt *= Riemann_out.Face_Density;
        
        double cmag = 0.0;
        double c_max = 0.0;
#ifdef MAGNETIC
        double B_interface[3];
        /* should use the solution in the appropriate face of the Riemann problem for interface values */
        for(k=0;k<3;k++)
        {
            //B_interface[k] = 0.5 * (local.BPred[k] + BPred_j[k]) * All.cf_a2inv;
            B_interface[k] = Riemann_out.Face_B[k];
        }
        double B_interface_mag = 0.0;
        double B_interface_dot_grad_T = 0.0;
        for(k=0;k<3;k++)
        {
            B_interface_dot_grad_T += B_interface[k] * (wt_i*local.Gradients.InternalEnergy[k]
                                                        + wt_j*SphP[j].Gradients.InternalEnergy[k]);
            B_interface_mag += B_interface[k] * B_interface[k];
        }
        if(B_interface_mag > 0)
        {
            for(k=0;k<3;k++)
            {
                c_max += Face_Area_Vec[k] * kernel.dp[k];
                cmag += B_interface[k] * Face_Area_Vec[k];
            }
            cmag *= B_interface_dot_grad_T / B_interface_mag;
        } else {
            /* no magnetic field; use isotropic conduction equation */
            for(k=0;k<3;k++)
            {
                c_max+=Face_Area_Vec[k] * kernel.dp[k];
                cmag += Face_Area_Vec[k] * (wt_i*local.Gradients.InternalEnergy[k]
                                            + wt_j*SphP[j].Gradients.InternalEnergy[k]);
            }
        }
#else
        for(k=0;k<3;k++)
        {
            c_max += Face_Area_Vec[k] * kernel.dp[k];
            cmag += Face_Area_Vec[k] * (wt_i*local.Gradients.InternalEnergy[k]
                                        + wt_j*SphP[j].Gradients.InternalEnergy[k]);
        }
#endif
        /* slope-limiter to ensure heat always flows from hot to cold */
        c_max *= rinv*rinv;
        double du_cond = local.InternalEnergyPred-SphP[j].InternalEnergyPred;
        cmag = MINMOD(MINMOD(MINMOD(cmag , c_max*du_cond), fabs(c_max)*du_cond) , Face_Area_Norm*du_cond*rinv);
        //double c_max = 1.0 * Face_Area_Norm * (local.InternalEnergyPred-SphP[j].InternalEnergyPred) * rinv; // inter-particle gradient times tolerance //
        //cmag = MINMOD(c_max,cmag);
        
        /* now multiply through the coefficient to get the actual flux */
        cmag *= -conduction_wt;
        
        /* follow that with a fluxlimiter as well */
        conduction_wt = dt_hydrostep * cmag; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a flux limiter for stability (to prevent overshoot) //
            double du_ij_cond = 0.5*DMIN(DMIN(0.5*fabs(DMIN(local.Mass,P[j].Mass)*(local.InternalEnergyPred-SphP[j].InternalEnergyPred)),
                                              local.Mass*local.InternalEnergyPred),
                                         P[j].Mass*SphP[j].InternalEnergyPred);
            if(fabs(conduction_wt)>du_ij_cond) {conduction_wt *= du_ij_cond/fabs(conduction_wt);}
            Fluxes.p += conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
        
#endif // end of SPH/NOT SPH check
        
    } // close check that kappa and particle masses are positive
}
