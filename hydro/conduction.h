/* --------------------------------------------------------------------------------- */
/* ... real conduction evaluation ... */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    if((local.Kappa_Conduction>0)&&(SphP[j].Kappa_Conduction>0))
    {
//#ifdef HYDRO_SPH_TEMPORARY (PFH: in progress)
        /* First, we have the usual SPH anisotropic conduction definition */
        kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);
        double conduction_wt = P[j].Mass * kernel.dwk_ij / (kernel.r * local.Density * SphP[j].Density) * All.cf_atime; // physical units //
        double du_ij_cond = kernel.spec_egy_u_i - Particle_Internal_energy_i(j);
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
/*
#else
        // NOT SPH: Now we use the more accurate finite-volume formulation, with the effective faces we have already calculated //

        Face_Area_Vec[k] = ;

        double conduction_wt = 0.0;
#ifdef MAGNETIC
        double conduction_wt_bT = 0.0, conduction_wt_bA = 0.0, conduction_b2 = 0.0;
        for(k=0;k<3;k++)
        {
            // use B from the final Riemann state (x,y,z); with face area, finally gradients //
            conduction_wt_bT += B_interface[k] * 0.5 * (local.Gradients.T[k] + SphP[j].Gradients.T[k]);
            conduction_wt_bA += B_interface[k] * Face_Area_Vec[k]; // ??? should be B_interface //
            conduction_b2 += B_interface[k] * B_interface[k];
        }
        conduction_wt = conduction_wt_bT * conduction_wt_bA / conduction_b2;
#else
        for(k=0;k<3;k++)
        {
            conduction_wt += Face_Area_Vec[k] * 0.5 * (local.Gradients.T[k] + SphP[j].Gradients.T[k]);
        }
        conduction_wt /= NUMDIMS;
#endif
        conduction_wt *= 2.0 * local.Kappa_Conduction*SphP[j].Kappa_Conduction/(local.Kappa_Conduction + SphP[j].Kappa_Conduction);
#endif // end of SPH/NOT SPH check
 */
        
        
        // compute actual change that will occur this timestep //
        conduction_wt *= 0.5 * dt_hydrostep; // all in physical units //
        if(fabs(conduction_wt) > 0)
        {
            // enforce a limiter for stability (to prevent artificial oscillations) //
            du_ij_cond = 0.005*DMIN(DMIN(fabs(du_ij_cond),kernel.spec_egy_u_i),Particle_Internal_energy_i(j));
            if(fabs(conduction_wt)>du_ij_cond) conduction_wt *= du_ij_cond/fabs(conduction_wt);
            // now apply time rate of change to particle 'i'
            Fluxes.p += local.Mass * conduction_wt / dt_hydrostep;
        } // if(conduction_wt > 0)
    }

}
