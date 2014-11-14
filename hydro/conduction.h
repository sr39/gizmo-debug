/* --------------------------------------------------------------------------------- */
/* ... real conduction evaluation ... */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    if((local.Kappa_Conduction>0)&&(SphP[j].Kappa_Conduction>0))
    {
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
            Bpro2_i = local.BPred[0]*kernel.dx + local.BPred[1]*kernel.dy + local.BPred[2]*kernel.dz;
            Bpro2_i *= Bpro2_i / (kernel.b2_i * kernel.r*kernel.r);
        }
        if(kernel.b2_j>0)
        {
            Bpro2_j = BPred_j[0]*kernel.dx + BPred_j[1]*kernel.dy + BPred_j[2]*kernel.dz;
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
        //    equally valid (slightly more accurate, but less stable) is to use arithmetic sum: = (local.Kappa_Conduction + SphP[j].Kappa_Conduction)
#endif
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
