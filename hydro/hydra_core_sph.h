/* --------------------------------------------------------------------------------- */
/* this is the sub-routine where we actually evaluate the SPH equations of motion */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    /* basic overhead variables and zero-ing fluxes for the computation */
    Fluxes.rho = Fluxes.p = Fluxes.v[0] = Fluxes.v[1] = Fluxes.v[2] = 0;
    double du_ij;
    kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);
    double vdotr2_phys = kernel.vdotr2;
    if(All.ComovingIntegrationOn) vdotr2_phys -= All.cf_hubble_a2 * r2;
    
    /* --------------------------------------------------------------------------------- */
    /* --------------------------------------------------------------------------------- */
    /* ... EQUATION OF MOTION (HYDRO) ... */
    /* --------------------------------------------------------------------------------- */
    /* --------------------------------------------------------------------------------- */
    double vi_dot_r,hfc,hfc_visc,hfc_i,hfc_j,hfc_dwk_i,hfc_dwk_j,hfc_egy=0;
    /* 'Standard' (Lagrangian) Density Formulation: the acceleration term is identical whether we use 'entropy' or 'energy' sph */
    /* (this step is the same in both 'Lagrangian' and 'traditional' SPH */

#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    /* Pressure-Energy and/or Pressure-Entropy form of SPH (using 'constant mass in kernel' h-constraint */
    /* -- note that, using appropriate definitions, both forms have an identical EOM in appearance here -- */
    double p_over_rho2_j = SphP[j].Pressure / (SphP[j].EgyWtDensity * SphP[j].EgyWtDensity);
    hfc_i = kernel.p_over_rho2_i * (SphP[j].InternalEnergyPred/local.InternalEnergyPred) *
        (1 + local.DhsmlHydroSumFactor / (P[j].Mass * SphP[j].InternalEnergyPred));
    hfc_j = p_over_rho2_j * (local.InternalEnergyPred/SphP[j].InternalEnergyPred) *
        (1 + SphP[j].DhsmlHydroSumFactor / (local.Mass * local.InternalEnergyPred));
#else
    /* Density-Entropy (or Density-Energy) formulation: x_tilde=1, x=mass */
    double p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
    hfc_i = kernel.p_over_rho2_i * (1 + local.DhsmlHydroSumFactor / P[j].Mass);
    hfc_j = p_over_rho2_j * (1 + SphP[j].DhsmlHydroSumFactor / local.Mass);
#endif
    
    hfc_egy = hfc_i; /* needed to follow the internal energy explicitly; note this is the same for any of the formulations above */
        
    /* use the traditional 'kernel derivative' (dwk) to compute derivatives */
    hfc_dwk_i = hfc_dwk_j = local.Mass * P[j].Mass / kernel.r;
    hfc_dwk_i *= kernel.dwk_i; /* grad-h terms have already been multiplied in here */
    hfc_dwk_j *= kernel.dwk_j;
    hfc = hfc_i*hfc_dwk_i + hfc_j*hfc_dwk_j;
        
    /* GASOLINE-like equation of motion: */
    /* hfc_dwk_i = 0.5 * (hfc_dwk_i + hfc_dwk_j); */
    /* hfc = (p_over_rho2_j + kernel.p_over_rho2_i) * hfc_dwk_i */
        
    /* RSPH equation of motion */
    /* hfc_dwk_i = 0.5 * (hfc_dwk_i + hfc_dwk_j); */
    /* hfc = (local.Pressure-SphP[j].Pressure)/(local.Density*SphP[j].Density) * hfc_dwk_i */
        
    Fluxes.v[0] += -hfc * kernel.dx; /* momentum */
    Fluxes.v[1] += -hfc * kernel.dy;
    Fluxes.v[2] += -hfc * kernel.dz;
    vi_dot_r = local.Vel[0]*kernel.dx + local.Vel[1]*kernel.dy + local.Vel[2]*kernel.dz;
    Fluxes.p += hfc_egy * hfc_dwk_i * vdotr2_phys - hfc * vi_dot_r; /* total energy */
    
    
    /* --------------------------------------------------------------------------------- */
    /* ... Magnetic evolution (ideal MHD terms all below) ... */
    /* --------------------------------------------------------------------------------- */
#ifdef MAGNETIC
    Fluxes.B[0] = Fluxes.B[1] = Fluxes.B[2] = 0;
    V_j = P[j].Mass / SphP[j].Density;
    double mj_r = P[j].Mass / kernel.r;
    double mf_i = kernel.mf_i * mj_r * kernel.dwk_i;
    double mf_j = kernel.mf_j * mj_r * kernel.dwk_j / (SphP[j].Density * SphP[j].Density);
    
    /* ---------------------------------------------------------------------------------
     * ... induction equation ...
     * (the SPH induction equation is inherently non-symmetric: the result for particle
     *  'a' and particle 'b' must each be separately computed. However, they do not rely on the 
     *   densities or smoothing lengths of the neighbor particles. Therefore they can and should 
     *   be computed in the density loop, not the hydro loop)
     * --------------------------------------------------------------------------------- */
    
#ifdef DIVBCLEANING_DEDNER
    /* --------------------------------------------------------------------------------- */
    /* ... update to the Dedner divergence-damping scalar field ... */
    /* --------------------------------------------------------------------------------- */
    /* PFH: this is the symmetric estimator of grad phi: used be used IFF div.dot.B is estimated by the 'direct difference' operator (recommended) */
    double phifac = -(mf_i*local.PhiPred + mf_j*PhiPred_j); // *All.cf_atime*All.cf_atime; // All.cf_atime NOT needed for convention [Phicode]=[Bcode][vcode]
    Fluxes.B[0] += phifac * kernel.dx;
    Fluxes.B[1] += phifac * kernel.dy;
    Fluxes.B[2] += phifac * kernel.dz;
    // GradPhi should have units of [Phicode]/[rcode] = [Bcode]*[vcode]/[rcode] = [DtB]=[Fluxes.B]
#endif // DIVBCLEANING_DEDNER
    
    /* --------------------------------------------------------------------------------- */
    /* ... magnetic acceleration (with correction/limiter term) ... */
    /* --------------------------------------------------------------------------------- */
    double dBx = local.BPred[0] - SphP[j].BPred[0];
    double dBy = local.BPred[1] - SphP[j].BPred[1];
    double dBz = local.BPred[2] - SphP[j].BPred[2];
    int k1;
    for(k = 0; k < 3; k++)
    {
        for(k1 = 0; k1 < 3; k1++)
            mm_j[k][k1] = BPred_j[k] * BPred_j[k1];
    }
    for(k = 0; k < 3; k++)
        mm_j[k][k] -= 0.5 * kernel.b2_j;
    
    /* need to be able to subtract component of force which owes to the non-zero value of div.B */
    Fluxes.B_normal_corrected = (local.BPred[0] * mf_i + BPred_j[0] * mf_j) * kernel.dx +
                                (local.BPred[1] * mf_i + BPred_j[1] * mf_j) * kernel.dy +
                                (local.BPred[2] * mf_i + BPred_j[2] * mf_j) * kernel.dz;
    for(k = 0; k < 3; k++)
    {
        /* momentum flux from MHD forces */
        magfluxv[k] = (mm_i[k][0] * mf_i + mm_j[k][0] * mf_j) * kernel.dx +
                          (mm_i[k][1] * mf_i + mm_j[k][1] * mf_j) * kernel.dy +
                          (mm_i[k][2] * mf_i + mm_j[k][2] * mf_j) * kernel.dz;
        Fluxes.v[k] += magfluxv[k];
    }
    
        
        
    /* --------------------------------------------------------------------------------- */
    /* ... Magnetic dissipation/diffusion terms (artificial resitivity evaluation) ... */
    /* --------------------------------------------------------------------------------- */
#ifdef MAGNETIC_DISSIPATION
    double mf_dissInd = local.Mass * mj_r * kernel.dwk_ij * kernel.rho_ij_inv * kernel.rho_ij_inv / fac_mu;
    double vsigb = 0.5 * sqrt(kernel.alfven2_i + kernel.alfven2_j);
#if defined(TRICCO_RESISTIVITY_SWITCH)
    double eta = 0.5 * (local.Balpha + SphP[j].Balpha) * vsigb * kernel.r;
#else
    double eta = All.ArtMagDispConst * vsigb * kernel.r;
#endif
    mf_dissInd *= eta;
    Fluxes.B[0] += mf_dissInd * dBx;
    Fluxes.B[1] += mf_dissInd * dBy;
    Fluxes.B[2] += mf_dissInd * dBz;
    Fluxes.p -= 0.5 * fac_magnetic_pressure * mf_dissInd * (dBx * dBx + dBy * dBy + dBz * dBz);
#endif
#endif /* end MAGNETIC */
        
    

    
    
    /* --------------------------------------------------------------------------------- */
    /* ... artificial viscosity evaluation ... */
    /* --------------------------------------------------------------------------------- */
    if(kernel.vdotr2 < 0) // no viscosity applied if particles are moving away from each other //
    {
        double c_ij = 0.5 * (kernel.sound_i + kernel.sound_j);
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        double BulkVisc_ij = 0.5 * (local.alpha + SphP[j].alpha_limiter * SphP[j].alpha);
        double mu_ij = fac_mu * kernel.vdotr2 / kernel.r;
        double visc = -BulkVisc_ij * mu_ij * (c_ij - mu_ij) * kernel.rho_ij_inv; /* this method should use beta/alpha=1 */
#else
        double BulkVisc_ij = 0.5 * All.ArtBulkViscConst * (local.alpha + SphP[j].alpha_limiter);
        double h_ij = KERNEL_CORE_SIZE * 0.5 * (kernel.h_i + kernel.h_j);
        double mu_ij = fac_mu * h_ij * kernel.vdotr2 / (r2 + 0.0001 * h_ij * h_ij); /* note: this is negative! */
        double visc = -BulkVisc_ij * mu_ij * (c_ij - 2*mu_ij) * kernel.rho_ij_inv; /* this method should use beta/alpha=2 */
#endif
#ifndef NOVISCOSITYLIMITER
        double dt = 2 * IMAX(local.Timestep, (P[j].TimeBin ? (1 << P[j].TimeBin) : 0)) * All.Timebase_interval;
        if(dt > 0 && kernel.dwk_ij < 0)
            visc = DMIN(visc, 0.5 * fac_vsic_fix * kernel.vdotr2 / ((local.Mass + P[j].Mass) * kernel.dwk_ij * kernel.r * dt));
#endif
        hfc_visc = -local.Mass * P[j].Mass * visc * kernel.dwk_ij / kernel.r;
        Fluxes.v[0] += hfc_visc * kernel.dx; /* this is momentum */
        Fluxes.v[1] += hfc_visc * kernel.dy;
        Fluxes.v[2] += hfc_visc * kernel.dz;
        Fluxes.p += hfc_visc * (vi_dot_r - 0.5*vdotr2_phys); /* remember, this is -total- energy now */
    } // kernel.vdotr2 < 0 -- triggers artificial viscosity
    
    
    
    /* --------------------------------------------------------------------------------- */
    /* ... artificial conductivity (thermal diffusion) evaluation ... */
    /* --------------------------------------------------------------------------------- */
#ifdef SPHAV_ARTIFICIAL_CONDUCTIVITY
#ifdef BP_REAL_CRs
    double vsigu = sqrt(fabs(local.Pressure - SphP[j].Pressure - local.CRpPressure + SphP[j].CRpPressure) * kernel.rho_ij_inv);
    double u_i = (local.Pressure - local.CRpPressure) / (GAMMA_MINUS1 * local.Density);
    double u_j = (SphP[j].Pressure - SphP[j].CRpPressure) / (GAMMA_MINUS1 * SphP[j].Density);
    Fluxes.p += local.Mass * P[j].Mass * All.ArtCondConstant * vsigu * (u_i - u_j) * kernel.rho_ij_inv * kernel.dwk_ij;
#else
    double vsigu = (kernel.sound_i + kernel.sound_j - 3 * fac_mu * kernel.vdotr2 / kernel.r) / fac_mu; // want in code velocity units
    if(vsigu > 0) // implicitly sets vsig=0 if 3*w_ij > (c_i+c_j)
    {
        vsigu *= fabs(local.Pressure - SphP[j].Pressure)/(local.Pressure + SphP[j].Pressure);
        du_ij = kernel.spec_egy_u_i - Particle_Internal_energy_i(j);
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        du_ij *= 0.5 * (local.alpha + SphP[j].alpha_limiter * SphP[j].alpha); // in this case, All.ArtCondConstant is just a multiplier -relative- to art. visc.
#endif
        Fluxes.p += local.Mass * All.ArtCondConstant * P[j].Mass * kernel.rho_ij_inv * vsigu * du_ij * kernel.dwk_ij;
    }
#endif
#endif // SPHAV_ARTIFICIAL_CONDUCTIVITY

    
    
    /* --------------------------------------------------------------------------------- */
    /* ... cosmic ray conductivity ... */
    /* --------------------------------------------------------------------------------- */
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
    int Nbin;
    double vsigu_cr = sqrt(fabs(local.CRpPressure - SphP[j].CRpPressure) * kernel.rho_ij_inv);
    for( Nbin = 0; Nbin < BP_REAL_CRs; Nbin++ )
    {
        out.DtCRpE[Nbin] += P[j].Mass * All.CRsArtCondConstant * vsigu_cr * (local.CRpE[Nbin] - SphP[j].CRpE[Nbin]) * kernel.rho_ij_inv * kernel.dwk_ij;
        out.DtCRpN[Nbin] += P[j].Mass * All.CRsArtCondConstant * vsigu_cr * (local.CRpN[Nbin] - SphP[j].CRpN[Nbin]) * kernel.rho_ij_inv * kernel.dwk_ij;
    }
#endif
    
    
    /* --------------------------------------------------------------------------------- */
    /* convert everything to PHYSICAL units! */
    /* --------------------------------------------------------------------------------- */
    Fluxes.p *= All.cf_afac2 / All.cf_atime;
    for(k=0;k<3;k++)
        Fluxes.v[k] *= All.cf_afac2;
}
