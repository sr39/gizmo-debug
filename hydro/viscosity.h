/* --------------------------------------------------------------------------------- */
/* ... real anisotropic viscosity evaluation ...
 *
 * For SPH, this relys on the (noisy and zeroth-order inconsistent) SPH second-derivative
 *  operator. So a large kernel is especially useful to minimize the systematic errors.
 *  For MFM/MFV methods, the consistent finite-volume formulation is used.
 *  In either case, since we solve the viscous equations explicitly, a stronger timestep
 *  restriction is necessary (since the equations are not strictly hyperbolic); this is in timestep.c
 *
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
    if( (((local.Eta_ShearViscosity>0)&&(SphP[j].Eta_ShearViscosity>0)) ||
         ((local.Zeta_BulkViscosity>0)&&(SphP[j].Zeta_BulkViscosity>0))) &&
       ((local.Mass>0)&&(P[j].Mass>0)) )
    {
        int k_v;
        double v_interface[3];
        double cmag[3];
        double wt_i,wt_j;
        wt_i = wt_j = 0.5;
        for(k=0;k<3;k++) {v_interface[k] = wt_i*local.Vel[k] + wt_j*VelPred_j[k];} // should use interface solution from Riemann problem //
        
        // use a geometric average, since we want to weight the smaller of the two coefficients //
        double eta = 0.5 * (local.Eta_ShearViscosity + SphP[j].Eta_ShearViscosity);
        if(eta > 0) {eta = local.Eta_ShearViscosity * SphP[j].Eta_ShearViscosity / eta * All.cf_a2inv;} else {eta = 0;} // also converts to physical units
        double zeta = 0.5 * (local.Zeta_BulkViscosity + SphP[j].Zeta_BulkViscosity);
        if(zeta > 0) {zeta = local.Zeta_BulkViscosity * SphP[j].Zeta_BulkViscosity / zeta * All.cf_a2inv;} else {zeta = 0;} // also converts to physical units
        double viscous_wt_physical = DMAX(eta,zeta) / All.cf_a2inv;
        // we need a minus sign at some point; its handy to just include it now in the weights //
        wt_i *= -1.; wt_j *= -1.;
        int do_non_mhd = 1;
        
#ifdef HYDRO_SPH
        /* SPH anisotropic viscosity equations */
        double Face_Area_Norm = local.Mass * P[j].Mass * fabs(kernel.dwk_i+kernel.dwk_j) / (local.Density * SphP[j].Density);
        
#ifdef MAGNETIC
        double BPred[3],Bproj=0, BMag=0;
        for(k=0;k<3;k++)
        {
            BPred[k] = local.BPred[k] + BPred_j[k];
            Bproj += BPred[k] * kernel.dp[k];
            BMag += BPred[k] * BPred[k];
        }
        if(BMag > 0)
        {
            /* Braginskii formulation of leading-order anisotropic viscosity */
            do_non_mhd = 0; // only do the MHD loop here! //
            BMag = sqrt(BMag);
            for(k=0;k<3;k++) {BPred[k] /= BMag;}
            Bproj /= BMag;
            double rhs = 0;
            for(k_v=0;k_v<3;k_v++)
            {
                double tmp;
                double tmp_kv = BPred[k_v]*BPred[k_v] - 1./3.;
                double b_tensor_dot_face = 0;
                // note that because of the symmetry of the tensors here, the order of k,k_v doesn't matter //
                for(k=0;k<3;k++)
                {
                    if(k==k_v) {tmp=tmp_kv;} else {tmp=BPred[k_v]*BPred[k];}
                    rhs += tmp * kernel.dp[k_v] * kernel.dv[k];
                }
            }
            rhs /= (kernel.r * kernel.r * kernel.r);
            for(k_v=0;k_v<3;k_v++) {cmag[k_v] = Face_Area_Norm * eta * (3.*BPred[k_v]*Bproj - kernel.dp[k_v]) * rhs;}
        }
#endif
        /* standard Navier-Stokes equations: viscosity is decomposed into the shear and bulk viscosity terms */
        if(do_non_mhd==1)
        {
            double divv=0;
            for(k=0;k<3;k++) {divv += kernel.dv[k] * kernel.dp[k];}
            divv *= (zeta - eta*(2./3.));
            
            double Pxx = 2.*eta*kernel.dv[0]*kernel.dp[0] + divv;
            double Pyy = 2.*eta*kernel.dv[1]*kernel.dp[1] + divv;
            double Pzz = 2.*eta*kernel.dv[2]*kernel.dp[2] + divv;
            double Pxy = eta * (kernel.dv[0]*kernel.dp[1] + kernel.dv[1]*kernel.dp[0]);
            double Pxz = eta * (kernel.dv[0]*kernel.dp[2] + kernel.dv[2]*kernel.dp[0]);
            double Pyz = eta * (kernel.dv[1]*kernel.dp[2] + kernel.dv[2]*kernel.dp[1]);
            
            double r3inv = Face_Area_Norm / (kernel.r * kernel.r * kernel.r);
            cmag[0] = (Pxx*kernel.dp[0] + Pxy*kernel.dp[1] + Pxz*kernel.dp[2]) * r3inv;
            cmag[1] = (Pxy*kernel.dp[0] + Pyy*kernel.dp[1] + Pyz*kernel.dp[2]) * r3inv;
            cmag[2] = (Pxz*kernel.dp[0] + Pyz*kernel.dp[1] + Pzz*kernel.dp[2]) * r3inv;
        }
        for(k_v=0;k_v<3;k_v++) {cmag[k_v]*=-1;} // sign of the above is flipped
        
        
#else // not SPH
        
        
        
        
#ifdef MAGNETIC
        /* should use the solution in the appropriate face of the Riemann problem for interface values */
        double B_interface[3],B_interface_mag2=0;
        for(k=0;k<3;k++)
        {
            B_interface[k] = Riemann_out.Face_B[k];
            B_interface_mag2 += B_interface[k]*B_interface[k];
        }
        double bhat_dot_gradvhat = 1;
        double bhat_dot_gradvhat_direct = 1;
        double grad_v_mag = 0;
        double Bi_proj = 1;
        if(B_interface_mag2 >= 0)
        {
            /* Braginskii formulation of leading-order anisotropic viscosity */
            do_non_mhd = 0;
            double rhs = 0, rhs_direct = 0;
            double one_third = 1./3. * B_interface_mag2;
            Bi_proj = (kernel.dp[0]*B_interface[0]+kernel.dp[1]*B_interface[1]+kernel.dp[2]*B_interface[2]) / B_interface_mag2;
            for(k_v=0;k_v<3;k_v++)
            {
                double tmp;
                double tmp_kv = B_interface[k_v]*B_interface[k_v] - one_third;
                double b_tensor_dot_face = 0;
                // note that because of the symmetry of the tensors here, the order of k,k_v doesn't matter //
                for(k=0;k<3;k++)
                {
                    double grad_v = wt_i*local.Gradients.Velocity[k_v][k]+wt_j*SphP[j].Gradients.Velocity[k_v][k];
                    grad_v = MINMOD_G( grad_v, -kernel.dv[k_v] * kernel.dp[k] * rinv*rinv);
                    
                    if(k==k_v) {tmp=tmp_kv;} else {tmp=B_interface[k_v]*B_interface[k];}
                    rhs += tmp * grad_v;
                    rhs_direct += tmp * kernel.dv[k] * kernel.dp[k_v];
                    grad_v_mag += grad_v*grad_v;
                    b_tensor_dot_face += tmp * Face_Area_Vec[k];
                }
                cmag[k_v] = b_tensor_dot_face / B_interface_mag2;
            }
            rhs /= B_interface_mag2;
            grad_v_mag = sqrt(grad_v_mag);
            bhat_dot_gradvhat = rhs / grad_v_mag;
            bhat_dot_gradvhat_direct = 3.*rinv*rinv * (rhs_direct / B_interface_mag2);
            /* ok now just multipy the scalar contraction of the B tensor and shear tensor to get the fluxes */
            for(k_v=0;k_v<3;k_v++) {cmag[k_v] *= 3.*eta*rhs;}
        }
#endif
        
        /* standard Navier-Stokes equations: viscosity is decomposed into the shear and bulk viscosity terms */
        double cmag_dir[3];
        if(do_non_mhd==1)
        {
            double divv_i=0,divv_j=0;
            for(k=0;k<3;k++)
            {
                divv_i += local.Gradients.Velocity[k][k];
                divv_j += SphP[j].Gradients.Velocity[k][k];
            }
            double divv = (wt_i*divv_i + wt_j*divv_j) * (zeta - eta*(2./3.));
            wt_i*=eta; wt_j*=eta;
            
            double Pxx = 2*(wt_i*local.Gradients.Velocity[0][0]+wt_j*SphP[j].Gradients.Velocity[0][0]) + divv;
            double Pyy = 2*(wt_i*local.Gradients.Velocity[1][1]+wt_j*SphP[j].Gradients.Velocity[1][1]) + divv;
            double Pzz = 2*(wt_i*local.Gradients.Velocity[2][2]+wt_j*SphP[j].Gradients.Velocity[2][2]) + divv;
            double Pxy = wt_i*(local.Gradients.Velocity[0][1]+local.Gradients.Velocity[1][0]) +
            wt_j*(SphP[j].Gradients.Velocity[0][1]+SphP[j].Gradients.Velocity[1][0]);
            double Pxz = wt_i*(local.Gradients.Velocity[0][2]+local.Gradients.Velocity[2][0]) +
            wt_j*(SphP[j].Gradients.Velocity[0][2]+SphP[j].Gradients.Velocity[2][0]);
            double Pyz = wt_i*(local.Gradients.Velocity[1][2]+local.Gradients.Velocity[2][1]) +
            wt_j*(SphP[j].Gradients.Velocity[1][2]+SphP[j].Gradients.Velocity[2][1]);
            
            cmag[0] = Pxx*Face_Area_Vec[0] + Pxy*Face_Area_Vec[1] + Pxz*Face_Area_Vec[2];
            cmag[1] = Pxy*Face_Area_Vec[0] + Pyy*Face_Area_Vec[1] + Pyz*Face_Area_Vec[2];
            cmag[2] = Pxz*Face_Area_Vec[0] + Pyz*Face_Area_Vec[1] + Pzz*Face_Area_Vec[2];
            
            double dv_dir = (zeta-eta*2./3.)*(kernel.dv[0]*kernel.dp[0]+kernel.dv[1]*kernel.dp[1]+kernel.dv[2]*kernel.dp[2]);
            double Pxx_direct = eta*2.*kernel.dv[0]*kernel.dp[0] + dv_dir;
            double Pyy_direct = eta*2.*kernel.dv[1]*kernel.dp[1] + dv_dir;
            double Pzz_direct = eta*2.*kernel.dv[2]*kernel.dp[2] + dv_dir;
            double Pxy_direct = eta*(kernel.dv[0]*kernel.dp[1]+kernel.dv[1]*kernel.dp[0]);
            double Pxz_direct = eta*(kernel.dv[0]*kernel.dp[2]+kernel.dv[2]*kernel.dp[0]);
            double Pyz_direct = eta*(kernel.dv[2]*kernel.dp[1]+kernel.dv[1]*kernel.dp[2]);
            cmag_dir[0] = Pxx_direct*kernel.dp[0] + Pxy_direct*kernel.dp[1] + Pxz_direct*kernel.dp[2];
            cmag_dir[1] = Pxy_direct*kernel.dp[0] + Pyy_direct*kernel.dp[1] + Pyz_direct*kernel.dp[2];
            cmag_dir[2] = Pxz_direct*kernel.dp[0] + Pyz_direct*kernel.dp[1] + Pzz_direct*kernel.dp[2];
            
            cmag_dir[0] = (Pxx_direct*Face_Area_Vec[0] + Pxy_direct*Face_Area_Vec[1] + Pxz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
            cmag_dir[1] = (Pxy_direct*Face_Area_Vec[0] + Pyy_direct*Face_Area_Vec[1] + Pyz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
            cmag_dir[2] = (Pxz_direct*Face_Area_Vec[0] + Pyz_direct*Face_Area_Vec[1] + Pzz_direct*Face_Area_Vec[2])/Face_Area_Norm*kernel.r;
        }
        
        
        /* slope-limit this to be sure that viscosity always acts in the proper direction when there is local noise */
        double rho_i = local.Density*All.cf_a3inv, rho_j = SphP[j].Density*All.cf_a3inv, rho_ij=0.5*(rho_i+rho_j);
        for(k_v=0;k_v<3;k_v++)
        {
#ifdef MAGNETIC
            double dv_visc = 3. * (B_interface[k_v]*Bi_proj - kernel.dp[k_v]/3.) * bhat_dot_gradvhat_direct;
            double b_hll_eff = DMAX(DMIN(1. , 3.*bhat_dot_gradvhat*bhat_dot_gradvhat) , 0.01);
            double thold_ptot_hll = 0.1 * exp(-2. * grad_v_mag * kernel.r / (1.e-30 + fabs(kernel.dv[k_v])));
#else
            double b_hll_eff=1;
            double dv_visc = 0.5*kernel.dv[k_v];
            dv_visc = cmag_dir[k_v] * rinv*rinv / DMAX(eta,zeta);
            double thold_ptot_hll = 0.03;
#endif
            /* obtain HLL correction terms for Reimann problem solution */
            double hll_tmp = rho_ij * HLL_correction(dv_visc,-dv_visc,rho_ij,viscous_wt_physical) / All.cf_atime;
            double fluxlimiter_absnorm = -DMAX(eta,zeta) * sqrt(b_hll_eff) * Face_Area_Norm * dv_visc*rinv;
            double ptot = DMIN(local.Mass,P[j].Mass)*sqrt(kernel.dv[0]*kernel.dv[0]+
                                                          kernel.dv[1]*kernel.dv[1]+
                                                          kernel.dv[2]*kernel.dv[2]) / (1.e-37 + dt_hydrostep);
            double thold_hll = 0.1*ptot;
            if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll/fabs(hll_tmp);}
            hll_tmp *= b_hll_eff;
            if(cmag[k_v] * hll_tmp <= 0)
            {
                thold_hll = b_hll_eff * thold_ptot_hll * ptot;
                if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll / fabs(hll_tmp);}
            } else {
                thold_hll = b_hll_eff * DMAX(fabs(0.5*cmag[k_v]) , 0.3*thold_ptot_hll*ptot);
                if(fabs(hll_tmp) > thold_hll) {hll_tmp *= thold_hll / fabs(hll_tmp);}
            }
            double cmag_corr = cmag[k_v] + hll_tmp;
            cmag[k_v] = MINMOD(cmag[k_v], cmag_corr);
            if((fluxlimiter_absnorm*cmag[k_v] < 0) && (fabs(fluxlimiter_absnorm) > fabs(cmag[k_v]))) {cmag[k_v] = 0;}
        }
        
        double v_dot_dv=0; for(k=0;k<3;k++) {v_dot_dv += kernel.dv[k] * cmag[k];}
        if(v_dot_dv>0)
        {
            for(k=0;k<3;k++) {cmag[k] = 0;}
        } else {
            double KE_com=0; for(k=0;k<3;k++) {KE_com += kernel.dv[k]*kernel.dv[k];}
            KE_com *= 0.25 * (local.Mass + P[j].Mass);
            double dKE_q = fabs(v_dot_dv) * dt_hydrostep / (1.e-40 + KE_com); // COSMO UNITS!???????
            double threshold_tmp = 1.0;
            double lim_corr=1; if(dKE_q > threshold_tmp) {lim_corr = threshold_tmp/dKE_q;}
            for(k=0;k<3;k++) {cmag[k] *= lim_corr;}
        }
        
#endif // end of SPH/NOT SPH check
        
        /* now add a flux-limiter to prevent overshoot (even when the directions are correct) */
        double cmag_E = cmag[0]*v_interface[0] + cmag[1]*v_interface[1] + cmag[2]*v_interface[2];
        if(dt_hydrostep > 0)
        {
            double cmag_lim = 0.5 * (v_interface[0]*v_interface[0]+v_interface[1]*v_interface[1]+v_interface[2]*v_interface[2]) / dt_hydrostep;
            if(fabs(cmag_E) > cmag_lim)
            {
                double corr_visc = cmag_lim / fabs(cmag_E);
                cmag[0]*=corr_visc; cmag[1]*=corr_visc; cmag[2]*=corr_visc; cmag_E*=corr_visc;
            }
        }
        /* ok now we can finally add this to the numerical fluxes */
        for(k=0;k<3;k++) {Fluxes.v[k] += cmag[k];}
        Fluxes.p += cmag_E;
        
    } // close check that kappa and particle masses are positive
}
