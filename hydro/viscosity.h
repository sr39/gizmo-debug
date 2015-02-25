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
#ifdef HYDRO_SPH
        /* fill in SPH viscosity equations; see the cautions above, though! */
        
        
#else
        int k_v;
        double v_interface[3];
        double cmag[3];
        double wt_i,wt_j;
        wt_i = wt_j = 0.5;
        //wt_i = PPP[j].Hsml / (PPP[j].Hsml + local.Hsml); wt_j = 1.-wt_i; // this is consistent with our second-order face location //
        for(k=0;k<3;k++) {v_interface[k] = wt_i*local.Vel[k] + wt_j*VelPred_j[k];} // should use interface solution from Riemann problem //
        
        // estimate the interface coefficients with a simple arithmetic average //
        double eta  = wt_i*local.Eta_ShearViscosity + wt_j*SphP[j].Eta_ShearViscosity;
        double zeta = wt_i*local.Zeta_BulkViscosity + wt_j*SphP[j].Zeta_BulkViscosity;
        // we need a minus sign at some point; its handy to just include it now in the weights //
        wt_i *= -1.; wt_j *= -1.;
        
#ifdef MAGNETIC
        /* should use the solution in the appropriate face of the Riemann problem for interface values */
        double B_interface[3],B_interface_mag=0;
        for(k=0;k<3;k++)
        {
            B_interface[k] = Riemann_out.Face_B[k];
            B_interface_mag += B_interface[k]*B_interface[k];
        }
        
        if(B_interface_mag <= 0)
        {
            /* no magnetic field: use isotropic viscosity */
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
            
        } else {
            /* Braginskii formulation of leading-order anisotropic viscosity */
            double rhs = 0;
            double one_third = 1./3. * B_interface_mag;
            for(k_v=0;k_v<3;k_v++)
            {
                double tmp;
                double tmp_kv = B_interface[k_v]*B_interface[k_v] - one_third;
                double b_tensor_dot_face = 0;
                // note that because of the symmetry of the tensors here, the order of k,k_v doesn't matter //
                for(k=0;k<3;k++)
                {
                    if(k==k_v) {tmp=tmp_kv;} else {tmp=B_interface[k_v]*B_interface[k];}
                    rhs += tmp * (wt_i*local.Gradients.Velocity[k_v][k]+wt_j*SphP[j].Gradients.Velocity[k_v][k]);
                    b_tensor_dot_face += tmp * Face_Area_Vec[k];
                }
                cmag[k_v] = b_tensor_dot_face;
            }
            /* ok now just multipy the scalar contraction of the B tensor and shear tensor to get the fluxes */
            rhs *= 3.0 * eta / (B_interface_mag*B_interface_mag);
            for(k_v=0;k_v<3;k_v++) {cmag[k_v] *= rhs;}
        }
#else
        /* standard Navier-Stokes equations: viscosity is decomposed into the shear and bulk viscosity terms */
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
#endif
        
        /* slope-limit this to be sure that viscosity always acts in the proper direction when there is local noise */
        double c_max_norm = -2.0 * Face_Area_Norm * DMAX(eta,zeta) * rinv;
        for(k_v=0;k_v<3;k_v++)
        {
            double c_max = c_max_norm * (local.Vel[k_v]-VelPred_j[k_v]); // inter-particle gradient times tolerance //
            cmag[k_v] = MINMOD(c_max,cmag[k_v]);
        }
        
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
        
#endif // end of SPH/NOT SPH check
        
    } // close check that kappa and particle masses are positive
}
