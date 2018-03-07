/* quantum pressure-tensor computation to calculate the relevant fluxes between faces of particles, within the gravity routine */
#ifdef DM_FUZZY
if((ptype==1) && (ptype_sec==1)) // only acts between DM particles of type 1 (can easily change if desired)
{
    /* since this isn't a super-expensive calculation, and we need to walk the gravity tree for both 'sides' anyways, we will effectively do this twice each timestep */
    int k,m; double Face_Area_Vec[3],dp[3],wk_i=0,wk_j=0,dummy,rho_i,rho_j,flux[3],vface_i=0,vface_j=0,Face_Area_Norm=0,dv[3]; h3_inv=h_inv*h_inv*h_inv; h_p3_inv=h_p_inv*h_p_inv*h_p_inv; u_p=r*h_p_inv; double h_eff = 0.5*(1./h_inv + 1./h_p_inv);
    dp[0]=pos_x-P[j0_sec_for_ags].Pos[0]; dp[1]=pos_y-P[j0_sec_for_ags].Pos[1]; dp[2]=pos_z-P[j0_sec_for_ags].Pos[2]; // recalculate positions (to make sure have in correct format needed below //
#ifdef BOX_PERIODIC  /* find the closest image in the given box size  */
    NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
    // calculate element effective volumes (for weights below) //
    double V_i=local_V_i, V_j = get_particle_volume_ags(j0_sec_for_ags);
    rho_i = pmass/V_i*All.cf_a3inv; rho_j = mass/V_j*All.cf_a3inv; // calculate densities (in physical units)
    // calculate kernel weight functions (for faces below) //
    if(u<1) {kernel_main(u, h3_inv, h3_inv*h_inv, &wk_i, &dummy, -1);} // get wk [these variables all pre-defined if here in loop] //
    if(u_p<1) {kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wk_j, &dummy, -1);} // likewise these should all be pre-calculated above, already [riding on adaptive gravsoft routines here] //
    // calculate effective faces and face velocity between elements //
    for(k=0;k<3;k++)
    {
        Face_Area_Vec[k] = (wk_i*V_i * (local_NV_T[k][0]*dp[0] + local_NV_T[k][1]*dp[1] + local_NV_T[k][2]*dp[2]) + wk_j*V_j * (P[j0_sec_for_ags].NV_T[k][0]*dp[0] + P[j0_sec_for_ags].NV_T[k][1]*dp[1] + P[j0_sec_for_ags].NV_T[k][2]*dp[2])) * All.cf_atime*All.cf_atime; // physical units
        Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k]; // physical units
        double vi = targetVel[k] / All.cf_atime; // physical units
        double vj = P[j0_sec_for_ags].Vel[k] / All.cf_atime; // physical units
        dv[k] = vi - vj; // physical units
        if(All.ComovingIntegrationOn) {dv[k] += All.cf_hubble_a * (dp[k]*All.cf_atime);} // add hubble flow in cosmological runs
        vi = vj + dv[k]; // reset to capture velocity difference correctly
        vface_i += vi * Face_Area_Vec[k]; // physical units
        vface_j += vj * Face_Area_Vec[k]; // physical units
    }
    Face_Area_Norm = sqrt(Face_Area_Norm); vface_i /= Face_Area_Norm; vface_j /= Face_Area_Norm;
    // convert everything needed below into physical units //
    double igrad[3], jgrad[3], i2grad[3][3], j2grad[3][3], fac_g = All.cf_a3inv/All.cf_atime, fac_g2 = All.cf_a3inv*All.cf_a2inv;
    for(k=0;k<3;k++)
    {
        dp[k] *= All.cf_atime;
        igrad[k] = fac_g * local_AGS_Gradients_Density[k];
        jgrad[k] = fac_g * P[j0_sec_for_ags].AGS_Gradients_Density[k];
        for(m=0;m<3;m++)
        {
            i2grad[k][m] = fac_g2 * local_AGS_Gradients2_Density[k][m];
            j2grad[k][m] = fac_g2 * P[j0_sec_for_ags].AGS_Gradients2_Density[k][m];
        }
    }
    double dt = targetdt_step * All.Timebase_interval / All.cf_hubble_a, m_mean = 0.5 * ( mass+pmass );
    double prev_acc = All.G*All.cf_a2inv * 0.5*(pmass*aold/All.ErrTolForceAcc + mass*P[j0_sec_for_ags].OldAcc);
    double HLLwt = (0.5*(wk_i/h3_inv + wk_j/h_p3_inv)) * (h_eff/r); HLLwt = 10.*HLLwt*HLLwt; // strong dissipation terms allowed for very-close particles, where second-derivative diverges, otherwise weak (no diffusion) //
    // actually compute the fluxes now, this is the key routine, below //
    do_dm_fuzzy_flux_computation(HLLwt, dt, m_mean, prev_acc, dp, dv, jgrad,  igrad, j2grad, i2grad, rho_j, rho_i, vface_j, vface_i, Face_Area_Vec, flux);
    double fac = 1. / (pmass * All.G * All.cf_a2inv); // 'flux' now holds dmomentum/dt in physical units; need to convert back to grav-acc routine acceleration units ~ m/r^2 (will be multipled by 'G' later) //
    acc_x += fac*flux[0]; acc_y += fac*flux[1]; acc_z += fac*flux[2]; // assign back to particles
} // master bracket (for variable protection
#endif
