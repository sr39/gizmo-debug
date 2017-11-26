#ifdef CBE_INTEGRATOR

/* --------------------------------------------------------------------------------- */
/* this is the kernel where the core of the CBE flux computation is performed,
 *  to calculate the relevant fluxes between faces of partciels, within the gravity routine
/* --------------------------------------------------------------------------------- */
{
    int k, j=no; // secondary has index 'no' for gravity loop
    double Face_Area_Norm=0, Face_Area_Vec[3]={0}, dp[3]={0}, wk_i=0, wk_j=0, dummy=0;
    double rho_i=0, rho_j=0, psi_i=0, psi_j=0, vface[3]={0};
    double V_i, V_j, L_j = Get_Particle_Size(j);
    // recalculate positions (to make sure have in correct format needed below //
    dp[0] = pos_x - P[j].Pos[0];
    dp[1] = pos_y - P[j].Pos[1];
    dp[2] = pos_z - P[j].Pos[2];
#ifdef BOX_PERIODIC  /* find the closest image in the given box size  */
    NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
    // calculate element effective volumes (for weights below) //
#if (NUMDIMS==1)
    V_j = L_j;
#elif (NUMDIMS==2)
    V_j = L_j*L_j;
#else
    V_j = L_j*L_j*L_j;
#endif
    rho_i = pmass/V_i; rho_j = mass/V_j; // calculate densities
    // calculate position between particles (psi) term for interpolation
    psi_i = 1. / (1. + h_p_inv/h_inv); // fraction of weight to i
    psi_j = 1-psi_i; // fraction of weight to j
    // rho never appears without a weight factor prefacing it, so just multiply here to get them
    rho_i *= psi_i; rho_j *= psi_j;
    // calculate kernel weight functions (for weights below) //
    if(u<1) {kernel_main(u, h3_inv, h3_inv*h_inv, &wk_i, &dummy, -1);} // get wk [these variables all pre-defined if here in loop] //
    if(u_p<1) {kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wk_j, &dummy, 1);} // likewise these should all be pre-calculated above, already [riding on adaptive gravsoft routines here] //
    // calculate effective faces and face velocity between elements //
    for(k=0;k<3;k++)
    {
        Face_Area_Vec[k] = wk_i*V_i * (local.NV_T[k][0]*dp[0] + local.NV_T[k][1]*dp[1] + local.NV_T[k][2]*dp[2])
                         + wk_j*V_j * (P[j].NV_T[k][0]*dp[0] + P[j].NV_T[k][1]*dp[1] + P[j].NV_T[k][2]*dp[2]);
        Face_Area_Vec[k] *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
        Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k];
        vface[k] = (rho_i*local.Vel[k] + rho_j*P[i].Vel[k]) / (rho_i+rho_j);
    }
    
    // OK, now we've done the easy bit -- we're ready for the actual flux computations //
    
    // first loop over pairs to determine closest a-to-b, closest b-to-a //
    for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
    {
    }

    // then loop over each basis for each particle and compute fluxes //
    for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
    {
        double dm = -
        for(k=1;k<10;k++) {}
    }

    
}
#endif
