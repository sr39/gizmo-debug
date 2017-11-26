#ifdef CBE_INTEGRATOR

/* --------------------------------------------------------------------------------- */
/* this is the kernel where the core of the CBE flux computation is performed,
 *  to calculate the relevant fluxes between faces of partciels, within the gravity routine
/* --------------------------------------------------------------------------------- */
{
/* since we compute both 'sides' fluxes together (and it is an expensive operation), we
    only need to do it once -- we need to figure out which particle is active on the smaller timestep, and
    compute for it only */
    int do_cbe_calculation = 1; // default to doing the calculation
    int k, j=no; // secondary has index 'no' for gravity loop
    int TimeStep_J = (P[j].TimeBin ? (1 << P[j].TimeBin) : 0);
    double pos_i[3]={pos_x,pos_y,pos_z};
    if(local.Timestep > TimeStep_J) do_cbe_calculation = 0; /* compute from particle with smaller timestep */
    if(local.Timestep == TimeStep_J) // same timestep, randomly break degeneracy with positions
    {
        k=0; if(pos_i[k] == P[j].Pos[k]) {k++; if(pos_i[k] == P[j].Pos[k]) {k++;}}
        if(pos_i[k] < P[j].Pos[k]) {do_cbe_calculation = 0;}
    }
    if(do_cbe_calculation == 1) // ok, go forward with the calculation //
    {
        double Face_Area_Norm=0, Face_Area_Vec[3]={0}, dp[3]={0}, wk_i=0, wk_j=0, dummy=0;
        double rho_i=0, rho_j=0, psi_i=0, psi_j=0, vface[3]={0};
        double V_i, V_j, L_j = Get_Particle_Size(j);
        // recalculate positions (to make sure have in correct format needed below //
        for(k=0;k<3;k++) {dp[k] = pos_i[k] - P[j].Pos[k];}
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
        } // for(k=0;k<3;k++)
        
        // OK, now we've done the easy bit -- we're ready for the actual flux computations //
        
        // first loop over pairs to determine closest a-to-b, closest b-to-a //
        int matching_basis_j_for_basis_in_i[CBE_INTEGRATOR_NBASIS], matching_basis_i_for_basis_in_j[CBE_INTEGRATOR_NBASIS];
        for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        {
            int m_j;
            for(m_j=m;m_j<CBE_INTEGRATOR_NBASIS;m_j++)
            {
                
            }
        } // for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        
        // then loop over each basis for each particle and compute fluxes //
        double wt_prefac_i = -rho_i / pmass; // normalized the fluxes, no need to re-compute below
        double wt_prefac_j = -rho_j / mass;
        double vface_dot_A = vface[0]*Face_Area_Vec[0] + vface[1]*Face_Area_Vec[1] + vface[2]*Face_Area_Vec[2]; // v_face . A_face
        double flux[10]={0};
        for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        {
            int j_m = matching_basis_j_for_basis_in_i[m], i_m = matching_basis_i_for_basis_in_j[m]; // id matched bases for fluxes below //
            // fluxes from "i" side
            do_cbe_flux_computation(local.CBE_basis_moments[m] , vface_dot_A, Face_Area_Vec, flux);
            for(k=0;k<10;k++)
            {
                flux[k] *= wt_prefac_i; // normalize appropriately
                out.CBE_basis_moments_dt[m][k] += flux[k]; // flux out of "i"
                if(TimeBinActive[P[j].TimeBin]) {P[j].CBE_basis_moments_dt[j_m][k] -= flux[k];} // flux into "j" (if j is active)
            }
            // fluxes from "j" side
            do_cbe_flux_computation(P[j].CBE_basis_moments[m] , vface_dot_A, Face_Area_Vec, flux);
            for(k=0;k<10;k++)
            {
                flux[k] *= wt_prefac_j; // normalize appropriately
                out.CBE_basis_moments_dt[i_m][k] += flux[k]; // flux out of "i"
                if(TimeBinActive[P[j].TimeBin]) {P[j].CBE_basis_moments_dt[j][k] -= flux[k];} // flux into "j" (if j is active)
            } // normalize appropriately
        } // for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
    } // if(do_cbe_calculation == 1)
} // master bracket (for variable protection
#endif

