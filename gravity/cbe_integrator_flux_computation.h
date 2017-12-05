#ifdef CBE_INTEGRATOR

/* ---------------------------------------------------------------------------------
 this is the kernel where the core of the CBE flux computation is performed,
 to calculate the relevant fluxes between faces of partciels, within the gravity routine
 --------------------------------------------------------------------------------- */
{
/* since we compute both 'sides' fluxes together (and it is an expensive operation), we
    only need to do it once -- we need to figure out which particle is active on the smaller timestep, and
    compute for it only */
    int do_cbe_calculation = 1; // default to doing the calculation
    int m, k, j = j0_sec_for_ags; // secondary has index 'j0_sec_for_ags' for gravity loop (cant use 'no' because this can get opened to the next neighbor before these operations
    double pos_i[3]={pos_x,pos_y,pos_z};
    if(targetdt_step > P[j].dt_step) do_cbe_calculation = 0; /* compute from particle with smaller timestep */
    if(targetdt_step == P[j].dt_step) // same timestep, randomly break degeneracy with positions
    {
        k=0; if(pos_i[k] == P[j].Pos[k]) {k++; if(pos_i[k] == P[j].Pos[k]) {k++;}}
        if(pos_i[k] < P[j].Pos[k]) {do_cbe_calculation = 0;}
    }
    if((All.Time <= All.TimeBegin) || (j < 0)) {do_cbe_calculation=0;}
    if(do_cbe_calculation == 1) // ok, go forward with the calculation //
    {
        double Face_Area_Vec[3]={0}, dp[3]={0}, wk_i=0, wk_j=0, dummy=0;
        double rho_i=0, rho_j=0, psi_i=0, psi_j=0, vface[3]={0};
        double V_i = local_V_i, V_j = 0, L_j = Get_Particle_Size_AGS(j);
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
        // calculate position between particles (psi) term for interpolation
        psi_i = 1. / (1. + h_p_inv/h_inv); // fraction of weight to i
        psi_j = 1-psi_i; // fraction of weight to j
        // calculate densities. rho never appears without a weight factor prefacing it, so just multiply here to get them; also convert to physical units //
        rho_i = psi_i * pmass/V_i * All.cf_a3inv; rho_j = psi_j * mass/V_j * All.cf_a3inv;
        // calculate kernel weight functions (for weights below) //
        if(u<1) {kernel_main(u, h3_inv, h3_inv*h_inv, &wk_i, &dummy, -1);} // get wk [these variables all pre-defined if here in loop] //
        if(u_p<1) {kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wk_j, &dummy, -1);} // likewise these should all be pre-calculated above, already [riding on adaptive gravsoft routines here] //
        // calculate effective faces and face velocity between elements //
        for(k=0;k<3;k++)
        {
            Face_Area_Vec[k] = wk_i*V_i * (local_NV_T[k][0]*dp[0] + local_NV_T[k][1]*dp[1] + local_NV_T[k][2]*dp[2])
                             + wk_j*V_j * (P[j].NV_T[k][0]*dp[0] + P[j].NV_T[k][1]*dp[1] + P[j].NV_T[k][2]*dp[2]);
            Face_Area_Vec[k] *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
            vface[k] = (rho_i*targetVel[k] + rho_j*P[j].Vel[k]) / (All.cf_atime * (rho_i+rho_j)); // normalized face velocity needed for zero mass flux (in physical units)
        } // for(k=0;k<3;k++)
        
        // OK, now we've done the easy bit -- we're ready for the actual computations //
        
        // first loop over pairs to determine closest a-to-b, closest b-to-a //
        int matching_basis_j_for_basis_in_i[CBE_INTEGRATOR_NBASIS], matching_basis_i_for_basis_in_j[CBE_INTEGRATOR_NBASIS], m_j;
        double wt_i[CBE_INTEGRATOR_NBASIS], wt_j[CBE_INTEGRATOR_NBASIS], imag_i[CBE_INTEGRATOR_NBASIS], imag_j[CBE_INTEGRATOR_NBASIS], cos_ij;
        for(m=0;m<CBE_INTEGRATOR_NBASIS;m++) // first pass to normalize and initialize quantities for both sides
        {
            wt_i[m] = -1000.; // large negative value (<-1)
            wt_j[m] = -1000.; // large negative value (<-1)
            double norm_tmp;
            norm_tmp=0; for(k=0;k<3;k++) {norm_tmp += local_CBE_basis_moments[m][k+1]*local_CBE_basis_moments[m][k+1];} // squared weight of momentum (vectors we'll use for matching below //
            if(norm_tmp > 0) {norm_tmp = 1./sqrt(norm_tmp);} // compute inverse-weight, for use below
            imag_i[m] = norm_tmp; // assign
            norm_tmp=0; for(k=0;k<3;k++) {norm_tmp += P[j].CBE_basis_moments[m][k+1]*P[j].CBE_basis_moments[m][k+1];} // squared weight of momentum (vectors we'll use for matching below //
            if(norm_tmp > 0) {norm_tmp = 1./sqrt(norm_tmp);} // compute inverse-weight, for use below
            imag_j[m] = norm_tmp; // assign
        }
        for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        {
            for(m_j=0;m_j<CBE_INTEGRATOR_NBASIS;m_j++)
            {
                cos_ij=0; for(k=1;k<4;k++) {cos_ij += local_CBE_basis_moments[m][k]*P[j].CBE_basis_moments[m_j][k];} // product of vector momenta
                cos_ij *= imag_i[m] * imag_j[m_j]; // ok this is now the dot product p_hat_i_alpha . p_hat_j_beta
                if(cos_ij > wt_i[m]) // better match found for i
                {
                    wt_i[m] = cos_ij; // note the best match so far
                    matching_basis_j_for_basis_in_i[m] = m_j; // assign the basis function to the list
                }
                if(cos_ij > wt_j[m_j]) // better match found for j
                {
                    wt_j[m_j] = cos_ij; // note the best match so far
                    matching_basis_i_for_basis_in_j[m_j] = m; // assign the basis function to the list
                }
            }
        } // for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        
        // now loop over each basis for each particle and actually compute fluxes //
        double wt_prefac_i = -rho_i / pmass; // normalized the fluxes, no need to re-compute below
        double wt_prefac_j = -rho_j / mass;
        double vface_dot_A = vface[0]*Face_Area_Vec[0] + vface[1]*Face_Area_Vec[1] + vface[2]*Face_Area_Vec[2]; // v_face . A_face
        double flux[10]={0};
        for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
        {
            int j_m = matching_basis_j_for_basis_in_i[m], i_m = matching_basis_i_for_basis_in_j[m]; // id matched bases for fluxes below //
            // fluxes from "i" side
            do_cbe_flux_computation(local_CBE_basis_moments[m] , vface_dot_A, Face_Area_Vec, flux); // moments are physical, these are as well //
            for(k=0;k<10;k++)
            {
                flux[k] *= wt_prefac_i; // normalize appropriately
                out_CBE_basis_moments_dt[m][k] += flux[k]; // flux out of "i"
                if(TimeBinActive[P[j].TimeBin]) {P[j].CBE_basis_moments_dt[j_m][k] -= flux[k];} // flux into "j" (if j is active)
            }
            // fluxes from "j" side
            do_cbe_flux_computation(P[j].CBE_basis_moments[m] , vface_dot_A, Face_Area_Vec, flux); // moments are physical, these are as well //
            for(k=0;k<10;k++)
            {
                flux[k] *= wt_prefac_j; // normalize appropriately
                out_CBE_basis_moments_dt[i_m][k] += flux[k]; // flux out of "i"
                if(TimeBinActive[P[j].TimeBin]) {P[j].CBE_basis_moments_dt[m][k] -= flux[k];} // flux into "j" (if j is active)
            } // normalize appropriately
        } // for(m=0;m<CBE_INTEGRATOR_NBASIS;m++)
    } // if(do_cbe_calculation == 1)
} // master bracket (for variable protection
#endif

