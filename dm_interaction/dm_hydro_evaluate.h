#ifdef DM_BARYON_INTERACTION
int dm_hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, k, n, m, l, startnode, numngb, kernel_mode, listindex;
    double hinv_i,hinv3_i,hinv4_i,hinv_j,hinv3_j,hinv4_j,V_i,V_j,r2,rinv,rinv_soft,u;
    double v_hll,k_hll,b_hll; v_hll=k_hll=0,b_hll=1;
    struct dm_kernel_hydra kernel;
    struct dm_hydrodata_in local;
    struct dm_hydrodata_out out;
    struct dm_Conserved_var_Riemann Fluxes;
    double dt_dmhydrostep;
    double LENGTH_physical = All.UnitLength_in_cm * All.cf_atime / All.HubbleParam;
    double MASS_physical = All.UnitMass_in_g / All.HubbleParam; 
   
    double n_i, n_j, N_i, N_j;
    double sigma_u, sigma_dm_proton, rate; 

    listindex = 0;
    memset(&out, 0, sizeof(struct dm_hydrodata_out));
    memset(&kernel, 0, sizeof(struct dm_kernel_hydra));
    memset(&Fluxes, 0, sizeof(struct dm_Conserved_var_Riemann));

    if(mode == 0)
    {
        dm_particle2in_hydra(&local, target); // this setup allows for all the fields we need to define (don't hard-code here)
    }
    else
    {
        local = DM_HydroDataGet[target]; // this setup allows for all the fields we need to define (don't hard-code here)
    }

    /* --------------------------------------------------------------------------------- */
    /* pre-define Particle-i based variables (so we save time in the loop below) */
    /* --------------------------------------------------------------------------------- */
    kernel.sound_i = local.SoundSpeed;
    kernel.spec_egy_u_i = local.InternalEnergyPred;

if(local.dm_Hsml > local.Hsml)
    {
      kernel.h_i = local.dm_Hsml;
    }else{
      kernel.h_i = local.Hsml;
    }
    kernel_hinv(kernel.h_i, &hinv_i, &hinv3_i, &hinv4_i);
    hinv_j=hinv3_j=hinv4_j=0;
    V_i = local.Mass / local.Density;
    kernel_mode = 0;

if(mode == 0)
    {   
        startnode = All.MaxPart;        /* root node */
    }
    else
    {
#ifndef DONOTUSENODELIST
        startnode = DM_HydroDataGet[target].NodeList[0];   
        startnode = Nodes[startnode].u.d.nextnode;      /* open it */
#else   
        startnode = All.MaxPart;        /* root node */
#endif
    }

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            /* --------------------------------------------------------------------------------- */
            /* get the neighbor list */
            /* --------------------------------------------------------------------------------- */
           // numngb = ngb_treefind_variable_threads_targeted(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
           //                            exportnodecount, exportindex, ngblist, TARGET_BITMASK);
            numngb = ngb_treefind_variable_threads_targeted(local.Pos, kernel.h_i, target, &startnode,
                                  mode, exportflag, exportnodecount, exportindex, ngblist, 3);

            if(numngb < 0) return -1;
            if(local.dm_count == 0) continue;
            if(local.baryon_count == 0) continue;

            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n];
                   if(P[j].Mass <= 0) continue;
       
                    kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                    kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                    kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC  /* find the closest image in the given box size  */
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                    r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];

          
              if(r2 >= kernel.h_i * kernel.h_i) continue;
              if(r2 <= 0) continue;

               kernel.r = sqrt(r2);
               
                  if(kernel.r < kernel.h_i)
                   {
                  kernel.vrel[0] = ((local.vx_dm / local.dm_count) - (local.vx_baryon / local.baryon_count)) / All.cf_atime;
                  kernel.vrel[1] = ((local.vy_dm / local.dm_count) - (local.vy_baryon / local.baryon_count)) / All.cf_atime; 
                  kernel.vrel[2] = ((local.vz_dm / local.dm_count) - (local.vz_baryon / local.baryon_count)) / All.cf_atime;

                  kernel.v_rel = sqrt(kernel.vrel[0] * kernel.vrel[0] + kernel.vrel[1] * kernel.vrel[1] + kernel.vrel[2] * kernel.vrel[2]);                      
  
                  u = kernel.r * hinv_i;
                  kernel_main(u, hinv3_i, hinv4_i, &kernel.wk_i, &kernel.dwk_i, kernel_mode);
         

           sigma_dm_proton = (CROSS_SECTION_XP / LENGTH_physical / LENGTH_physical) * pow(kernel.v_rel / REF_VELOCITY_XP, P_XP-4);//(cm^2)
           N_j = (P[j].Mass * All.UnitMass_in_g / All.HubbleParam) / DARKMATTERMASS;
           n_j = (local.dm_density * All.cf_a3inv * All.UnitDensity_in_cgs / All.HubbleParam / All.HubbleParam) / DARKMATTERMASS;//ダークマター粒子としての数密度
           N_i = (P[target].Mass * All.UnitMass_in_g / All.HubbleParam) / PROTONMASS;
           n_i = (local.Density * All.cf_a3inv * All.UnitDensity_in_cgs / All.HubbleParam / All.HubbleParam) / PROTONMASS;//プロトン粒子としての数密度
           sigma_u = sqrt((GAMMA - 1) * local.InternalEnergyPred);

//        double interaction_probability = dm_dv * N_i * kernel.wk_i * sigma_dm_proton * dt_dmhydrostep / All.HubbleParam;

//        if(interaction_probability <= 0)continue;

//        init_genrand((unsigned)time(NULL));

//        if(genrand_real1() < interaction_probability)
//        {  
         rate = energy_transfer_rate_in_m5_over_s3(kernel.v_rel, sigma_u); 
         
        if(P[j].Type == 0)
         {
           if(kernel.r < local.Hsml)
            {
              SphP[j].dm_DtInternalEnergy += (PROTONMASS * DARKMATTERMASS / (PROTONMASS + DARKMATTERMASS) / MASS_physical) * n_i * n_j * pow(LENGTH_physical, 6) * rate * V_i;//dE/dt(baryon)

              for(k = 0; k < 3; k++)
              {
               SphP[j].baryon_dtVel[k] += ( DARKMATTERMASS / (PROTONMASS + DARKMATTERMASS)) * n_j * pow(LENGTH_physical, 3) * rate * kernel.vrel[k]  / (kernel.v_rel * kernel.v_rel);//dv/dt(baryon)
              }
           
           }
         }

        if(P[j].Type == 1)
         {   
           if(kernel.r < local.dm_Hsml)
            {
              for(k = 0; k < 3; k++)
              {
                P[j].dm_dtVel[k] -= ( PROTONMASS / (PROTONMASS + DARKMATTERMASS)) * n_i * pow(LENGTH_physical, 3) * rate * kernel.vrel[k] / (kernel.v_rel * kernel.v_rel); // dv/dt(DM)
              }
           }
         }
        }//if(kernel.dm_r < kernel.dm_h_i)
       } // for(n = 0; n < numngb; n++) //
      } // while(startnode >= 0) /
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DM_HydroDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1) //
#endif

    } // while(startnode >= 0) //
  
   
 
    /* Now collect the result at the right place */
    if(mode == 0)
        dm_out2particle_hydra(&out, target, 0);
    else
        DM_HydroDataResult[target] = out;
    
    return 0;
}
#endif                                                                                                    
