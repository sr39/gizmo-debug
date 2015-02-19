/* --------------------------------------------------------------------------------- */
/* this is the sub-routine where we actually extrapolate quantities to the faces 
    and set up, then solve the pair-wise Riemann problem for the method */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
//#define DO_HALFSTEP_FOR_MESHLESS_METHODS 1
    
    double s_star_ij,s_i,s_j,v_frame[3],n_unit[3];
    double distance_from_i[3],distance_from_j[3];
    face_vel_i=face_vel_j=Face_Area_Norm=0;
#ifdef COSMIC_RAYS
    Fluxes.CosmicRayPressure = 0;
#endif
    
    /* --------------------------------------------------------------------------------- */
    /* define volume elements and interface position */
    /* --------------------------------------------------------------------------------- */
    V_j = P[j].Mass / SphP[j].Density;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    s_star_ij = 0;
#else
    s_star_ij = 0.5 * kernel.r * (PPP[j].Hsml - local.Hsml) / (local.Hsml + PPP[j].Hsml);
#endif
    /* ------------------------------------------------------------------------------------------------------------------- */
    /* now we're ready to compute the volume integral of the fluxes (or equivalently an 'effective area'/face orientation) */
    /* ------------------------------------------------------------------------------------------------------------------- */
    if(SphP[j].ConditionNumber*SphP[j].ConditionNumber > cnumcrit2)
    {
        /* the effective gradient matrix is ill-conditioned: for stability, we revert to the "RSPH" EOM */
        /* we need to evaluate dwk, since it is not done by default */
        if(kernel.r < kernel.h_i)
        {
            u = kernel.r * hinv_i;
            kernel_main(u, hinv3_i, hinv4_i, &kernel.wk_i, &kernel.dwk_i, 1);
        }
        if(kernel.r < kernel.h_j)
        {
            u = kernel.r * hinv_j;
            kernel_main(u, hinv3_j, hinv4_j, &kernel.wk_j, &kernel.dwk_j, 1);
        }
        Face_Area_Norm = -(V_i*V_i*kernel.dwk_i + V_j*V_j*kernel.dwk_j) / kernel.r;
        Face_Area_Norm *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
        Face_Area_Vec[0] = Face_Area_Norm * kernel.dp[0];
        Face_Area_Vec[1] = Face_Area_Norm * kernel.dp[1];
        Face_Area_Vec[2] = Face_Area_Norm * kernel.dp[2];
        Face_Area_Norm = Face_Area_Norm * Face_Area_Norm * r2;
    } else {
        /* the effective gradient matrix is well-conditioned: we can safely use the consistent EOM */
        for(k=0;k<3;k++)
        {
            Face_Area_Vec[k] = kernel.wk_i * V_i * (local.NV_T[k][0]*kernel.dp[0] + local.NV_T[k][1]*kernel.dp[1] + local.NV_T[k][2]*kernel.dp[2])
                       + kernel.wk_j * V_j * (SphP[j].NV_T[k][0]*kernel.dp[0] + SphP[j].NV_T[k][1]*kernel.dp[1] + SphP[j].NV_T[k][2]*kernel.dp[2]);
            Face_Area_Vec[k] *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
            Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k];
        }
    }
    if(Face_Area_Norm == 0)
    {
        memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
        Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
    } else {
        
        if((Face_Area_Norm<=0)||(isnan(Face_Area_Norm)))
        {
            printf("PANIC! Face_Area_Norm=%g Mij=%g/%g wk_ij=%g/%g Vij=%g/%g dx/dy/dz=%g/%g/%g NVT=%g/%g/%g NVT_j=%g/%g/%g \n",Face_Area_Norm,local.Mass,P[j].Mass,kernel.wk_i,
                   kernel.wk_j,V_i,V_j,kernel.dp[0],kernel.dp[1],kernel.dp[2],local.NV_T[0][0],local.NV_T[0][1],local.NV_T[0][2],SphP[j].NV_T[0][0],SphP[j].NV_T[0][1],
                   SphP[j].NV_T[0][2]);fflush(stdout);
        }
        Face_Area_Norm = sqrt(Face_Area_Norm);
        for(k=0;k<3;k++) {n_unit[k] = Face_Area_Vec[k] / Face_Area_Norm;}
        

        /* --------------------------------------------------------------------------------- */
        /* extrapolate the conserved quantities to the interaction face between the particles */
        /* first we define some useful variables for the extrapolation */
        /* --------------------------------------------------------------------------------- */
        s_i =  0.5 * kernel.r;
        s_j = -0.5 * kernel.r;
        //(simple up-winding formulation: use if desired instead of time-centered fluxes)//
        //delta_halfstep_i = kernel.sound_i*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_i>s_i) {delta_halfstep_i=s_i;}
        //delta_halfstep_j = kernel.sound_j*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_j>-s_j) {delta_halfstep_j=-s_j;}
#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing
            certain types of noise and oscillations (but not always!) */
        s_i += 0.5 * dt_hydrostep * (local.Vel[0]*kernel.dp[0] + local.Vel[1]*kernel.dp[1] + local.Vel[2]*kernel.dp[2]) * rinv;
        s_j += 0.5 * dt_hydrostep * (VelPred_j[0]*kernel.dp[0] + VelPred_j[1]*kernel.dp[1] + VelPred_j[2]*kernel.dp[2]) * rinv;
#endif
        s_i = s_star_ij - s_i; //+ delta_halfstep_i; /* projection element for gradients */
        s_j = s_star_ij - s_j; //- delta_halfstep_j;
        distance_from_i[0]=kernel.dp[0]*rinv; distance_from_i[1]=kernel.dp[1]*rinv; distance_from_i[2]=kernel.dp[2]*rinv;
        for(k=0;k<3;k++) {distance_from_j[k] = distance_from_i[k] * s_j; distance_from_i[k] *= s_i;}
        //for(k=0;k<3;k++) {v_frame[k] = 0.5 * (VelPred_j[k] + local.Vel[k]);}
        for(k=0;k<3;k++) {v_frame[k] = rinv * (-s_i*VelPred_j[k] + s_j*local.Vel[k]);} // allows for face to be off-center (to second-order)
        // (note that in the above, the s_i/s_j terms are crossed with the opposing velocity terms: this is because the face is closer to the
        //   particle with the smaller smoothing length; so it's values are slightly up-weighted //
        
        
        /* now we do the reconstruction (second-order reconstruction at the face) */
        reconstruct_face_states(local.Density, local.Gradients.Density, SphP[j].Density, SphP[j].Gradients.Density,
                                distance_from_i, distance_from_j, &Riemann_vec.L.rho, &Riemann_vec.R.rho, 1);
        reconstruct_face_states(local.Pressure, local.Gradients.Pressure, SphP[j].Pressure, SphP[j].Gradients.Pressure,
                                distance_from_i, distance_from_j, &Riemann_vec.L.p, &Riemann_vec.R.p, 1);
#ifdef NON_IDEAL_EOS
        reconstruct_face_states(local.InternalEnergyPred, local.Gradients.InternalEnergy, SphP[j].InternalEnergyPred, SphP[j].Gradients.InternalEnergy,
                                distance_from_i, distance_from_j, &Riemann_vec.L.u, &Riemann_vec.R.u, 1);
        reconstruct_face_states(kernel.sound_i, local.Gradients.SoundSpeed, kernel.sound_j, SphP[j].Gradients.SoundSpeed,
                                distance_from_i, distance_from_j, &Riemann_vec.L.cs, &Riemann_vec.R.cs, 1);
#endif
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.Vel[k]-v_frame[k], local.Gradients.Velocity[k], VelPred_j[k]-v_frame[k], SphP[j].Gradients.Velocity[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.v[k], &Riemann_vec.R.v[k], 1);
        }
#ifdef MAGNETIC
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.BPred[k], local.Gradients.B[k], BPred_j[k], SphP[j].Gradients.B[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.B[k], &Riemann_vec.R.B[k], 1);
        }
#ifdef DIVBCLEANING_DEDNER
        reconstruct_face_states(local.PhiPred, local.Gradients.Phi, PhiPred_j, SphP[j].Gradients.Phi,
                                distance_from_i, distance_from_j, &Riemann_vec.L.phi, &Riemann_vec.R.phi, 1);
#endif
#endif

#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing 
            certain types of noise and oscillations */
        double dt_half = 0.5*dt_hydrostep;
        for(k=0;k<3;k++)
        {
            Riemann_vec.R.rho -= dt_half * local.Density * local.Gradients.Velocity[k][k];
            Riemann_vec.L.rho -= dt_half * SphP[j].Density * SphP[j].Gradients.Velocity[k][k];
            Riemann_vec.R.p -= dt_half * GAMMA * local.Pressure * local.Gradients.Velocity[k][k];
            Riemann_vec.L.p -= dt_half * GAMMA * SphP[j].Pressure * SphP[j].Gradients.Velocity[k][k];
            double dv_l_half = -dt_half * local.Gradients.Pressure[k] / local.Density;
            double dv_r_half = -dt_half * SphP[j].Gradients.Pressure[k] / SphP[j].Density;
            /* // this part doesn't need to be done, because our derivatives are co-moving; but if they are not, use it //
             Riemann_vec.R.rho -= dt_half * (local.Vel[k]-v_frame[k]) * local.Gradients.Density[k];
             Riemann_vec.L.rho -= dt_half * (VelPred_j[k]-v_frame[k]) * SphP[j].Gradients.Density[k];
             Riemann_vec.R.p -= dt_half * (local.Vel[k]-v_frame[k]) * local.Gradients.Pressure[k];
             Riemann_vec.L.p -= dt_half * (VelPred_j[k]-v_frame[k]) * SphP[j].Gradients.Pressure[k];
             for(int kx=0;kx<3;kx++)
             {
                dv_r_half -= dt_half * (local.Vel[kx]-v_frame[kx]) * local.Gradients.Velocity[k][kx];
                dv_l_half -= dt_half * (VelPred_j[kx]-v_frame[kx]) * SphP[j].Gradients.Velocity[k][kx];
             }
            */
            Riemann_vec.R.v[k] += 0.5 * (dv_l_half - dv_r_half);
            Riemann_vec.L.v[k] += 0.5 * (dv_r_half - dv_l_half);
            v_frame[k] += 0.5*(dv_l_half + dv_l_half);
        }
#endif
        
        /* --------------------------------------------------------------------------------- */
        /* Alright! Now we're actually ready to solve the Riemann problem at the particle interface */
        /* --------------------------------------------------------------------------------- */
        Riemann_solver(Riemann_vec, &Riemann_out, n_unit);
        /* before going on, check to make sure we have a valid Riemann solution */
        if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
        {
            /* go to a linear reconstruction of P, rho, and v, and re-try */
            Riemann_vec.R.p = local.Pressure; Riemann_vec.L.p = SphP[j].Pressure;
            Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
            for(k=0;k<3;k++) {Riemann_vec.R.v[k]=local.Vel[k]-v_frame[k]; Riemann_vec.L.v[k]=VelPred_j[k]-v_frame[k];}
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
            Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef NON_IDEAL_EOS
            Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
            Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
            Riemann_solver(Riemann_vec, &Riemann_out, n_unit);
            if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
            {
                /* ignore any velocity difference between the particles: this should gaurantee we have a positive pressure! */
                Riemann_vec.R.p = local.Pressure; Riemann_vec.L.p = SphP[j].Pressure;
                Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
                for(k=0;k<3;k++) {Riemann_vec.R.v[k]=0; Riemann_vec.L.v[k]=0;}
#ifdef MAGNETIC
                for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
                Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef NON_IDEAL_EOS
                Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
                Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
                Riemann_solver(Riemann_vec, &Riemann_out, n_unit);
                if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
                {
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
                    printf("Riemann Solver Failed to Find Positive Pressure!: PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g H_ij=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g BL=%g/%g/%g BR=%g/%g/%g phiL/R=%g/%g \n",
                           Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,local.Hsml,PPP[j].Hsml,
                           local.Vel[0]-v_frame[0],local.Vel[1]-v_frame[1],local.Vel[2]-v_frame[2],
                           VelPred_j[0]-v_frame[0],VelPred_j[1]-v_frame[1],VelPred_j[2]-v_frame[2],
                           n_unit[0],n_unit[1],n_unit[2],
                           Riemann_vec.L.B[0],Riemann_vec.L.B[1],Riemann_vec.L.B[2],
                           Riemann_vec.R.B[0],Riemann_vec.R.B[1],Riemann_vec.R.B[2],
                           Riemann_vec.L.phi,Riemann_vec.R.phi);
#else
                    printf("Riemann Solver Failed to Find Positive Pressure!: PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g \n",
                           Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,
                           Riemann_vec.L.v[0],Riemann_vec.L.v[1],Riemann_vec.L.v[2],
                           Riemann_vec.R.v[0],Riemann_vec.R.v[1],Riemann_vec.R.v[2],n_unit[0],n_unit[1],n_unit[2]);
#endif
                    exit(1234);
                }
            }
        } // closes loop of alternative reconstructions if invalid pressures are found //
        
        
        /* --------------------------------------------------------------------------------- */
        /* Calculate the fluxes (EQUATION OF MOTION) -- all in physical units -- */
        /* --------------------------------------------------------------------------------- */
        if((Riemann_out.P_M>0)&&(!isnan(Riemann_out.P_M)))
        {
            if(All.ComovingIntegrationOn) {for(k=0;k<3;k++) v_frame[k] /= All.cf_atime;}
#if !defined(HYDRO_MESHLESS_FINITE_VOLUME) && !defined(MAGNETIC)
            double facenorm_pm = Face_Area_Norm * Riemann_out.P_M;
            Fluxes.p = 0;
            for(k=0;k<3;k++)
            {
                Fluxes.v[k] = facenorm_pm * n_unit[k]; /* total momentum flux */
                Fluxes.p += facenorm_pm * (Riemann_out.S_M*n_unit[k] + v_frame[k]) * n_unit[k]; /* total energy flux = v_frame.dot.mom_flux */
            }
#else
            /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
             which we do following Pakmor et al. 2011 */
            for(k=0;k<3;k++)
            {
                /* Riemann_out->Fluxes.rho is un-modified */
                Riemann_out.Fluxes.p += (0.5*v_frame[k]*v_frame[k])*Riemann_out.Fluxes.rho + v_frame[k]*Riemann_out.Fluxes.v[k];
                Riemann_out.Fluxes.v[k] += v_frame[k] * Riemann_out.Fluxes.rho; /* just boost by frame vel (as we would in non-moving frame) */
            }
#ifdef MAGNETIC
            for(k=0;k<3;k++)
            {
                Riemann_out.Fluxes.B[k] += -v_frame[k] * Riemann_out.B_normal_corrected; /* v dotted into B along the normal to the face (careful of sign here) */
                face_vel_i += local.Vel[k] * n_unit[k] / All.cf_atime;
                face_vel_j += VelPred_j[k] * n_unit[k] / All.cf_atime;
            }
            face_area_dot_vel = -(-s_i*face_vel_j + s_j*face_vel_i) * rinv;
#endif
            
            /* ok now we can actually apply this to the EOM */
            Fluxes.rho = Face_Area_Norm * Riemann_out.Fluxes.rho;
            Fluxes.p = Face_Area_Norm * Riemann_out.Fluxes.p; // this is really Dt of --total-- energy, need to subtract KE component for e */
            for(k=0;k<3;k++)
            {
                Fluxes.v[k] = Face_Area_Norm * Riemann_out.Fluxes.v[k]; // momentum flux (need to divide by mass) //
            }
#if defined(COSMIC_RAYS) && defined(HYDRO_MESHLESS_FINITE_VOLUME)
            /* here we simply assume that if there is mass flux, the cosmic ray fluid is advected -with the mass flux-, taking an
             implicit constant (zeroth-order) reconstruction of the CR energy density at the face (we could reconstruct the CR
             properties at the face, and calculate a more accurate advection term; however at that stage we should actually be
             including them self-consistently in the Riemann problem */
            if(Fluxes.rho < 0)
            {
                Fluxes.CosmicRayPressure = Fluxes.rho * (local.CosmicRayPressure*V_i/(GAMMA_COSMICRAY_MINUS1*local.Mass));
            } else {
                Fluxes.CosmicRayPressure = Fluxes.rho * (CosmicRayPressure_j*V_j/(GAMMA_COSMICRAY_MINUS1*P[j].Mass));
            }
#endif // cosmic_rays
            
#ifdef MAGNETIC
            for(k=0;k<3;k++)
            {
                Fluxes.B[k] = Face_Area_Norm * Riemann_out.Fluxes.B[k]; // flux of magnetic flux (B*V) //
            }
            Fluxes.B_normal_corrected = -Riemann_out.B_normal_corrected * Face_Area_Norm;
#ifdef DIVBCLEANING_DEDNER
            Fluxes.phi = Riemann_out.Fluxes.phi * Face_Area_Norm; // much more accurate than mass-based flux //
#endif // DIVBCLEANING_DEDNER
#endif // MAGNETIC
#endif
        } else {
            /* nothing but bad riemann solutions found! */
            memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
            Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
        }
    } // Face_Area_Norm != 0
}
