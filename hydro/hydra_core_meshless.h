/* --------------------------------------------------------------------------------- */
/* this is the sub-routine where we actually extrapolate quantities to the faces 
    and set up, then solve the pair-wise Riemann problem for the method */
/* --------------------------------------------------------------------------------- */
{
    double s_star_ij,s_i,s_j,v_frame[3],n_unit[3];
    double distance_from_i[3],distance_from_j[3];
    Fluxes.rho = Fluxes.p = Fluxes.v[0] = Fluxes.v[1] = Fluxes.v[2] = 0;
    
    /* --------------------------------------------------------------------------------- */
    /* define volume elements and interface position */
    /* --------------------------------------------------------------------------------- */
    V_j = P[j].Mass / SphP[j].Density;
    s_star_ij = 0;

    /* ------------------------------------------------------------------------------------------------------------------- */
    /* now we're ready to compute the volume integral of the fluxes (or equivalently an 'effective area'/face orientation) */
    /* ------------------------------------------------------------------------------------------------------------------- */
    double dwk_tmp[3],dwk_norm=0;
    if(SphP[j].ConditionNumber*SphP[j].ConditionNumber > cnumcrit2)
    {
        /* the effective gradient matrix is ill-conditioned: for stability, we revert to the "RSPH" EOM */
#if !(defined(HYDRO_SPH) || defined(CONDUCTION_EXPLICIT) || defined(TURB_DIFFUSION))
        /* we may need to evaluate dwk, if it is not done by default */
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
#endif
        dwk_norm = -(V_i*V_i*kernel.dwk_i + V_j*V_j*kernel.dwk_j) / kernel.r;
        dwk_tmp[0] = dwk_norm * kernel.dx;
        dwk_tmp[1] = dwk_norm * kernel.dy;
        dwk_tmp[2] = dwk_norm * kernel.dz;
        dwk_norm = dwk_norm * dwk_norm * r2;
    } else {
        /* the effective gradient matrix is well-conditioned: we can safely use the consistent EOM */
        for(k=0;k<3;k++)
        {
            dwk_tmp[k] = kernel.wk_i * V_i * (local.NV_T[k][0]*kernel.dx + local.NV_T[k][1]*kernel.dy + local.NV_T[k][2]*kernel.dz)
                       + kernel.wk_j * V_j * (SphP[j].NV_T[k][0]*kernel.dx + SphP[j].NV_T[k][1]*kernel.dy + SphP[j].NV_T[k][2]*kernel.dz);
            dwk_norm += dwk_tmp[k]*dwk_tmp[k];
        }
    }
    if(dwk_norm == 0)
    {
    } else {
        
        if((dwk_norm<=0)||(isnan(dwk_norm)))
        {
            printf("PANIC! dwk_norm=%g Mij=%g/%g wk_ij=%g/%g Vij=%g/%g dx/dy/dz=%g/%g/%g NVT=%g/%g/%g NVT_j=%g/%g/%g \n",dwk_norm,local.Mass,P[j].Mass,kernel.wk_i,
                   kernel.wk_j,V_i,V_j,kernel.dx,kernel.dy,kernel.dz,local.NV_T[0][0],local.NV_T[0][1],local.NV_T[0][2],SphP[j].NV_T[0][0],SphP[j].NV_T[0][1],
                   SphP[j].NV_T[0][2]);fflush(stdout);
        }
        dwk_norm = sqrt(dwk_norm);
        for(k=0;k<3;k++) {n_unit[k] = dwk_tmp[k] / dwk_norm;}
        dwk_norm *= All.cf_atime*All.cf_atime; /* dwk_norm has units of area, need to convert to physical */
            

        /* --------------------------------------------------------------------------------- */
        /* extrapolate the conserved quantities to the interaction face between the particles */
        /* first we define some useful variables for the extrapolation */
        /* --------------------------------------------------------------------------------- */
        s_i =  0.5 * kernel.r;
        s_j = -0.5 * kernel.r;
        delta_halfstep_i = kernel.sound_i*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_i>s_i) {delta_halfstep_i=s_i;}
        delta_halfstep_j = kernel.sound_j*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_j>-s_j) {delta_halfstep_j=-s_j;}
        s_i = s_star_ij - s_i + delta_halfstep_i; /* projection element for gradients */
        s_j = s_star_ij - s_j - delta_halfstep_j;
        distance_from_i[0]=kernel.dx*rinv; distance_from_i[1]=kernel.dy*rinv; distance_from_i[2]=kernel.dz*rinv;
        for(k=0;k<3;k++) {distance_from_j[k] = distance_from_i[k] * s_j; distance_from_i[k] *= s_i;}
        for(k=0;k<3;k++) {v_frame[k] = 0.5 * (local.Vel[k] + SphP[j].VelPred[k]);}
        
        
        /* now we do the reconstruction (second-order reconstruction at the face) */
        reconstruct_face_states(local.Density, local.Gradients.Density, SphP[j].Density, SphP[j].Gradients.Density,
                                distance_from_i, distance_from_j, &Riemann_vec.L.rho, &Riemann_vec.R.rho, 1);
        reconstruct_face_states(local.Pressure, local.Gradients.Pressure, SphP[j].Pressure, SphP[j].Gradients.Pressure,
                                distance_from_i, distance_from_j, &Riemann_vec.L.p, &Riemann_vec.R.p, 1);
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.Vel[k]-v_frame[k], local.Gradients.Velocity[k], SphP[j].VelPred[k]-v_frame[k], SphP[j].Gradients.Velocity[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.v[k], &Riemann_vec.R.v[k], 1);
        }
#ifdef MAGNETIC
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.Bpred[k], local.Gradients.B[k], SphP[j].BPred[k], SphP[j].Gradients.B[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.B[k], &Riemann_vec.R.B[k], 1);
        }
#ifdef DIVBCLEANING_DEDNER
        reconstruct_face_states(local.PhiPred, local.Gradients.Phi, SphP[j].PhiPred, SphP[j].Gradients.Phi,
                                distance_from_i, distance_from_j, &Riemann_vec.L.phi, &Riemann_vec.R.phi, 1);
#endif
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
            for(k=0;k<3;k++) {Riemann_vec.R.v[k]=local.Vel[k]-v_frame[k]; Riemann_vec.L.v[k]=SphP[j].VelPred[k]-v_frame[k];}
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.Bpred[k]; Riemann_vec.L.B[k]=SphP[j].Bpred[k];}
#ifdef DIVBCLEANING_DEDNER
            Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = SphP[j].PhiPred;
#endif
#endif
            Riemann_solver(Riemann_vec, &Riemann_out, n_unit);
            if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
            {
                /* ignore any velocity difference between the particles: this should gaurantee we have a positive pressure! */
                Riemann_vec.R.p = local.Pressure; Riemann_vec.L.p = SphP[j].Pressure;
                Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
                for(k=0;k<3;k++) {Riemann_vec.R.v[k]=0; Riemann_vec.L.v[k]=0;}
#ifdef MAGNETIC
                for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.Bpred[k]; Riemann_vec.L.B[k]=SphP[j].Bpred[k];}
#ifdef DIVBCLEANING_DEDNER
                Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = SphP[j].PhiPred;
#endif
#endif
                Riemann_solver(Riemann_vec, &Riemann_out, n_unit);
                if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
                {
                    printf("Riemann Solver Failed to Find Positive Pressure!: PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g \n",
                           Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,
                           Riemann_vec.L.v[0],Riemann_vec.L.v[1],Riemann_vec.L.v[2],
                           Riemann_vec.R.v[0],Riemann_vec.R.v[1],Riemann_vec.R.v[2],n_unit[0],n_unit[1],n_unit[2]);
                    exit(1234);
                }
            }
        }
        
        
        /* --------------------------------------------------------------------------------- */
        /* Calculate the fluxes (EQUATION OF MOTION) -- all in physical units -- */
        /* --------------------------------------------------------------------------------- */
        if((Riemann_out.P_M>0)&&(!isnan(Riemann_out.P_M)))
        {
            if(All.ComovingIntegrationOn) {for(k=0;k<3;k++) v_frame[k] /= All.cf_atime;}
#if !defined(HYDRO_MESHLESS_FINITE_VOLUME) && !defined(MAGNETIC)
            dwk_norm *= Riemann_out.P_M;
            for(k=0;k<3;k++)
            {
                Fluxes.v[k] = dwk_norm * n_unit[k]; /* total momentum flux */
                Fluxes.p += dwk_norm * (Riemann_out.S_M*n_unit[k] + v_frame[k]) * n_unit[k]; /* total energy flux = v_frame.dot.mom_flux */
            }
#else
            /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
             which we do following Pakmor et al. 2011 */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME 
            /* limiter to prevent crazy mass fluxes per-particle per-timestep */
            if(dt_hydrostep > 0)
            {
                double dmass_min = 0.01 * DMIN(DMIN(local.Mass,SphP[j].MassTrue),P[j].Mass) / dt_hydrostep;
                if(fabs(Riemann_out.Fluxes.rho) > dmass_min) {Riemann_out.Fluxes.rho *= dmass_min / fabs(Riemann_out.Fluxes.rho);}
            }
#endif
            for(k=0;k<3;k++)
            {
                /* Riemann_out->Fluxes.rho is un-modified */
                Riemann_out.Fluxes.p += (0.5*v_frame[k]*v_frame[k])*Riemann_out.Fluxes.rho + v_frame[k]*Riemann_out.Fluxes.v[k];
                Riemann_out.Fluxes.v[k] += v_frame[k] * Riemann_out.Fluxes.rho; /* just boost by frame vel (as we would in non-moving frame) */
            }
            /* ok now we can actually apply this to the EOM */
            Fluxes.rho = dwk_norm * Riemann_out.Fluxes.rho;
            Fluxes.p = dwk_norm * Riemann_out.Fluxes.p; // this is really Dt of --total-- energy, need to subtract KE component for e */
            for(k=0;k<3;k++)
            {
                Fluxes.v[k] = dwk_norm * Riemann_out.Fluxes.v[k]; // momentum flux (need to divide by mass) //
            }
#endif
        }
    } // dwk_norm != 0
}
