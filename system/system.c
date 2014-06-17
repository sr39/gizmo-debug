#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>


#include "../allvars.h"
#include "../proto.h"


void merge_and_split_particles(void)
{
    int target_for_merger,dummy=0,numngb_inbox,startnode,i,j,n;
    double threshold_val;
    int n_particles_merged,n_particles_split,MPI_n_particles_merged,MPI_n_particles_split;
    Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    Gas_split=0;
    n_particles_merged=0;
    n_particles_split=0;
    MPI_n_particles_merged=0;
    MPI_n_particles_split=0;
    /* loop over active particles */
    for(i=0; i<NumPart; i++)
    {
        /* check if we're a gas particle */
        if((P[i].Type==0)&&(TimeBinActive[P[i].TimeBin]))
        {
            /* we have a gas particle, ask if it needs to be merged */
            if(does_particle_need_to_be_merged(i))
            {
                /* if merging: do a neighbor loop ON THE SAME DOMAIN to determine the target */
                startnode=All.MaxPart;
                numngb_inbox = ngb_treefind_variable_threads(P[i].Pos,PPP[i].Hsml,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
                if(numngb_inbox>0)
                {
                    target_for_merger = -1;
                    threshold_val = MAX_REAL_NUMBER;
                    /* loop over neighbors */
                    for(n=0; n<numngb_inbox; n++)
                    {
                        j = Ngblist[n];
                        /* make sure we're not taking the same particle (and that its available to be merged into)! */
                        if((j>=0)&&(j!=i)&&(P[j].Type==0)&&(P[j].Mass > P[i].Mass)&&(P[i].Mass+P[j].Mass < All.MaxMassForParticleSplit))
                        {
                            if(P[j].Mass<threshold_val) {threshold_val=P[j].Mass; target_for_merger=j;} // mass-based //
                        }
                    } // for(n=0; n<numngb_inbox; n++)
                    if(target_for_merger >= 0)
                    {
                        /* we have a valid target! now we can actually do the merger operation! */
                        merge_particles_ij(i,target_for_merger);
                        n_particles_merged++;
                    }
                } // if(numngb_inbox>0)
            } // if(does_particle_need_to_be_merged(i))
            /* alright, the particle merger operations are complete! */
            
            /* now ask if the particle needs to be split */
            if(does_particle_need_to_be_split(i))
            {
                /* if splitting: do a neighbor loop ON THE SAME DOMAIN to determine the nearest particle (so dont overshoot it) */
                startnode=All.MaxPart;
                numngb_inbox = ngb_treefind_variable_threads(P[i].Pos,PPP[i].Hsml,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
                if(numngb_inbox>0)
                {
                    target_for_merger = -1;
                    threshold_val = MAX_REAL_NUMBER;
                    /* loop over neighbors */
                    for(n=0; n<numngb_inbox; n++)
                    {
                        j = Ngblist[n];
                        /* make sure we're not taking the same particle */
                        if((j>=0)&&(j!=i))
                        {
                            int k; double r2=0; for(k=0;k<3;k++) {r2+=(P[i].Pos[k]-P[j].Pos[k])*(P[i].Pos[k]-P[j].Pos[k]);}
                            if(r2<threshold_val) {threshold_val=r2; target_for_merger=j;} // position-based //
                            if(P[j].Mass<threshold_val) {threshold_val=P[j].Mass; target_for_merger=j;} // mass-based //
                        }
                    } // for(n=0; n<numngb_inbox; n++)
                    if(target_for_merger>=0)
                    {
                        /* some neighbors were found, we can true we're not going to crash the tree by splitting */
                        split_particle_i(i, n_particles_split,target_for_merger,threshold_val);
                        n_particles_split++;
                    }
                } // if(numngb_inbox>0)
            }
            /* alright, particle splitting operations are complete! */
        } // P[i].Type==0
    } // for(i = 0; i < NumPart; i++)
    myfree(Ngblist);
    MPI_Allreduce(&n_particles_merged, &MPI_n_particles_merged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("Particle split/merge check: %d particles merged, %d particles split \n",
               MPI_n_particles_merged,MPI_n_particles_split);
        fflush(stdout);
    }
    /* the reduction or increase of n_part by MPI_n_particles_merged will occur in rearrange_particle_sequence, which -must- 
        be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas += (long long)MPI_n_particles_split;
    Gas_split = n_particles_split; // specific to the local processor //
}


int does_particle_need_to_be_merged(MyIDType i)
{
    /* here we can insert any desired criteria for particle mergers: by default, this will occur 
        when particles fall below some minimum mass threshold */
    if(P[i].Mass <= 0) return 0;
    if(P[i].Mass <= All.MinMassForParticleMerger) return 1;
    return 0;
}

int does_particle_need_to_be_split(MyIDType i)
{
    /* here we can insert any desired criteria for particle splitting: by default, this will occur
        when particles become too massive, but it could also be done when Hsml gets very large, densities are high, etc */
    if(P[i].Mass >= All.MaxMassForParticleSplit) return 1;
    return 0;
}


void split_particle_i(MyIDType i, int n_particles_split, MyIDType i_nearest, double r2_nearest)
{
    double mass_of_new_particle;
    if(NumPart + n_particles_split >= All.MaxPart)
    {
        printf ("On Task=%d with NumPart=%d we try to split a particle. Sorry, no space left...(All.MaxPart=%d)\n", ThisTask, NumPart, All.MaxPart);
        fflush(stdout);
        endrun(8888);
    }
    
    /* here is where the details of the split are coded, the rest is bookkeeping */
    mass_of_new_particle = 0.5;
    
    int k;
    double phi = 2.0*M_PI*get_random_number(i+1+ThisTask); // random from 0 to 2pi //
    double cos_theta = 2.0*(get_random_number(i+3+2*ThisTask)-0.5); // random between 1 to -1 //
    double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    d_r = DMIN(d_r , 0.35 * sqrt(r2_nearest)); // use a 'buffer' to limit to some multiple of the distance to the nearest particle //
    
    /* find the first non-gas particle and move it to the end of the particle list */
    long j = NumPart + n_particles_split;
    /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */
    P[j] = P[i];
    SphP[j] = SphP[i];
    /* the particle needs to be 'born active' and added to the active set */
    NextActiveParticle[j] = FirstActiveParticle;
    FirstActiveParticle = j;
    NumForceUpdate++;
    /* likewise add it to the counters that register how many particles are in each timebin */
    TimeBinCount[P[j].TimeBin]++;
    PrevInTimeBin[j] = i;
    NextInTimeBin[j] = NextInTimeBin[i];
    if(NextInTimeBin[i] >= 0)
        PrevInTimeBin[NextInTimeBin[i]] = j;
    NextInTimeBin[i] = j;
    if(LastInTimeBin[P[i].TimeBin] == i)
        LastInTimeBin[P[i].TimeBin] = j;
    /* the particle needs an ID: we give it a bit-flip from the original particle to signify the split */
    unsigned int bits;
    int SPLIT_GENERATIONS = 10;
    for(bits = 0; SPLIT_GENERATIONS > (1 << bits); bits++);
    P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));
    /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
    SphP[i].ConditionNumber *= 10.0;
    SphP[j].ConditionNumber = SphP[i].ConditionNumber;
    /* assign masses to both particles (so they sum correctly) */
    P[j].Mass = mass_of_new_particle * P[i].Mass;
    P[i].Mass -= P[j].Mass;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    double dmass = mass_of_new_particle * SphP[i].DtMass;
    SphP[j].DtMass = dmass;
    SphP[i].DtMass -= dmass;
    dmass = mass_of_new_particle * SphP[i].dMass;
    SphP[j].dMass = dmass;
    SphP[i].dMass -= dmass;
    for(k=0;k<3;k++)
    {
        SphP[j].GravWorkTerm[k] = mass_of_new_particle * SphP[i].GravWorkTerm[k];
        SphP[i].GravWorkTerm[k] -= SphP[j].GravWorkTerm[k];
    }
    SphP[j].MassTrue = mass_of_new_particle * SphP[i].MassTrue;
    SphP[i].MassTrue -= SphP[j].MassTrue;
#endif
    
    
    /* shift the particle locations according to the random number we drew above */
    double dx, dy, dz;
#ifdef ONEDIM 
    dy=dz=0; dx=d_r; // here the split direction is trivial //
#else
    /* in 2D and 3D its not so trivial how to split the directions */
    double sin_theta = sqrt(1 - cos_theta*cos_theta);
    dx = d_r * sin_theta * cos(phi);
    dy = d_r * sin_theta * sin(phi);
    dz = d_r * cos_theta;
#ifdef TWODIMS
    dz=0; dx=d_r*cos(phi); dy=d_r*sin(phi);
#endif
    double norm=0, dp[3]; int m; dp[0]=dp[1]=dp[2]=0;
    for(k = 0; k < NUMDIMS; k++)
    {
        for(m = 0; m < NUMDIMS; m++) dp[k] += SphP[i].NV_T[k][m];
        dp[k] = SphP[i].Gradients.Density[k]; // ???
        norm += dp[k] * dp[k];
    }
    if(norm > 0)
    {
        norm = 1/sqrt(norm);
        for(k=0;k<NUMDIMS;k++) dp[k] *= norm;
        dx=d_r*dp[0]; dy=d_r*dp[1]; dz=d_r*dp[2];
        
        /* rotate to 90-degree offset from above orientation */
        if(dp[2]==1)
        {
            dx=d_r; dy=0; dz=0;
        } else {
            dz = sqrt(dp[1]*dp[1] + dp[0]*dp[0]);
            dx = -d_r * dp[1]/dz;
            dy = d_r * dp[0]/dz;
            dz = 0.0;
        }
    }
#endif
    
    P[i].Pos[0] += dx;
    P[j].Pos[0] -= dx;
    P[i].Pos[1] += dy;
    P[j].Pos[1] -= dy;
    P[i].Pos[2] += dz;
    P[j].Pos[2] -= dz;
    /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    force_add_star_to_tree(i, j);// (buggy)
    /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
}



void merge_particles_ij(MyIDType i, MyIDType j)
{
    /* routine to merge particle 'i' into particle 'j' */
    /* volumetric quantities (density, pressure, gradients, smoothing lengths)
     can be re-estimated after the merger operation. but we need to make sure
     all conserved quantities are appropriately dealt with */
    int k;
    if(P[i].Mass <= 0)
    {
        P[i].Mass = 0;
        return;
    }
    if(P[j].Mass <= 0)
    {
        P[j].Mass = 0;
        return;
    }
    double mtot = P[j].Mass + P[i].Mass;
    double wt_i = P[i].Mass / mtot;
    double wt_j = P[j].Mass / mtot;
    if(P[i].TimeBin < P[j].TimeBin)
    {
#ifdef WAKEUP
        SphP[j].wakeup = 1;
#endif
    }
    double dm_i=0,dm_j=0,de_i=0,de_j=0,dp_i[3],dp_j[3],dm_ij,de_ij,dp_ij[3];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    dm_i += SphP[i].DtMass;
    dm_j += SphP[j].DtMass;
#endif
    de_i = P[i].Mass * SphP[i].DtInternalEnergy + dm_i*SphP[i].InternalEnergy;
    de_j = P[j].Mass * SphP[j].DtInternalEnergy + dm_j*SphP[j].InternalEnergy;
    for(k=0;k<3;k++)
    {
        dp_i[k] = P[i].Mass * SphP[i].HydroAccel[k] + dm_i * SphP[i].VelPred[k] / All.cf_atime;
        dp_j[k] = P[j].Mass * SphP[j].HydroAccel[k] + dm_j * SphP[j].VelPred[k] / All.cf_atime;
        de_i += dp_i[k] * SphP[i].VelPred[k] / All.cf_atime - 0.5 * dm_i * SphP[i].VelPred[k] * SphP[i].VelPred[k] * All.cf_a2inv;
        de_j += dp_j[k] * SphP[j].VelPred[k] / All.cf_atime - 0.5 * dm_j * SphP[j].VelPred[k] * SphP[j].VelPred[k] * All.cf_a2inv;
        dp_ij[k] = dp_i[k] + dp_j[k];
    }
    dm_ij = dm_i+dm_j;
    de_ij = de_i+de_j;
    
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    SphP[j].MassTrue += SphP[i].MassTrue;
    SphP[i].MassTrue = 0;
    SphP[j].DtMass=dm_ij;
    SphP[i].DtMass=0;
    SphP[j].dMass = SphP[i].dMass + SphP[j].dMass;
    SphP[i].dMass = 0;
#endif

    /* make sure to update the conserved variables correctly: mass and momentum are easy, energy is non-trivial */
    double egy_old = 0, pos_new_xyz = 0;
    egy_old += mtot * (wt_j*SphP[j].InternalEnergy + wt_i*SphP[i].InternalEnergy); // internal energy //
    for(k=0;k<3;k++)
    {
        egy_old += mtot*wt_j * 0.5 * P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv; // kinetic energy (j) //
        egy_old += mtot*wt_i * 0.5 * P[i].Vel[k]*P[i].Vel[k]*All.cf_a2inv; // kinetic energy (i) //
        // gravitational energy terms need to be added (including work for moving particles 'together') //
        pos_new_xyz = wt_j*P[j].Pos[k] + wt_i*P[i].Pos[k];
        // Egrav = m*g*h = m * (-grav_acc) * (position relative to zero point) //
        egy_old += mtot*wt_j * (P[j].Pos[k] - pos_new_xyz)*All.cf_atime * (-P[j].GravAccel[k])*All.cf_a2inv; // work (j) //
        egy_old += mtot*wt_i * (P[i].Pos[k] - pos_new_xyz)*All.cf_atime * (-P[i].GravAccel[k])*All.cf_a2inv; // work (i) //
#ifdef PMGRID
        egy_old += mtot*wt_j * (P[j].Pos[k] - pos_new_xyz)*All.cf_atime * (-P[j].GravPM[k])*All.cf_a2inv; // work (j) [PMGRID] //
        egy_old += mtot*wt_i * (P[i].Pos[k] - pos_new_xyz)*All.cf_atime * (-P[i].GravPM[k])*All.cf_a2inv; // work (i) [PMGRID] //
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].GravWorkTerm[k] = 0; // since we're accounting for the work above and dont want to accidentally double-count //
#endif
    }
    
    
    SphP[j].InternalEnergy = wt_j*SphP[j].InternalEnergy + wt_i*SphP[i].InternalEnergy;
    SphP[j].InternalEnergyPred = wt_j*SphP[j].InternalEnergyPred + wt_i*SphP[i].InternalEnergyPred;
    double p_old_i[3],p_old_j[3];
    for(k=0;k<3;k++)
    {
        p_old_i[k] = P[i].Mass * P[i].Vel[k];
        p_old_j[k] = P[j].Mass * P[j].Vel[k];
    }
    for(k=0;k<3;k++)
    {
        P[j].Pos[k] = wt_j*P[j].Pos[k] + wt_i*P[i].Pos[k]; // center-of-mass conserving //
        P[j].Vel[k] = wt_j*P[j].Vel[k] + wt_i*P[i].Vel[k]; // momentum-conserving //
        SphP[j].VelPred[k] = wt_j*SphP[j].VelPred[k] + wt_i*SphP[i].VelPred[k]; // momentum-conserving //
        P[j].GravAccel[k] = wt_j*P[j].GravAccel[k] + wt_i*P[i].GravAccel[k]; // force-conserving //
#ifdef PMGRID
        P[j].GravPM[k] = wt_j*P[j].GravPM[k] + wt_i*P[i].GravPM[k]; // force-conserving //
#endif
    }

    /* correct our 'guess' for the internal energy with the residual from exact energy conservation */
    double egy_new = mtot * SphP[j].InternalEnergy;
    for(k=0;k<3;k++) {egy_new += mtot * 0.5*P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv;}
    egy_new = (egy_old - egy_new) / mtot; /* this residual needs to be put into the thermal energy */
    if(egy_new < -0.5*SphP[j].InternalEnergy) egy_new = -0.5 * SphP[j].InternalEnergy;
    SphP[j].InternalEnergy += egy_new;
    SphP[j].InternalEnergyPred += egy_new;
    if(SphP[j].InternalEnergyPred<0.5*SphP[j].InternalEnergy) SphP[j].InternalEnergyPred=0.5*SphP[j].InternalEnergy;
    
    
    // now use the conserved variables to correct the derivatives to primitive variables //
    de_ij -= dm_ij * SphP[j].InternalEnergyPred;
    for(k=0;k<3;k++)
    {
        SphP[j].HydroAccel[k] = (dp_ij[k] - dm_ij * SphP[j].VelPred[k]) / mtot;
        de_ij -= mtot * SphP[j].VelPred[k]/All.cf_atime * SphP[j].HydroAccel[k] + 0.5 * dm_ij * SphP[j].VelPred[k]*SphP[j].VelPred[k]*All.cf_a2inv;
    }
    SphP[j].DtInternalEnergy = de_ij;
    // to be conservative adopt the maximum signal velocity and smoothing length //
    SphP[j].MaxSignalVel = sqrt(SphP[j].MaxSignalVel*SphP[j].MaxSignalVel + SphP[i].MaxSignalVel*SphP[i].MaxSignalVel); /* need to be conservative */
    PPP[j].Hsml = pow(pow(PPP[j].Hsml,NUMDIMS)+pow(PPP[i].Hsml,NUMDIMS),1.0/NUMDIMS); /* sum the volume of the two particles */
    SphP[j].ConditionNumber = SphP[j].ConditionNumber + SphP[i].ConditionNumber; /* sum to be conservative */
    SphP[j].MaxKineticEnergyNgb = DMAX(SphP[j].MaxKineticEnergyNgb,SphP[i].MaxKineticEnergyNgb); /* for the entropy/energy switch condition */

    // below, we need to take care of additional physics //
#ifdef MAGNETIC
    for(k=0;k<3;k++)
    {
        SphP[j].B[k] = wt_j*SphP[j].B[k] + wt_i*SphP[i].B[k];
        SphP[j].BPred[k] = wt_j*SphP[j].BPred[k] + wt_i*SphP[i].BPred[k];
        SphP[j].DtB[k] = wt_j*SphP[j].DtB[k] + wt_i*SphP[i].DtB[k];
    }
#ifdef DIVBCLEANING_DEDNER
    SphP[j].Phi = wt_j*SphP[j].Phi + wt_i*SphP[i].Phi;
    SphP[j].PhiPred = wt_j*SphP[j].PhiPred + wt_i*SphP[i].PhiPred;
    SphP[j].DtPhi = wt_j*SphP[j].DtPhi + wt_i*SphP[i].DtPhi;
#endif
#endif
#ifdef EOS_DEGENERATE
    SphP[j].temp = wt_j*SphP[j].temp + wt_i*SphP[i].temp;
    SphP[j].dp_drho = wt_j*SphP[j].dp_drho + wt_i*SphP[i].dp_drho;
    for(k=0;k<EOS_NSPECIES;k++)
    {
        SphP[j].xnuc[k] = wt_j*SphP[j].xnuc[k] + wt_i*SphP[i].xnuc[k];
        SphP[j].dxnuc[k] = wt_j*SphP[j].dxnuc[k] + wt_i*SphP[i].dxnuc[k];
        SphP[j].xnucPred[k] = wt_j*SphP[j].xnucPred[k] + wt_i*SphP[i].xnucPred[k];
    }
#endif
#if defined(RADTRANSFER)
    for(k=0;k<6;k++) SphP[j].ET[k] = wt_j*SphP[j].ET[k] + wt_i*SphP[i].ET[k];
    for(k=0;k<N_BINS;k++) SphP[j].n_gamma[k] = wt_j*SphP[j].n_gamma[k] + wt_i*SphP[i].n_gamma[k];
#endif
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++)
        P[j].Metallicity[k] = wt_j*P[j].Metallicity[k] + wt_i*P[i].Metallicity[k]; /* metal-mass conserving */
#endif
    /* finally zero out the particle mass so it will be deleted */
    P[i].Mass = 0;
    P[j].Mass = mtot;
    for(k=0;k<3;k++)
    {
        /* momentum shift for passing to tree (so we know how to move it) */
        P[i].dp[k] += P[i].Mass*P[i].Vel[k] - p_old_i[k];
        P[j].dp[k] += P[j].Mass*P[j].Vel[k] - p_old_j[k];
    }
    return;
}



void rearrange_particle_sequence(void)
{
    int i, j, flag = 0, flag_sum;
    int count_elim, count_gaselim, count_bhelim, tot_elim, tot_gaselim, tot_bhelim;
    struct particle_data psave;
    struct sph_particle_data sphsave;
    
    int do_loop_check = 0;
    if(Gas_split>0)
    {
        N_gas += Gas_split;
        NumPart += Gas_split;
        Gas_split = 0;
        do_loop_check = 1;
    }
#ifdef GALSF
    if(Stars_converted)
    {
        N_gas -= Stars_converted;
        Stars_converted = 0;
        do_loop_check = 1;
    }
#endif
    if(NumPart <= N_gas) do_loop_check=0;
    if(N_gas <= 0) do_loop_check=0;
    
    /* if more gas than stars, need to be sure the block ordering is correct (gas first, then stars) */
    if(do_loop_check)
    {
        for(i = 0; i < N_gas; i++) /* loop over the gas block */
            if(P[i].Type != 0) /* and look for a particle converted to non-gas */
            {
                /* ok found a non-gas particle: */
                for(j = N_gas; j < NumPart; j++) /* loop from N_gas to Numpart, to find first labeled as gas */
                    if(P[j].Type == 0) break; /* break on that to record the j of interest */
                if(j >= NumPart) endrun(181170); /* if that j is too large, exit with error */
                
                psave = P[i]; /* otherwise, save the old pointer */
                P[i] = P[j]; /* now set the pointer equal to this new P[j] */
                P[j] = psave; /* now set the P[j] equal to the old, saved pointer */
                /* so we've swapped the two P[i] and P[j] */
                sphsave = SphP[i];
                SphP[i] = SphP[j];
                SphP[j] = sphsave;  /* have the gas particle take its sph pointer with it */
                /* ok we've now swapped the ordering so the gas particle is still inside the block */
                flag = 1;
            }
    }
    
    count_elim = 0;
    count_gaselim = 0;
    count_bhelim = 0;
    /* loop over entire block looking for things with zero mass, which need to be eliminated */
    for(i = 0; i < NumPart; i++)
        if(P[i].Mass <= 0)
        {
            P[i].Mass = 0;
            TimeBinCount[P[i].TimeBin]--;
            
            if(TimeBinActive[P[i].TimeBin])
                NumForceUpdate--;
            
            if(P[i].Type == 0)
            {
                TimeBinCountSph[P[i].TimeBin]--;
                
                P[i] = P[N_gas - 1];
                SphP[i] = SphP[N_gas - 1];
                /* swap with properties of last gas particle (i-- below will force a check of this so its ok) */
                
                P[N_gas - 1] = P[NumPart - 1]; /* redirect the final gas pointer to go to the final particle (BH) */
                N_gas--; /* shorten the total N_gas count */
                count_gaselim++; /* record that a BH was eliminated */
            }
            else
            {
                if(P[i].Type == 5) {count_bhelim++;} /* record elimination if BH */
                P[i] = P[NumPart - 1]; /* re-directs pointer for this particle to pointer at final particle -- so we
                                        swap the two; note that ordering -does not- matter among the non-SPH particles
                                        so its fine if this mixes up the list ordering of different particle types */
            }
            NumPart--;
            i--;
            count_elim++;
        }
    
    MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_bhelim, &tot_bhelim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if(count_elim)
        flag = 1;
    
    if(ThisTask == 0)
    {
        printf("Rearrange: Eliminated %d/%d gas/star particles and merged away %d black holes.\n",
               tot_gaselim, tot_elim - tot_gaselim - tot_bhelim, tot_bhelim);
        fflush(stdout);
    }
    
    All.TotNumPart -= tot_elim;
    All.TotN_gas -= tot_gaselim;
#ifdef BLACK_HOLES
    All.TotBHs -= tot_bhelim;
#endif
    
    MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(flag_sum)
        reconstruct_timebins();
}



/*  This function aborts the simulations. If a single processors
 *  wants an immediate termination,  the function needs to be
 *  called with ierr>0. A bunch of MPI-error messages will also
 *  appear in this case.
 *  For ierr=0, MPI is gracefully cleaned up, but this requires
 *  that all processors call endrun().
 */
void endrun(int ierr)
{
    if(ierr)
    {
        printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, ierr);
        exit(0);
    }
    
    MPI_Finalize();
    exit(0);
}


#ifdef DEBUG
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
  struct rlimit rlim;
  extern int feenableexcept(int __excepts);

  /* enable floating point exceptions */

  /*
     feenableexcept(FE_DIVBYZERO | FE_INVALID);
   */

  /* Note: FPU exceptions appear not to work properly
   * when the Intel C-Compiler for Linux is used
   */

  /* set core-dump size to infinity */
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.
   */
  /*
     signal(SIGSEGV, catch_fatal);
     signal(SIGBUS, catch_fatal);
     signal(SIGFPE, catch_fatal);
     signal(SIGINT, catch_fatal);
   */

  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);

  /* Establish a handler for SIGABRT signals. */
  signal(SIGABRT, catch_abort);
}


void catch_abort(int sig)
{
  MPI_Finalize();
  exit(0);
}

void catch_fatal(int sig)
{
  terminate_processes();
  MPI_Finalize();

  signal(sig, SIG_DFL);
  raise(sig);
}


void terminate_processes(void)
{
  pid_t my_pid;
  char buf[500], hostname[500], *cp;
  char commandbuf[500];
  FILE *fd;
  int i, pid;

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  my_pid = getpid();

  if((fd = fopen(buf, "r")))
    {
      for(i = 0; i < NTask; i++)
	{
	  int ret;

	  ret = fscanf(fd, "%s %d", hostname, &pid);

	  cp = hostname;
	  while(*cp)
	    {
	      if(*cp == '.')
		*cp = 0;
	      else
		cp++;
	    }

	  if(my_pid != pid)
	    {
	      sprintf(commandbuf, "ssh %s kill -ABRT %d", hostname, pid);
	      printf("--> %s\n", commandbuf);
	      fflush(stdout);
#ifndef NOCALLSOFSYSTEM
	      ret = system(commandbuf);
#endif
	    }
	}

      fclose(fd);
    }
}

void write_pid_file(void)
{
  pid_t my_pid;
  char mode[8], buf[500];
  FILE *fd;
  int i;

  my_pid = getpid();

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
	{
	  if(ThisTask == 0)
	    sprintf(mode, "w");
	  else
	    sprintf(mode, "a");

	  if((fd = fopen(buf, mode)))
	    {
	      fprintf(fd, "%s %d\n", getenv("HOST"), (int) my_pid);
	      fclose(fd);
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
}
#endif

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger()
{
  int continue_run = 0, my_pid;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  my_pid = getpid();
  printf("PID %d on %s ready for debugger to attach\n", my_pid, hostname);
  fflush(stdout);
  
  int *all_pids;
  if (ThisTask == 0) all_pids = (int *) malloc(sizeof(int)*NTask);
  
  MPI_Gather(&my_pid, 1, MPI_INT, all_pids, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (ThisTask == 0)
  {
    FILE *fd = fopen("pid_list_for_debugger.txt","w");
    int j;
    for (j=0; j<NTask; ++j) fprintf(fd, "%i \n", all_pids[j]);
    fclose(fd);   
    free(all_pids);
    printf("PID file written\n");
    fflush(stdout);
  }
  
  while (continue_run == 0) /* continue_run needsto be set to 1 ("set var continue_run = 1") by debugger to continue */
    sleep(1);
}
#endif

/*
double get_random_number(unsigned int id)
{
  return RndTable[(id % RNDTABLE)];
}
*/

double get_random_number(MyIDType id)
{
  return RndTable[(int) (id % RNDTABLE)];
}

void set_random_numbers(void)
{
  int i;

  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
}


/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double my_second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double) clock()) / CLOCKS_PER_SEC;
#endif

  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time(void)	/* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = my_second();
  dt = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

double report_time(void)       /* strategy: call this to measure sub-times of functions*/
{
  double t, dt;

  t = my_second();
  dt = t - WallclockTime;

  return dt;
}


/* returns the time difference between two measurements
 * obtained with my_second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}




#ifdef X86FIX

#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#define _FPU_EXTENDED 0x0300
#define _FPU_DOUBLE   0x0200

void x86_fix(void)
{
  unsigned short dummy, new_cw;
  unsigned short *old_cw;

  old_cw = &dummy;

  _FPU_GETCW(*old_cw);
  new_cw = (*old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
}

#endif


void minimum_large_ints(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
                MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = src[j];

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      if(res[j] > numlist[i * n + j])
        res[j] = numlist[i * n + j];

  myfree(numlist);
}

void sumup_large_ints(int n, int *src, long long *res)
{
  int i, j, *numlist;

  numlist = (int *) mymalloc("numlist", NTask * n * sizeof(int));
  MPI_Allgather(src, n, MPI_INT, numlist, n, MPI_INT, MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

void sumup_longs(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

size_t sizemax(size_t a, size_t b)
{
  if(a < b)
    return b;
  else
    return a;
}


void report_VmRSS(void)
{
  pid_t my_pid;
  FILE *fd;
  char buf[1024];

  my_pid = getpid();

  sprintf(buf, "/proc/%d/status", my_pid);

  if((fd = fopen(buf, "r")))
    {
      while(1)
	{
	  if(fgets(buf, 500, fd) != buf)
	    break;

	  if(strncmp(buf, "VmRSS", 5) == 0)
	    {
	      printf("ThisTask=%d: %s", ThisTask, buf);
	    }
	  if(strncmp(buf, "VmSize", 6) == 0)
	    {
	      printf("ThisTask=%d: %s", ThisTask, buf);
	    }
	}
      fclose(fd);
    }
}

long long report_comittable_memory(long long *MemTotal,
				   long long *Committed_AS, long long *SwapTotal, long long *SwapFree)
{
  FILE *fd;
  char buf[1024];

  if((fd = fopen("/proc/meminfo", "r")))
    {
      while(1)
	{
	  if(fgets(buf, 500, fd) != buf)
	    break;

	  if(bcmp(buf, "MemTotal", 8) == 0)
	    {
	      *MemTotal = atoll(buf + 10);
	    }
	  if(strncmp(buf, "Committed_AS", 12) == 0)
	    {
	      *Committed_AS = atoll(buf + 14);
	    }
	  if(strncmp(buf, "SwapTotal", 9) == 0)
	    {
	      *SwapTotal = atoll(buf + 11);
	    }
	  if(strncmp(buf, "SwapFree", 8) == 0)
	    {
	      *SwapFree = atoll(buf + 10);
	    }
	}
      fclose(fd);
    }

  return (*MemTotal - *Committed_AS);
}

void mpi_report_comittable_memory(long long BaseMem)
{
  long long *sizelist, maxsize[6], minsize[6];
  double avgsize[6];
  int i, imem, mintask[6], maxtask[6];
  long long Mem[6];
  char label[512];

  Mem[0] = report_comittable_memory(&Mem[1], &Mem[2], &Mem[3], &Mem[4]);
  Mem[5] = Mem[1] - Mem[0];

  for(imem = 0; imem < 6; imem++)
    {
      sizelist = (long long *) malloc(NTask * sizeof(long long));
      MPI_Allgather(&Mem[imem], sizeof(long long), MPI_BYTE, sizelist, sizeof(long long), MPI_BYTE,
		    MPI_COMM_WORLD);

      for(i = 1, mintask[imem] = 0, maxtask[imem] = 0, maxsize[imem] = minsize[imem] =
	  sizelist[0], avgsize[imem] = sizelist[0]; i < NTask; i++)
	{
	  if(sizelist[i] > maxsize[imem])
	    {
	      maxsize[imem] = sizelist[i];
	      maxtask[imem] = i;
	    }
	  if(sizelist[i] < minsize[imem])
	    {
	      minsize[imem] = sizelist[i];
	      mintask[imem] = i;
	    }
	  avgsize[imem] += sizelist[i];
	}

      free(sizelist);
    }

  if(ThisTask == 0)
    {
      printf("-------------------------------------------------------------------------------------------\n");
      for(imem = 0; imem < 6; imem++)
	{
	  switch (imem)
	    {
	    case 0:
	      sprintf(label, "AvailMem");
	      break;
	    case 1:
	      sprintf(label, "Total Mem");
	      break;
	    case 2:
	      sprintf(label, "Committed_AS");
	      break;
	    case 3:
	      sprintf(label, "SwapTotal");
	      break;
	    case 4:
	      sprintf(label, "SwapFree");
	      break;
	    case 5:
	      sprintf(label, "AllocMem");
	      break;
	    }
	  printf
	    ("%s:\t Largest = %10.2f Mb (on task=%d), Smallest = %10.2f Mb (on task=%d), Average = %10.2f Mb\n",
	     label, maxsize[imem] / (1024.), maxtask[imem], minsize[imem] / (1024.), mintask[imem],
	     avgsize[imem] / (1024. * NTask));
	}
      printf("-------------------------------------------------------------------------------------------\n");
    }
  if(ThisTask == maxtask[2])
    {
      printf("Task with the maximum commited memory");
      system("echo $HOST");
    }


  fflush(stdout);
}
