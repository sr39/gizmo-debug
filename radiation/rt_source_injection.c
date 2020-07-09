#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file rt_source_injection.c
 *  \brief inject luminosity from point sources to neighboring gas particles
 *
 *  This file contains a loop modeled on the gas density computation which will
 *    share luminosity from non-gas particles to the surrounding gas particles,
 *    so that it can be treated within e.g. the flux-limited diffusion or other
 *    radiation-hydrodynamic approximations. Basically the same concept as
 *    injecting the radiation 'in a cell' surrounding the particle
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef RT_SOURCE_INJECTION


#define MASTER_FUNCTION_NAME rt_sourceinjection_evaluate /*! name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int MASTER_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(rt_sourceinjection_active_check(i)) /*! function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /*! pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/*! Structure for communication during the kernel computation. Holds data that is sent to other processors  */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3]; MyFloat Hsml, KernelSum_Around_RT_Source, Luminosity[N_RT_FREQ_BINS], Vel[3];
    int NodeList[NODELISTLENGTH];
#ifdef RT_AREAWEIGHT_INJECTION
    MyDouble V_i;
#endif
#if defined(RT_REPROCESS_INJECTED_PHOTONS) && defined(RT_CHEM_PHOTOION)
    MyDouble RHII;
#endif
}
*DATAIN_NAME, *DATAGET_NAME;

/*! subroutine to insert the data needed to be passed to other processors: here for convenience, match to structure above  */
void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    for(k=0; k<3; k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = PPP[i].Hsml;
    //if(P[i].Type==0) {in->KernelSum_Around_RT_Source = SphP[i].Density;} else {in->KernelSum_Around_RT_Source = P[i].DensAroundStar;}
    in->KernelSum_Around_RT_Source = P[i].KernelSum_Around_RT_Source;
#ifdef RT_AREAWEIGHT_INJECTION
    in->V_i = 4.*M_PI/3. * PPP[i].Hsml * PPP[i].Hsml * PPP[i].Hsml / (All.DesNumNgb);
#ifdef BLACK_HOLES
    if(P[i].Type == 5) in->V_i /= All.BlackHoleNgbFactor;
#endif
#endif
    /* luminosity is set to zero here for gas particles because their self-illumination is handled trivially in a single loop, earlier */
    double lum[N_RT_FREQ_BINS];
    int active_check = rt_get_source_luminosity(i,0,lum);
    double dt = 1; // make this do nothing unless flags below are set:
#if defined(RT_INJECT_PHOTONS_DISCRETELY)
    dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#if defined(RT_EVOLVE_FLUX)
    for(k=0; k<3; k++) {if(P[i].Type==0) {in->Vel[k] = SphP[i].VelPred[k];} else {in->Vel[k] = P[i].Vel[k];}}
#endif
#endif
    for(k=0; k<N_RT_FREQ_BINS; k++) {if(P[i].Type==0 || active_check==0) {in->Luminosity[k]=0;} else {in->Luminosity[k] = lum[k] * dt;}}
#if defined(RT_REPROCESS_INJECTED_PHOTONS) && defined(RT_CHEM_PHOTOION)
    double stellum = lum[RT_FREQ_BIN_H0] * UNIT_LUM_IN_CGS;
    double rho = P[i].DensAroundStar;
    in->RHII = 4.78e-9*pow(stellum,0.333)*pow(rho*All.cf_a3inv*UNIT_DENSITY_IN_CGS,-0.66667);
    in->RHII /= All.cf_atime*UNIT_LENGTH_IN_CGS;
#endif
}


/*! this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME
{ /* define variables below as e.g. "double X;" */
}
*DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/*! this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{  /* "i" is the particle to which data from structure "out" will be assigned. mode=0 for local communication,
    =1 for data sent back from other processors. you must account for this. */
    /* example: ASSIGN_ADD(P[i].X,out->X,mode); which is short for: if(mode==0) {P[i].X=out->X;} else {P[i].X+=out->X;} */
}



/*! determine if an element is active as a source */
int rt_sourceinjection_active_check(int i);
int rt_sourceinjection_active_check(int i)
{
    if(PPP[i].NumNgb <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(P[i].Mass <= 0) return 0;
    double lum[N_RT_FREQ_BINS];
    return rt_get_source_luminosity(i,-1,lum);
}


/*! operations that need to be performed before entering the main loop */
void rt_source_injection_initial_operations_preloop(void);
void rt_source_injection_initial_operations_preloop(void)
{
    /* first, we do a loop over the gas particles themselves. these are trivial -- they don't need to share any information,
     they just determine their own source functions. so we don't need to do any loops. and we can zero everything before the loop below. */
    int j;
    for(j=0;j<NumPart;j++) {
        if(P[j].Type==0) {
            double lum[N_RT_FREQ_BINS]; int k;
            for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[j].Rad_Je[k]=0;} // need to zero -before- calling injection //
            int active_check = rt_get_source_luminosity(j,0,lum);
            /* here is where we would need to code some source luminosity for the gas */
            for(k=0;k<N_RT_FREQ_BINS;k++) if(active_check) {SphP[j].Rad_Je[k]=lum[k];}
        }
    }
}




/*! subroutine that actually distributes the luminosity as desired to neighbor particles in the kernel */
int rt_sourceinjection_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    /* Load the data for the particle */
    int j, k, n, startnode, numngb_inbox, listindex = 0;
    struct INPUT_STRUCT_NAME local;
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    /* basic calculations */
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    double hinv, hinv3, hinv4, h2=local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &hinv, &hinv3, &hinv4);
    
    /* Now start the actual operations for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;/* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -1;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                double dp[3]; for(k=0; k<3; k++) {dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); /* find the closest image in the given box size  */
                double r2=0,r,c_light_eff; for(k=0;k<3;k++) {r2 += dp[k]*dp[k];}
                if(r2<=0) continue; // same particle //
                if(r2>=h2) continue; // outside kernel //
                // calculate kernel quantities //
                r = sqrt(r2); c_light_eff = C_LIGHT_CODE_REDUCED;
#ifdef RT_AREAWEIGHT_INJECTION
		double u=r*hinv, wk, dwk;
		if(u>0) kernel_main(u, hinv3, hinv4, &wk, &dwk, 1);
		double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j; /* note these lines and many below assume 3D sims! */
		double wk_j = 0, dwk_j = 0, u_j = r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = P[j].Mass / SphP[j].Density;

		if(u_j<1) {kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);} else {wk_j=dwk_j=0;}
		if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}           	       
		if(V_j<0 || isnan(V_j)) {V_j=0;}
		double sph_area = fabs(local.V_i*local.V_i*dwk + V_j*V_j*dwk_j); // effective face area //
		wk = (1 - 1/sqrt(1 + sph_area / (M_PI*r2))) / local.KernelSum_Around_RT_Source; // corresponding geometric weight //
#else
                double wk = (1 - r2*hinv*hinv) / local.KernelSum_Around_RT_Source;
#endif
#if defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
                double dv0 = -1. / (c_light_eff * r) * All.cf_atime;
                double lmax_0 = DMAX(local.Hsml, r);
#ifdef RT_EVOLVE_INTENSITIES
                int kx; double angle_wt_Inu_sum=0, angle_wt_Inu[N_RT_INTENSITY_BINS];
                // pre-compute a set of weights based on the projection of the particle position along the radial direction for the radiation direction //
                for(kx=0;kx<N_RT_INTENSITY_BINS;kx++)
                {
                    double cos_t=0; int kq; for(kq=0;kq<3;kq++) {cos_t+=All.Rad_Intensity_Direction[kx][kq]*dp[kq]/r;}
                    double wt_function = cos_t*cos_t*cos_t*cos_t; if(cos_t < 0) {wt_function=0;}
                    angle_wt_Inu[kx] = wt_function; angle_wt_Inu_sum += angle_wt_Inu[kx];
                }
#endif
#endif
                // now actually apply the kernel distribution
                for(k=0;k<N_RT_FREQ_BINS;k++) 
                {
                    double dE = wk * local.Luminosity[k];
#if defined(RT_INJECT_PHOTONS_DISCRETELY)
#if defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
                    // add discrete photon momentum from un-resolved absorption //
                    double x_abs = 2. * SphP[j].Rad_Kappa[k] * (SphP[j].Density*All.cf_a3inv) * (DMAX(2.*Get_Particle_Size(j),lmax_0)*All.cf_atime); // effective optical depth through particle
                    double slabfac_x = x_abs * slab_averaging_function(x_abs); // 1-exp(-x)
                    if(isnan(slabfac_x)||(slabfac_x<=0)) {slabfac_x=0;}
                    if(slabfac_x>1) {slabfac_x=1;}
                    double dv = slabfac_x * dv0 * dE / P[j].Mass; // total absorbed momentum (needs multiplication by dp[kv] for directionality)
                    int kv; for(kv=0;kv<3;kv++) {P[j].Vel[kv] += dv*dp[kv]; SphP[j].VelPred[kv] += dv*dp[kv];}

#ifdef RT_REPROCESS_INJECTED_PHOTONS //conserving photon energy, put only the un-absorbed component of the current band into that band, putting the rest in its "donation" bin (ionizing->optical, all others->IR). This would happen anyway during the routine for resolved absorption, but this may more realistically handle situations where e.g. your dust destruction front is at totally unresolved scales and you don't want to spuriously ionize stuff on larger scales. Assume isotropic re-radiation, so inject only energy for the donated bin and not net flux/momentum.
		    int donation_bin;
		    double dE_donation;
		    donation_bin = rt_get_donation_target_bin(k);
#ifdef RT_CHEM_PHOTOION
		    if((k!=RT_FREQ_BIN_H0) || (r < local.RHII)) // don't inject ionizing photons outside the Stromgren radius
#endif
		    {
		        dE_donation = fabs(1-slabfac_x)*dE;
                        dE *= slabfac_x;
		    }
#ifdef RT_CHEM_PHOTOION
		    else {
                        dE_donation = 0;
		    }    
#endif
#endif
#if defined(RT_EVOLVE_FLUX)
                    double dflux = -dE * c_light_eff / r;
                    for(kv=0;kv<3;kv++) {SphP[j].Rad_Flux[k][kv] += dflux*dp[kv]; SphP[j].Rad_Flux_Pred[k][kv] += dflux*dp[kv];}
#endif
#ifdef RT_EVOLVE_INTENSITIES
                    double dflux = dE / angle_wt_Inu_sum;
                    for(kv=0;kv<N_RT_INTENSITY_BINS;kv++) {SphP[j].Rad_Intensity[k][kv] += dflux * angle_wt_Inu[N_RT_INTENSITY_BINS]; SphP[j].Rad_Intensity_Pred[k][kv] += dflux * angle_wt_Inu[N_RT_INTENSITY_BINS];}
#endif
#endif // local extinction-corrected version gets the 'full' thin flux above: more general formulation allows these to build up self-consistently, since we don't know what the flux 'should' be in fact
                    SphP[j].Rad_E_gamma[k] += dE;
#ifdef RT_REPROCESS_INJECTED_PHOTONS
		    if(donation_bin > -1) {SphP[j].Rad_E_gamma[donation_bin] += dE_donation;}
#endif
#ifdef RT_EVOLVE_ENERGY
                    SphP[j].Rad_E_gamma_Pred[k] += dE; // dump discreetly (noisier, but works smoothly with large timebin hierarchy)
#ifdef RT_REPROCESS_INJECTED_PHOTONS
		    if(donation_bin > -1) {SphP[j].Rad_E_gamma_Pred[donation_bin] += dE_donation;}
#endif
#endif                    
#if defined(RT_EVOLVE_FLUX) // add relativistic corrections here, which should be there in general. however we will ignore [here] the 'back-reaction' term, since we're assuming the source is a star or something like that, where this would be negligible. gas self gain/loss is handled separately.
                    {int kv; for(kv=0;kv<3;kv++) {SphP[j].Rad_Flux[k][kv] += dE*local.Vel[kv]/All.cf_atime; SphP[j].Rad_Flux_Pred[k][kv] += dE*local.Vel[kv]/All.cf_atime;}}
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                    {double dflux=dE*C_LIGHT_CODE_REDUCED; SphP[j].Rad_Flux_Pred[k][2]+=dflux; SphP[j].Rad_Flux[k][2]+=dflux;}
#endif
#endif
#else // end discrete injection clause
                    SphP[j].Rad_Je[k] += dE; // treat continuously
#endif
                }
#ifdef SINGLE_STAR_FB_RAD                
                SphP[j].wakeup = 1; NeedToWakeupParticles_local = 1;
#endif                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DATAGET_NAME[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1)
#endif
    } // while(startnode >= 0)
    return 0;
} 



/*! routine to do the master loop over particles, for the source injection (photons put into surrounding gas) */
void rt_source_injection(void)
{
    PRINT_STATUS(" ..injecting radiation onto grid for RHD steps");
    rt_source_injection_initial_operations_preloop(); /* operations before the main loop */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    CPU_Step[CPU_RTNONFLUXOPS] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */





#endif
