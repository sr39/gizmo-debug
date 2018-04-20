#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"

/* Routines for gas equation-of-state terms (collects things like calculation of gas pressure)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* this pair of functions: 'return_user_desired_target_density' and 'return_user_desired_target_pressure' should be used
 together with 'HYDRO_GENERATE_TARGET_MESH'. This will attempt to move the mesh and mass
 towards the 'target' pressure profile. Use this to build your ICs.
 The 'desired' pressure and density as a function of particle properties (most commonly, position) should be provided in the function below */
double return_user_desired_target_density(int i)
{
    return 1; // uniform density everywhere -- will try to generate a glass //
    /*
     // this example would initialize a constant-density (density=rho_0) spherical cloud (radius=r_cloud) with a smooth density 'edge' (width=interp_width) surrounded by an ambient medium of density =rho_0/rho_contrast //
     double dx=P[i].Pos[0]-boxHalf_X, dy=P[i].Pos[1]-boxHalf_Y, dz=P[i].Pos[2]-boxHalf_Z, r=sqrt(dx*dx+dy*dy+dz*dz);
     double rho_0=1, r_cloud=0.5*boxHalf_X, interp_width=0.1*r_cloud, rho_contrast=10.;
     return rho_0 * ((1.-1./rho_contrast)*0.5*erfc(2.*(r-r_cloud)/interp_width) + 1./rho_contrast);
     */
}
double return_user_desired_target_pressure(int i)
{
    return 1; // uniform pressure everywhere -- will try to generate a constant-pressure medium //
    /*
     // this example would initialize a radial pressure gradient corresponding to a self-gravitating, spherically-symmetric, infinite power-law
     //   density profile rho ~ r^(-b) -- note to do this right, you need to actually set that power-law for density, too, in 'return_user_desired_target_density' above
     double dx=P[i].Pos[0]-boxHalf_X, dy=P[i].Pos[1]-boxHalf_Y, dz=P[i].Pos[2]-boxHalf_Z, r=sqrt(dx*dx+dy*dy+dz*dz);
     double b = 2.; return 2.*M_PI/fabs((3.-b)*(1.-b)) * pow(return_user_desired_target_density(i),2) * r*r;
     */
}




/* return the pressure of particle i */
double get_pressure(int i)
{
    MyFloat press = GAMMA_MINUS1 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); /* ideal gas EOS (will get over-written it more complex EOS assumed) */
    
#ifdef GALSF_EFFECTIVE_EQS
    /* modify pressure to 'interpolate' between effective EOS and isothermal, with the Springel & Hernquist 2003 'effective' EOS */
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh)
        press = All.FactorForSofterEQS * press + (1 - All.FactorForSofterEQS) * All.cf_afac1 * GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;
#endif
    
    
#ifdef EOS_HELMHOLTZ
    /* pass the necessary quantities to wrappers for the Timms EOS */
    struct eos_input eos_in;
    struct eos_output eos_out;
    eos_in.rho  = SphP[i].Density;
    eos_in.eps  = SphP[i].InternalEnergyPred;
    eos_in.Ye   = SphP[i].Ye;
    eos_in.Abar = SphP[i].Abar;
    eos_in.temp = SphP[i].Temperature;
    int ierr = eos_compute(&eos_in, &eos_out);
    assert(!ierr);
    press              = eos_out.press;
    SphP[i].SoundSpeed = eos_out.csound;
    SphP[i].Temperature= eos_out.temp;
#endif

    
#ifdef EOS_TILLOTSON
    press = calculate_eos_tillotson(i);
#endif
    
#ifdef EOS_ENFORCE_ADIABAT
    press = EOS_ENFORCE_ADIABAT * pow(SphP[i].Density, GAMMA);
#endif
    
    
#ifdef COSMIC_RAYS
    press += Get_Particle_CosmicRayPressure(i);
    /* we will also compute the CR contribution to the effective soundspeed here */
    if((P[i].Mass > 0) && (SphP[i].Density>0) && (SphP[i].CosmicRayEnergyPred > 0))
    {
        SphP[i].SoundSpeed = sqrt(GAMMA*GAMMA_MINUS1 * SphP[i].InternalEnergyPred + GAMMA_COSMICRAY*GAMMA_COSMICRAY_MINUS1 * SphP[i].CosmicRayEnergyPred / P[i].Mass);
    } else {
        SphP[i].SoundSpeed = sqrt(GAMMA*GAMMA_MINUS1 * SphP[i].InternalEnergyPred);
    }
#endif
    
    
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
    /* add an artificial pressure term to suppress fragmentation at/below the explicit resolution scale */
    double h_eff = DMAX(Get_Particle_Size(i), All.ForceSoftening[0]/2.8); /* need to include latter to account for inter-particle spacing << grav soft cases */
    /* standard finite-volume formulation of this (note there is some geometric ambiguity about whether there should be a "pi" in the equation below, but this 
        can be completely folded into the (already arbitrary) definition of NJeans, so we just use the latter parameter */
    double NJeans = 4; // set so that resolution = lambda_Jeans/NJeans -- fragmentation with Jeans/Toomre scales below this will be artificially suppressed now
    double xJeans = (NJeans * NJeans / GAMMA) * All.G * h_eff*h_eff * SphP[i].Density * SphP[i].Density * All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
    SphP[i].SoundSpeed = sqrt(GAMMA * press / Particle_density_for_energy_i(i));
#endif
    
    
#if defined(HYDRO_GENERATE_TARGET_MESH)
    press = return_user_desired_target_pressure(i) * (SphP[i].Density / return_user_desired_target_density(i)); // define pressure by reference to 'desired' fluid quantities //
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = press / (GAMMA_MINUS1 * SphP[i].Density);
#endif
    
    return press;
}




/* trivial function to check if particle falls below the minimum allowed temperature */
void check_particle_for_temperature_minimum(int i)
{
    if(All.MinEgySpec)
    {
        if(SphP[i].InternalEnergy < All.MinEgySpec)
        {
            SphP[i].InternalEnergy = All.MinEgySpec;
            SphP[i].DtInternalEnergy = 0;
        }
    }
}



double INLINE_FUNC Particle_density_for_energy_i(int i)
{
#ifdef HYDRO_PRESSURE_SPH
    return SphP[i].EgyWtDensity;
#endif
    return SphP[i].Density;
}




double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{
#ifdef EOS_GENERAL
    return SphP[i].SoundSpeed;
#endif
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
}



#ifdef COSMIC_RAYS
double INLINE_FUNC Get_Particle_CosmicRayPressure(int i)
{
    if((P[i].Mass > 0) && (SphP[i].Density>0) && (SphP[i].CosmicRayEnergyPred > 0))
    {
        return GAMMA_COSMICRAY_MINUS1 * (SphP[i].CosmicRayEnergyPred * SphP[i].Density) / P[i].Mass; // cosmic ray pressure = (4/3-1) * e_cr = 1/3 * (E_cr/Vol) //
    } else {
        return 0;
    }
}


double Get_CosmicRayGradientLength(int i)
{
    /* now we need the cosmic ray pressure or energy density scale length, defined as :
        L = (e_cr + p_cr) / |gradient_p_cr| = cr_enthalpy / |gradient(p_cr)| */
    double CRPressureGradMag = 0.0;
    int k; for(k=0;k<3;k++) {CRPressureGradMag += SphP[i].Gradients.CosmicRayPressure[k]*SphP[i].Gradients.CosmicRayPressure[k];}
    /* limit the scale length: if too sharp, need a slope limiter at around the particle size */
    double L_gradient_min = Get_Particle_Size(i) * All.cf_atime;
    /* limit this scale length; if the gradient is too shallow, there is no information beyond a few smoothing lengths, so we can't let streaming go that far */
    double L_gradient_max = DMAX(200.*L_gradient_min, 100.0*PPP[i].Hsml*All.cf_atime);

    /* also, physically, cosmic rays cannot stream/diffuse with a faster coefficient than ~v_max*L_mean_free_path, where L_mean_free_path ~ 2.e20 * (cm^-3/n) */
    double nH_cgs = SphP[i].Density * All.cf_a3inv * ( All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam ) / PROTONMASS ;
    double L_mean_free_path = (3.e25 / nH_cgs) / (All.UnitLength_in_cm / All.HubbleParam);
    L_gradient_max = DMIN(L_gradient_max, L_mean_free_path);
    
    double CRPressureGradScaleLength = Get_Particle_CosmicRayPressure(i) / sqrt(1.0e-33 + CRPressureGradMag) * All.cf_atime;
    if(CRPressureGradScaleLength > 0) {CRPressureGradScaleLength = 1.0/(1.0/CRPressureGradScaleLength + 1.0/L_gradient_max);} else {CRPressureGradScaleLength=0;}
    CRPressureGradScaleLength = sqrt(L_gradient_min*L_gradient_min + CRPressureGradScaleLength*CRPressureGradScaleLength);
    return CRPressureGradScaleLength; /* this is returned in -physical- units */
}

double Get_CosmicRayStreamingVelocity(int i)
{
    /* in the weak-field (high-beta) case, the streaming velocity is approximately the sound speed */
    double v_streaming = sqrt(GAMMA*GAMMA_MINUS1 * SphP[i].InternalEnergyPred); // thermal ion sound speed //
#ifdef MAGNETIC
    /* in the strong-field (low-beta) case, it's actually the Alfven velocity: interpolate between these */
    double vA_2 = 0.0; double cs_stream = v_streaming;
    int k; for(k=0;k<3;k++) {vA_2 += Get_Particle_BField(i,k)*Get_Particle_BField(i,k);}
    vA_2 *= All.cf_afac1 / (All.cf_atime * SphP[i].Density);
    v_streaming = DMIN(1.0e6*cs_stream, sqrt(cs_stream*cs_stream + vA_2));
#endif
    v_streaming *= All.cf_afac3; // converts to physical units and rescales according to chosen coefficient //
    return v_streaming;
}


double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode)
{
    /* routine to do the drift/kick operations for CRs: mode=0 is kick, mode=1 is drift */
    if(dt_entr <= 0) {return 0;} // no update
    int k;
    double eCR, u0;
    if(mode==0) {eCR=SphP[i].CosmicRayEnergy; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred; u0=SphP[i].InternalEnergyPred;} // initial energy
    if(eCR < 0) {eCR=0;} // limit to physical values
    
#ifdef COSMIC_RAYS_M1
    // this is the exact solution for the CR flux-update equation over a finite timestep dt:
    //   it needs to be solved this way [implicitly] as opposed to explicitly for dt because
    //   in the limit of dt_cr_dimless being large, the problem exactly approaches the diffusive solution
    double DtCosmicRayFlux[3]={0}, flux[3]={0}, CR_veff[3]={0}, CR_vmag=0, q_cr = 0, cr_speed = COSMIC_RAYS_M1;// * (C/All.UnitVelocity_in_cm_per_s);
    cr_speed = DMAX( All.cf_afac3*SphP[i].MaxSignalVel , DMIN(COSMIC_RAYS_M1 , fabs(SphP[i].CosmicRayDiffusionCoeff)/(Get_Particle_Size(i)*All.cf_atime)));// * (C/All.UnitVelocity_in_cm_per_s);
    for(k=0;k<3;k++) {DtCosmicRayFlux[k] = -fabs(SphP[i].CosmicRayDiffusionCoeff) * (P[i].Mass/SphP[i].Density) * (SphP[i].Gradients.CosmicRayPressure[k]/GAMMA_COSMICRAY_MINUS1);}
#ifdef MAGNETIC // do projection onto field lines
    double B0[3]={0}, Bmag2=0, DtCRDotBhat=0;
    for(k=0;k<3;k++)
    {
        if(mode==0) {B0[k]=SphP[i].B[k];} else {B0[k]=SphP[i].BPred[k];}
        DtCRDotBhat += DtCosmicRayFlux[k] * B0[k]; Bmag2 += B0[k]*B0[k];
    }
    if(Bmag2 > 0) {for(k=0;k<3;k++) {DtCosmicRayFlux[k] = DtCRDotBhat * B0[k] / Bmag2;}}
#endif
    double dt_cr_dimless = dt_entr * cr_speed*cr_speed * GAMMA_COSMICRAY_MINUS1 / (MIN_REAL_NUMBER + fabs(SphP[i].CosmicRayDiffusionCoeff));
    if((dt_cr_dimless > 0)&&(dt_cr_dimless < 20.)) {q_cr = exp(-dt_cr_dimless);} // factor for CR interpolation
    if(mode==0) {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFlux[k];}} else {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFluxPred[k];}}
    for(k=0;k<3;k++) {flux[k] = q_cr*flux[k] + (1.-q_cr)*DtCosmicRayFlux[k];} // updated flux
    for(k=0;k<3;k++) {CR_veff[k]=flux[k]/(eCR+MIN_REAL_NUMBER); CR_vmag+=CR_veff[k]*CR_veff[k];} // effective streaming speed
    if((CR_vmag <= 0) || (isnan(CR_vmag))) // check for valid numbers
    {
        for(k=0;k<3;k++) {flux[k]=0; for(k=0;k<3;k++) {CR_veff[k]=0;}} // zero if invalid
    } else {
        CR_vmag = sqrt(CR_vmag);
        if(CR_vmag > cr_speed) {for(k=0;k<3;k++) {flux[k]*=cr_speed/CR_vmag; CR_veff[k]*=cr_speed/CR_vmag;}} // limit flux to free-streaming speed [as with RT]
    }
    if(mode==0) {for(k=0;k<3;k++) {SphP[i].CosmicRayFlux[k]=flux[k];}} else {for(k=0;k<3;k++) {SphP[i].CosmicRayFluxPred[k]=flux[k];}}
#endif
    
    // now we update the CR energies. since this is positive-definite, some additional care is needed //
    double dCR = SphP[i].DtCosmicRayEnergy*dt_entr, dCRmax = 1.e10*(eCR+MIN_REAL_NUMBER);
#ifdef GALSF
    dCRmax = DMAX(2.0*eCR , 0.1*u0*P[i].Mass);
#endif
    if(dCR > dCRmax) {dCR=dCRmax;} // don't allow excessively large values
    if(dCR < -eCR) {dCR=-eCR;} // don't allow it to go negative
    eCR += dCR; if((eCR<0)||(isnan(eCR))) {eCR=0;}
    if(mode==0) {SphP[i].CosmicRayEnergy=eCR;} else {SphP[i].CosmicRayEnergyPred=eCR;} // updated energy
    double eCR_0 = eCR; // save this value for below
    
    /* now need to account for the adiabatic heating/cooling of the cosmic ray fluid, here: its an ultra-relativistic fluid with gamma=4/3 */
    double d_div = (-GAMMA_COSMICRAY_MINUS1 * P[i].Particle_DivVel*All.cf_a2inv) * dt_entr;
    /* adiabatic term from Hubble expansion (needed for cosmological integrations */
    if(All.ComovingIntegrationOn) {d_div += (-3.*GAMMA_COSMICRAY_MINUS1 * All.cf_hubble_a) * dt_entr;}
    
#if 0 //def COSMIC_RAYS_M1 // per recent discussions, this is not matched to standard diffusive approaches, although question of whether this transport does trigger streaming-like instabilities remains (physically) open
    // need to get pressure gradient scale-lengths to slope-limit pressure gradient to resovable values //
    double Pmag=0; for(k=0;k<3;k++) {Pmag += SphP[i].Gradients.CosmicRayPressure[k]*SphP[i].Gradients.CosmicRayPressure[k];} // magnitude of pressure gradient
    if((Pmag<=0)||(isnan(Pmag))) {Pmag=0;} else {Pmag=sqrt(Pmag);} // check for unphysical values
    double CR_P = Get_Particle_CosmicRayPressure(i); // CR pressure to normalize gradient (since only 1/scale actually matters)
    double Pscale=CR_P/Pmag, Pscale_min=Get_Particle_Size(i), Pgrad_hat[3]={0}, fluxterm=0; // values for slope-limiter
    for(k=0;k<3;k++) {Pgrad_hat[k] = GAMMA_COSMICRAY_MINUS1 * SphP[i].Gradients.CosmicRayPressure[k] / CR_P;} // pressure gradient / pressure
    if(Pscale_min > Pscale) {for(k=0;k<3;k++) {Pgrad_hat[k] *= Pscale/Pscale_min;}} // slope-limited gradient
    for(k=0;k<3;k++) {fluxterm += dt_entr * CR_veff[k] * (Pgrad_hat[k]/All.cf_atime);} // calculate the adiabatic term [physical units]
    double fluxterm_isotropic = -DMIN(cr_speed,CR_vmag) / (DMAX(Pscale,Pscale_min)*All.cf_atime) * dt_entr; // dissipative flux if the streaming were exactly along the CR pressure gradient [isotropic diffusion limit]
    if(fluxterm > 0) {fluxterm = -DMIN(fabs(fluxterm), fabs(fluxterm_isotropic));} // limit anisotropic compression term since this should come only from integration error
    if(isnan(fluxterm)||(eCR<=0)||(isnan(eCR))) {fluxterm=0;} // check against nans
    d_div += fluxterm;
#endif
    
    double dCR_div = DMIN(eCR*d_div , 0.5*u0*P[i].Mass); // limit so don't take away all the gas internal energy [to negative values]
    if(dCR_div + eCR < 0) {dCR_div = -eCR;}
    eCR += dCR_div; if((eCR<0)||(isnan(eCR))) {eCR=0;}
    dCR_div = eCR - eCR_0; // actual change that is going to be applied
    if(mode==0)
    {
        SphP[i].CosmicRayEnergy += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;
    } else {
        SphP[i].CosmicRayEnergyPred += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;
    }
    return 1;
}
#endif


