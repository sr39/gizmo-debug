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
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh) {press = All.FactorForSofterEQS * press + (1 - All.FactorForSofterEQS) * All.cf_afac1 * GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;}
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
#ifdef COSMIC_RAYS_ALFVEN
    press += (GAMMA_ALFVEN_CRS-1) * SphP[i].Density * (SphP[i].CosmicRayAlfvenEnergy[0]+SphP[i].CosmicRayAlfvenEnergy[1]);
    SphP[i].SoundSpeed = sqrt(SphP[i].SoundSpeed*SphP[i].SoundSpeed + GAMMA_ALFVEN_CRS*(GAMMA_ALFVEN_CRS-1)*(SphP[i].CosmicRayAlfvenEnergy[0]+SphP[i].CosmicRayAlfvenEnergy[1])/P[i].Mass);
#endif
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
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = return_user_desired_target_pressure(i) / (GAMMA_MINUS1 * SphP[i].Density);
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
    /* now we need the -parallel- cosmic ray pressure or energy density scale length */
    double CRPressureGradMag = 0.0;
    int k; for(k=0;k<3;k++) {CRPressureGradMag += SphP[i].Gradients.CosmicRayPressure[k]*SphP[i].Gradients.CosmicRayPressure[k];}
    CRPressureGradMag = sqrt(1.e-46 + CRPressureGradMag); // sqrt to make absolute value
#ifdef MAGNETIC /* with anisotropic transport, we really want the -parallel- gradient scale-length, so need another factor here */
    double B2_tot=0.0, b=0; CRPressureGradMag=0; for(k=0;k<3;k++) {b=Get_Particle_BField(i,k); B2_tot+=b*b; CRPressureGradMag+=b*SphP[i].Gradients.CosmicRayPressure[k];} // note, this is signed!
    CRPressureGradMag = sqrt((1.e-40 + CRPressureGradMag*CRPressureGradMag) / (1.e-46 + B2_tot)); // divide B-magnitude to get scalar magnitude, and take sqrt[(G.P)^2] to get absolute value
#endif
    
    /* limit the scale length: if too sharp, need a slope limiter at around the particle size */
    double L_gradient_min = Get_Particle_Size(i) * All.cf_atime;
    /* limit this scale length; if the gradient is too shallow, there is no information beyond a few smoothing lengths, so we can't let streaming go that far */
    double L_gradient_max = DMAX(1000.*L_gradient_min, 500.0*PPP[i].Hsml*All.cf_atime);

    /* also, physically, cosmic rays cannot stream/diffuse with a faster coefficient than ~v_max*L_mean_free_path, where L_mean_free_path ~ 2.e20 * (cm^-3/n) */
    double nH_cgs = SphP[i].Density * All.cf_a3inv * ( All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam ) / PROTONMASS ;
    double L_mean_free_path = (3.e25 / nH_cgs) / (All.UnitLength_in_cm / All.HubbleParam);
    L_gradient_max = DMIN(L_gradient_max, L_mean_free_path);
    
    double CRPressureGradScaleLength = Get_Particle_CosmicRayPressure(i) / CRPressureGradMag * All.cf_atime;
    if(CRPressureGradScaleLength > 0) {CRPressureGradScaleLength = 1.0/(1.0/CRPressureGradScaleLength + 1.0/L_gradient_max);} else {CRPressureGradScaleLength=0;}
    CRPressureGradScaleLength = sqrt(L_gradient_min*L_gradient_min + CRPressureGradScaleLength*CRPressureGradScaleLength);
    return CRPressureGradScaleLength; /* this is returned in -physical- units */
}


double Get_Gas_Ionized_Fraction(int i)
{
#ifdef COOLING
    double ne=SphP[i].Ne, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    if(SphP[i].DelayTimeHII>0) {nh0=0;} // account for our effective ionization model here
#endif
    double f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
    return f_ion;
#endif
    return 1;
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
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_2 /= Get_Gas_Ionized_Fraction(i); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    v_streaming = DMIN(1.0e6*cs_stream, sqrt(MIN_REAL_NUMBER*cs_stream*cs_stream + vA_2)); // limit to Alfven speed //
#endif
#ifdef COSMIC_RAYS_M1
    v_streaming = DMIN(v_streaming , COSMIC_RAYS_M1); // limit to maximum transport speed //
#endif
    v_streaming *= All.cf_afac3; // converts to physical units and rescales according to chosen coefficient //
    return v_streaming;
}


void CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(int i)
{
    double CRPressureGradScaleLength = Get_CosmicRayGradientLength(i), CR_kappa_streaming = 0, Z_charge_CR, E_CRs_GeV, E_CRs_GeV_over_Z; SphP[i].CosmicRayDiffusionCoeff=0; int k;
    Z_charge_CR = 1, E_CRs_GeV = 1, E_CRs_GeV_over_Z = E_CRs_GeV/Z_charge_CR; // charge and energy and resonant Alfven wavenumber (in gyro units) of the CR population we're evolving
#ifndef COSMIC_RAYS_DISABLE_STREAMING /* self-consistently calculate the diffusion coefficients for cosmic ray fluids; first the streaming part of this (kappa~v_stream*L_CR_grad) following e.g. Wentzel 1968, Skilling 1971, 1975, Holman 1979, as updated in Kulsrud 2005, Yan & Lazarian 2008, Ensslin 2011 */
    double v_streaming = Get_CosmicRayStreamingVelocity(i);
    CR_kappa_streaming = GAMMA_COSMICRAY * v_streaming * CRPressureGradScaleLength; /* the diffusivity is now just the product of these two coefficients (all physical units) */
#endif
        
    /* calculate a bunch of different properties for various diffusion models used below */
    double p_scale=0, unit_kappa_code=0, b_muG=0, f_ion=1, temperature=0; unit_kappa_code=All.UnitVelocity_in_cm_per_s*All.UnitLength_in_cm/All.HubbleParam; b_muG=sqrt( SphP[i].Pressure*All.cf_a3inv*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam / 4.0e-14 ); // this assumes beta=1, if no MHD enabled
#ifdef MAGNETIC /* get actual B-field */
    double gizmo2gauss = sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam), b2_mag = 0.0;
    for(k=0;k<3;k++) {b2_mag += Get_Particle_BField(i,k) * Get_Particle_BField(i,k);}
    b_muG=sqrt(DMAX(b2_mag,0)) * All.cf_a2inv * gizmo2gauss / 1.0e-6; b_muG = sqrt(b_muG*b_muG + 1.e-6); /* B-field in units of physical microGauss */
#endif
    /* pressure gradient scale length [with various limiters] for e.g. estimate of cascade driving length */
    p_scale = 0.0; for(k=0;k<3;k++) {p_scale += SphP[i].Gradients.Pressure[k]*SphP[i].Gradients.Pressure[k];}
    p_scale = SphP[i].Pressure / (1.e-33 + sqrt(p_scale)); double p_scale_min = 0.5 * Get_Particle_Size(i); // sets a 'floor' at some multiple of the particle size (unresolved below this) //
    p_scale = sqrt(p_scale_min*p_scale_min + p_scale*p_scale); double p_scale_max = 1000.*PPP[i].Hsml; // sets a maximum beyond which we have no meaningful information //
    double codelength_to_kpc = (All.UnitLength_in_cm / All.HubbleParam) / (3.086e21); /* unit conversion to kpc [for physical units] */
    p_scale = 1./(1./p_scale + 1./p_scale_max); p_scale *= codelength_to_kpc / All.cf_atime; /* physical pressure scale length in units of kpc */
    if(p_scale > 1000.) {p_scale=1000.;} // limit at 1000 kpc
#ifdef COOLING
    double ne=SphP[i].Ne, nh0=0, nHe0=0, nHepp, nhp, nHeII, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
    rho_cgs  = rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; // used by several of the models below
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    if(SphP[i].DelayTimeHII>0) {nh0=DMIN(1.e-4,nh0); nHe0=DMIN(1.e-5,nHe0);} // account for our effective ionization model here
#endif
    f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL < 0) /* disable CR diffusion, specifically */
    SphP[i].CosmicRayDiffusionCoeff = 0; // no diffusion (but -can- allow streaming
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 0) /* set diffusivity to a constant  */
    SphP[i].CosmicRayDiffusionCoeff = All.CosmicRayDiffusionCoeff; //  this is the input value of the diffusivity, for constant-kappa models
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 3) /* Farber et al. 2018 -- higher coeff in neutral gas, lower in ionized gas */
    SphP[i].CosmicRayDiffusionCoeff = (3.e29/unit_kappa_code) * (1.-f_ion + f_ion/30.); // 30x lower in neutral (note use f_ion directly here, not temperature as they do)
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 4) /* Wiener et al. 2017 style pure-streaming but with larger streaming speeds and limited losses */
    double ni_m3=f_ion*(rho_cgs/PROTONMASS)/1.e-3, T6=temperature/1.e6, Lturbkpc=p_scale, Lgradkpc=CRPressureGradScaleLength*(All.UnitLength_in_cm/All.HubbleParam)/3.086e21, h0_fac=Get_Particle_Size(i)*All.cf_atime*All.cf_a2inv*All.UnitVelocity_in_cm_per_s/1.e6;
    double dv2_10=0; for(k=0;k<3;k++) {int j; for(j=0;j<3;j++) {dv2_10 += SphP[i].Gradients.Velocity[j][k]*SphP[i].Gradients.Velocity[j][k]*h0_fac*h0_fac;}}
    double ecr_14 = SphP[i].CosmicRayEnergyPred * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * ((All.UnitEnergy_in_cgs/All.HubbleParam)/pow(All.UnitLength_in_cm/All.HubbleParam,3)) / 1.0e-14; // CR energy density in CGS units //
    v_streaming = Get_CosmicRayStreamingVelocity(i) + (1.e5/All.UnitVelocity_in_cm_per_s)*(4.1*pow(MIN_REAL_NUMBER+ni_m3*T6,0.25)/pow(MIN_REAL_NUMBER+ecr_14*Lgradkpc,0.5) + 1.2*pow(MIN_REAL_NUMBER+dv2_10*ni_m3,0.75)/(MIN_REAL_NUMBER+ecr_14*sqrt(Lturbkpc)));
    CR_kappa_streaming = GAMMA_COSMICRAY * v_streaming * CRPressureGradScaleLength; // convert to effective diffusivity
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 5) /* streaming at fast MHD wavespeed */
    v_streaming = sqrt(v_streaming*v_streaming + GAMMA*GAMMA_MINUS1*SphP[i].InternalEnergyPred);
    CR_kappa_streaming = GAMMA_COSMICRAY * v_streaming * CRPressureGradScaleLength;
#endif
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 1) || (COSMIC_RAYS_DIFFUSION_MODEL == 2) || (COSMIC_RAYS_DIFFUSION_MODEL == 6) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* extrinsic turbulence OR streaming/diffusion at self-confinement equilibrium */
    double EPSILON_SMALL=1.e-50, cs_thermal=sqrt(GAMMA*GAMMA_MINUS1*u0), B[3]={0}, Bmag=0, cos_Bgrad=0, dPmag=0, Bmag_Gauss=0, r_turb_driving=0; // internal energy,  thermal sound speed, B-field properties
    for(k=0;k<3;k++)
    {
        B[k]=SphP[i].BPred[k]*SphP[i].Density/P[i].Mass*All.cf_a2inv; Bmag+=B[k]*B[k]; // B magnitude in code units (physical)
        cos_Bgrad+=B[k]*SphP[i].Gradients.CosmicRayPressure[k]; // dot product of B_hat and gradient_P_hat
        dPmag+=SphP[i].Gradients.CosmicRayPressure[k]*SphP[i].Gradients.CosmicRayPressure[k]; // magnitude of CR pressure gradient
        r_turb_driving += SphP[i].Gradients.Pressure[k]*SphP[i].Gradients.Pressure[k]; // pressure-scale length defines initial guess for turb driving scale
    }
    Bmag=sqrt(Bmag); dPmag=sqrt(dPmag); cos_Bgrad/=(Bmag*dPmag + EPSILON_SMALL); Bmag_Gauss=Bmag*gizmo2gauss; // compute proper magnitudes
    r_turb_driving = DMAX(SphP[i].Pressure/(EPSILON_SMALL+sqrt(r_turb_driving)),Get_Particle_Size(i))*All.cf_atime; double k_turb=1./r_turb_driving; // maximum of gradient scale length or resolution scale
    double vA_code=sqrt(Bmag*Bmag/(SphP[i].Density*All.cf_a3inv)), vA_noion=vA_code; // Alfven speed^2 in code units [recall B units such that there is no 4pi here]
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_code /= sqrt(f_ion); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    double clight_code=C/All.UnitVelocity_in_cm_per_s, Omega_gyro=8987.34*(Bmag_Gauss/E_CRs_GeV_over_Z)*(All.UnitTime_in_s/All.HubbleParam), r_L=clight_code/Omega_gyro, kappa_0=r_L*clight_code, k_L=1./r_L; // gyro frequency of the CR population we're evolving, CR gyro radius, reference diffusivity
    double E_B=0.5*Bmag*Bmag*(P[i].Mass/(SphP[i].Density*All.cf_a3inv)), E_CR=SphP[i].CosmicRayEnergyPred, x_EB_ECR=(E_B+EPSILON_SMALL)/(E_CR+EPSILON_SMALL); // B-field energy (energy density times volume, for ratios with energies above), CR-energy, ratio
    
    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,x_LL,M_A,h0,fturb_multiplier=1; // factor which will represent which cascade model we are going to use
    for(i1=0;i1<3;i1++)
    {
        v2_t += SphP[i].VelPred[i1]*SphP[i].VelPred[i1]; b2_t += Get_Particle_BField(i,i1) * Get_Particle_BField(i,i1);
        for(i2=0;i2<3;i2++) {dv2_t += SphP[i].Gradients.Velocity[i1][i2]*SphP[i].Gradients.Velocity[i1][i2]; db2_t += SphP[i].Gradients.B[i1][i2]*SphP[i].Gradients.B[i1][i2];}
    }
    v2_t=sqrt(v2_t); b2_t=sqrt(b2_t); dv2_t=sqrt(dv2_t); db2_t=sqrt(db2_t); dv2_t/=All.cf_atime; db2_t/=All.cf_atime; b2_t*=All.cf_a2inv; db2_t*=All.cf_a2inv; v2_t/=All.cf_atime; dv2_t/=All.cf_atime; h0=Get_Particle_Size(i)*All.cf_atime; // physical units

    int use_shear_corrected_vturb = 1;
    if(use_shear_corrected_vturb==1) 
    {
        double dv2_t = sqrt((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) *
            (SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) + (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) *
            (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) + (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) * (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
            (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] + SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) - (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] + SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv;
        double db2_t = sqrt((1./2.)*((SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) * (SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) +
            (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) * (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) +
            (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2]) * (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2])) +
            (2./3.)*((SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[0][0] + SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[1][1] +
            SphP[i].Gradients.B[2][2]*SphP[i].Gradients.B[2][2]) - (SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[2][2] +
            SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[1][1] + SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[2][2]))) * All.cf_a3inv;

        double db_v_equiv = h0*db2_t * vA_noion / (EPSILON_SMALL + b2_t);
        double dv_e = sqrt(h0*dv2_t*h0*dv2_t + db_v_equiv*db_v_equiv + EPSILON_SMALL);
        double vA_eff = sqrt(vA_noion*vA_noion + db_v_equiv*db_v_equiv + EPSILON_SMALL);
        M_A = (EPSILON_SMALL + dv_e) / (EPSILON_SMALL + vA_eff);
    } else {
        M_A = h0*(EPSILON_SMALL + dv2_t) / (EPSILON_SMALL + vA_noion); M_A = DMAX(M_A , h0*(EPSILON_SMALL + db2_t) / (EPSILON_SMALL + b2_t)); 
    }
    M_A = DMAX( EPSILON_SMALL , M_A ); // proper calculation of the local Alfven Mach number
    x_LL = clight_code / (Omega_gyro * h0); x_LL=DMAX(x_LL,EPSILON_SMALL);
    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 6) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* terms below only needed if we try to calculate the self-confinement-based diffusivity */
    k_turb = 1./h0; // scale at which turbulence is being measured here //
    fturb_multiplier = pow(M_A,3./2.); // corrects to Alfven scale, for correct estimate according to Farmer and Goldreich, Lazarian, etc.
    if(M_A<1.) {fturb_multiplier*=DMIN(sqrt(M_A),pow(M_A,7./6.)/pow(x_LL,1./6.));} else {fturb_multiplier*=DMIN(1.,1./(pow(M_A,1./2.)*pow(x_LL,1./6.)));} /* Lazarian+ 2016 multi-part model for where the resolved scales lie on the cascade */
    
    //fturb_multiplier = pow(M_A,3./2.) / pow(x_LL,1./10.); // GS anisotropic but perp cascade is IK
    //fturb_multiplier = pow(M_A,3./2.) * 1./(pow(M_A,1./2.)*pow(x_LL,1./6.)); // pure-Kolmogorov 
    //fturb_multiplier = pow(M_A,3./2.) * 100.; // arbitrary multiplier
    //fturb_multiplier = pow(M_A,3./2.) * 10.; // arbitrary multiplier

    /* ok now we finally have all the terms needed to calculate the various damping rates that determine the equilibrium diffusivity */
    double G_ion_neutral = 5.77e-11 * (rho_cgs/PROTONMASS) * (0.97*nh0 + 0.03*nHe0) * sqrt(temperature) * (All.UnitTime_in_s/All.HubbleParam); if(Z_charge_CR > 1) {G_ion_neutral /= sqrt(2.*Z_charge_CR);} // ion-neutral damping: need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]
    double G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI)*cs_thermal/4.) * sqrt(k_turb*k_L) * fturb_multiplier,  G0 = G_ion_neutral + G_turb_plus_linear_landau, Gamma_effective = G0; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)
    double phi_0 = (sqrt(M_PI)/6.)*(fabs(cos_Bgrad))*(1./(x_EB_ECR+EPSILON_SMALL))*(cs_thermal*vA_code*k_L/(CRPressureGradScaleLength*G0*G0 + EPSILON_SMALL)); // parameter which determines whether NLL dominates
    if(isfinite(phi_0) && (phi_0>0.01)) {Gamma_effective *= phi_0/(2.*(sqrt(1.+phi_0)-1.));} // this accounts exactly for the steady-state solution for the Thomas+Pfrommer formulation, including both the linear [Landau,turbulent,ion-neutral] + non-linear terms. can estimate (G_nonlinear_landau_effective = Gamma_effective - G0)
    
    /* with damping rates above, equilibrium transport is equivalent to pure streaming, with v_stream = vA + (diffusive equilibrium part) give by the solution below, proportional to Gamma_effective and valid to O(v^2/c^2) */
    double v_st_eff = vA_code * (1. + 4. * kappa_0 * Gamma_effective * x_EB_ECR * (1. + 2.*vA_code*vA_code/(clight_code*clight_code)) / (M_PI*vA_code*vA_code + EPSILON_SMALL)); // effective equilibrium streaming speed for all terms accounted
    CR_kappa_streaming = GAMMA_COSMICRAY*v_st_eff*CRPressureGradScaleLength; // convert to effective diffusivity from 'streaming'
    CR_kappa_streaming = DMAX(kappa_0,DMAX(1.e20/unit_kappa_code,DMIN(1.e35/unit_kappa_code,CR_kappa_streaming))); // and limit (can get extreme) to prevent numerical overflow errors
#endif
#endif

    
#if (COSMIC_RAYS_DIFFUSION_MODEL == 1) || (COSMIC_RAYS_DIFFUSION_MODEL == 2) || (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* textbook extrinsic turbulence model: kappa~v_CR*r_gyro * B_bulk^2/(B_random[scale~r_gyro]^2) v_CR~c, r_gyro~p*c/(Z*e*B)~1e12 cm * RGV *(3 muG/B)  (RGV~1 is the magnetic rigidity). assuming a Kolmogorov spectrum */
    // double l_alfven = h0/(EPSILON_SMALL + M_A*M_A*M_A), l_alfven_kpc = l_alfven*codelength_to_kpc; /* L_alfven defined as scale for Kolmogorov cascade to extrapolate to dv ~ v_Alfven; second term just a conversion to convenient units below */
    // double f_cas_ET = 1; /* pure Goldreich-Shridhar cascade, ignoring anisotropic effects per Chandran-00 */
    // if(M_A < 1.) {f_cas_ET = pow(1./DMAX(M_A*M_A , x_LL), 1./3.);} /* Lazarian '16 modification for weak cascade in sub-Alfvenic turbulence */
    // SphP[i].CosmicRayDiffusionCoeff = (9.e28/unit_kappa_code) * pow(E_CRs_GeV_over_Z*l_alfven_kpc*l_alfven_kpc/b_muG,1./3.) * f_cas_ET; /* kappa~9e28 * (l_alfven/kpc)^(2/3) * RGV^(1/3) * (B/muG)^(-1/3) * f_cas_ET,  follows Jokipii 1966, with our corrections for spectral shape */
    double f_cas_ET = 0.007 * clight_code / vA_noion; /* damping in Alfvenic turbulence following an anisotropic Goldreich-Shridar cascade, per Chandran 2000 */

#if (COSMIC_RAYS_DIFFUSION_MODEL == 2) || (COSMIC_RAYS_DIFFUSION_MODEL == 7)
    /* Snodin et al. 2016 -- different expression for extrinsic MHD-turb diffusivity, using the proper definition of L_Alfven for scaling to the correct limit */
    // SphP[i].CosmicRayDiffusionCoeff = (3.e29/unit_kappa_code) * (l_alfven_kpc + 0.1*pow(E_CRs_GeV_over_Z*l_alfven_kpc*l_alfven_kpc/b_muG,1./3.) + 2.4e-7*E_CRs_GeV_over_Z/b_muG);
    double n1=rho_cgs/PROTONMASS, h0_kpc=h0*codelength_to_kpc, T4=temperature/1.e4, gL=E_CRs_GeV_over_Z, fcasET_colless = 0.04*cs_thermal/vA_noion; /* collisionless [Landau] damping of fast modes */
    double fcasET_viscBrg = 0.03*pow(EPSILON_SMALL + M_A,4./3.)*T4/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*gL*T4,1./6.); /* Spitzer/Braginski viscous damping of fast modes */
    double fcasET_viscMol = 0.41*pow(EPSILON_SMALL + M_A,4./3.)*nh0/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*gL/(EPSILON_SMALL + T4),1./6.); /* atomic/molecular collisional damping of fast modes */
    double f_cas_ET_fast = fcasET_colless + fcasET_viscBrg + fcasET_viscMol; /* fast modes, accounting for damping, following Yan+Lazarian 2005 */
    f_cas_ET = 1./(EPSILON_SMALL + 1./(EPSILON_SMALL+f_cas_ET) + 1./(EPSILON_SMALL+f_cas_ET_fast)); /* combine fast-mode and Alfvenic scattering */
#endif

    SphP[i].CosmicRayDiffusionCoeff = (1.e32/unit_kappa_code) * h0_kpc / (EPSILON_SMALL + M_A*M_A) * f_cas_ET;
#endif
    

#if (COSMIC_RAYS_DIFFUSION_MODEL == 7) /* 'combined' extrinsic turbulence + self-confinement model */
    double kappa_diff_extrinsicturb = SphP[i].CosmicRayDiffusionCoeff; // expression which should be calculated above for turbulent part //
    SphP[i].CosmicRayDiffusionCoeff = 0; // re-zero because we will calculate the more appropriate rate below
    CR_kappa_streaming = 1. / (EPSILON_SMALL +  1./(CR_kappa_streaming+EPSILON_SMALL) + 1./(kappa_diff_extrinsicturb+EPSILON_SMALL) ); // if scattering rates add linearly, this is a rough approximation to the total transport (essentially, smaller of the two dominates)

    double kappa_max=1.e34/unit_kappa_code, kappa_min=1.e26/unit_kappa_code; CR_kappa_streaming=DMIN(DMAX(CR_kappa_streaming,kappa_min),kappa_max);    
#endif
    
    SphP[i].CosmicRayDiffusionCoeff += CR_kappa_streaming; //  add effective streaming coefficient
    
#ifndef COSMIC_RAYS_M1 /* now we apply a limiter to prevent the coefficient from becoming too large: cosmic rays cannot stream/diffuse with v_diff > c */
    // [all of this only applies if we are using the pure-diffusion description: the M1-type description should -not- use a limiter here, or negative kappa]
    double diffusion_velocity_limit = 1.0 * C; /* maximum diffusion velocity (set <C if desired) */
    double Lscale = DMIN(20.*Get_Particle_Size(i)*All.cf_atime , CRPressureGradScaleLength);
    double kappa_diff_vel = SphP[i].CosmicRayDiffusionCoeff * GAMMA_COSMICRAY_MINUS1 / Lscale * All.UnitVelocity_in_cm_per_s;
    SphP[i].CosmicRayDiffusionCoeff *= 1 / (1 + kappa_diff_vel/diffusion_velocity_limit); /* caps maximum here */
#ifdef GALSF /* for multi-physics problems, we suppress diffusion where it is irrelevant */
    SphP[i].CosmicRayDiffusionCoeff *= 1 / (1 + kappa_diff_vel/(0.01 * C)); /* caps maximum here */
    double P_cr_Ratio = Get_Particle_CosmicRayPressure(i) / (MIN_REAL_NUMBER + SphP[i].Pressure);
    double P_min = 1.0e-4; if(P_cr_Ratio < P_min) {SphP[i].CosmicRayDiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
    P_min = 1.0e-6; if(P_cr_Ratio < P_min) {SphP[i].CosmicRayDiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
#endif
    SphP[i].CosmicRayDiffusionCoeff /= GAMMA_COSMICRAY_MINUS1; // ensure correct units for subsequent operations //
#endif
    
    if((SphP[i].CosmicRayDiffusionCoeff<=0)||(isnan(SphP[i].CosmicRayDiffusionCoeff))) {SphP[i].CosmicRayDiffusionCoeff=0;} /* nan check! */
}



double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode)
{
    /* routine to do the drift/kick operations for CRs: mode=0 is kick, mode=1 is drift */
    if(dt_entr <= 0) {return 0;} // no update
    int k; double eCR, u0; k=0;
    if(mode==0) {eCR=SphP[i].CosmicRayEnergy; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred; u0=SphP[i].InternalEnergyPred;} // initial energy
    if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforced throughout code
    if(eCR < 0) {eCR=0;} // limit to physical values
    
#if defined(COSMIC_RAYS_M1) && !defined(COSMIC_RAYS_ALFVEN)
    // this is the exact solution for the CR flux-update equation over a finite timestep dt: it needs to be solved this way [implicitly] as opposed to explicitly for dt because in the limit of dt_cr_dimless being large, the problem exactly approaches the diffusive solution
    double DtCosmicRayFlux[3]={0}, flux[3]={0}, CR_veff[3]={0}, CR_vmag=0, q_cr = 0, cr_speed = COSMIC_RAYS_M1;// * (C/All.UnitVelocity_in_cm_per_s);
    cr_speed = DMAX( All.cf_afac3*SphP[i].MaxSignalVel , DMIN(COSMIC_RAYS_M1 , 10.*fabs(SphP[i].CosmicRayDiffusionCoeff)/(Get_Particle_Size(i)*All.cf_atime)));// * (C/All.UnitVelocity_in_cm_per_s);
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
#ifdef MAGNETIC // do projection onto field lines
    double fluxmag=0, fluxdot=0, Bmag=0; for(k=0;k<3;k++) {fluxmag+=flux[k]*flux[k]; fluxdot+=flux[k]*B0[k];}
    if(fluxmag>0) {fluxmag=sqrt(fluxmag);} else {fluxmag=0;}
    if(fluxdot<0) {fluxmag*=-1;} // points down-field
    if(Bmag2>0) {for(k=0;k<3;k++) {flux[k] = fluxmag * B0[k] / sqrt(Bmag2);}} // re-assign to be along field
#endif
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
    
    
    // now update all scalar fields (CR energies and Alfvenic energies, if those are followed) from fluxes and adiabatic terms //
    int q_whichupdate, q_N_updates = 1;
#ifdef COSMIC_RAYS_ALFVEN
    q_N_updates = 3; // update the Alfvenic energy terms from (0) advection with gas [already solved], (1) their fluxes, and (2) their adiabatic terms. this should be basically identical to the CR term.
#endif
    for(q_whichupdate=0; q_whichupdate<q_N_updates; q_whichupdate++)
    {
        // first update the CR energies from fluxes. since this is positive-definite, some additional care is needed //
        double dCR_dt = SphP[i].DtCosmicRayEnergy, gamma_eff = GAMMA_COSMICRAY, eCR_tmp = eCR;
#ifdef COSMIC_RAYS_ALFVEN
        if(q_whichupdate>0) {dCR_dt=SphP[i].DtCosmicRayAlfvenEnergy[q_whichupdate-1]; gamma_eff=GAMMA_ALFVEN_CRS; if(mode==0) {eCR_tmp=SphP[i].CosmicRayAlfvenEnergy[q_whichupdate-1];} else {eCR_tmp=SphP[i].CosmicRayAlfvenEnergyPred[q_whichupdate-1];}}
#endif
        double dCR = dCR_dt*dt_entr, dCRmax = 1.e10*(eCR_tmp+MIN_REAL_NUMBER);
#if defined(GALSF) && !defined(COSMIC_RAYS_ALFVEN)
        dCRmax = DMAX(2.0*eCR_tmp , 0.1*u0*P[i].Mass);
#endif
        if(dCR > dCRmax) {dCR=dCRmax;} // don't allow excessively large values
        if(dCR < -eCR_tmp) {dCR=-eCR_tmp;} // don't allow it to go negative
	    double eCR_00 = eCR_tmp;
        eCR_tmp += dCR; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        if(q_whichupdate==0) {if(mode==0) {SphP[i].CosmicRayEnergy=eCR_tmp;} else {SphP[i].CosmicRayEnergyPred=eCR_tmp;}} // updated energy
#ifdef COSMIC_RAYS_ALFVEN
        if(q_whichupdate>0) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[q_whichupdate-1]=eCR_tmp;} else {SphP[i].CosmicRayAlfvenEnergyPred[q_whichupdate-1]=eCR_tmp;}} // updated energy
#endif
        double eCR_0 = eCR_tmp; // save this value for below
        // now need to account for the adiabatic heating/cooling of the 'fluid', here, with gamma=gamma_eff //
        double d_div = (-(gamma_eff-1.) * P[i].Particle_DivVel*All.cf_a2inv) * dt_entr;
        if(All.ComovingIntegrationOn) {d_div += (-3.*(gamma_eff-1.) * All.cf_hubble_a) * dt_entr;} /* adiabatic term from Hubble expansion (needed for cosmological integrations */
        double dCR_div = DMIN(eCR_tmp*d_div , 0.5*u0*P[i].Mass); // limit so don't take away all the gas internal energy [to negative values]
        if(dCR_div + eCR_tmp < 0) {dCR_div = -eCR_tmp;} // check against energy going negative
        eCR_tmp += dCR_div; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        dCR_div = eCR_tmp - eCR_0; // actual change that is going to be applied
        if(dCR_div < -0.5*P[i].Mass*u0) {dCR_div=-0.5*P[i].Mass*u0;} // before re-coupling, ensure this will not cause negative energies
        if(dCR_div < -0.9*eCR_00) {dCR_div=-0.9*eCR_00;} // before re-coupling, ensure this will not cause negative energies 
        if(q_whichupdate==0) {if(mode==0) {SphP[i].CosmicRayEnergy += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;} else {SphP[i].CosmicRayEnergyPred += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;}}
#ifdef COSMIC_RAYS_ALFVEN
        if(q_whichupdate>0) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[q_whichupdate-1] += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;} else {SphP[i].CosmicRayAlfvenEnergyPred[q_whichupdate-1] += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;}}
#endif
    }
    
    
#ifdef COSMIC_RAYS_ALFVEN
    double Z_charge_CR = 1, E_CRs_Gev = 1; // charge and energy and resonant Alfven wavenumber (in gyro units) of the CR population we're evolving

    // ok, the updates from [0] advection w gas, [1] fluxes, [2] adiabatic, [-] catastrophic (in cooling.c) are all set, just need exchange terms b/t CR and Alfven //
    double EPSILON_SMALL = 1.e-77; // want a very small number here 
    double bhat[3], Bmag=0, Bmag_Gauss, clight_code=C/All.UnitVelocity_in_cm_per_s, Omega_gyro, eA[2], vA_code, vA2_c2, E_B, fac, flux_G, fac_Omega, flux[3], f_CR, f_CR_dot_B, cs_thermal, r_turb_driving, G_ion_neutral=0, G_turb_plus_linear_landau=0, G_nonlinear_landau_prefix=0;
    double ne=0, f_ion=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs=rho*All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
#ifdef COOLING
    ne=SphP[i].Ne, temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
#ifdef GALSF_FB_FIRE_RT_HIIHEATING
    if(SphP[i].DelayTimeHII>0) {nh0=0;} // account for our effective ionization model here
#endif
    f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
#endif
    for(k=0;k<3;k++) {if(mode==0) {bhat[k]=SphP[i].B[k];} else {bhat[k]=SphP[i].BPred[k];}} // grab whichever B field we need for our mode
    if(mode==0) {eCR=SphP[i].CosmicRayEnergy; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred; u0=SphP[i].InternalEnergyPred;} // initial energy
    if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforce the usual minimum thermal energy the code requires
    for(k=0;k<2;k++) {if(mode==0) {eA[k]=SphP[i].CosmicRayAlfvenEnergy[k];} else {eA[k]=SphP[i].CosmicRayAlfvenEnergyPred[k];}} // Alfven energy
    if(mode==0) {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFlux[k];}} else {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFluxPred[k];}} // load flux
    f_CR=0; f_CR_dot_B=0; for(k=0;k<3;k++) {f_CR+=flux[k]*flux[k]; f_CR_dot_B+=bhat[k]*flux[k];} // compute the magnitude of the flux density
    f_CR=sqrt(f_CR); if(f_CR_dot_B<0) {f_CR*=-1;} // initialize the flux density variable from the previous timestep, appropriately signed with respect to the b-field
    for(k=0;k<3;k++) {Bmag+=bhat[k]*bhat[k];} // compute magnitude
    Bmag = sqrt(Bmag); for(k=0;k<3;k++) {bhat[k]/=(EPSILON_SMALL+Bmag);} // now it's bhat we have here
    Bmag *= SphP[i].Density/P[i].Mass * All.cf_a2inv; // convert to actual B in physical units
    E_B = 0.5*Bmag*Bmag * (P[i].Mass/(SphP[i].Density*All.cf_a3inv)); // B-field energy (energy density times volume, for ratios with energies above)
    double Eth_0 = EPSILON_SMALL + 1.e-8 * P[i].Mass*u0; // set minimum magnetic energy relative to thermal (maximum plasma beta ~ 1e8) to prevent nasty divergences
    if(E_B < Eth_0) {Bmag = sqrt(2.*Eth_0/((P[i].Mass/(SphP[i].Density*All.cf_a3inv))));} // enforce this maximum beta for purposes of "B" to insert below 
    E_B = 0.5*Bmag*Bmag * (P[i].Mass/(SphP[i].Density*All.cf_a3inv)); // B-field energy (energy density times volume, for ratios with energies above)
    Bmag_Gauss = Bmag * sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam); // turn it into Gauss
    Omega_gyro = 8987.34 * Bmag_Gauss * (Z_charge_CR/E_CRs_Gev) * (All.UnitTime_in_s/All.HubbleParam); // gyro frequency of the CR population we're evolving
    vA_code = sqrt( Bmag*Bmag / (SphP[i].Density*All.cf_a3inv) ); double vA_noion=vA_code; // Alfven speed^2 in code units [recall B units such that there is no 4pi here]
#ifdef COSMIC_RAYS_ION_ALFVEN_SPEED
    vA_code /= sqrt(f_ion); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    cs_thermal = sqrt(GAMMA*(GAMMA-1.) * u0); // thermal sound speed at appropriate drift-time [in code units, physical]
    vA2_c2 = vA_code*vA_code / (clight_code*clight_code); // Alfven speed vs c_light
    fac_Omega = (3.*M_PI/16.) * Omega_gyro * (1.+2.*vA2_c2); // factor which will be used heavily below
    /* for turbulent (anisotropic and linear landau) damping terms: need to know the turbulent driving scale: assume a cascade with a driving length equal to the pressure gradient scale length */
    r_turb_driving = 0; for(k=0;k<3;k++) {r_turb_driving += SphP[i].Gradients.Pressure[k]*SphP[i].Gradients.Pressure[k];} // compute gradient magnitude
    r_turb_driving = DMAX( SphP[i].Pressure / (EPSILON_SMALL + sqrt(r_turb_driving)) , Get_Particle_Size(i) ) * All.cf_atime; // maximum of gradient scale length or resolution scale
    double k_turb = 1./r_turb_driving, k_L = Omega_gyro / clight_code;
    
    // before acting on the 'stiff' sub-system, account for the 'extra' advection term that accounts for 'twisting' of B:
    fac=0; for(k=0;k<3;k++) {fac += bhat[k] * (bhat[0]*SphP[i].Gradients.Velocity[k][0] + bhat[1]*SphP[i].Gradients.Velocity[k][1] + bhat[2]*SphP[i].Gradients.Velocity[k][2]);}
    if(All.ComovingIntegrationOn) {fac += All.cf_hubble_a;} // adds cosmological/hubble flow term here [not included in peculiar velocity gradient]
    fac *= -All.cf_a2inv*dt_entr; if(!isfinite(fac)) {fac=0;} else {if(fac>2.) {fac=2.;} else {if(fac<-2.) {fac=-2.;}}} // limit factor for change here, should be small given Courant factor
    f_CR *= exp(fac); // update flux term accordingly, before next step //

    // because the equations below will very much try to take things to far-too-small values for numerical precision, we need to define a bunch of sensible bounds for values to allow, to prevent divergences, but also enforce conservation
    // calculate minimum eA,eCR to enforce; needed because if eA is identically zero, nothing can get amplified, and it will always be zero. but for large enough seed to amplify, results should not depend on seed //
    eA[0]=DMAX(eA[0],0); eA[1]=DMAX(eA[1],0); eCR=DMAX(eCR,0); // enforce non-negative energies 
    double Min_Egy=0, e_tot=0, e_tot_new=0, fmax=0; e_tot = eCR + eA[0] + eA[1] + EPSILON_SMALL; // sum total energy, enforce positive-definite: will use this to ensure total energy conservation when enforcing minima below
    { 
        double h=Get_Particle_Size(i)*All.cf_atime; int k2; for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {Min_Egy+=SphP[i].Gradients.B[k][k2]*SphP[i].Gradients.B[k][k2];}}
        Min_Egy=h*sqrt(Min_Egy/9.)*All.cf_a2inv; Min_Egy=DMIN(Min_Egy,Bmag); r_turb_driving=DMAX(h,r_turb_driving); Min_Egy=DMIN(Min_Egy,Bmag*pow(h/r_turb_driving,1./3.)); Min_Egy=Min_Egy*pow(DMIN(clight_code/Omega_gyro,DMIN(h,r_turb_driving))/h,1./3.); // Min_Egy is now magnetic field extrap to r_gyro
        Min_Egy = 0.5 * (Min_Egy*Min_Egy) * P[i].Mass/(SphP[i].Density*All.cf_a3inv); // magnetic energy at this scale, from the above //
        double epsilon = 1.e-15; Min_Egy *= epsilon; // minimum energy is a tiny fraction of B at the dissipation scale
        if(Min_Egy <= 0 || !isfinite(Min_Egy)) {Min_Egy = 1.e-15*eCR;} // if this minimum-energy calculation failed, enforce a tiny fraction of the CR energy
        if(Min_Egy <= 0 || !isfinite(Min_Egy)) {Min_Egy = 1.e-15*P[i].Mass*u0;} // if this minimum-energy calculation failed, enforce a tiny fraction of the thermal energy
        if(Min_Egy <= 0) {Min_Egy = EPSILON_SMALL;} // if this still failed, simply enforce a tiny positive-definite value
    }
    eCR=DMAX(eCR,Min_Egy); eA[0]=DMAX(eA[0],Min_Egy); eA[1]=DMAX(eA[1],Min_Egy); // enforce

    // ok, now all the advection and adiabatic operations should be complete. they are split above. 
    //  what remains is the stiff, coupled subsystem of wave growth+damping, which needs to be treated 
    //  more carefully or else we get very large over/under-shoots

    // first define some convenient units and dimensionless quantities, and enforce limits on values of input quantities
    double eCR_0 = 1.e-6*(E_B + P[i].Mass*u0) + eCR + eA[0] + eA[1] + fabs(f_CR/COSMIC_RAYS_ALFVEN); // this can be anything, just need a normalization for the characteristic energy scale of the problem //
    double ceff2_va2=(COSMIC_RAYS_ALFVEN*COSMIC_RAYS_ALFVEN)/(vA_code*vA_code), t0=1./(fac_Omega*(eCR_0/E_B)*vA2_c2), gammCR=GAMMA_COSMICRAY, f_unit=vA_code*eCR_0, volume=P[i].Mass/(SphP[i].Density*All.cf_a3inv); // factors used below , and for units
    double x_e=eCR/eCR_0, x_f=f_CR/f_unit, x_up=eA[0]/eCR_0, x_um=eA[1]/eCR_0, dtau=dt_entr/t0; e_tot/=eCR_0; Min_Egy/=eCR_0; // initial values in relevant units
    Min_Egy=DMAX(DMIN(Min_Egy,DMIN(x_e,DMIN(x_up,x_um))),EPSILON_SMALL); if(!isfinite(Min_Egy)) {Min_Egy=EPSILON_SMALL;} // enforce positive-definite-ness
    // we can more robustly define a minimum and maximum e_A by reference to a minimum and maximum 'effective diffusivity' over which it is physically meaningful, and numerically possible to evolve them
    double ref_diffusivity = 4.4e26 / (All.UnitVelocity_in_cm_per_s * All.UnitLength_in_cm / All.HubbleParam); // define a unit diffusivity in code units for reference below
    double xkappa_min = DMAX(vA_code*vA_code*t0/(3.e8*ref_diffusivity) , EPSILON_SMALL); // maximum diffusivity ~1e35, but be non-zero
    double xkappa_max = DMAX(DMIN(vA_code*vA_code*t0/(3.e-8*ref_diffusivity) , 0.5*E_B/eCR_0), xkappa_min); // minimum diffusivity at ~1e19, but cannot have more energy in eAp+eAm than total magnetic energy! (equations below assume -small- fraction of E_B in eA!, or growth rates non-linearly modified)
    if(e_tot < Min_Egy || !isfinite(e_tot)) {e_tot = Min_Egy;} // enforce minima/maxima
    if(x_e   < Min_Egy || !isfinite(x_e)  ) {x_e   = Min_Egy;} // enforce minima/maxima
    if(x_um<EPSILON_SMALL || !isfinite(x_um)) {x_um=EPSILON_SMALL;} else {if(x_um>xkappa_max) {x_um=xkappa_max;}} // enforce minima/maxima
    if(x_up<EPSILON_SMALL || !isfinite(x_up)) {x_up=EPSILON_SMALL;} else {if(x_up>xkappa_max) {x_up=xkappa_max;}} // enforce minima/maxima
    if(x_um+x_up<xkappa_min) {fac=xkappa_min/(x_um+x_up); x_um*=fac; x_up*=fac;} // only want to enforce -sum- having effective diffusivity, not both
    e_tot_new=x_e+x_um+x_up; x_e*=e_tot/e_tot_new; x_up*=e_tot/e_tot_new; x_um*=e_tot/e_tot_new; // check energy after limit-enforcement
    fmax = x_e*sqrt(ceff2_va2); if(!isfinite(x_f)) {x_f=0;} else {if(x_f>fmax) {x_f=fmax;} else {if(x_f<-fmax) {x_f=-fmax;}}} // check for flux maximum/minimum

    // calculate the dimensionless flux source term for the stiff part of the equations
    flux_G=0; for(k=0;k<3;k++) {flux_G += bhat[k] * SphP[i].Gradients.CosmicRayPressure[k];} // b.gradient[P] -- flux source term
    double psifac = flux_G * (vA_code*t0) / (eCR/volume); // this gives the strength of the gradient source term, should remain fixed over stiff part of loop

    // calculate the wave-damping rates (again in appropriate dimensionless units)
    /* ion-neutral damping: need thermodynamic information (neutral fractions, etc) to compute self-consistently */
    G_ion_neutral = 5.77e-11 * (rho_cgs/PROTONMASS) * nh0 * sqrt(temperature) * (All.UnitTime_in_s/All.HubbleParam); // need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]
    if(Z_charge_CR > 1) {G_ion_neutral /= sqrt(2.*Z_charge_CR);}

    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,x_LL,M_A,h0,fturb_multiplier=1; // factor which will represent which cascade model we are going to use
    for(i1=0;i1<3;i1++)
    {
        v2_t += SphP[i].VelPred[i1]*SphP[i].VelPred[i1]; b2_t += Get_Particle_BField(i,i1) * Get_Particle_BField(i,i1);
        for(i2=0;i2<3;i2++) {dv2_t += SphP[i].Gradients.Velocity[i1][i2]*SphP[i].Gradients.Velocity[i1][i2]; db2_t += SphP[i].Gradients.B[i1][i2]*SphP[i].Gradients.B[i1][i2];}
    }
    v2_t=sqrt(v2_t); b2_t=sqrt(b2_t); dv2_t=sqrt(dv2_t); db2_t=sqrt(db2_t); dv2_t/=All.cf_atime; db2_t/=All.cf_atime; b2_t*=All.cf_a2inv; db2_t*=All.cf_a2inv; v2_t/=All.cf_atime; dv2_t/=All.cf_atime; h0=Get_Particle_Size(i)*All.cf_atime; // physical units
    M_A = h0*(EPSILON_SMALL + dv2_t) / (EPSILON_SMALL + vA_noion); M_A = DMAX(M_A , h0*(EPSILON_SMALL + db2_t) / (EPSILON_SMALL + b2_t)); M_A = DMAX( EPSILON_SMALL , M_A ); // proper calculation of the local Alfven Mach number
    x_LL = clight_code / (Omega_gyro * h0); x_LL=DMAX(x_LL,EPSILON_SMALL); k_turb = 1./h0; // scale at which turbulence is being measured here //
    fturb_multiplier = pow(M_A,3./2.); // corrects to Alfven scale, for correct estimate according to Farmer and Goldreich, Lazarian, etc.
    if(M_A<1.) {fturb_multiplier*=DMIN(sqrt(M_A),pow(M_A,7./6.)/pow(x_LL,1./6.));} else {fturb_multiplier*=DMIN(1.,1./(pow(M_A,1./2.)*pow(x_LL,1./6.)));} /* Lazarian+ 2016 multi-part model for where the resolved scales lie on the cascade */
    G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI/16.)*cs_thermal) * sqrt(k_turb*k_L) * fturb_multiplier; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)

    G_nonlinear_landau_prefix = (sqrt(M_PI)/8.) * (1./E_B) * (cs_thermal*k_L); // non-linear Landau damping (will be multiplied by eA)
    double gamma_in_t_ll = (G_ion_neutral + G_turb_plus_linear_landau) * t0; // dimensionless now and appropriate code units
    double gamma_nll = G_nonlinear_landau_prefix * eCR_0 * t0; // dimensionless now and appropriate code units


    // now we are ready to actually integrate these equations, in a numerically-stable manner, with protection from over/under-shooting
    double dtau_max = 1.e-5;
    double dx_e,dx_f,dx_up,dx_um,x_e_0=x_e,x_f_0=x_f,x_up_0=x_up,x_um_0=x_um,dtaux=0.,efmax=50.,expfac;
    double x_e_prev,x_f_prev,x_up_prev,x_um_prev,n_eqm_loops=1.; fmax=1./EPSILON_SMALL; // (need to set initial fmax to large value)
    long n_iter=0, n_iter_max=100000; // sets the maximum number of sub-cycles which we will allow below for any sub-process
    while(1)
    {
        /* here's the actual set of remaining stiff equations to be solved
            dx_e  = gammCR*(x_up+x_um)*x_e + (x_um-x_up)*x_f;                     // deCR_dt
            dx_f  = -ceff2_va2*(psifac + (x_um-x_up)*x_e + (x_up+x_um)*x_f);      // df_dt
            dx_up = -x_up*(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_up - x_f);    // deAp_dt
            dx_um = -x_um*(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_um + x_f);    // deAm_dt
        */
        x_e_prev=x_e; x_f_prev=x_f; x_up_prev=x_up; x_um_prev=x_um; // reset values at the beginning of the loop (these will be cycled multiple times below)

        // for eqm: if psi>0, f<0, um->grows, up->pure-damping //
        double f_eqm, up_eqm, um_eqm, tinv_u, tinv_f;
        double q_tmp = (gammCR-1.)*x_e + gamma_in_t_ll, q_inner = (4.*gamma_nll*fabs(psifac)) / (q_tmp*q_tmp);
        if(q_inner < 1.e-4) {q_inner=q_inner/2.;} else {q_inner=sqrt(1.+q_inner)-1.;}
        double x_nonzero = q_tmp * q_inner / (2.*gamma_nll); if(fabs(gamma_nll) < EPSILON_SMALL) {x_nonzero = fabs(psifac)/q_tmp;}    
        double x_f_magnitude = x_e + fabs(psifac) / x_nonzero;
        if(psifac > 0)
        {
            up_eqm=xkappa_min; um_eqm=x_nonzero; f_eqm=-x_f_magnitude;
            tinv_u = EPSILON_SMALL + fabs(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_um + x_f);
        } else {
            um_eqm=xkappa_min; up_eqm=x_nonzero; f_eqm=+x_f_magnitude;
            tinv_u = EPSILON_SMALL + fabs(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_up - x_f);
        }
        tinv_f = fabs( ceff2_va2*(psifac + (x_um-x_up)*x_e + (x_up+x_um)*x_f) ) * (1./(EPSILON_SMALL + fabs(f_eqm)) + 1./(EPSILON_SMALL+fabs(x_f)));
        double t_eqm = 1./(tinv_u + tinv_f); // timescale to approach equilibrium solution
        
        // set timestep (steadily  growing from initial conservative value ) //
        dtaux = dtau; if(dtaux > dtau) {dtaux=dtau;}
        if(dtaux > dtau_max) {dtaux = dtau_max;}
        dtau_max *= 2.; if(dtaux > 10.) {dtaux=10.;}
        
        double jump_fac = 0.5; // fraction towards equilibrium to 'jump' each time
        //if(dtaux >= 0.33*jump_fac*t_eqm)
        if(dtau >= jump_fac*t_eqm)
        {
            // timestep is larger than the timescale to approach the equilibrium solution, 
            //  so move the systems towards equilibrium, strictly
            //
            if(dtaux > jump_fac*t_eqm) {dtaux = jump_fac*t_eqm;} else {jump_fac = dtaux/t_eqm;} // initial 'step' is small fraction of equilibrium
            dtaux = n_eqm_loops * t_eqm; jump_fac = 1. + (jump_fac-1.)/n_eqm_loops; n_eqm_loops*=1.1; // each sub-cycle consecutively in eqm, allow longer step
            if((x_f<=-fmax && f_eqm<=-fmax) || (x_f>=+fmax && f_eqm>=+fmax)) {t_eqm=1./EPSILON_SMALL; jump_fac=1.; dtaux=dtau;} // if slamming into limits, terminate cycle with big step
            if(f_eqm > 0)
            {
                x_up = exp( log(x_up)*(1.-jump_fac) + log(up_eqm)*jump_fac );
                if(x_f > 0) {x_f = +exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );} else {
                    if(fabs(x_f)<10.*fabs(f_eqm)) {x_f=x_f*(1.-jump_fac)+f_eqm*jump_fac;} else {
                        x_f = -exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );}}
            } else {
                x_um = exp( log(x_um)*(1.-jump_fac) + log(um_eqm)*jump_fac );
                if(x_f < 0) {x_f = -exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );} else {
                    if(fabs(x_f)<10.*fabs(f_eqm)) {x_f=x_f*(1.-jump_fac)+f_eqm*jump_fac;} else {
                        x_f = +exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );}}
            }

        } else {

            // timestep is smaller than the timescale to approach equilibrium, so integrate directly, 
            //  but we will still use a fully implicit backwards-Euler type scheme for the two 'stiffest' 
            //  components of the system (namely, the flux and eA term corresponding to the multiplicative direction)
            //  [the other terms, e.g. the damped energy change and the CR energy change, can be dealt with after]
            //
            double x_dum=0, x_out=0; n_eqm_loops=1.; // (if we enter this,  we need to terminate the parent loop above)
            if(f_eqm>0) {x_dum=x_um;} else {x_dum=x_up;}
            double q0 = 1.+ceff2_va2*dtaux*x_dum, g00 = gammCR*x_e + gamma_in_t_ll, psi00 = psifac + x_dum*x_e;
            double a_m1 = -x_up_prev/dtaux, a_0 = g00 + 1./dtaux , a_1 = gamma_nll, c2dt = ceff2_va2*dtaux;
            if(f_eqm<0) {psi00=psifac-x_dum*x_e; a_1=-gamma_nll; a_0=-(g00 + 1./dtaux); a_m1=x_um_prev/dtaux;}
            double d0 = -a_m1*q0, c0 = x_f_prev - a_0*q0 - c2dt*(a_m1 + psi00), b0 = -a_1*q0 - c2dt*(a_0-x_e), a0 = -a_1*c2dt; 
            if(f_eqm<0) {b0 = -a_1*q0 - c2dt*(a_0+x_e);;}
            if(fabs(a0) < EPSILON_SMALL)
            {
                if(fabs(b0) < EPSILON_SMALL)
                {
                    x_out = fabs(d0/c0); // linear solve
                } else {
                    d0/=c0; b0/=c0; c0=fabs(4.*b0*d0); if(c0<1.e-4) {c0*=0.5;} else {c0=sqrt(1.+c0)-1.;}
                    x_out = c0/(2.*fabs(b0)); // quadratic solve
                }    
            } else {
                // cubic solve
                double p0=-b0/(3.*a0), q0=p0*p0*p0 + (b0*c0-3.*a0*d0)/(6.*a0*a0), r0=c0/(3.*a0), f0=r0-p0*p0, g0=q0*q0 + f0*f0*f0; 
                if(g0 >= 0.)
                {
                    g0=sqrt(g0); a0=q0+g0; b0=q0-g0; 
                    q0=pow(fabs(a0),1./3.); if(a0<0) {q0*=-1.;}
                    r0=pow(fabs(b0),1./3.); if(b0<0) {r0*=-1.;}
                } else {
                    g0=sqrt(-g0); a0=sqrt(q0*q0+g0*g0); b0=atan(g0/q0); r0=0.; q0=2.*a0*cos(b0/3.);
                }
                x_out = fabs(p0 + q0 + r0);
            }
            x_f = a_m1/x_out + a_0 + a_1*x_out;
            if(f_eqm>0) {x_up=x_out;} else {x_um=x_out;}

        }

        // now deal with the non-stiff part of the equations, namely the other Alfven-energy component + CR energy
        if(f_eqm > 0) // do the evolution for the eA term -not- involved in the stiff part of the equations //
        {
            double g0 = gammCR*x_e + gamma_in_t_ll + x_f; // pure-damping for um
            if(g0 > 0) {
                expfac=g0*dtaux; if(expfac>efmax) {expfac=efmax;} 
                if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
                x_um /= (1. + (1.+gamma_nll*x_um/g0)*expfac);
            } else {x_um -= x_um*(g0 + gamma_nll*x_up)*dtaux;} // (linear if x_f hasn't behaved yet)
        } else {
            double g0 = gammCR*x_e + gamma_in_t_ll - x_f; // pure-damping for up
            if(g0 > 0) {
                expfac=g0*dtaux; if(expfac>efmax) {expfac=efmax;} 
                if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
                x_up /= (1. + (1.+gamma_nll*x_up/g0)*expfac);
            } else {x_up -= x_up*(g0 + gamma_nll*x_up)*dtaux;} // (linear if x_f hasn't behaved yet)
        }
        
        // calculate total-energy damping (needed for deriving change in e_cr, which is then given by energy conservation) //
        double x_um_eff=0.5*(x_um+x_um_prev), x_up_eff=0.5*(x_up+x_up_prev), x_f_eff=0.5*(x_f+x_f_prev); // effective values for use in damping rates below
        expfac=gamma_in_t_ll*dtaux; if(expfac>efmax) {expfac=efmax;} 
        if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
        double de_damp = x_um_eff/(1.+1./(expfac*(1.+gamma_nll*x_um_eff/gamma_in_t_ll))) + 
                         x_up_eff/(1.+1./(expfac*(1.+gamma_nll*x_up_eff/gamma_in_t_ll))); // energy loss to thermalized wave-damping
        if(!isfinite(de_damp)) {de_damp=0;}
        double e_tot = DMAX(x_um_prev,0) + DMAX(x_up_prev,0) + DMAX(x_e_prev,0) - de_damp; // total energy (less damping) before step
        double x_e_egycon = DMAX(e_tot-(x_up+x_um), Min_Egy);
        expfac = gammCR*(x_up_eff+x_um_eff)*dtaux; double x_numer = x_e_prev + dtaux*(x_um_eff-x_up_eff)*x_f_eff;
        if(expfac<0.9 && x_numer>0.) {x_e=x_numer/(1.-expfac);} else {
            if(expfac*x_e_prev+x_numer > 0.01*x_e) {x_e=expfac*x_e_prev+x_numer;} else {x_e*=0.01;}}
        if(fabs(x_e_egycon-x_e_prev) < fabs(x_e-x_e_prev)) {x_e=x_e_egycon;}
        if(e_tot < Min_Egy || !isfinite(e_tot)) {e_tot = Min_Egy;} // enforce minima/maxima
        if(x_e   < Min_Egy || !isfinite(x_e)  ) {x_e   = Min_Egy;} // enforce minima/maxima
        if(x_um<EPSILON_SMALL || !isfinite(x_um)) {x_um=EPSILON_SMALL;} else {if(x_um>xkappa_max) {x_um=xkappa_max;}} // enforce minima/maxima
        if(x_up<EPSILON_SMALL || !isfinite(x_up)) {x_up=EPSILON_SMALL;} else {if(x_up>xkappa_max) {x_up=xkappa_max;}} // enforce minima/maxima
        if(x_um+x_up<xkappa_min) {expfac=xkappa_min/(x_um+x_up); x_um*=expfac; x_up*=expfac;} // only want to enforce -sum- having effective diffusivity, not both
        double e_tot_new=x_e+x_um+x_up; x_e*=e_tot/e_tot_new; x_up*=e_tot/e_tot_new; x_um*=e_tot/e_tot_new; // check energy after limit-enforcement
	    fmax = x_e*sqrt(ceff2_va2); if(!isfinite(x_f)) {x_f=0;} else {if(x_f>fmax) {x_f=fmax;} else {if(x_f<-fmax) {x_f=-fmax;}}} // check for flux maximum/minimum
        
        // calculate change in parameters to potentially break the cycle here
        dx_e  = (x_e - x_e_prev) / (EPSILON_SMALL + x_e + x_e_prev);
        dx_up = (x_up - x_up_prev) / (EPSILON_SMALL + x_up + x_up_prev);
        dx_um = (x_um - x_um_prev) / (EPSILON_SMALL + x_um + x_um_prev);
        dx_f  = (x_f - x_f_prev) / (EPSILON_SMALL + fabs(x_f) + fabs(x_f_prev));
        double dx_max = sqrt(dx_e*dx_e + dx_up*dx_up + dx_um*dx_um + dx_f*dx_f); // sum in quadrature
        if(!isfinite(dx_max)) {dx_max=1;} // enforce validity for check below 

        if((n_iter > 0) && (n_iter % 10000 == 0)) // print diagnostics if the convergence is happening slowly
        {
            printf("niter/max=%ld/%ld dtau/step=%g/%g d_params=%g (init/previous/now) xeCR=%g/%g/%g xeAp=%g/%g/%g xeAm=%g/%g/%g xflux=%g/%g/%g ceff2_va2=%g damp_g_intll=%g damp_g_nll=%g psi_gradientfac=%g heat_term=%g xkappa_min/max=%g/%g egy_min=%g flux_max=%g \n",
                n_iter,n_iter_max,dtau,dtaux,dx_max,x_e_0,x_e_prev,x_e,x_up_0,x_up_prev,x_up,x_um_0,x_um_prev,x_um,x_f_0,x_f_prev,x_f,ceff2_va2,gamma_in_t_ll,gamma_nll,psifac,(x_e_0+x_up_0+x_um_0)-(x_e+x_up+x_um),xkappa_min,xkappa_max,Min_Egy,fmax);
            fflush(stdout);
        }
        dtau -= dtaux; // subtract the time we've already integrated from the total timestep
        n_iter++; // count the number of iterations
        if(dtau <= 0.) break; // we have reached the end of the integration time for our sub-stepping. we are done!
        if(dx_max <= 1.e-3*dtaux/dtau) break; // the values of -all- the parameters are changing by much less than the floating-point errors. we are done!
        if(n_iter > n_iter_max) break; // we have reached the maximum allowed number of iterations. we give up!
    }

    // ok! done with the main integration/sub-cycle loop, now just do various clean-up operations
    double thermal_heating = eCR_0 * ((x_e_0+x_up_0+x_um_0)-(x_e+x_up+x_um)); // net thermalized energy from damping terms
    if(thermal_heating < 0 || !isfinite(thermal_heating)) // if this is less than zero (from residual floating-point error), then the energy goes up, which shouldnt happen: set to zero
    {
        e_tot=x_e_0+x_up_0+x_um_0; e_tot_new = x_e+x_up+x_um; // initial and final energies should be equal in this case
        x_e *= e_tot/e_tot_new; x_up *= e_tot/e_tot_new; x_um *= e_tot/e_tot_new; // enforce that equality
        thermal_heating=0; // set thermal change to nil
    }
    eCR=eCR_0*x_e; eA[0]=eCR_0*x_up; eA[1]=eCR_0*x_um; f_CR=f_unit*x_f; Min_Egy*=eCR_0; xkappa_min*=eCR_0; xkappa_max*=eCR_0; // re-assign dimensional quantities

    // assign the updated values back to the resolution elements, finally!
    if(mode==0) {SphP[i].CosmicRayEnergy=eCR;} else {SphP[i].CosmicRayEnergyPred=eCR;} // CR energy
    for(k=0;k<2;k++) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[k]=eA[k];} else {SphP[i].CosmicRayAlfvenEnergyPred[k]=eA[k];}} // Alfven energy
    if(mode==0) {for(k=0;k<3;k++) {SphP[i].CosmicRayFlux[k]=f_CR*bhat[k];}} else {for(k=0;k<3;k++) {SphP[i].CosmicRayFluxPred[k]=f_CR*bhat[k];}} // assign to flux vector
    if(mode==0) {SphP[i].InternalEnergy+=thermal_heating/P[i].Mass;} else {SphP[i].InternalEnergyPred+=thermal_heating/P[i].Mass;} // heating term from damping
    SphP[i].CosmicRayDiffusionCoeff = 1. / (fac_Omega*((eA[0]+eA[1])/E_B)/(clight_code*clight_code)); // effective diffusion coefficient in code units

#endif

    return 1;
}
#endif


