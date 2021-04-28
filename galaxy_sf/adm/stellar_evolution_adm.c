#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif

/* Routines for models that require stellar evolution: luminosities, mass loss, SNe rates, etc.
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifdef GALSF
#ifdef ADM



/* this routine tells the feedback algorithms what to 'inject' when a stellar feedback event occurs.
    you must define the mass, velocity (which defines the momentum and energy), and metal content (yields)
    of the ejecta for the event[s] of interest. Mass [Msne] and velocity [SNe_v_ejecta] should
    be in code units. yields[k] should be defined for all metal species [k], and in dimensionless units
    (mass fraction of the ejecta in that species). */
void particle2in_addFB_fromstars_adm(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration)
{
//#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
//#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
//    if(fb_loop_iteration == 0) {particle2in_addFB_SNe(in,i);}
//    if(fb_loop_iteration == 1) {particle2in_addFB_winds(in,i);}
//    if(fb_loop_iteration == 2) {particle2in_addFB_Rprocess(in,i);}
//    if(fb_loop_iteration == 3) {particle2in_addFB_ageTracer(in,i);}
//    return;
//#endif
    if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    else {in->Msne=0; return;} // For now, I don't want any ADM feedback. Subject to change later.

    // 'dummy' example model assumes all SNe are identical with IMF-averaged properties from the AGORA model (Kim et al., 2016 ApJ, 833, 202)
//    in->Msne = P[i].SNe_ThisTimeStep * (14.8/UNIT_MASS_IN_SOLAR); // assume every SNe carries 14.8 solar masses (IMF-average)
//    in->SNe_v_ejecta = 2607. / UNIT_VEL_IN_KMS; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun], which is IMF-averaged
//#ifdef SINGLE_STAR_SINK_DYNAMICS // if single-star exploding or returning mass, use its actual mass & assumed energy to obtain the velocity
//    in->Msne = DMIN(1.,P[i].SNe_ThisTimeStep) * P[i].Mass; // mass fraction of star being returned this timestep
//    in->SNe_v_ejecta = sqrt(2.*(1.e51/UNIT_ENERGY_IN_CGS)/P[i].Mass); // for SNe [total return], simple v=sqrt(2E/m)should be fine without relativistic corrections
//    if(P[i].SNe_ThisTimeStep<1) {double m_msun=P[i].Mass*UNIT_MASS_IN_SOLAR; in->SNe_v_ejecta = (616. * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / UNIT_VEL_IN_KMS;} // scaling from size-mass relation+eddington factor, assuming line-driven winds //
//#endif
//#ifdef METALS
//    int k; for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
//    if(NUM_LIVE_SPECIES_FOR_COOLTABLES>=10) {in->yields[1] = 0.4;} // (catch for Helium, which the above scaling would give bad values for)
//#endif
//#endif
}


/* this routine calculates the event rates for different types of mechanical/thermal feedback
    algorithms. things like SNe rates, which determine when energy/momentum/mass are injected, should go here.
    you can easily modify this to accomodate any form of thermal or mechanical feedback/injection of various
    quantities from stars. */
double mechanical_fb_calculate_eventrates_adm(int i, double dt)
{
//#if defined(GALSF_FB_MECHANICAL) && defined(GALSF_FB_FIRE_STELLAREVOLUTION) // FIRE-specific stellar population version: separate calculation for SNe, stellar mass loss, R-process injection //
//    double SNe_rate = mechanical_fb_calculate_eventrates_SNe(i,dt);
//    mechanical_fb_calculate_eventrates_Winds(i,dt);
//    mechanical_fb_calculate_eventrates_Rprocess(i,dt);
//    mechanical_fb_calculate_eventrates_Agetracers(i,dt);
//    return SNe_rate;
//#endif

//#if defined(SINGLE_STAR_SINK_DYNAMICS) && !defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION) /* SINGLE-STAR version: simple implementation of single-star wind mass-loss and SNe rates */
//    double m_sol,l_sol; m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR; l_sol=bh_lum_bol(0,P[i].Mass,i)*UNIT_LUM_IN_SOLAR;
//#ifdef SINGLE_STAR_FB_WINDS
//    double gam=DMIN(0.5,3.2e-5*l_sol/m_sol), alpha=0.5+0.4/(1.+16./m_sol), q0=(1.-alpha)*gam/(1.-gam), k0=1./30.; // Eddington factor (~L/Ledd for winds), capped at 1/2 for sanity reasons, approximate scaling for alpha factor with stellar type (weak effect)
//    P[i].SNe_ThisTimeStep = DMIN(0.5, (2.338 * alpha * pow(l_sol,7./8.) * pow(m_sol,0.1845) * (1./q0) * pow(q0*k0,1./alpha) / m_sol) * (dt*UNIT_TIME_IN_GYR)); // Castor, Abbot, & Klein scaling
//#endif
//#ifdef SINGLE_STAR_FB_SNE
//    double t_lifetime_Gyr = 10.*(m_sol/l_sol) + 0.003; /* crude estimate of main-sequence lifetime, capped at 3 Myr*/
//    if(evaluate_stellar_age_Gyr(P[i].StellarAge) >= t_lifetime_Gyr) {P[i].SNe_ThisTimeStep=1;}
//#endif
//    return 1;
//#endif

//#ifdef GALSF_FB_THERMAL /* STELLAR-POPULATION version: pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly */
//    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
//    if(evaluate_stellar_age_Gyr(P[i].StellarAge) < 0.005) {return 0;} // enforce age limit of 5 Myr
//    P[i].SNe_ThisTimeStep = P[i].Mass*UNIT_MASS_IN_SOLAR / 91.; // 1 event per 91 solar masses
//    return 1;
//#endif

//#ifdef GALSF_FB_MECHANICAL /* STELLAR-POPULATION version: mechanical feedback: 'dummy' example model below assumes a constant SNe rate for t < 30 Myr, then nothing. experiment! */
//    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
//    if(star_age < 0.03)
//    {
//        double RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
//        double p = RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
//        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
//        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
//        return RSNe;
//    }
//#endif

    return 0;
}







#endif // ADM
#endif /* GALSF */
