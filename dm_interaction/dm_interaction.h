#include <stdio.h>
#include "../allvars.h"

int init_fyp_xp(int output_fyp);
int fin_fyp_xp(void);
double get_fyp_xp(double y);
double energy_transfer_rate_in_m5_over_s3(double v_rel, double sigma_u);
void do_sph_interaction_kick(int i, integertime tstart, integertime tend, double dt_entr);

int dm_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *dm_density_evaluate_primary(void *p);
void *dm_density_evaluate_secondary(void *p);
void *dm_hydro_evaluate_primary(void *p);
void *dm_hydro_evaluate_secondary(void *p);
int dm_density_isactive(int n);



// Following prameters specify the momentum transfer cross section between baryon and dark matter.
// // Cross section here is dependent of relative velocity as sigma_xp(v):=sigma_xp(v_ref)*(v/v_ref)^(p-4)
 #define CROSS_SECTION_XP 1.e-30/// cross section in (cm^2) at reference relative velocity
 #define REF_VELOCITY_XP 1.e3 // reference relative velocity in (km/sec)
 #define P_XP 4 // slope of velocity-dependence; P_XP >= 0 is required; 4 for velocity-independent cross section
//
// // Maximum and minimum of y:=(relative gas-DM velocity)/(baryon velocity dispersion) between 
// // which energy transfer rate is calculated. Above and below this value, analytic formulae are used.
// // Probably there's no need to change these values.
 #define Y_MAX_XP 1.e2
 #define Y_MIN_XP 1.e-2
 #define NUM_Y_XP 100 // number of y-bins with logarithmically uniform spacing
//

//
