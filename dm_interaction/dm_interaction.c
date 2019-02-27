#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include "dm_interaction.h"
#include "../allvars.h"
#include "../proto.h"
#ifdef DM_BARYON_INTERACTION

static double integrate_fyp_xp(double y);
static double integrant_fyp_xp(double z, void *params);


static gsl_interp_accel *acc_fyp_xp;
static gsl_spline *spline_fyp_xp;

double energy_transfer_rate_in_m5_over_s3(double v_rel, double sigma_u)
{
  double km_in_m = 1.e3, cm_in_m = 1.e-2;
  double y =v_rel / sigma_u;
  double fyp, rate;
  double LENGTH_physical = All.UnitLength_in_cm * All.cf_atime / All.HubbleParam;

	fyp = get_fyp_xp(y);
	rate = pow(sigma_u,3) * (CROSS_SECTION_XP / LENGTH_physical / LENGTH_physical) * pow(v_rel / REF_VELOCITY_XP, P_XP-4) * fyp;
	return rate;

}

double integrate_fyp_xp(double y)
{
  double fyp, err;
  gsl_function f;
  gsl_integration_workspace *w;

  size_t n = 1000;
  double z_max = 6.e0/y, z_min = 0.e0;
  double tol = 1.e-7;

  if(P_XP == 0) z_max = z_max<1 ? z_max: 1.e0;

  w= gsl_integration_workspace_alloc(n);
  f.function = &integrant_fyp_xp;
  f.params = &y;

  gsl_integration_qags(&f, z_min, z_max, 0, tol, n, w, &fyp, &err);
  fyp *= -pow(y,P_XP+2)/sqrt(2*M_PI);

  gsl_integration_workspace_free(w);

  return fyp;

}

double integrant_fyp_xp(double z, void *params)
{
  double gyp;
  double y, azm1;

  y = *(double *) params;

  if(P_XP == 0){

    gyp = 0.e0;
    if(z <= 1) gyp = -2*z;
  
  }else if(P_XP == 1){
 
    gyp = -z;
    azm1 = fabs(z-1);
    if(azm1 > 1.e-7) gyp += (z*z-1) * log((z+1)/azm1) /2;

  }else{

    gyp = (z-P_XP) * pow(z+1,P_XP);
    if(z>1){
      gyp -= (z+P_XP) * pow(z-1,P_XP);
    }else{
      gyp +=(z+P_XP) * pow(1-z,P_XP);
    }

   gyp /= P_XP*P_XP-1; 
  
  }

  gyp *= z*exp(-y*y*z*z/2);

  return gyp;

}

int init_fyp_xp(int output_fyp)
{

  double y, fyp, lny[NUM_Y_XP], lnfyp[NUM_Y_XP], dlny, lny_max, lny_min;
  int i;

  if(P_XP < 0){
    printf("error: P_XP should not be negative");
    return 1;
  }

  lny_max = log(Y_MAX_XP);
  lny_min = log(Y_MIN_XP);
  dlny = (lny_max-lny_min)/(NUM_Y_XP-1);

  for(i = 1; i <= NUM_Y_XP; i++){

    y = lny_min + dlny * (i-1);
    lny[i-1] = y;
    y = exp(y);
    fyp = integrate_fyp_xp(y);
    lnfyp[i-1] = log(fyp);

    if(output_fyp) printf("%e %e %e\n", P_XP, y, fyp);

}

  acc_fyp_xp = gsl_interp_accel_alloc();
  spline_fyp_xp = gsl_spline_alloc(gsl_interp_cspline,NUM_Y_XP);
  gsl_spline_init(spline_fyp_xp,lny,lnfyp,NUM_Y_XP);

  return 0;
}

int fin_fyp_xp(void)
{
  gsl_spline_free(spline_fyp_xp);
  gsl_interp_accel_free(acc_fyp_xp);

  return 0;

}

double get_fyp_xp(double y)
{

  double fyp, lny;

  if(spline_fyp_xp == NULL) init_fyp_xp(0);

  if(y>Y_MAX_XP/2){

    fyp = pow(y,P_XP-1)*(1+P_XP*(P_XP-3)/2/y/y);

  }else if(y<Y_MIN_XP*2){


    fyp = y*y*pow(2,P_XP/2+1)/3/sqrt(2*M_PI) * gsl_sf_gamma(P_XP/2+1) * (1+(P_XP-3)*y*y/10);

  }else{

    lny =log(y);
    fyp = gsl_spline_eval(spline_fyp_xp,lny,acc_fyp_xp);
    fyp = exp(fyp);

  }

  return fyp;

}
#endif
