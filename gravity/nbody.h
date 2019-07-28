/* File nbody.h
Contains function prototypes for routines used to do sub-cycled few-N-body evolution in star cluster simulations
 */

#ifndef gizmo_nbody_h
#define gizmo_nbody_h
#endif

void kepler_timestep(int i, double dt, double kick_dv[3], double drift_dx[3], int mode);
double gravfac(double r, double mass);
double gravfac2(double r, double mass);
void grav_accel_jerk(double mass, double dx[3], double dv[3], double accel[3], double jerk[3]);
double eccentric_anomaly(double mean_anomaly, double ecc);
