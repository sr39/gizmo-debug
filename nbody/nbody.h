/* File nbody.h
Contains function prototypes for routines used to do sub-cycled few-N-body evolution in star cluster simulations
 */

#ifndef gizmo_nbody_h
#define gizmo_nbody_h
#endif

void kepler_timestep(int i, double dt, double kick_dv[3], double drift_dx[3], int mode);
double eccentric_anomaly(double mean_anomaly, double ecc);
