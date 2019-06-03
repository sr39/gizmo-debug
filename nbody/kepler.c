#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#ifdef SINGLE_STAR_SUPERTIMESTEPPING
// Solve Kepler's equation to convert mean anomaly into eccentric anomaly
double eccentric_anomaly(double mean_anomaly, double ecc){
    double x0 = mean_anomaly;
    double err = 1e100;
    int iterations = 0;
    double twopi = 2*M_PI;
    while(fabs(err/twopi) > 1e-14 && iterations < 20){ // do Newton iterations
	err = (x0 - ecc*sin(x0) - mean_anomaly)/(1 - ecc*cos(x0));
        x0 -= err;
	x0 = fmod(x0, twopi);
        iterations += 1;
    }
    return x0;
}


/* 
Advances the binary by timestep dt 
mode 0 - Just fill out the particle's kick and drift for the timestep, without doing the update
mode 1 - Actually update the binary separation and relative velocity. This should be done on the full-step drift.
*/

void kepler_timestep(int i, double dt, double kick_dv[3], double drift_dx[3], int mode){
    double h[3]; //Specific angular momentum vector
    double dr = sqrt(P[i].comp_dx[0]*P[i].comp_dx[0] + P[i].comp_dx[1]*P[i].comp_dx[1] + P[i].comp_dx[2]*P[i].comp_dx[2]);
    double dv = sqrt(P[i].comp_dv[0]*P[i].comp_dv[0] + P[i].comp_dv[1]*P[i].comp_dv[1] + P[i].comp_dv[2]*P[i].comp_dv[2]);

    double dx_normalized[3] = {P[i].comp_dx[0]/dr, P[i].comp_dx[1]/dr, P[i].comp_dx[2]/dr};
    double dx_new[3], dv_new[3];
    double n_x[3]; // normalized Laplace-Runge-Lenz vector, just to get the unit vector along the major axis of the binary
    double n_y[3]; // normalized unit vector along the minor axis of the binary
    double norm, h2, true_anomaly, mean_anomaly, ecc_anomaly;
    double x = 0, y =0, vx =0, vy = 0; // Coordinates in the frame aligned with the binary
    int k,l,m;
    double Mtot = P[i].Mass + P[i].comp_Mass;

    double specific_energy = .5*dv*dv - All.G * Mtot / dr;
    double semimajor_axis = -All.G * Mtot / (2*specific_energy);

    for(k=0; k<3; k++, l=(k+1)%3, m=(k+2)%3){ // dx cross dv to get specific angular momentum vector
        h[k] = P[i].comp_dx[l]*P[i].comp_dv[m] - P[i].comp_dx[m]*P[i].comp_dv[l];
    }

    double hSqr = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    double ecc = sqrt(1 + 2 * specific_energy * hSqr / (All.G*All.G*Mtot*Mtot)); 

    for(k=0; k<3; k++, l=(k+1)%3, m=(k+2)%3){ // Get the LRL vector dv x h - GM dx/r
        l = (k+1)%3; m = (k+2)%3;
        n_x[k] = P[i].comp_dv[l]*h[m] - P[i].comp_dv[m]*h[l] - All.G * Mtot * dx_normalized[k]; // Worry about cancellation error for low eccentricity?
    }

    norm = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);    
    for(k=0; k<3; k++) n_x[k] /= -norm; // take the opposite direction of the LRL vector, so x points from periapsis to apoapsis

    for(k=0; k<3; k++, l=(k+1)%3, m=(k+2)%3){ // cross product of n_x with angular momentum to get a vector along the minor axis
        l = (k+1)%3; m = (k+2)%3;
        n_y[k] = n_x[l] * h[m] - n_x[m] * h[l];
        n_y[k] /= sqrt(h2);
    }

    // Transform to coordinates in the plane of the ellipse
    for(k=0; k<3; k++){
        x += P[i].comp_dx[k]*n_x[k];
        y += P[i].comp_dx[k]*n_y[k];
    }
    
    true_anomaly = atan2(y,x);
    mean_anomaly = atan2(sqrt(1 - ecc*ecc) * sin(true_anomaly), ecc + cos(true_anomaly));
    mean_anomaly += dt/P[i].min_bh_t_orbital * 2 * M_PI;
    mean_anomaly = fmod(mean_anomaly, 2*M_PI);
    ecc_anomaly = eccentric_anomaly(mean_anomaly, ecc);

    x = semimajor_axis * (cos(ecc_anomaly) - ecc);
    y = semimajor_axis * sqrt(1 - ecc*ecc) * sin(ecc_anomaly);
    
    dr = sqrt(x*x + y*y);
    dv = sqrt(All.G * Mtot * (2/dr - 1/semimajor_axis)); // We conserve energy exactly

    double v_phi = sqrt(h2) / dr; // conserving angular momentum
    double v_r = sqrt(DMAX(0, dv*dv - v_phi*v_phi));
    if(ecc_anomaly < M_PI) v_r = -v_r; // if radius is decreasing, make sure v_r is negative
    
    //relative velocities in the frame aligned with the ellipse:
    vx = v_phi * (-y/dr) + v_r * x/dr;
    vy = v_phi * (x/dr) + v_r * y/dr;

    // transform back to global coordinates
    for(k=0; k<3; k++){
    dx_new[k] = x * n_x[k] + y * n_y[k];
    dv_new[k] = vx * n_x[k] + vy * n_y[k];
    drift_dx[k] = dx_new[k] - P[i].comp_dx[k];
    kick_dv[k] = dv_new[k] - P[i].comp_dv[k];

    if(mode==1){ // if we want to do the actual self-consistent binary update
        P[i].comp_dx[k] = dx_new[k];
        P[i].comp_dv[k] = dv_new[k];
    }
    }
}

#endif
