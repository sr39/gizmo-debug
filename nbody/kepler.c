#ifdef SINGLE_STAR_SUPERTIMESTEPPING
/* Advances the binary by timestep dt */

void kepler_timestep(int i, dt){
  //  MyDouble comp_Pos[3]; //position of binary companion                        
  //  MyDouble comp_Vel[3]; //velocity of binary companion                        
  //  MyDouble comp_Mass; //mass of binary companion
    double h[3]; //Specific angular momentum vector
    double dr = sqrt(P[i].comp_dx[0]*P[i].comp_dx[0] + P[i].comp_dx[1]*P[i].comp_dx[1] + P[i].comp_dx[2]*P[i].comp_dx[2]);
    double dx_normalized[3] = {P[i].comp_dx[0]/dr, P[i].comp_dx[1]/dr, P[i].comp_dx[2]/dr};
    double n_x[3]; // normalized Laplace-Runge-Lenz vector, just to get the unit vector along the major axis of the binary
    double n_y[3]; // normalized unit vector along the minor axis of the binary
    double norm, h2;
    double x = 0, y =0, vx =0, vy = 0; // Coordinates in the frame aligned with the binary
    int k,l,m;
    double Mtot = P[i].Mass + P[i].comp_mass;
//    double mu = P[i].Mass * P[i].comp_mass / Mtot;
    

    for(k=0; k<3; l=(k+1)%3; m=(k+2)%3; k++){ // dx cross dv to get specific angular momentum vector
	h[k] = P[i].comp_dx[l]*P[i].comp_dv[m] - P[i].comp_dx[m]*P[i].comp_dv[l]
    }

    h2 = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    for(k=0; k<3; l=(k+1)%3; m=(k+2)%3; k++){ // dx cross dv to get specific angular momentum vector
	l = (k+1)%3; m = (k+2)%3;
	n_x[k] = P[i].comp_dv[l]*h[m] - P[i].comp_dv[m]*h[l] - All.G * Mtot * dx_normalized[k]; // Worry about cancellation error for low eccentricity?
    }

    norm = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);    
    for(k=0; k<3; k++) n_x[k] /= -norm; // take the opposite direction of the LRL vector, so x points from periapsis to apoapsis

    for(k=0; k<3; l=(k+1)%3; m=(k+2)%3; k++){ // cross product of n_x with angular momentum to get a vector along the minor axis
	l = (k+1)%3; m = (k+2)%3;
	n_y[k] = n_x[l] * h[m] - n_x[m] * h[l];
	n_y[k] /= sqrt(h2);
    }

    for(k=0; k<3; k++){
	x += dx[k]*n_x[k];
	y += dx[k]*n_y[k];
	vx += dv[k]*n_x[k];
	vy += dv[k]*n_y[k];
    }
    
    

}

// Solve Kepler's equation to turn eccentric anomaly into mean anomaly
double eccentric_anomaly(double mean_anomaly){

}

#endif
