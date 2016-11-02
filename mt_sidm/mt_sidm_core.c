#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#define GSLWORKSIZE 100000

/*! This routine sets the kicks for each particle after it has been decided that they will
 *  interact. It uses an algorithm tha conserves energy and momentum but picks a random direction
 *  so it does not conserves angular momentum.
*/
void mt_calculate_interact_kick(Mass, double Vtarget[3], double Vno[3], double new_Vtarget[3], double new_Vno[3])
{
	double dV, theta, phi, dvx, dvy, dvz;
	double M, Vcm[3]

 	// calculate the center of mass velocity
	M = 2*Mass;
	Vcm[0] = Vtarget[0] + Vno[0];
	Vcm[1] = Vtarget[1] + Vno[1];
	Vcm[2] = Vtarget[1] + Vno[2];

	//
	dvx = Vno[0]-Vtarget[0];
	dvy = Vno[1]-Vtarget[1];
	dvz = Vno[2]-Vtarget[2];
	dV = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

	// determin unit vector randomly
	theta = acos(2*gsl_rng_uniform(random_generator)-1.0);
	phi = gsl_rng_uniform(random_generator)*2.0*M_PI;

	// 
	new_Vtarget[0] = Vcm[0] + dvx * 0.5 * sin(theta)*cos(phi);
	new_Vtarget[1] = Vcm[1] + dvx * 0.5 * sin(theta)*sin(phi);
	new_Vtarget[2] = Vcm[2] + dvx * 0.5 * cos(theta);
}

void mt_calculate_cross_section(double Vtarget[3], double Vno[3], double *CrossSection)
{
	double dV, dvx, dvy, dvz;

	dvx = Vno[0]-Vtarget[0];
	dvy = Vno[1]-Vtarget[1];
	dvz = Vno[2]-Vtarget[2];
	dV = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

	CrossSection = All.mtSIDMparameterA * pow(dV, All.mtSIDMparameterarpha); 
}
