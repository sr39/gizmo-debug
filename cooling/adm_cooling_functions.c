#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "./adm_cooling_functions.h"


//////////////////////////////
//  Collisional Excitation  //
//////////////////////////////


// General functions to be inserted in the numerical integrator

double g0(double x) {
	//return (1 / x)* (1 / pow(1 + 4.0 * pow(x, 2) / 9.0, 6)); 
	return 9.0 * (299619.0 + 374220.0 * pow(x, 2.0) + 203040.0 * pow(x, 4.0) + 51840.0 * pow(x, 6.0) + 5120.0 * pow(x, 8.0)) / (40.0 * pow(9.0 + 4.0 * pow(x, 2.0), 5.0)) + log(x) - 0.5 * log(9.0 + 4.0 * pow(x, 2.0));
}

double g(double u, void* params) {
	double f;
	double y2 = *(double*)params;
	double xplus = (u / sqrt(y2)) * (1.0 + sqrt(1.0 - 3.0 * y2 * pow(1.0 / u, 2.0) / 4.0));
	double xminus = (u / sqrt(y2)) * (1.0 - sqrt(1.0 - 3.0 * y2 * pow(1.0 / u, 2.0) / 4.0));
	if (xminus == 0.0) {
		double x = 0.75 * y2 * pow(1.0 / u, 2.0);
		xminus = (u / sqrt(y2)) * (0.5 * x + (1.0 / 8.0) * pow(x, 2.0) + (1.0 / 16.0) * pow(x, 3.0));
	}
	f = (u * exp(-pow(u, 2.0)) / (1.0 + 7.0 * y2 / (4 * pow(u, 2)))) * (g0(xplus) - g0(xminus));
	return f;
}

// Numerical Integrator

double g_integral(double y2) {
	int GSLWORKSIZE = 1000;
	double result, error;
	gsl_function F; gsl_integration_workspace* workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
	F.function = &g; F.params = &y2;
	gsl_integration_qagiu(&F, (sqrt(3.0 * y2) / 2), 0, 1.0e-5, GSLWORKSIZE, workspace, &result, &error);
	gsl_integration_workspace_free(workspace);
	return result;
}

////////////////////////////////
////  Collisional Ionisation  //
////////////////////////////////

// Function to be integrated

double f(double u, void* params) {
	double y2 = *(double*)params;
	double u2 = pow(u, 2.0);
	return (u * exp(-u2) / (1 + 2 * y2 / u2)) * (1 - y2 / u2 - 0.5 * (1 - pow(y2 / u2, 2.0)) * log(y2 / u2) + (y2 / u2) * log(y2 / u2) / (1 + y2 / u2));
}


// Numerical Integrator

double f_integral(double y2) {
	int GSLWORKSIZE = 1000;
	double result, error;
	gsl_function F; gsl_integration_workspace* workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
	F.function = &f; F.params = &y2;
	gsl_integration_qagiu(&F, sqrt(y2), 0, 1.0e-6, GSLWORKSIZE, workspace, &result, &error);
	gsl_integration_workspace_free(workspace);
	return result;
}


