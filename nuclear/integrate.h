#include "./network.h"
/*
 *  This code is place-holder, inherited from GADGET3,
 *   to be replaced by David Radice's version (written completely independently)
 */

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw);
int network_integrate( double temp, double rho, const double *x, double *dx, double dt, double *dedt, double *drhodt, const struct network_data *nd, struct network_workspace *nw );
