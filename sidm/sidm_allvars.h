/*! \file sidm_allvars.h
 *  \brief declares global variables for rochas_sidm module.
 */

#define GEOFACTOR_TABLE_LENGTH 1000    /*!< length of the table used for the geometric factor spline */
#define INTERACTION_TABLE_LENGTH 5000  /*!< This should be about the maximum number of interactions expected at each timestep */
#define PARTICLE_MAX_INTERACTIONS 1000 /*!< Maximum number of interactions a particle can have at each time step */

extern MyDouble GeoFactorTable[GEOFACTOR_TABLE_LENGTH];
extern MyIDType** InteractionTable;

