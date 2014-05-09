/*! \file sidm_proto.h
 *  \brief this file contains all function prototypes for sidm module
 */

double prob_of_interaction(double r, double h_si,  double Vtarget[3], double Vno[3], int dt_step);
double g_geo(double r);
void calculate_interact_kick(double Vtarget[3], double Vno[3], double kick_target[3], double kick_no[3]);

void init_geofactor_table(void);
double geofactor_integ(double x, void * params);
double geofactor_angle_integ(double u, void * params);
double kernel(double u);

void update_interaction_table(MyIDType id1, MyIDType id2);
int  check_interaction_table(MyIDType id1, MyIDType id2);
void AllocateInteractionTable(int x, int y);
void init_self_interactions();
void log_self_interactions();
