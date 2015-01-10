//
//  blackhole.h
//  gizmo
//
//  Created by Paul Torrey on 10/21/14.
//
//

#ifndef gizmo_blackhole_h
#define gizmo_blackhole_h


#ifndef BH_CSND_FRAC_BH_MERGE
/* Relative velocity fraction (in units of soundspeed) for merging blackholes, default=1.0 */
#define BH_CSND_FRAC_BH_MERGE 1.0
#endif



#if defined(BLACK_HOLES)
void blackhole_start(void);
void blackhole_end(void);
void blackhole_environment_loop(void);
void blackhole_properties_loop(void);
void blackhole_feedback_loop(void);
void blackhole_final_loop(void);
void blackhole_feed_loop(void);
int blackhole_environment_evaluate(int target, int mode, int *nexport, int *nSend_local);
void out2particle_blackhole(struct blackhole_temp_particle_data *out, int target, int mode);

//void check_for_bh_merger(int j, MyIDType id);
void normalize_temp_info_struct(int i);
void set_blackhole_mdot(int i, int n, double dt);
void set_blackhole_new_mass(int i, int n, double dt);
#if defined(BH_DRAG) || defined(BH_DYNFRICTION)
void set_blackhole_drag(int i, int n, double dt);
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
void set_blackhole_long_range_rp(int i, int n);
#endif


#endif


#endif
