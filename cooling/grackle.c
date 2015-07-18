#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#include "./cooling.h"

#ifdef GRACKLE

#include <grackle.h>

#define ENDRUNVAL 12345

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return pressure
//     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0)
//
double CallGrackle(double u_old, double rho, double dt, double *ne_guess, int target, int mode)
{
  gr_float returnval = 0.0;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];
  int i;
  for (i = 0;i < 3;i++) {
    grid_dimension[i] = 1; // the active dimension not including ghost zones.
    grid_start[i]     = 0;
    grid_end[i]       = 0;
  }
  grid_dimension[0] = field_size;
  grid_end[0]       = field_size - 1;
  
  gr_float density, metal_density, energy, vel;
  gr_float cooling_time, temperature, pressure, gamma;
  vel           = 0.0;
  density       = rho;
  energy        = u_old;
  metal_density = density * P[target].Metallicity[0];
  
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
  gr_float ne_density;
  gr_float HI_density, HII_density, HM_density;
  gr_float HeI_density, HeII_density, HeIII_density;
  gr_float H2I_density, H2II_density;
  gr_float DI_density, DII_density, HDI_density;
  gr_float tiny = 1.0e-20;
  
  // Atomic
  ne_density    = density * *ne_guess;
  
  HI_density    = density * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
  HII_density   = density * SphP[target].grHII;
  HM_density    = density * SphP[target].grHM;
  
  HeI_density   = density * SphP[target].grHeI;
  HeII_density  = density * SphP[target].grHeII;
  HeIII_density = density * SphP[target].grHeIII;
  
  H2I_density  = density * tiny;
  H2II_density = density * tiny;
  DI_density   = density * tiny;
  DII_density  = density * tiny;
  HDI_density  = density * tiny;        
  
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
  H2I_density  = density * SphP[target].grH2I;
  H2II_density = density * SphP[target].grH2II;
#endif
  
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
  DI_density   = density * SphP[target].grDI;
  DII_density  = density * SphP[target].grDII;
  HDI_density  = density * SphP[target].grHDI;
#endif
  
  switch(mode) {
  case 0:  //solve chemistry & update values
    if(solve_chemistry(&All.GrackleUnits,
		       All.cf_atime, dt,
		       grid_rank, grid_dimension,
		       grid_start, grid_end,
		       &density, &energy,
		       &vel, &vel, &vel,
		       &HI_density, &HII_density, &HM_density,
		       &HeI_density, &HeII_density, &HeIII_density,
		       &H2I_density, &H2II_density,
		       &DI_density, &DII_density, &HDI_density,
		       &ne_density, &metal_density) == 0) {
      fprintf(stderr, "Error in solve_chemistry.\n");
      endrun(ENDRUNVAL);
    }
    
    // Assign variables back
    *ne_guess            = ne_density    / density;
    
    SphP[target].grHI    = HI_density    / density;
    SphP[target].grHII   = HII_density   / density;
    SphP[target].grHM    = HM_density    / density;
    
    SphP[target].grHeI   = HeI_density   / density;
    SphP[target].grHeII  = HeII_density  / density;
    SphP[target].grHeIII = HeIII_density / density;    
    
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    SphP[target].grH2I   = H2I_density   / density;
    SphP[target].grH2II  = H2II_density  / density;
#endif
    
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    SphP[target].grDI    = DI_density    / density;
    SphP[target].grDII   = DII_density   / density;
    SphP[target].grHDI   = HDI_density   / density;
#endif
    returnval = energy;
    break;
  
  case 1:  //cooling time
    if(calculate_cooling_time(&All.GrackleUnits,
			      All.cf_atime,
			      grid_rank, grid_dimension,
			      grid_start, grid_end,
			      &density, &energy,
			      &vel, &vel, &vel,
			      &HI_density, &HII_density, &HM_density,
			      &HeI_density, &HeII_density, &HeIII_density,
			      &H2I_density, &H2II_density,
			      &DI_density, &DII_density, &HDI_density,
			      &ne_density, &metal_density,
			      &cooling_time) == 0) {
      fprintf(stderr, "Error in calculate_cooling_time.\n");
      endrun(ENDRUNVAL);
    }
    returnval = cooling_time;
    break;
  case 2:  //calculate temperature
    if(calculate_temperature(&All.GrackleUnits,
			     grid_rank, grid_dimension,
			     &density, &energy,
			     &HI_density, &HII_density, &HM_density,
			     &HeI_density, &HeII_density, &HeIII_density,
			     &H2I_density, &H2II_density,
			     &DI_density, &DII_density, &HDI_density,
			     &ne_density, &metal_density,
			     &temperature) == 0) {
      fprintf(stderr, "Error in calculate_temperature.\n");
      endrun(ENDRUNVAL);
    }
    returnval = temperature;
    break;
  case 3:  //calculate pressure
    if(calculate_pressure(&All.GrackleUnits,
			  grid_rank, grid_dimension,
			  &density, &energy,
			  &HI_density, &HII_density, &HM_density,
			  &HeI_density, &HeII_density, &HeIII_density,
			  &H2I_density, &H2II_density,
			  &DI_density, &DII_density, &HDI_density,
			  &ne_density, &metal_density,
			  &pressure) == 0) {
      fprintf(stderr, "Error in calculate_temperature.\n");
      endrun(ENDRUNVAL);
    }
    returnval = pressure;
    break;
  case 4:  //calculate gamma
    if(calculate_gamma(&All.GrackleUnits,
		       grid_rank, grid_dimension,
		       &density, &energy,
		       &HI_density, &HII_density, &HM_density,
		       &HeI_density, &HeII_density, &HeIII_density,
		       &H2I_density, &H2II_density,
		       &DI_density, &DII_density, &HDI_density,
		       &ne_density, &metal_density,
		       &gamma) == 0) {
      fprintf(stderr, "Error in calculate_gamma.\n");
      endrun(ENDRUNVAL);
    }
    returnval = gamma;
    break;
  } //end switch

#else // tabular
  
  switch(mode){
  case 0:  //solve chemistry & update values (table)
    if(solve_chemistry_table(&All.GrackleUnits,
			     All.cf_atime, dt,
			     grid_rank, grid_dimension,
			     grid_start, grid_end,
			     &density, &energy,
			     &vel, &vel, &vel,
			     &metal_density) == 0){
      fprintf(stderr, "Error in solve_chemistry_table.\n");
      endrun(ENDRUNVAL);
    }    
    convert_u_to_temp(energy, rho, ne_guess, target); //need to update *ne_guess for tabular!!, this may be wrong
    returnval = energy;
    break;
  case 1:  //cooling time (table)
    if(calculate_cooling_time_table(&All.GrackleUnits,
				    All.cf_atime, dt,
				    grid_rank, grid_dimension,
				    grid_start, grid_end,
				    &density, &energy,
				    &vel, &vel, &vel,
				    &metal_density,
				    &cooling_time) == 0){
      fprintf(stderr, "Error in calculate_cooling_time.\n");
      endrun(ENDRUNVAL);
    }  
    returnval = cooling_time;
    break;
  case 2:  //calculate temperature (table)
    if(calculate_temperature_table(&All.GrackleUnits,
				   All.cf_atime, dt,
				   grid_rank, grid_dimension,
				   &density, &energy,
				   &metal_density,
				   &temperature) == 0){
      fprintf(stderr, "Error in calculate_temperature.\n");
      endrun(ENDRUNVAL);
    }  
    returnval = temperature;
    break;
  case 3:  //calculate pressure (table)
    if(calculate_pressure_table(&All.GrackleUnits,
				All.cf_atime, dt,
				grid_rank, grid_dimension,
				&density, &energy,
				&pressure) == 0){
      fprintf(stderr, "Error in calculate_pressure.\n");
      endrun(ENDRUNVAL);
    }  
    returnval = pressure;
    break;
  } //end switch

#endif // GRACKLE_CHEMISTRY

  return returnval;
}

//Initialize Grackle
void InitGrackle(void)
{
  grackle_verbose = 0;
  // Enable output
  if(ThisTask == 0)
    grackle_verbose = 1;
  
  // First, set up the units system.
  // These are conversions from code units to cgs.
  All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
  All.GrackleUnits.density_units        = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  All.GrackleUnits.length_units         = All.UnitLength_in_cm / All.HubbleParam;
  All.GrackleUnits.time_units           = All.UnitTime_in_s / All.HubbleParam;
  All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
  All.GrackleUnits.a_units              = 1.0; // units for the expansion factor
  
  // Second, create a chemistry object for parameters and rate data.
  if (set_default_chemistry_parameters() == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    exit(ENDRUNVAL);
  }
  // Set parameter values for chemistry.
  grackle_data.use_grackle            = 1;                   // chemistry on
  grackle_data.with_radiative_cooling = 1;                   // cooling on
  grackle_data.metal_cooling          = 1;                   // metal cooling on
  grackle_data.UVbackground           = 1;                   // UV background on
  grackle_data.grackle_data_file      = All.GrackleDataFile; // data file
  
  // optional stuff
  /*
  grackle_data.h2_on_dust                       = 0;
  grackle_data.cmb_temperature_floor            = 1;
  grackle_data.Gamma                            = GAMMA;
  grackle_data.three_body_rate                  = 0;
  grackle_data.cie_cooling                      = 0;
  grackle_data.h2_optical_depth_approximation   = 0;
  grackle_data.photoelectric_heating            = 0;
  grackle_data.photoelectric_heating_rate       = 8.5e-26;
  grackle_data.Compton_xray_heating             = 0;
  grackle_data.LWbackground_intensity           = 0;
  grackle_data.LWbackground_sawtooth_supression = 0;
  */

#ifndef GRACKLE_CHEMISTRY
  grackle_data.primordial_chemistry = 0;                     // fully tabulated cooling
#else
  grackle_data.primordial_chemistry = GRACKLE_CHEMISTRY;
#endif 

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  double a_value = 1.0;
  if(All.ComovingIntegrationOn)
      a_value = All.TimeBegin;
  
  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&All.GrackleUnits, a_value) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    exit(ENDRUNVAL);
  }
  
  if(ThisTask == 0)
    printf("GRACKLE INITIALIZED\n");
}

#endif  //GRACKLE