#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;
#ifdef _OPENMP
  int tid;
#endif
  if(ThisTask == 0)
    {
      printf("\nThis is GIZMO, version %s.\n", GIZMO_VERSION);
      printf("\nRunning on %d MPI tasks.\n", NTask);
#ifdef _OPENMP
#pragma omp parallel private(tid)
      {
#pragma omp master
	printf("\nUsing %d OpenMP threads\n", omp_get_num_threads());

	tid = omp_get_thread_num();
	/*
	   printf("Hello from thread = %d\n", tid);
	 */
      }
#endif
      printf("\nCode was compiled with settings:\n\n");

      output_compile_time_options();

      printf("Size of particle structure       %d  [bytes]\n", (int) sizeof(struct particle_data));
      printf("\nSize of sph particle structure   %d  [bytes]\n", (int) sizeof(struct sph_particle_data));
#ifdef DETACH_BLACK_HOLES
      printf("\nSize of BH particle structure    %d  [bytes]\n\n", (int) sizeof(struct bh_particle_data));
#endif

    }

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

  mymalloc_init();

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

  set_units();

#ifdef SUB_TURB_DRIVING
  sub_turb_read_table();
#endif

#ifdef COOLING
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  InitCool();
#endif

#ifdef GALSF_EFFECTIVE_EQS
  init_clouds();
#endif


#ifdef PERIODIC
  ewald_init();
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
  inverse_boxSize = 1. / boxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
  inverse_boxSize_X = 1. / boxSize_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
  inverse_boxSize_Y = 1. / boxSize_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
  inverse_boxSize_Z = 1. / boxSize_Z;
#endif
#endif

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

  set_random_numbers();

#ifdef PMGRID
#ifndef ADAPTIVE_GRAVSOFT_FORALL
  if(RestartFlag != 3 && RestartFlag != 4)
#endif
    long_range_init();
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
  long_range_init();
#endif
#endif

#ifdef SUBFIND
  GrNr = -1;
#endif

#ifdef EOS_DEGENERATE
  eos_init(All.EosTable, All.EosSpecies);
#endif

#ifdef NUCLEAR_NETWORK
  network_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses,
	       All.NetworkWeakrates, &All.nd);
  network_workspace_init(&All.nd, &All.nw);
#endif

#if defined (VS_TURB) || defined (AB_TURB) || defined(TURB_DRIVING)
  init_turb();
#endif

#ifdef SIDM
    AllocateInteractionTable(INTERACTION_TABLE_LENGTH, PARTICLE_MAX_INTERACTIONS + 1);
    init_geofactor_table();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 3 || RestartFlag == 4 || RestartFlag == 5 || RestartFlag == 6)
    {

mpi_printf("needs to be in init if ...\n");

      init();			/* ... read in initial model */

mpi_printf("... this caught the error \n");
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets
				   all variables in the struct `All'.
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.
				 */

      set_random_numbers();

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MinGasHsmlFractional = all.MinGasHsmlFractional;
      All.MinGasTemp = all.MinGasTemp;
        
        /* allow softenings to be modified during the run */
        if(All.ComovingIntegrationOn)
        {
        All.SofteningGasMaxPhys = all.SofteningGasMaxPhys;
        All.SofteningHaloMaxPhys = all.SofteningHaloMaxPhys;
        All.SofteningDiskMaxPhys = all.SofteningDiskMaxPhys;
        All.SofteningBulgeMaxPhys = all.SofteningBulgeMaxPhys;
        All.SofteningStarsMaxPhys = all.SofteningStarsMaxPhys;
        All.SofteningBndryMaxPhys = all.SofteningBndryMaxPhys;
        }
        All.SofteningGas = all.SofteningGas;
        All.SofteningHalo = all.SofteningHalo;
        All.SofteningDisk = all.SofteningDisk;
        All.SofteningBulge = all.SofteningBulge;
        All.SofteningStars = all.SofteningStars;
        All.SofteningBndry = all.SofteningBndry;
#ifdef SINKS
        All.SinkHsml = all.SinkHsml;
#endif

      All.MaxHsml = all.MaxHsml;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);
      memcpy(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength);

#ifdef GALSF
      All.CritPhysDensity = all.CritPhysDensity;
      All.MaxSfrTimescale = all.MaxSfrTimescale;
#endif
        
#ifdef SIDM
        All.SIDMSmoothingFactor = all.SIDMSmoothingFactor;
#endif


#ifdef AV_CD10_VISCOSITY_SWITCH
      All.ArtBulkViscConst = all.ArtBulkViscConst;
      All.ViscosityAMin = all.ViscosityAMin;
      All.ViscosityAMax = all.ViscosityAMax;
#endif
#ifdef TURB_DIFFUSION
      All.TurbDiffusion_Coefficient = all.TurbDiffusion_Coefficient;
#endif
#ifdef AV_ARTIFICIAL_CONDUCTIVITY
      All.ArtCondConstant = all.ArtCondConstant;
#endif

#if defined(MAGNETIC_DISSIPATION)
      All.ArtMagDispConst = all.ArtMagDispConst;
#endif

#ifdef DIVBCLEANING_DEDNER
      All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
      All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
      All.DivBcleanQ = all.DivBcleanQ;
#endif
#ifdef BLACK_HOLES
      All.BlackHoleMaxAccretionRadius = all.BlackHoleMaxAccretionRadius;
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORALL
      /* Allow the tolerance over the number of neighbours to vary during the run:
       * If it was initially set to a very strict value, convergence in ngb-iteration may at some point fail */
      All.AGS_MaxNumNgbDeviation = all.AGS_MaxNumNgbDeviation;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.TimebinFile, all.TimebinFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

#ifdef EOS_DEGENERATE
      strcpy(All.EosTable, all.EosTable);
      strcpy(All.EosSpecies, all.EosSpecies);
#endif

#ifdef RELAXOBJECT
      All.RelaxBaseFac = all.RelaxBaseFac;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(All.NetworkRates, all.NetworkRates);
      strcpy(All.NetworkPartFunc, all.NetworkPartFunc);
      strcpy(All.NetworkMasses, all.NetworkMasses);
      strcpy(All.NetworkWeakrates, all.NetworkWeakrates);
      All.nd = all.nd;
      All.nw = all.nw;
      All.NetworkTempThreshold = all.NetworkTempThreshold;
#endif

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);
    }

  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  open_outputfiles();

#ifdef PMGRID
  long_range_init_regionsize();
#endif
  reconstruct_timebins();

#ifdef COSMIC_RAYS
  int CRpop;

  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
    CR_initialize_beta_tabs(All.CR_Alpha[CRpop], CRpop);
  CR_Tab_Initialize();
#ifdef COSMIC_RAY_TEST
  CR_test_routine();
#endif

#endif

#ifdef TWODIMS
  int i;

  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;
      P[i].Vel[2] = 0;

      P[i].GravAccel[2] = 0;

      if(P[i].Type == 0)
	{
	  SphP[i].VelPred[2] = 0;
	  SphP[i].HydroAccel[2] = 0;
	}
    }
#endif


#ifdef RADTRANSFER
  if(RestartFlag == 0)
    radtransfer_set_simple_inits();
  
  rt_get_sigma();

#ifdef RT_MULTI_FREQUENCY
#if defined(EDDINGTON_TENSOR_STARS) || defined(EDDINGTON_TENSOR_SFR)
    rt_get_lum_stars();
#endif
#endif

#endif

  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
  else if(RestartFlag == 1)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the cgs-system
 */
void set_units(void)
{
  double meanweight;

#if defined(CONDUCTION_EXPLICIT)
  double coulomb_log;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;
#ifdef TIMEDEPGRAV
  All.Gini = All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

#ifdef DISTORTIONTENSORPS
  /* 5.609589206e23 is the factor to convert from g to GeV/c^2, the rest comes from All.UnitDensity_in_cgs */
  All.UnitDensity_in_Gev_per_cm3 = 5.609589206e23 / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g;
#endif
  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
#ifdef DISTORTIONTENSORPS
      printf("Annihilation radiation units:\n");
      printf("UnitDensity_in_Gev_per_cm3 = %g\n", All.UnitDensity_in_Gev_per_cm3);
#endif

      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#if defined(GALSF)
  set_units_sfr();
#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)


#if defined(CONDUCTION_EXPLICIT)

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  /* assuming full ionization */

  coulomb_log = 37.8;
  /* accordin1g to Sarazin's book */

  All.ConductionCoeff *= (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif
#endif /* CONDUCTION */

#if defined(CR_DIFFUSION)
  if(All.CR_DiffusionDensZero == 0.0)
    {
      /* Set reference density for CR Diffusion to rhocrit at z=0 */
      All.CR_DiffusionDensZero = 3.0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
    }

  if(All.CR_DiffusionInternalEnergyZero == 0.0)
    {
      All.CR_DiffusionInternalEnergyZero = 1.0e4;
    }

  /* Set reference entropic function to correspond to
     Reference Temperature @ ReferenceDensity */
  if(ThisTask == 0)
    {
      printf("CR Diffusion: T0 = %g\n", All.CR_DiffusionInternalEnergyZero);
    }

  /* convert Temperature value in Kelvin to thermal energy per unit mass in internal units */
  All.CR_DiffusionInternalEnergyZero *=
    BOLTZMANN / (4.0 * PROTONMASS / (3.0 * HYDROGEN_MASSFRAC + 1.0)) *
    All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(ThisTask == 0)
    {
      printf("CR Diffusion: Rho0 = %g -- A0 = %g\n", All.CR_DiffusionDensZero, All.CR_DiffusionInternalEnergyZero);
    }

#endif /* CR_DIFFUSION */

}




/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  if(ThisTask == 0)
    mkdir(All.OutputDir, 02755);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_details", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(buf, "%sblackhole_details/blackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimebinFile);
  if(!(FdTimebin = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1],
	  CPU_SymbolImbalance[CPU_TREEWALK1]);
  fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2],
	  CPU_SymbolImbalance[CPU_TREEWALK2]);
  fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1],
	  CPU_SymbolImbalance[CPU_TREEWAIT1]);
  fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2],
	  CPU_SymbolImbalance[CPU_TREEWAIT2]);
  fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND],
	  CPU_SymbolImbalance[CPU_TREESEND]);
  fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV],
	  CPU_SymbolImbalance[CPU_TREERECV]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC],
	  CPU_SymbolImbalance[CPU_TREEMISC]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE],
	  CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
  fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT],
	  CPU_SymbolImbalance[CPU_DENSWAIT]);
  fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC],
	  CPU_SymbolImbalance[CPU_DENSMISC]);
  fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE],
	  CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
  fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT],
	  CPU_SymbolImbalance[CPU_HYDWAIT]);
  fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC],
	  CPU_SymbolImbalance[CPU_HYDMISC]);
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLES],
	  CPU_SymbolImbalance[CPU_BLACKHOLES]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
  fprintf(FdBalance, "\n");

#ifdef SCFPOTENTIAL
  sprintf(buf, "%s%s", All.OutputDir, "scf_coeff.txt");
  if(!(FdSCF = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef GALSF
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

    
#ifdef GALSF_FB_GASRETURN
    sprintf(buf, "%s%s", All.OutputDir, "GasReturn.txt");
    if(!(FdGasReturn = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
#endif
#ifdef GALSF_FB_RPWIND_LOCAL
    sprintf(buf, "%s%s", All.OutputDir, "MomWinds.txt");
    if(!(FdMomWinds = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
#endif
#ifdef GALSF_FB_HII_HEATING
    sprintf(buf, "%s%s", All.OutputDir, "HIIheating.txt");
    if(!(FdHIIHeating = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }
#endif
#ifdef GALSF_FB_SNE_HEATING
    sprintf(buf, "%s%s", All.OutputDir, "SNeIIheating.txt");
    if(!(FdSneIIHeating = fopen(buf, mode)))
    {
        printf("error in opening file '%s'\n", buf);
        endrun(1);
    }  
#endif
    
    
#ifdef RADTRANSFER
  sprintf(buf, "%s%s", All.OutputDir, "radtransfer.txt");
  if(!(FdRad = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef RT_MULTI_FREQUENCY
  sprintf(buf, "%s%s", All.OutputDir, "star_lum.txt");
  if(!(FdStar = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  sprintf(buf, "%s%s", All.OutputDir, "turb.txt");
  if(!(FdTurb = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif


#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "dH dG ");
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
    }
#endif

}




/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdTimebin);
  fclose(FdBalance);

#ifdef SCFPOTENTIAL
  fclose(FdSCF);
#endif

#ifdef GALSF
  fclose(FdSfr);
#endif

#ifdef GALSF_FB_GASRETURN
    fclose(FdGasReturn);
#endif
#ifdef GALSF_FB_RPWIND_LOCAL
    fclose(FdMomWinds);
#endif
#ifdef GALSF_FB_HII_HEATING
    fclose(FdHIIHeating);
#endif
#ifdef GALSF_FB_SNE_HEATING
    fclose(FdSneIIHeating);
#endif
    
#ifdef RADTRANSFER
  fclose(FdRad);
#ifdef RT_MULTI_FREQUENCY
  fclose(FdStar);
#endif
#endif

#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif

#ifdef DARKENERGY
  fclose(FdDE);
#endif

}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

  All.StarformationOn = 0;	/* defaults */

#ifdef COSMIC_RAYS
  char tempBuf[20];

  int CRpop, x;

  double tempAlpha;
#endif


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimebinFile");
      addr[nt] = All.TimebinFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

#ifdef GRAIN_FLUID
        strcpy(tag[nt],"Grain_Density_CGS");
        addr[nt] = &All.Grain_Internal_Density;
        id[nt++] = REAL;
        
        strcpy(tag[nt],"Grain_Size_Init");
        addr[nt] = &All.InitGrainSize;
        id[nt++] = REAL;
#endif
        
#if defined(GALSF_FB_GASRETURN) || defined(GALSF_FB_SNE_HEATING)
        strcpy(tag[nt],"GasReturnFraction");
        addr[nt] = &All.GasReturnFraction;
        id[nt++] = REAL;
#endif
        
#ifdef GALSF_FB_GASRETURN
        strcpy(tag[nt],"AGBGasTemp");
        addr[nt] = &All.AGBGasEnergy;
        id[nt++] = REAL;
#endif
        
#ifdef BH_BAL_WINDS
        strcpy(tag[nt],"BAL_f_accretion");
        addr[nt] = &All.BAL_f_accretion;
        id[nt++] = REAL;
        
        strcpy(tag[nt],"BAL_v_outflow");
        addr[nt] = &All.BAL_v_outflow;
        id[nt++] = REAL;
#endif
        
#ifdef BH_PHOTONMOMENTUM
        strcpy(tag[nt],"BH_FluxMomentumFactor");
        addr[nt] = &All.BH_FluxMomentumFactor;
        id[nt++] = REAL;
#endif

#if defined(GALSF_FB_GASRETURN) || defined(GALSF_FB_RPWIND_LOCAL) || defined(GALSF_FB_HII_HEATING) || defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_RT_PHOTONMOMENTUM)
        strcpy(tag[nt],"InitMetallicity");
        addr[nt] = &All.InitMetallicityinSolar;
        id[nt++] = REAL;
        
        strcpy(tag[nt],"InitStellarAge");
        addr[nt] = &All.InitStellarAgeinGyr;
        id[nt++] = REAL;
#endif
        
#ifdef GALSF_FB_SNE_HEATING
        strcpy(tag[nt], "SNeIIEnergyFrac");
        addr[nt] = &All.SNeIIEnergyFrac;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "SNeIIBW_Radius_Factor");
        addr[nt] = &All.SNeIIBW_Radius_Factor;
        id[nt++] = REAL;
#endif
        
#ifdef GALSF_FB_HII_HEATING
        strcpy(tag[nt], "HIIRegion_fLum_Coupled");
        addr[nt] = &All.HIIRegion_fLum_Coupled;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "HIIRegion_Temp");
        addr[nt] = &All.HIIRegion_Temp;
        id[nt++] = REAL;
#endif

#ifdef GALSF_FB_RT_PHOTONMOMENTUM
        strcpy(tag[nt], "PhotonMomentum_Coupled_Fraction");
        addr[nt] = &All.PhotonMomentum_Coupled_Fraction;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "PhotonMomentum_fUV");
        addr[nt] = &All.PhotonMomentum_fUV;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "PhotonMomentum_fOPT");
        addr[nt] = &All.PhotonMomentum_fOPT;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "PhotonMomentum_fIR");
        addr[nt] = &All.PhotonMomentum_fIR;
        id[nt++] = REAL;
#endif

        
        
#ifdef SIDM
        strcpy(tag[nt], "InteractionCrossSection");
        addr[nt] = &All.InteractionCrossSection;
        id[nt++] = REAL;

        strcpy(tag[nt], "SIDMSmoothingFactor");
        addr[nt] = &All.SIDMSmoothingFactor;
        id[nt++] = REAL;
#endif


#ifdef SUBFIND
      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = REAL;

#if defined(ISOTHERM_EQS) || defined(VS_TURB) || defined (AB_TURB)
      strcpy(tag[nt], "IsoSoundSpeed");
      addr[nt] = &All.IsoSoundSpeed;
      id[nt++] = REAL;
#endif


      strcpy(tag[nt], "MaxHsml");
      addr[nt] = &All.MaxHsml;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = REAL;

#ifdef AV_ARTIFICIAL_CONDUCTIVITY
      strcpy(tag[nt], "ArtCondConstant");
      addr[nt] = &All.ArtCondConstant;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;


#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;
#endif

#if (defined(AB_TURB) || defined(VS_TURB)) && defined(POWERSPEC_GRID)
      strcpy(tag[nt], "TimeBetTurbSpectrum");
      addr[nt] = &All.TimeBetTurbSpectrum;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

#ifdef DISTORTIONTENSORPS
      strcpy(tag[nt], "TidalCorrection");
      addr[nt] = &All.TidalCorrection;
      id[nt++] = REAL;

      strcpy(tag[nt], "DM_velocity_dispersion");
      addr[nt] = &All.DM_velocity_dispersion;
      id[nt++] = REAL;
#endif
#ifdef SCALARFIELD
      strcpy(tag[nt], "ScalarBeta");
      addr[nt] = &All.ScalarBeta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ScalarScreeningLength");
      addr[nt] = &All.ScalarScreeningLength;
      id[nt++] = REAL;
#endif

#ifdef OUTPUTLINEOFSIGHT
      strcpy(tag[nt], "TimeFirstLineOfSight");
      addr[nt] = &All.TimeFirstLineOfSight;
      id[nt++] = REAL;
#endif


#if defined(BUBBLES) || defined(MULTI_BUBBLES)
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleTimeInterval");
      addr[nt] = &All.BubbleTimeInterval;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "FirstBubbleRedshift");
      addr[nt] = &All.FirstBubbleRedshift;
      id[nt++] = REAL;
#endif

#ifdef MULTI_BUBBLES
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "ClusterMass200");
      addr[nt] = &All.ClusterMass200;
      id[nt++] = REAL;

      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;

#endif

#ifdef BH_BUBBLES
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadioTriggeringFactor");
      addr[nt] = &All.BlackHoleRadioTriggeringFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "DefaultICMDensity");
      addr[nt] = &All.DefaultICMDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioFeedbackFactor");
      addr[nt] = &All.RadioFeedbackFactor;
      id[nt++] = REAL;
#ifdef UNIFIED_FEEDBACK
      strcpy(tag[nt], "RadioThreshold");
      addr[nt] = &All.RadioThreshold;
      id[nt++] = REAL;
#endif
#endif

#ifdef COSMIC_RAYS
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  sprintf(tempBuf, "CR_SpectralIndex_%i", CRpop);
	  strcpy(tag[nt], tempBuf);
	  addr[nt] = &All.CR_Alpha[CRpop];
	  id[nt++] = REAL;
	}

      /* sort the All.CR_Alpha[CRpop] array in assending order */
      if(NUMCRPOP > 1)
	for(x = 0; x < NUMCRPOP - 1; x++)
	  for(CRpop = 0; CRpop < NUMCRPOP - x - 1; CRpop++)
	    if(All.CR_Alpha[CRpop] > All.CR_Alpha[CRpop + 1])
	      {
		tempAlpha = All.CR_Alpha[CRpop];
		All.CR_Alpha[CRpop] = All.CR_Alpha[CRpop + 1];
		All.CR_Alpha[CRpop + 1] = tempAlpha;
	      }

      strcpy(tag[nt], "CR_SupernovaEfficiency");
      addr[nt] = &All.CR_SNEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaSpectralIndex");
      addr[nt] = &All.CR_SNAlpha;
      id[nt++] = REAL;

#if defined(CR_DIFFUSION)
      strcpy(tag[nt], "CR_DiffusionCoeff");
      addr[nt] = &All.CR_DiffusionCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionDensityScaling");
      addr[nt] = &All.CR_DiffusionDensScaling;
      id[nt++] = REAL;

      /* CR Diffusion scaling: reference density rho_0.
         If value is 0, then rho_crit @ z=0 is used. */
      strcpy(tag[nt], "CR_DiffusionReferenceDensity");
      addr[nt] = &All.CR_DiffusionDensZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTemperatureScaling");
      addr[nt] = &All.CR_DiffusionInternalEnergyScaling;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionReferenceTemperature");
      addr[nt] = &All.CR_DiffusionInternalEnergyZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionMaxSizeTimestep");
      addr[nt] = &All.CR_DiffusionMaxSizeTimestep;
      id[nt++] = REAL;
#endif /* CR_DIFFUSION */

#ifdef CR_SHOCK
      strcpy(tag[nt], "CR_ShockEfficiency");
      addr[nt] = &All.CR_ShockEfficiency;
      id[nt++] = REAL;
#if ( CR_SHOCK == 1 )		/* Constant Spectral Index method */
      strcpy(tag[nt], "CR_ShockSpectralIndex");
      addr[nt] = &All.CR_ShockAlpha;
      id[nt++] = REAL;
#else /* Mach-Number - Dependent method */
      strcpy(tag[nt], "CR_ShockCutoffFac");
      addr[nt] = &All.CR_ShockCutoff;
      id[nt++] = REAL;
#endif
#endif /* CR_SHOCK */

#ifdef FIX_QINJ
      strcpy(tag[nt], "Shock_Fix_Qinj");
      addr[nt] = &All.Shock_Fix_Qinj;
      id[nt++] = REAL;
#endif

#ifdef CR_BUBBLES
      strcpy(tag[nt], "CR_AGNEff");
      addr[nt] = &All.CR_AGNEff;
      id[nt++] = REAL;
#endif
#endif /* COSMIC_RAYS */


#ifdef MACHNUM
      strcpy(tag[nt], "Shock_LengthScale");
      addr[nt] = &All.Shock_Length;
      id[nt++] = REAL;

      strcpy(tag[nt], "Shock_DeltaDecayTimeMax");
      addr[nt] = &All.Shock_DeltaDecayTimeMax;
      id[nt++] = REAL;
#endif

#if defined(BLACK_HOLES) || defined(GALSF_SUBGRID_VARIABLEVELOCITY)
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = REAL;
#endif

#ifdef BLACK_HOLES

#ifdef DETACH_BLACK_HOLES
      strcpy(tag[nt], "BHFormationFactor");
      addr[nt] = &All.BHfactor;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;


      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleNgbFactor");
      addr[nt] = &All.BlackHoleNgbFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
      addr[nt] = &All.BlackHoleMaxAccretionRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
      addr[nt] = &All.BlackHoleRadiativeEfficiency;
      id[nt++] = REAL;

#ifdef FOF
      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

#endif /* BLACK_HOLES */

#ifdef GALSF
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;
        
#ifdef GALSF_EFFECTIVE_EQS
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;
        
      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;
#endif
        
#ifdef GALSF_FB_RPWIND_LOCAL
        strcpy(tag[nt], "WindMomentumLoading");
        addr[nt] = &All.WindMomentumLoading;
        id[nt++] = REAL;
        
        strcpy(tag[nt], "WindInitialVelocityBoost");
        addr[nt] = &All.WindInitialVelocityBoost;
        id[nt++] = REAL;
#endif
        
#ifdef GALSF_SUBGRID_WINDS
      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;

#ifdef GALSF_SUBGRID_VARIABLEVELOCITY
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "VariableWindSpecMomentum");
      addr[nt] = &All.VariableWindSpecMomentum;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloConcentrationNorm");
      addr[nt] = &All.HaloConcentrationNorm;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloConcentrationSlope");
      addr[nt] = &All.HaloConcentrationSlope;
      id[nt++] = REAL;
#endif // GALSF_SUBGRID_VARIABLEVELOCITY
#endif // GALSF_SUBGRID_WINDS

#endif

#ifdef GALSF_EFFECTIVE_EQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;
#endif
#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef AV_CD10_VISCOSITY_SWITCH
      strcpy(tag[nt], "ViscosityAMin");
      addr[nt] = &All.ViscosityAMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityAMax");
      addr[nt] = &All.ViscosityAMax;
      id[nt++] = REAL;
#endif
#ifdef TURB_DIFFUSION
      strcpy(tag[nt], "TurbDiffusionCoefficient");
      addr[nt] = &All.TurbDiffusion_Coefficient;
      id[nt++] = REAL;
#endif

#if defined(MAGNETIC_DISSIPATION)
      strcpy(tag[nt], "ArtificialMagneticDissipationConstant");
      addr[nt] = &All.ArtMagDispConst;
      id[nt++] = REAL;
#endif

#ifdef DIVBCLEANING_DEDNER
      strcpy(tag[nt], "DivBcleaningParabolicSigma");
      addr[nt] = &All.DivBcleanParabolicSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
      addr[nt] = &All.DivBcleanHyperbolicSigma;
      id[nt++] = REAL;
#endif

#ifdef MAGNETIC
#ifdef BINISET
      strcpy(tag[nt], "BiniX");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif
#endif /* MAGNETIC */

#ifdef EOS_DEGENERATE
      strcpy(tag[nt], "EosTable");
      addr[nt] = All.EosTable;
      id[nt++] = STRING;

      strcpy(tag[nt], "EosSpecies");
      addr[nt] = All.EosSpecies;
      id[nt++] = STRING;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(tag[nt], "NetworkRates");
      addr[nt] = All.NetworkRates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkPartFunc");
      addr[nt] = All.NetworkPartFunc;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkMasses");
      addr[nt] = All.NetworkMasses;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkWeakrates");
      addr[nt] = All.NetworkWeakrates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkTempThreshold");
      addr[nt] = &All.NetworkTempThreshold;
      id[nt++] = REAL;
#endif

#ifdef RELAXOBJECT
      strcpy(tag[nt], "RelaxBaseFac");
      addr[nt] = &All.RelaxBaseFac;
      id[nt++] = REAL;
#endif

#ifdef RADTRANSFER
#if defined(EDDINGTON_TENSOR_STARS) && defined(RT_TEST_SST)
      strcpy(tag[nt], "IonizingLumPerSolarMass");
      addr[nt] = &All.IonizingLumPerSolarMass;
      id[nt++] = REAL;
#ifdef RT_POPIII
      strcpy(tag[nt], "IonizingLumPerSolarMassPopIII");
      addr[nt] = &All.IonizingLumPerSolarMassPopIII;
      id[nt++] = REAL;
#endif
#endif

#ifdef RT_MULTI_FREQUENCY
      strcpy(tag[nt], "star_Teff");
      addr[nt] = &All.star_Teff;
      id[nt++] = REAL;
#ifdef RT_POPIII
      strcpy(tag[nt], "star_Teff_popIII");
      addr[nt] = &All.star_Teff_popIII;
      id[nt++] = REAL;
#endif
#endif

#if defined(EDDINGTON_TENSOR_SFR) && defined(RT_TEST_SST)
      strcpy(tag[nt], "IonizingLumPerSFR");
      addr[nt] = &All.IonizingLumPerSFR;
      id[nt++] = REAL;
#endif

#endif //end radtransfer

#ifdef SINKS
      strcpy(tag[nt], "SinkHsml");
      addr[nt] = &All.SinkHsml;
      id[nt++] = REAL;

      strcpy(tag[nt], "SinkDensThresh");
      addr[nt] = &All.SinkDensThresh;
      id[nt++] = REAL;
#endif

#ifdef BP_REAL_CRs
      strcpy(tag[nt], "MinCREnergy");
      addr[nt] = &All.ecr_min;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxCREnergy");
      addr[nt] = &All.ecr_max;
      id[nt++] = REAL;
#ifdef BP_SEED_CRs
      strcpy(tag[nt], "CRpInitSlope");
      addr[nt] = &All.pSlope_init;
      id[nt++] = REAL;

      strcpy(tag[nt], "CReInitSlope");
      addr[nt] = &All.eSlope_init;
      id[nt++] = REAL;
#endif
#ifdef BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY
      strcpy(tag[nt], "CRsArtCondConstant");
      addr[nt] = &All.CRsArtCondConstant;
      id[nt++] = REAL;
#endif
#endif


#ifdef ADAPTIVE_GRAVSOFT_FORALL
      strcpy(tag[nt], "AGS_DesNumNgb");
      addr[nt] = &All.AGS_DesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "AGS_MaxNumNgbDeviation");
      addr[nt] = &All.AGS_MaxNumNgbDeviation;
      id[nt++] = REAL;
#endif

#if defined(AB_TURB) || defined(TURB_DRIVING)
      strcpy(tag[nt], "ST_decay");
      addr[nt] = &All.StDecay;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_energy");
      addr[nt] = &All.StEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_DtFreq");
      addr[nt] = &All.StDtFreq;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_Kmin");
      addr[nt] = &All.StKmin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_Kmax");
      addr[nt] = &All.StKmax;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_SolWeight");
      addr[nt] = &All.StSolWeight;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_AmplFac");
      addr[nt] = &All.StAmplFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_SpectForm");
      addr[nt] = &All.StSpectForm;
      id[nt++] = INT;
#endif
#ifdef AB_TURB
      strcpy(tag[nt], "ST_Seed");
      addr[nt] = &All.StSeed;
      id[nt++] = REAL;
#endif
#ifdef TURB_DRIVING
      strcpy(tag[nt], "ST_Seed");
      addr[nt] = &All.StSeed;
      id[nt++] = INT;
#endif

#ifdef ADJ_BOX_POWERSPEC
      strcpy(tag[nt], "BoxWidth");
      addr[nt] = &All.BoxWidth;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_x");
      addr[nt] = &All.BoxCenter_x;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_y");
      addr[nt] = &All.BoxCenter_y;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_z");
      addr[nt] = &All.BoxCenter_z;
      id[nt++] = REAL;

      strcpy(tag[nt], "TransformSize");
      addr[nt] = &All.FourierGrid;
      id[nt++] = INT;
#endif


      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      printf("Obtaining parameters from file '%s':\n", fname);
	      while(!feof(fd))
		{

		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  fprintf(stdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  fprintf(stdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  fprintf(stdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
#ifdef ALLOWEXTRAPARAMS
		      fprintf(stdout, "WARNING from file %s:   Tag '%s' ignored !\n", fname, buf1);
#else
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
#endif
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      printf("\n");

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}

      /* Counts number of CR populations in parameter file */
      /*
         #ifdef COSMIC_RAYS
         CRpop=0;
         while ((All.CR_Alpha[CRpop] != 0.0 ) && ( CRpop < MaxNumCRpop ))
         CRpop++;
         NUMCRPOP = CRpop;
         #else
         NUMCRPOP = 1;
         #endif
       */
#ifdef COSMIC_RAYS
      printf(" NUMCRPOP = %i \n", NUMCRPOP);
#endif

#ifdef BP_REAL_CRs
      printf(" BP_REAL_CRs = %i \n", BP_REAL_CRs);
#endif

      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#ifdef EDDINGTON_TENSOR_BH
#ifndef BLACK_HOLES
  if(ThisTask == 0)
    {
      printf("Code was compiled with EDDINGTON_TENSOR_BH, but not with BLACK_HOLES.\n");
      printf("Switch on BLACK_HOLES.\n");
    }
  endrun(0);
#endif
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef GALSF
  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif



#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif


#ifdef BH_BUBBLES
#ifndef BLACK_HOLES
  if(ThisTask == 0)
    {
      printf("Code was compiled with BH_BUBBLES, but not with BLACK_HOLES.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES) || defined(EBUB_PROPTO_BHAR)
  if(ThisTask == 0)
    {
      printf
	("If the code is compiled with BH_BUBBLES, then BUBBLES, MULTI_BUBBLES or EBUB_PROPTO_BHAR options cannot be used.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif



#ifdef OMP_NUM_THREADS
#ifdef _OPENMP
  if(ThisTask == 0)
    printf("OMP_NUM_THREADS is incompatible with enabling OpenMP in the compiler options \n");
  endrun(0);
#endif
#endif

#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS


#ifdef COSMIC_RAYS
  if(ThisTask == 0)
    {
      printf("CR SN Efficiency: %g\n", All.CR_SNEff);
    }
#endif
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
	break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
	flag = 1;

      if(count == 1 || count == 2)
	{
	  if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
	    {
	      if(ThisTask == 0)
		printf("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n",
		       (int) MAXLEN_OUTPUTLIST);
	      endrun(13);
	    }

	  All.OutputListFlag[All.OutputListLength] = flag;
	  All.OutputListLength++;
	}
    }

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif
#ifdef CR_DIFFUSION
      All.CR_Diffusion_Ti_begstep /= 2;
      All.CR_Diffusion_Ti_endstep /= 2;
#endif

#if defined(AB_TURB) || defined(TURB_DRIVING)
      StTPrev /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_current /= 2;

	  if(P[i].TimeBin > 0)
	    {
	      P[i].TimeBin--;
	      if(P[i].TimeBin <= 0)
		{
		  printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
		  endrun(8765);
		}
	    }
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
