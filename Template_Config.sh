#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
#  Enable/Disable compile-time options as needed: this is where you determine how the code will act
#  From the list below, please activate/deactivate the
#       options that apply to your run. If you modify any of these options,
#       make sure that you recompile the whole code by typing "make clean; make".
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (to add new modules and clean
#   up the naming conventions and changed many of them to match the new GIZMO conventions)
#
####################################################################################################



####################################################################################################
#--------------------------------------- Boundary Conditions & Dimensions
####################################################################################################
PERIODIC                        # Use this if periodic boundaries are needed (otherwise open boundaries are assumed)
#BND_PARTICLES                  # particles with ID=0 are forced in place (their accelerations are set =0):
                                # use for special boundary conditions where these particles represent fixed "walls"
#LONG_X=140                     # modify box dimensions (non-square periodic box): multiply X (PERIODIC and NOGRAVITY required)
#LONG_Y=1                       # modify box dimensions (non-square periodic box): multiply Y
#LONG_Z=1                       # modify box dimensions (non-square periodic box): multiply Z
#REFLECT_BND_X                  # make the x-boundary reflecting (assumes a box 0<x<1, unless PERIODIC is set)
#REFLECT_BND_Y                  # make the y-boundary reflecting (assumes a box 0<y<1, unless PERIODIC is set)
#REFLECT_BND_Z                  # make the z-boundary reflecting (assumes a box 0<z<1, unless PERIODIC is set)
#SHEARING_BOX=1                 # shearing box boundaries: 1=r-z sheet (r,z,phi coordinates), 2=r-phi sheet (r,phi,z), 3=r-phi-z box, 4=as 3, with vertical gravity
#SHEARING_BOX_Q=(3./2.)         # shearing box q=-dlnOmega/dlnr; will default to 3/2 (Keplerian) if not set
#ONEDIM                         # Switch for 1D test problems: code only follows the x-line. requires NOGRAVITY, and all y=z=0
#TWODIMS                        # Switch for 2D test problems: code only follows the xy-plane. requires NOGRAVITY, and all z=0.
####################################################################################################



####################################################################################################
#--------------------------------------- Hydro solver method
####################################################################################################
HYDRO_MESHLESS_FINITE_MASS      # Lagrangian (constant-mass) finite-volume Godunov method
#HYDRO_MESHLESS_FINITE_VOLUME   # Moving (quasi-Lagrangian) finite-volume Godunov method
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- SPH methods:
#SPHEQ_DENSITY_INDEPENDENT_SPH  # force SPH to use the 'pressure-sph' formulation ("modern" SPH)
#SPHEQ_TRADITIONAL_SPH          # force SPH to use the 'density-sph' (GADGET-2 & GASOLINE SPH)
#--------------------------------------- SPH artificial diffusion options (use with SPH; not relevant for Godunov/Mesh modes)
#SPHAV_DISABLE_CD10_VISCOSITY   # Disable Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks); just use Balsara switch
#SPHAV_DISABLE_PM_CONDUCTIVITY  # Disable mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches)
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Kernel Options
#KERNEL_QUARTIC                 # Implementation of the Morris 1996 quartic spline kernel
#KERNEL_QUINTIC                 # Implementation of the Morris 1996 quintic spline kernel, requires more neighbours !
####################################################################################################



####################################################################################################
#--------------------------------------- Additional Fluid Physics
####################################################################################################
#--------------------------------- Polytropic Index of Gas (for an ideal gas law)
#GAMMA=(5.0/3.0)                # if not set and no other (more complex) EOS set, defaults to GAMMA=5/3
##-----------------------------------------------------------------------------------------------------
#--------------------------------- Magneto-Hydrodynamics
#MAGNETIC                       # master switch for MHD, regardless of which Hydro solver is used
#B_SET_IN_PARAMS                # set initial fields (Bx,By,Bz) in parameter file
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Conduction
#CONDUCTION                     # Thermal conduction solved *explicitly*: isotropic if MAGNETIC off, otherwise anisotropic
#CONDUCTION_SPITZER             # Spitzer conductivity accounting for saturation: otherwise conduction coefficient is constant
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Viscosity
#VISCOSITY                      # Navier-stokes equations solved *explicitly*: isotropic coefficients if MAGNETIC off, otherwise anisotropic
#VISCOSITY_BRAGINSKII           # Braginskii viscosity tensor for ideal MHD
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Radiative Cooling physics (mostly geared towards galactic/extragalactic cooling)
#COOLING                        # enables radiative cooling and heating: if GALSF, also external UV background read from file "TREECOOL"
#COOL_LOW_TEMPERATURES          # allow fine-structure and molecular cooling to ~10 K
#COOL_METAL_LINES_BY_SPECIES    # use full multi-species-dependent cooling tables (https://dl.dropbox.com/u/16659252/spcool_tables.tgz)
#GRACKLE                        # enable GRACKLE: cooling+chemistry package (requires COOLING above; https://grackle.readthedocs.org/en/latest/)
#GRACKLE_CHEMISTRY=1            # choose GRACKLE cooling chemistry: (0)=tabular, (1)=Atomic, (2)=(1)+H2+H2I+H2II, (3)=(2)+DI+DII+HD
#TRUELOVE_CRITERION_PRESSURE    # adds artificial pressure floor force Jeans length above resolution scale (not necessarily better in meshless methods!)
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Smagorinsky Turbulent Eddy Diffusion Model
#TURB_DIFF_ENERGY               # turbulent diffusion of internal energy (energy-conserving formulation)
#TURB_DIFF_VELOCITY             # turbulent diffusion of velocity (effective turbulent viscosity)
#TURB_DIFF_METALS               # turbulent diffusion of metals (passive scalars)
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#---------------------------------------- Aerodynamic Particles
#GRAIN_FLUID                    # two-fluid medium with aerodynamically-coupled grains (particle type 3 are grains); default is Stokes drag
#GRAIN_EPSTEIN=1                # uses the cross section for molecular hydrogen (times this number) to calculate Epstein drag
#GRAIN_BACKREACTION             # account for momentum of grains pushing back on gas (from drag terms)
#GRAIN_LORENTZFORCE             # charged grains feel Lorentz forces (requires MAGNETIC)
#GRAIN_COLLISIONS               # model collisions between grains (super-particles; so this is stochastic)
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#---------------------------------------- Cosmic Rays
#COSMIC_RAYS                    # two-fluid medium with CRs as an ultrarelativistic fluid: heating/cooling, anisotropic diffusion, streaming, injection by SNe
##-----------------------------------------------------------------------------------------------------
####################################################################################################



####################################################################################################
##------------------------ Gravity & Cosmological Integration Options ---------------------------------
####################################################################################################
#--------------------------------------- TreePM Options (recommended for cosmological sims)
#PMGRID=512                     # COSMO enable: resolution of particle-mesh grid
#PM_PLACEHIGHRESREGION=1+2+16   # COSMO enable: particle types to place high-res PMGRID around
#PM_HIRES_REGION_CLIPPING=1000  # for stability: clips particles that escape the hires region in zoom/isolated sims
#MULTIPLEDOMAINS=64             # Multi-Domain option for the top-tree level: iso=16,COSMO=64-128
##-----------------------------------------------------------------------------------------------------
#---------------------------------------- Adaptive Grav. Softening (including Lagrangian conservation terms!)
#ADAPTIVE_GRAVSOFT_FORGAS       # allows variable softening length (=Hsml) for gas particles
#ADAPTIVE_GRAVSOFT_FORALL       # enable adaptive gravitational softening lengths for all particle types
                                # (ADAPTIVE_GRAVSOFT_FORGAS should be disabled). the softening is set to the distance
                                # enclosing a neighbor number set in the parameter file. baryons search for other baryons,
                                # dm for dm, sidm for sidm, etc.
##-----------------------------------------------------------------------------------------------------
#NOGRAVITY                      # turn off self-gravity (compatible with analytic_gravity)
#GRAVITY_NOT_PERIODIC           # self-gravity is not periodic, even though the rest of the box is periodic
##-----------------------------------------------------------------------------------------------------
#ANALYTIC_GRAVITY               # Specific analytic gravitational force to use instead of/with self-gravity
                                #  (edit these to assign specific parameters desired in "gravity/analytic_gravity.h")
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Self-Interacting DM (Rocha et al. 2012)
#SIDM=2                         # Self-interacting particle types (specify the particle types which are self-interacting DM
                                # with a bit mask, as for PM_PLACEHIGHRESREGION above (see description)
                                # (previous "DMDISK_INTERACTIONS" is identical to setting SIDM=2+4)
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- SCF (potential expansion in basis functions)
#SCFPOTENTIAL                   # turn SCF on/off
#SCF_HYBRID=1                   # =1:tree:stars<->stars,DM<->DM,stars->DM/SCF:stars<-DM =2:tree:stars<->stars,stars->DM/SCF:stars<-DM,DM<->DM
#SCF_HQ_MASS=95.2401            # mass of halo of expansion basis
#SCF_HQ_A=29.7754               # scale length of halo of expansion basis
#SCF_NEFF=5000000               # effective particle number in halo
#SCF_NMAX=1                     # maximum n for expansion cut-off
#SCF_LMAX=1                     # maximum l for expansion cut-off
#SCF_SCALEFAC                   # readin scale factors for coefficients from file scf_scalefac.dat
##-----------------------------------------------------------------------------------------------------
#SCALARFIELD                    # Gravity is mediated by a long-range scalar field instead of DE or DM
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Dark energy (flags to allow complicated equations of state)
#DARKENERGY                     # enables Dark Energy with non-trivial equations-of-state (master switch)
#TIMEDEPDE                      # read w(z) from a DE file (pre-tabulated)
#EXTERNALHUBBLE                 # reads the hubble function from the DE file (pre-tabulated)
#TIMEDEPGRAV                    # rescales H and G according to DE model (pre-tabulated)
##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Fine-grained phase space structure analysis (M. Vogelsberger)
#DISTORTIONTENSORPS             # main switch: integrate phase-space distortion tensor
#OUTPUT_DISTORTIONTENSORPS      # write phase-space distortion tensor to snapshot
#OUTPUT_TIDALTENSORPS           # write configuration-space tidal tensor to snapshot
#OUTPUT_LAST_CAUSTIC            # write info on last passed caustic to snapshot
#GDE_TYPES=2+4+8+16+32          # track GDE for these types
#GDE_READIC                     # read initial sheet orientation/initial density/initial caustic count from ICs
#GDE_LEAN                       # lean version of GDE
####################################################################################################



####################################################################################################
#------------------ Galaxy formation / Star formation / Supermassive BH Models (with feedback)
####################################################################################################
# ---- basic/master switches ---- #
#GALSF                          # master switch for galactic star formation model: enables SF, stellar ages, metals, generations, etc.
#METALS                         # enable metallicities (with multiple species optional) for gas and stars
#SINKS                          # add sink particles (SF studies)
##GALSF_GENERATIONS=1           # the number of stars a gas particle may spawn (defaults to 1, set otherwise)
##-----------------------------------------------------------------------------------------------------
# ---- old sub-grid models (for large-volume simulations) ---- #
#GALSF_EFFECTIVE_EQS            # 'effective equation of state' model for the ISM and star formation
#GALSF_SUBGRID_WINDS            # sub-grid winds ('kicks' as in Oppenheimer+Dave,Springel+Hernquist,Boothe+Schaye,etc)
#GALSF_SUBGRID_VARIABLEVELOCITY # winds with velocity scaling based on halo properties (Oppenheimer+Dave); req.GALSF_SUBGRID_WINDS
#GALSF_SUBGRID_DMDISPERSION     # wind velocity scaling based on MV 13 paper, as used in Illustris. req.GALSF_SUBGRID_WINDS
#GALSF_WINDS_ISOTROPIC          # forces winds to have a random orientation (works with both subgrid+explicit winds)
#GALSF_WINDS_POLAR              # forces winds to have polar orientation (works for sub-grid winds)
#GALSF_TURNOFF_COOLING_WINDS    # turn off cooling for SNe-heated particles (as Stinson+ GASOLINE model; requires GALSF_FB_SNE_HEATING)
#GALSF_GASOLINE_RADHEATING      # heat gas with luminosity from young stars (as Stinson+ 2013 GASOLINE model; requires GALSF_FB_SNE_HEATING)
#
##-----------------------------------------------------------------------------------------------------
## ------ PFH physical models for star formation and feedback ------ ##
# ---- star formation law ---- #
#GALSF_SFR_MOLECULAR_CRITERION	# estimates molecular fraction in SF-ing gas, only SF from that is allowed
#GALSF_SFR_VIRIAL_SF_CRITERION	# only allow star formation in virialized sub-regions (alpha<1)
#GALSF_SFR_IMF_VARIATION        # determines the stellar IMF for each particle from the Gusjenov/Hopkins/Hennebelle/Chabrier/Padoan theory
# ---- physical stellar feedback mechanisms ---- #
#GALSF_FB_GASRETURN             # Paul Torrey's addition for stochastic gas return (modified for continuous return)
#GALSF_FB_HII_HEATING           # gas within HII regions around young stars is photo-heated to 10^4 K
#GALSF_FB_SNE_HEATING           # time-dependent heating from SNe (I & II) in shockwave radii around stars
#GALSF_FB_RPROCESS_ENRICHMENT=8 # tracks a set of 'dummy' species from neutron-star mergers (set to number: 8=extended model)
#GALSF_FB_RT_PHOTONMOMENTUM     # continuous acceleration from starlight (uses luminosity tree)
#GALSF_FB_RT_PHOTON_LOCALATTEN  # incident SED for GALSF_FB_RT_PHOTONMOMENTUM calculated w local attenuation of stars
#GALSF_FB_LOCAL_UV_HEATING      # use local estimate of spectral information for photoionization and photoelectric heating
#GALSF_FB_RPWIND_LOCAL          # turn on local radiation pressure coupling to gas
#GALSF_FB_RPWIND_FROMSTARS      # drive radiation pressure with local young stars (otherwise uses the gas SFR)
##-----------------------------------------------------------------------------------------------------
# ----------- deprecated options (most have been combined or optimized into the functions above, here for legacy)
##GALSF_FB_SEPARATELY_TRACK_LUMPOS  # keep stellar vs. gas positions separate in tree (useful if running in tree-only mode)
##GALSF_FB_RPWIND_FROMCLUMPS	# clump-find to for wind angles (deprecated; now use density routine to find)
##GALSF_FB_RPWIND_CONTINUOUS	# wind accel term is continuous (more expensive and introduces more artificial dissipation)
##GALSF_FB_RPWIND_DO_IN_SFCALC	# do IR wind loop in SFR routine (allows fof clump-finding, useful for very IR-thick, but slow)

####################################################################################################

#-------------------------------------- SMBH/AGN stuff; also heavily expanded with PFH models
#BLACK_HOLES                    # enables Black-Holes (master switch)
# ----- seed models
#BH_POPIII_SEEDS                # BHs seeded on-the-fly from dense, low-metallicity gas
# ----- accretion models/options
#BH_SWALLOWGAS                  # enables stochastic accretion of gas particles consistent with growth rate of hole
#BH_ALPHADISK_ACCRETION         # gas accreted into 'virtual' alpha-disk, and from there onto the BH
#BH_GRAVCAPTURE_SWALLOWS        # accretion determined only by resolved gravitational capture by the BH
##BH_GRAVCAPTURE_NOGAS          # as BH_GRAVCAPTURE_SWALLOWS, but excludes gas (use with other mdot estimator for gas)
#BH_GRAVACCRETION               # Gravitational instability accretion estimator from Hopkins & Quataert 2010
##BH_BONDI                      # Bondi-Hoyle style accretion model
##BH_VARIABLE_ACCRETION_FACTOR  # variable-alpha model as in Booth&Schaye 2009
##BH_USE_GASVEL_IN_BONDI        # surrounding gas velocity used with sounds speed in the Bondi rate
#BH_SUBGRIDBHVARIABILITY        # model variability below resolved dynamical time for BH
# ----- feedback models/options
#BH_BAL_WINDS                   # accreted particles are launched back out as high-vel BAL winds
#BH_PHOTONMOMENTUM              # continuous long-range IR radiation pressure acceleration from BH (needs GALSF_FB_RT_PHOTONMOMENTUM)
#BH_HII_HEATING                 # photo-ionization feedback from BH (needs GALSF_FB_HII_HEATING)
#BH_COMPTON_HEATING             # enable Compton heating in cooling function (needs BH_PHOTONMOMENTUM)
##BH_THERMALFEEDBACK            # couple a fraction of the BH luminosity into surrounding gas as thermal energy
# ----------- use the BH_DRAG options only in cosmological cases where M_BH is not >> other particle masses
#BH_DYNFRICTION                 # apply dynamical friction force to the BHs when m_bh not >> other particle mass
##BH_DRAG                       # Drag on black-holes due to accretion (w real mdot)
##BH_STRONG_DRAG                # Drag rate boosted as if BH is accreting at eddington (requires BH_DRAG)
##-----------------------------------------------------------------------------------------------------
# ----------- deprecated options (most have been combined or optimized into the functions above, here for legacy)
##EVALPOTENTIAL                 # computes gravitational potential
##BH_REPOSITION_ON_POTMIN       # repositions hole on potential minimum (requires EVALPOTENTIAL)
##DETACH_BLACK_HOLES            # Insert an independent data structure for BHs (currently exlicitly depends on SEPARATE_STELLARDOMAINDECOMP)
##BH_FOLLOW_ACCRETED_GAS_MOMENTUM # Follow momentum for each swallowed gas parcel (add to BH; this ignores dissipation on small scales)
##BH_SEED_ON_POTMIN             # Seed in minimal potential instead of max density
##BH_SEED_STAR_MASS_FRACTION=0.02 # minimum star mass fraction for BH seeding
##-----------------------------------------------------------------------------------------------------
#-------------------------------------- AGN-Bubble feedback (D. Sijacki)
#BUBBLES                        # generation of hot bubbles in an isolated halo or the the biggest halo in the run
#MULTI_BUBBLES                  # hot bubbles in all haloes above certain mass threshold (works only with FOF and without BUBBLES)
#EBUB_PROPTO_BHAR               # Energy content of the bubbles with cosmic time evolves as an integrated BHAR(z) over a Salpeter time (Di Matteo 2003 eq. [11])
#BH_BUBBLES                     # calculate bubble energy directly from the black hole accretion rate
#UNIFIED_FEEDBACK               # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot
####################################################################################################




####################################################################################################
#-------------------------------------- Driven turbulence (for turbulence tests, large-eddy sims)
####################################################################################################
#TURB_DRIVING                   # turns on turbulent driving/stirring. see begrun for parameters that must be set
#POWERSPEC_GRID=128             # activates on-the-fly calculation of the turbulent velocity, vorticity, and smoothed-velocity power spectra
#ADJ_BOX_POWERSPEC              # compiles in a code module that allows via restart-flag 6 the calculation of a gas velocity power spectrum of a snapshot with an adjustable box (user defined center and size)
####################################################################################################






####################################################################################################
#--------------------------------------- Multi-Threading (parallelization) options
####################################################################################################
#OPENMP=2                       # Masterswitch for explicit OpenMP implementation
#OMP_NUM_THREADS=4              # custom PTHREADs implementation (don't enable with OPENMP)
####################################################################################################



####################################################################################################
#--------------------------------------- Output/Input options
####################################################################################################
HAVE_HDF5						# needed when HDF5 I/O support is desired
#OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION       # input files assumed to be in double precision (otherwise float is assumed)
#OUTPUTPOTENTIAL                # forces code to compute+output potentials in snapshots
#OUTPUTACCELERATION             # output physical acceleration of each particle in snapshots
#OUTPUTCHANGEOFENERGY           # outputs rate-of-change of internal energy of gas particles in snapshots
#OUTPUT_VORTICITY				# outputs the vorticity vector
#OUTPUTTIMESTEP                 # outputs timesteps for each particle
#OUTPUTCOOLRATE					# outputs cooling rate, and conduction rate if enabled
#OUTPUTLINEOFSIGHT				# enables on-the-fly output of Ly-alpha absorption spectra
#OUTPUTLINEOFSIGHT_SPECTRUM
#OUTPUTLINEOFSIGHT_PARTICLES
#POWERSPEC_ON_OUTPUT            # compute and output power spectra (not used)
#RECOMPUTE_POTENTIAL_ON_OUTPUT	# update potential every output even it EVALPOTENTIAL is set
####################################################################################################



####################################################################################################
#-------------------------------------------- De-Bugging & special (usually test-problem only) behaviors
####################################################################################################
#DEVELOPER_MODE                 # allows you to modify various numerical parameters (courant factor, etc) at run-time
#GAMMA_ENFORCE_ADIABAT=(1.0)    # if set, this forces gas to lie -exactly- along the adiabat P=GAMMA_ENFORCE_ADIABAT*(rho^GAMMA)
#TEST_FOR_IDUNIQUENESS          # explicitly check if particles have unique id numbers (only use for special behaviors)
#LONGIDS                        # use long ints for IDs (needed for super-large simulations)
#ASSIGN_NEW_IDS                 # assign IDs on startup instead of reading from ICs
#READ_HSML                      # reads hsml from IC file

#USE_MPI_IN_PLACE               # MPI debugging: makes AllGatherV compatible with MPI_IN_PLACE definitions in some MPI libraries
#NO_ISEND_IRECV_IN_DOMAIN       # MPI debugging
#FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG # MPI debugging
#MPISENDRECV_SIZELIMIT=100      # MPI debugging
#MPISENDRECV_CHECKSUM           # MPI debugging
#DONOTUSENODELIST               # MPI debugging
#NOTYPEPREFIX_FFTW              # FFTW debugging (fftw-header/libraries accessed without type prefix, adopting whatever was
                                #chosen as default at compile of fftw). Otherwise, the type prefix 'd' for double is used.
#DOUBLEPRECISION_FFTW           # FFTW in double precision to match libraries
#DEBUG                          # enables core-dumps and FPU exceptions
#STOP_WHEN_BELOW_MINTIMESTEP    # forces code to quit when stepsize wants to go below MinSizeTimestep specified in the parameterfile
#SEPARATE_STELLARDOMAINDECOMP   # separate stars (ptype=4) and other non-gas particles in domain decomposition (may help load-balancing)
#DISABLE_SPH_PARTICLE_WAKEUP    # don't let gas particles move to lower timesteps based on neighbor activity (use for debugging)
####################################################################################################





####################################################################################################
####################################################################################################
##
## LEGACY CODE & PARTIALLY-IMPLEMENTED FEATURES: BEWARE EVERYTHING BELOW THIS LINE !!!
##
####################################################################################################
####################################################################################################


####################################################################################################
#---------------------------------------- On the fly FOF groupfinder
####################################################################################################
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_GROUP_MIN_LEN=32               # default is 32
#SUBFIND                            # enables substructure finder
#DENSITY_SPLIT_BY_TYPE=1+2+16+32    # 2^type for whch the densities should be calculated seperately
#MAX_NGB_CHECK=3                    # Max numbers of neighbours for sattlepoint detection (default = 2)
#SAVE_MASS_TAB                      # Saves the an additional array with the masses of the different components
#SUBFIND_SAVE_PARTICLELISTS         # Saves also phase-space and type variables parallel to IDs
#SO_VEL_DISPERSIONS                 # computes velocity dispersions for as part of FOF SO-properties
#ORDER_SNAPSHOTS_BY_ID
#SAVE_HSML_IN_IC_ORDER              # will store the hsml-values in the order of the particles in the IC file
#ONLY_PRODUCE_HSML_FILES            # only carries out density estimate
#KEEP_HSML_AS_GUESS                 # keep using hsml for gas particles in subfind_density
#LINKLENGTH=0.16                    # Linkinglength for FoF (default=0.2)
#NO_GAS_CLOUDS                      # Do not accept pure gaseous substructures
#WRITE_SUB_IN_SNAP_FORMAT           # Save subfind results in snap format
#DUSTATT=11                         # Includes dust attenuation into the luminosity calculation (using 11 radial bins)
#OBSERVER_FRAME                     # If defined, use CB07 Observer Frame Luminosities, otherwise CB07 Rest Frame Luminosities
#SO_BAR_INFO                        # Adds temperature, Lx, bfrac, etc to Groups
#SUBFIND_COUNT_BIG_HALOS=1e4        # Adds extra blocks for Halos with M_TopHat > SUBFIND_COUNT_BIG_HALOS
#KD_CHOOSE_PSUBFIND_LIMIT           # Increases the limit for the parallel subfind to the maximum possible
#KD_ALTERNATIVE_GROUP_SORT          # Alternative way to sort the Groups/SubGroupe before writing
#KD_CHOOSE_LINKING_LENGTH           # Special way to estimate the linking length
#SUBFIND_READ_FOF
#SUBFIND_COLLECTIVE_STAGE1
#SUBFIND_COLLECTIVE_STAGE2
#SUBFIND_ALTERNATIVE_COLLECTIVE
#SUBFIND_RESHUFFLE_CATALOGUE
#SUBFIND_RESHUFFLE_AND_POTENTIAL    #needs -DSUBFIND_RESHUFFLE_CATALOGUE and COMPUTE_POTENTIAL_ENERGY
#SUBFIND_DENSITY_AND_POTENTIAL      #only calculated density and potential and write them into snapshot
####################################################################################################



####################################################################################################
#--------------------------------------- Radiative Transfer
####################################################################################################
#RADTRANSFER                            # main switch; RT equation solved using the Conjugate Gradient iterative method
#RADTRANSFER_FLUXLIMITER                # introduces a flux limiter, so that the Ifront speed does not exceed c
#RADTRANSFER_MODIFY_EDDINGTON_TENSOR    # modifies the Eddington tensor to the fully anisotropic version
#RT_COOLING_PHOTOHEATING                # includes photoheating and cooling (using RT information)
#RT_RAD_PRESSURE                        # includes radiation pressure
#EDDINGTON_TENSOR_GAS                   # includes gas, works only together with RT_MULTI_FREQUENCY
#EDDINGTON_TENSOR_STARS                 # includes stars as sources of ionising photons
#EDDINGTON_TENSOR_SFR                   # uses sf partices (gas above certain density) as sources of ionizing photons
#EDDINGTON_TENSOR_BH                    # includes BH as source of ionising photons (not working)
#RT_OUTPUT_ET                           # outputs the eddington tensor (used for diagnostics)
#HYDROGEN_ONLY                          # sets hydrogen fraction to 1.0 (simplifies the chemistry)
#RT_INCLUDE_HE                          # includes helium cooling and collisional ionisation
#RT_MULTI_FREQUENCY                     # enables multi-frequency radiation transport. Here the integration
                                        # variable is the ionising intensity J_nu
####################################################################################################




####################################################################################################
#--------------------------------------- degenerate equation of state
####################################################################################################
#EOS_DEGENERATE
#EOS_IDEAL
#EOS_COULOMB_CORRECTIONS
#EOS_NSPECIES=3
#RELAXOBJECT
####################################################################################################



####################################################################################################
#--------------------------------------- nuclear network
####################################################################################################
#NUCLEAR_NETWORK
#FIXED_TEMPERATURE
#NEGLECT_DTDY_TERMS
#NETWORK_OUTPUT
#NETWORK_OUTPUT_BINARY
####################################################################################################







