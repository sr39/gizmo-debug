
CHIMES 
====== 

The CHIMES module tracks the time-dependent evolution of 157 ions and molecules, and uses the resulting ion and molecule abundances to compute the radiative cooling rates of the gas. This replaces the standard cooling routines in Gizmo. 

The details of the CHIMES model are described in `Richings+ 2014a` and `Richings+ 2014b`.  

.. _Richings+ 2014a: http://adsabs.harvard.edu/abs/2014MNRAS.440.3349R 
.. _Richings+ 2014b: http://adsabs.harvard.edu/abs/2014MNRAS.442.2780R 

If you would like to use this module for a new project, please contact the authors: 

* Alex Richings 
* Joop Schaye 
* Ben Oppenheimer 

Usage 
----- 

The CHIMES module can be switched on using the ``CHIMES`` Config option. There are several additional Config options that control the behaviour of CHIMES, which are described further below. CHIMES can also be used together with the standard FIRE sub-grid physics models. Most of the FIRE options can be used as normal, but there are a few exceptions (see the CHIMES Config.sh Options section below). 

To run CHIMES, you will also need the CHIMES data files that contain the various reaction rate coefficients, photoionisation cross sections etc. These can be downloaded from the following bitbucket repository: 

https://bitbucket.org/richings/chimes-data 

The paths to the various data files can then be specified via the parameter file, as described in the CHIMES Parameters section below. 

CHIMES also requires the CVode and KINsol libraries, from the Sundials package which can be downloaded here: 

https://computation.llnl.gov/projects/sundials 

Note that CHIMES hasn't been tested with the latest version of CVode and KINsol. It is safest to use version 2.4.0 of these libraries. The Gizmo Makefile contains examples for how to link to these libraries (see under the Quest-intel SYSTYPE). 

OPENMP
  The CHIMES module is compatible with multithreading using the OPENMP option. When this option is switched on, then the cooling loop, where the code loops over active gas particles and performs the chemistry and cooling, will be split between OPENMP threads. This option greatly improves the work-load balancing of the chemistry solver, because particles are dynamically distributed between threads as their chemistry is being integrated, which helps to avoid the scenario where one CPU has to integrate several particles that are much slower (e.g. if they are far from chemical equilibrium) while the remaining CPUs have already finished. In many cases you will find that the most efficient approach is to run with just one or two MPI tasks on each node, and use OPENMP to then fill the CPUs on each node with threads. 

CHIMES Config.sh Options
------------------------

CHIMES
  Main switch to enable the CHIMES module. 

CHIMES_HYDROGEN_ONLY
  Option to run with only hydrogen. This option is ignored if METALS is set. 

CHIMES_SOBOLEV_SHIELDING
  The shielding length is set by a Sobolev-like approximation, with L_shield = h_inter + rho / grad(rho), multiplied by the ``Shielding_length_factor`` parameter (see the next section). Column densities are then set to the cell density times this shielding length. 

CHIMES_INITIALISE_IN_EQM
  By default, Gizmo will read in the initial ion and molecule abundances from the ICs file. However, when this option is enabled, Gizmo will instead compute the initial abundances in chemical equilibrium, by evolving the chemistry of each gas cell for a long period of time. 

CHIMES_HII_REGIONS 
  Self-shielding is disabled for gas particles that are flagged as HII regions. Note that this should be used instead of FIRE's standard GALSF_FB_FIRE_RT_HIIHEATING option. 

CHIMES_STELLAR_FLUXES 
  Tracks stellar fluxes in different age bins through the luminosity tree, using the LEBRON method from FIRE, and couples these to the CHIMES chemistry routines. Note that this should be used instead of FIRE's standard GALSF_FB_FIRE_RT_UVHEATING routine. 

CHIMES_SFR_MOLECULAR_CRITERION 
  As GALSF_SFR_MOLECULAR_CRITERION, but using the H2 fractions from CHIMES. 

CHIMES_REDUCED_OUTPUT 
  By default, the full array of chemical abundances (with 157 species) is written out in every snapshot. When this option is turned on, the full abundance array is only written out in a subset of snapshots, with a frequency given by the ``N_chimes_full_output_freq`` parameter (see the next section). For all other snapshots, it just writes out a smaller array containing the abundances of electrons, HI, H2 and CO. 

CHIMES_NH_OUTPUT 
  Writes out the column densities of gas particles, as used within the chemistry solver, in the snapshots. This is useful if you want to compute equilibrium abundances from the snapshot in post-processing, as it ensures that you use the same column densities that went into the non-equilibrium calculation. 

CHIMES_OUTPUT_DENS_AROUND_STAR 
  Writes out the gas density computed at the position of each star particle in the snapshot. If both CHIMES_OUTPUT_DENS_AROUND_STAR and CHIMES_NH_OUTPUT are switched on, it will also write out the gas column density around each star particle, which can be used for computing fluxes from each star particle including local shielding around the star in post-processing. 

CHIMES_OUTPUT_DELAY_TIME_HII 
  Writes out the DelayTimeHII of each gas particle. Requires CHIMES_HII_REGIONS or GALSF_FB_FIRE_RT_HIIHEATING. This is useful for identifying gas particles that were flagged as being within an HII region. 

CHIMES_TURB_DIFF_IONS 
  Switches on the turbulent diffusion of individual ions and molecules. This is treated in the same way as FIRE's standard metal diffusion option. Requires TURB_DIFF_METALS and TURB_DIFF_METALS_LOWORDER. 

CHIMES_METAL_DEPLETION 
  Reduces the abundance of metals in the gas-phase according to observed density-dependent metal depletion factors from Jenkins (2009) and De Cia et al. (2016). It also computes a density-dependent dust to gas ratio that is consistent with these depletion factors. 

CHIMES Parameters
---------------------

* ``Chimes_data_path`` Path to the directory containing the CHIMES data files, which you downloaded from the above bitbucket repository. The directory structure of this repository should not be changed - CHIMES will then be able to find the various individual files that it needs within this directory. 

* ``PhotoIon_table_path`` Path to a text file that contains the path(s) to the photoionisation cross sections HDF5 file(s) that will be used for the given UV radiation field, with one path per line. The HDF5 file is specified in this way, via a text file, as this makes it easier to specify a variable number of radiation fields. If you are using the CHIMES_STELLAR_FLUXES option, you will need to first specify the cross sections file for the extragalactic UV background, and then those for the Starburst99 spectra from each of the 8 stellar age bins that are tracked separately. 

* ``Thermal_Evolution_On`` Flag to switch on temperature evolution in CHIMES. When this flag is 0, the chemical abundances will be evolved at constant temperature in CHIMES. When it is 1, CHIMES will evolve the temperature along with the chemical abundances. *Typical value: 1*. 

* ``Chemistry_eqm`` Flag to determine whether to folow the time-dependent chemical evolution, or set the abundances to chemical equilibrium. When this flag is 0, the full time-dependent evolution is used. When it is 1, the abundances are set to equilibrium using pre-computed equilibrium abundance tables (although currently, the path to this equilibrium table is hard-coded to a Dummy Table, so it is difficult to use this option at the moment). When it is 2, the equilibrium abundances are computed on the fly by setting the chemical rate equations to zero. *Typical value: 0*. 

* ``StaticMolCooling`` The molecular cooling from CO and H2O depends on line broadening. When this flag is set to 0, we take into account the divergence of the velocity here. When this flag is set to 1, we only include thermal broadening for the CO and H2O cooling. *Typical value: 0*. 

* ``CellSelfShieldingOn`` Flag to determine whether self shielding is included. When this flag is 0, self shielding is not included. When it is 1, it is included, and the column densities of individual species (e.g. HI, H2 etc.) are updated throughout the course of the integration of the rate equations. However, updating the column densities in this way is very slow. When this flag is 2, self shielding is included, but column densities are only updated at the beginning of the hydrodynamic time-step, and not throughout the course of the CHIMES integration. Note that, to switch on self shielding, you will also need to enable the ``CHIMES_SOBOLEV_SHIELDING`` Config option, which will define the shielding length. *Typical value: 2*. 

* ``Shielding_length_factor`` The shielding length is multiplied by this factor. This allows you to control the normalisation of the shielding length. *Typical value: 0.5 (for Sobolev shielding)*. 

* ``Grain_Temperature`` The temperature of dust grains in Kelvin, as used when computing the formation rate of H2 on dust grains. *Typical value: 10*. 

* ``CrRate`` Cosmic ray ionisation rate of HI. The cosmic ray ionisation rate of all other species are then scaled relative to this parameter. *Typical value: 1.8e-16*. 

* ``max_mol_temperature`` Molecules are excluded above this temperature. *Typical value: 1.0e5*. 

* ``Isotropic_photon_density`` This sets the normalisation of the UV radiation field. It is defined as the number density of hydrogen-ionising photons in cgs units. The flux of hydrogen-ionising photons is then the Isotropic_photon_density multiplied by the speed of light. If the ``CHIMES_STELLAR_FLUXES`` Config option is switched on, this should just be set to the normalisation of the extragalactic UV background component. The normalisations of the other radiation fields will then be set by the stellar fluxes computed within the simulation. 

* ``relativeTolerance`` This controls the accuracy of the thermo-chemistry integration. CHIMES will sub-cycle each hydro time-step, aiming to achieve a relative error as given by this parameter. *Typical value: 1.0e-4*. 

* ``absoluteTolerance`` This controls the accuracy of the chemistry integration. Species with an abundance much below the absolute tolerance are not taken into account when determining the sub-steps needed to achieve a given relative error. *Typical value: 1.0e-10*. 

* ``thermalAbsoluteTolerance`` As the ``absoluteTolerance`` parameter, but for the thermal energy. *Typical value: 1.0e-40*. 

* ``scale_metal_tolerances`` In cosmological simulations, the metal element abundances can be arbitrarily small. If an element abundance is zero, its ions will not be included in the chemical network. However, if it is non-zero but lower than the absolute tolerance, this can cause problems for the chemical integration, because none of that element's ions are taken into account when determining the sub-steps. To avoid this problem, set the ``scale_metal_tolerances`` parameter to 1. Then the absolute tolerance of each individual ion and molecule species will be set to the ``absoluteTolerance`` parameter multiplied by that species' corresponding element abundance. If this parameter is set to 0, CHIMES will use a constant absolute tolerance. *Typical value: 1 (cosmological), 0 (non-cosmological)*. 

* ``IncludeCarbon`` Set this flag to 1 to include carbon in the CHIMES network. Set it to 0 to exclude carbon. *Typical value: 1*. 

* ``IncludeNitrogen`` Set this flag to 1 to include nitrogen in the CHIMES network. Set it to 0 to exclude nitrogen. *Typical value: 1*. 

* ``IncludeOxygen`` Set this flag to 1 to include oxygen in the CHIMES network. Set it to 0 to exclude oxygen. *Typical value: 1*. 

* ``IncludeNeon`` Set this flag to 1 to include neon in the CHIMES network. Set it to 0 to exclude neon. *Typical value: 1*. 

* ``IncludeMagnesium`` Set this flag to 1 to include magnesium in the CHIMES network. Set it to 0 to exclude magnesium. *Typical value: 1*. 

* ``IncludeSilicon`` Set this flag to 1 to include silicon in the CHIMES network. Set it to 0 to exclude silicon. *Typical value: 1*. 

* ``IncludeSulphur`` Set this flag to 1 to include sulphur in the CHIMES network. Set it to 0 to exclude sulphur. *Typical value: 1*. 

* ``IncludeCalcium`` Set this flag to 1 to include calcium in the CHIMES network. Set it to 0 to exclude calcium. *Typical value: 1*. 

* ``IncludeIron`` Set this flag to 1 to include iron in the CHIMES network. Set it to 0 to exclude iron. *Typical value: 1*. 

* ``N_chimes_full_output_freq`` If the CHIMES_REDUCED_OUTPUT Config option is switched on, this parameter determines the frequency with which the full chemistry array is written out to the snapshot. *Typical value: 10*. 

* ``Chimes_f_esc_ion`` If the CHIMES_STELLAR_FLUXES Config option is switched on, this parameter sets the escape fraction of ionising photons (>13.6 eV) from HII regions. 

* ``Chimes_f_esc_G0`` If the CHIMES_STELLAR_FLUXES Config option is switched on, this parameter sets the escape fraction of photons in the 6-13.6 eV band from HII regions. 

The following parameters are related to an obsolete option in CHIMES that is no longer used. These will be removed from CHIMES in the future, but for now you should just set these to their typical values, so that this option is turned off. 

* ``Reduction_On`` *Typical value: 0*. 

* ``Reduction_N_Ions_Low`` *Typical value: 3*. 

* ``Reduction_N_Ions_Med`` *Typical value: 4*. 

* ``Reduction_N_Ions_High`` *Typical value: 5*. 

* ``macrostep_tolerance`` *Typical value: 1.0*. 

* ``min_macrostep`` *Typical value: 1.0e2*. 

