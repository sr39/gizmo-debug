#!/bin/bash

## emit Config.sh
echo "HYDRO_MESHLESS_FINITE_MASS
COOLING
MULTIPLEDOMAINS=64
ADAPTIVE_GRAVSOFT_FORGAS
METALS
GALSF
GALSF_SFR_MOLECULAR_CRITERION
GALSF_SFR_VIRIAL_SF_CRITERION=0" > Config.sh

## emit Makefile.systype
echo "SYSTYPE=\"bitbucket\"" > Makefile.systype

## re-build
make clean && make

## delete comments from TREECOOL file
sed '/^#/d' < cooling/TREECOOL > TREECOOL

## download initial conditions
wget "http://www.tapir.caltech.edu/~phopkins/sims/isodisk_ics.hdf5"

## run test problem
./GIZMO scripts/pipelines/isodisk.params
