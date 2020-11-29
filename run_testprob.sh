#!/bin/bash

## delete comments from TREECOOL file
sed '/^#/d' < cooling/TREECOOL > TREECOOL

## download initial conditions
wget "http://www.tapir.caltech.edu/~phopkins/sims/isodisk_ics.hdf5"

## run test problem
mpirun ./GIZMO scripts/test_problems/isodisk.params
