#!/bin/bash

sed '/^#/d' < cooling/TREECOOL > TREECOOL
wget "http://www.tapir.caltech.edu/~phopkins/sims/isodisk_ics.hdf5"
mpirun ./GIZMO scripts/test_problems/isodisk.params
