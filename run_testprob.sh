#!/bin/bash

## delete comments from TREECOOL file
sed '/^#/d' < cooling/TREECOOL > TREECOOL

## download initial conditions
wget "http://www.tapir.caltech.edu/~phopkins/sims/isodisk_ics.hdf5"

## Bitbucket CI runs as root inside a Docker container, so it's ok
mpirun --allow-run-as-root ./GIZMO scripts/test_problems/isodisk.params
