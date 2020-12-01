#!/bin/bash

## download reference solution
wget -O snapshot_001.ref.hdf5 https://cloudstor.aarnet.edu.au/plus/s/5vIwL81CEm8uAlR/download

## compare at machine precision
h5diff -r --count=50 --relative=1.0E-07 output/snapshot_001.hdf5 snapshot_001.ref.hdf5
