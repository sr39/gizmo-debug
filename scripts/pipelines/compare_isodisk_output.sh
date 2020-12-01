#!/bin/bash

## download reference solution
wget -O snapshot_001.ref.hdf5 https://cloudstor.aarnet.edu.au/plus/s/5vIwL81CEm8uAlR/download

## compare at machine precision
h5diff -v --use-system-epsilon output/snapshot_001.hdf5 snapshot_001.ref.hdf5
