#!/bin/bash

set -x

#pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.0_b3.dat \
#    -f 1.6 3.0 -b 3 -s 8 -p -n 1
#pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.0_b3.dat \
#    -f 1.6 3.0 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.0_b3.dat
#pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.0_b3.dat \
#    -f 1.6 3.0 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.0_b3.dat
#pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.0_b3.dat \
#    -f 1.6 3.0 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.0_b3.dat

pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.1_b3.dat \
    -f 1.6 3.1 -b 3 -s 8 -p -n 1
pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.1_b3.dat \
    -f 1.6 3.1 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.1_b3.dat
pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.1_b3.dat \
    -f 1.6 3.1 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.1_b3.dat
pysca data/scdata/LC_6443122.fits -o iau_6443122_sc_f1.6-3.1_b3.dat \
    -f 1.6 3.1 -b 3 -s 8 -p -n 1 -i iau_6443122_sc_f1.6-3.1_b3.dat
