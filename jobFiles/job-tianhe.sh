#!/bin/bash

rm *.out log.err
mkdir -p restart_dir1
mkdir -p restart_dir2

mkdir -p phi_dir
mkdir -p trackp_dir
#export OMP_STACKSIZE=4000M
export OMP_NUM_THREADS=6
yhrun  -N 16 -n 32 -c 6 -p TH_SR ./gtc
#yhrun -N 1 -n 1 -p TH_NEW ./gtc
#yhrun -N 64 -n 64 -p TH_NEW ./gtc
#yhrun -N 32 -n 32 -p TH_NEW ./gtc
