#!/bin/bash


rm -r trajectory
mpirun -np 1 python equil-sdr.py > output0.dat
mpirun -np 1 python sdr-ti.py > output.dat
