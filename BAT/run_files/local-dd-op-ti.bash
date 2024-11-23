#!/bin/bash


rm -r trajectory
mpirun -np 1 python equil-dd.py > output0.dat
mpirun -np 1 python dd-ti.py > output.dat
