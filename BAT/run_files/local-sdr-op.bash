#!/bin/bash


rm -r trajectory
mpirun -np 1 python sdr.py > output.dat
