#!/bin/bash


rm -r trajectory
mpirun -np 1 python rest.py > output.dat
