#!/bin/bash


rm -r trajectory
mpirun -np 1 python dd.py > output.dat
