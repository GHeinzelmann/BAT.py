#!/bin/bash


rm -r trajectory
python equil-dd.py > output0.dat
python dd-ti.py > output.dat
