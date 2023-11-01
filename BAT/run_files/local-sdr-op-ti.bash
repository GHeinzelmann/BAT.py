#!/bin/bash


rm -r trajectory
python equil-sdr.py > output0.dat
python sdr-ti.py > output.dat
