#!/bin/bash


rm -r trajectory
python equil-sdr-ee.py > output0-ee.dat
python sdr-ti-ee.py > output-ee.dat
cp complex_trajectory_0.dcd complex_trajectory_0-ee.dcd
cp restart.chk restart-ee.chk
rm -r trajectory
python equil-sdr-1v.py > output0-1v.dat
python sdr-ti-1v.py > output-1v.dat
cp complex_trajectory_0.dcd complex_trajectory_0-1v.dcd
rm -r trajectory
python equil-sdr-2v.py > output0-2v.dat
python sdr-ti-2v.py > output-2v.dat
cp complex_trajectory_0.dcd complex_trajectory_0-2v.dcd

