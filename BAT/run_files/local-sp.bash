#!/bin/bash


rm -r trajectory
python sdr-ee.py > output-ee.dat
rm -r trajectory
python sdr-1v.py > output-1v.dat
rm -r trajectory
python sdr-2v.py > output-2v.dat

