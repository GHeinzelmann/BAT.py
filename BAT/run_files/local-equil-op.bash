#!/bin/bash


echo "python equil-00.py > output-00.dat"
python equil-00.py > output-00.dat
i=1
while [ $i -le RANGE ] ;
do
x=`printf "%02.0f" $i`
echo "python equil-$x.py > output-$x.dat"
python equil-$x.py > output-$x.dat
let i=$i+1
done
