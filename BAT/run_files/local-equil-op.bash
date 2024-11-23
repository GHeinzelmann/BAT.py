#!/bin/bash


echo "mpirun -np 1 python equil-00.py > output-00.dat"
mpirun -np 1 python equil-00.py > output-00.dat
i=1
while [ $i -le RANGE ] ;
do
x=`printf "%02.0f" $i`
echo "mpirun -np 1 python equil-$x.py > output-$x.dat"
mpirun -np 1 python equil-$x.py > output-$x.dat
let i=$i+1
done
