#!/bin/bash

pmemd.cuda -O -i mdin-01 -p full.hmr.prmtop -c md00.rst7 -o md-01.out -r md01.rst7 -x md01.nc -ref full.inpcrd
pmemd.cuda -O -i mdin-02 -p full.hmr.prmtop -c md01.rst7 -o md-02.out -r md02.rst7 -x md02.nc -ref full.inpcrd
