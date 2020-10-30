#!/bin/bash

pmemd.cuda -O -i heat.in -p full.hmr.prmtop -c full.inpcrd -o md-heat.out -r md-heat.rst7 -x md-heat.nc -ref full.inpcrd
pmemd.cuda -O -i eqnpt.in -p full.hmr.prmtop -c md-heat.rst7 -o eqnpt.out -r eqnpt.rst7 -x traj000.nc -ref full.inpcrd
pmemd.cuda -O -i eqnpt.in -p full.hmr.prmtop -c eqnpt.rst7 -o eqnpt01.out -r eqnpt01.rst7 -x traj001.nc -ref full.inpcrd
pmemd.cuda -O -i eqnpt.in -p full.hmr.prmtop -c eqnpt01.rst7 -o eqnpt02.out -r eqnpt02.rst7 -x traj002.nc -ref full.inpcrd
pmemd.cuda -O -i eqnpt.in -p full.hmr.prmtop -c eqnpt02.rst7 -o eqnpt03.out -r eqnpt03.rst7 -x traj003.nc -ref full.inpcrd
pmemd.cuda -O -i eqnpt.in -p full.hmr.prmtop -c eqnpt03.rst7 -o eqnpt04.out -r eqnpt04.rst7 -x traj004.nc -ref full.inpcrd
pmemd.cuda -O -i mdin-00 -p full.hmr.prmtop -c eqnpt04.rst7 -o md-00.out -r md00.rst7 -x md00.nc -ref full.inpcrd
pmemd.cuda -O -i mdin-01 -p full.hmr.prmtop -c md00.rst7 -o md-01.out -r md01.rst7 -x md01.nc -ref full.inpcrd
pmemd.cuda -O -i mdin-02 -p full.hmr.prmtop -c md01.rst7 -o md-02.out -r md02.rst7 -x md02.nc -ref full.inpcrd
