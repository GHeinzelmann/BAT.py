


x=0
while [  $x -lt 10 ]; do
cd t0$x
sbatch SLURMM-run
cd ../
cd l0$x
sbatch SLURMM-run
cd ../
cd a0$x
sbatch SLURMM-run
cd ../
cd r0$x
sbatch SLURMM-run
cd ../
let x=x+1
done

x=0
while [  $x -lt 10 ]; do
cd c0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
