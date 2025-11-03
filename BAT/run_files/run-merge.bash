


x=0
while [  $x -lt 10 ]; do
cd m0$x
sbatch SLURMM-run
cd ../
cd c0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
