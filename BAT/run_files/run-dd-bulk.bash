

x=0
while [  $x -lt 10 ]; do
cd f0$x
sbatch SLURMM-run
cd ../
cd w0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
