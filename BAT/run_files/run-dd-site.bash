

x=0
while [  $x -lt 7 ]; do
cd e0$x
sbatch SLURMM-run
cd ../
cd v0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
