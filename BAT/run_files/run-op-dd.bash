cd rest
cd a-comp
sbatch SLURMM-run
cd ../
cd l-comp
sbatch SLURMM-run
cd ../
cd t-comp
sbatch SLURMM-run
cd ../
cd c-comp
sbatch SLURMM-run
cd ../
cd r-comp
sbatch SLURMM-run
cd ../
cd ../

cd dd

cd e-comp
x=0
while [  $x -lt 10 ]; do
cd e0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
cd ../

cd v-comp
x=0
while [  $x -lt 10 ]; do
cd v0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
cd ../

cd f-comp
x=0
while [  $x -lt 10 ]; do
cd f0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
cd ../

cd w-comp
x=0
while [  $x -lt 10 ]; do
cd w0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
cd ../

cd ../
