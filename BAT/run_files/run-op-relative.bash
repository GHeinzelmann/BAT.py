cd rest
cd m-comp
sbatch SLURMM-run
cd ../
cd n-comp
sbatch SLURMM-run
cd ../
cd ../

cd sdr

cd e-comp
x=0
while [  $x -lt 10 ]; do
cd e0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
if [ $x -ge 10 ]; then
while [  $x -lt 12 ]; do
cd e$x
sbatch SLURMM-run
cd ../
let x=x+1
done
fi
cd ../

cd x-comp
y=0
while [  $y -lt 10 ]; do
cd x0$y
sbatch SLURMM-run
cd ../
let y=y+1
done
if [ $y -ge 10 ]; then
while [  $y -lt 12 ]; do
cd x$y
sbatch SLURMM-run
cd ../
let y=y+1
done
fi
cd ../

cd ../
