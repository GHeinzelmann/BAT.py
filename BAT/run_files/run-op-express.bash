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
sbatch SLURMM-run
cd ../
cd v-comp
sbatch SLURMM-run
cd ../
cd ../
