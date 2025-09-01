#!/bin/bash

#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --tasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=e05-biosoft-blu
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos

module load cray-python

export OMP_NUM_THREADS=1

for i in $(seq 0 0)
do
   cd run-fssh-$i
       /work/e05/e05/fivanovic/flavoured-cptk-X-SH-single-phase/cp2k/exe/archer2/cp2k.sopt run.inp > run.log &
   cd ../
done

wait

