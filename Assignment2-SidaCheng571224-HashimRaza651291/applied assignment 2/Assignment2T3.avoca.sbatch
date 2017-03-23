#!/bin/sh

#SBATCH --account=VR0084
#SBATCH --time=0:10:0
#SBATCH --nodes=1
#SBATCH --job-name=Assignment2T3
#SBATCH --output=Assignment2T3.out
#SBATCH --error=Assignment2T3.err

export OMP_NUM_THREADS=1
./Assignment2T3
export OMP_NUM_THREADS=2
./Assignment2T3
export OMP_NUM_THREADS=4
./Assignment2T3
export OMP_NUM_THREADS=8
./Assignment2T3


srun -n 8 --nodes 1 --ntasks-per-node 8 Assignment2T3 

