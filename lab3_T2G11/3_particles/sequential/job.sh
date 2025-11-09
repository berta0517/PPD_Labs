#!/bin/bash

#SBATCH --job-name=ex3.1
#SBATCH -p std
#SBATCH --output=out_sequential_%j.out
#SBATCH --error=out_sequential_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=00:05:00
module purge
module load gcc/13.3.0 openmpi/5.0.3

make >> make.out || exit 1      # Exit if make fails

mpirun -np 1 ./partis_seq 1000 0
