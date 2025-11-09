#!/bin/bash

#SBATCH --job-name=ex1
#SBATCH -p std
#SBATCH --output=out_montecarlo_%j.out
#SBATCH --error=out_montecarlo_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=192
#SBATCH --nodes=1
#SBATCH --time=00:05:00
module purge
module load gcc/13.3.0 openmpi/5.0.3

make >> make.out || exit 1      # Exit if make fails

mpirun -np 1 ./montecarlo 10 100000000 10
mpirun -np 5 ./montecarlo 10 500000000 10
mpirun -np 10 ./montecarlo 10 1000000000 10
mpirun -np 20 ./montecarlo 10 2000000000 10
mpirun -np 40 ./montecarlo 10 4000000000 10
mpirun -np 60 ./montecarlo 10 6000000000 10
mpirun -np 80 ./montecarlo 10 8000000000 10
mpirun -np 100 ./montecarlo 10 10000000000 10
mpirun -np 150 ./montecarlo 10 15000000000 10
mpirun -np 192 ./montecarlo 10 19200000000 10