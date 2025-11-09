#!/bin/bash

#SBATCH --job-name=ex2
#SBATCH -p std
#SBATCH --output=out_fc_mpi_%j.out
#SBATCH --error=out_fc_mpi_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=80
#SBATCH --nodes=1
#SBATCH --time=00:08:00
module purge
module load gcc/13.3.0 openmpi/5.0.3

make >> make.out || exit 1      # Exit if make fails

mpirun -np 20 ./fc_mpi input_planes_10kk.txt 2 1 0
mpirun -np 40 ./fc_mpi input_planes_10kk.txt 2 1 0
mpirun -np 60 ./fc_mpi input_planes_10kk.txt 2 1 0
mpirun -np 80 ./fc_mpi input_planes_10kk.txt 2 1 0