#!/bin/bash

#SBATCH --job-name=ex2.1
#SBATCH -p gpu               # Usa la partición GPU
#SBATCH --output=out_matmul_%j.out
#SBATCH --error=out_matmul_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1         # Solicita 1 GPU
#SBATCH --time=00:05:00

module purge
module load nvhpc/24.9     # Ajusta según tu sistema

make >> make.out || exit 1 

./matmul 100 1
