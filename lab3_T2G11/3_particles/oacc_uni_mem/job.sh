#!/bin/bash

#SBATCH --job-name=ex3.2
#SBATCH -p gpu               # Usa la partición GPU
#SBATCH --output=out_partis_oacc_uni_mem_%j.out
#SBATCH --error=out_partis_oacc_uni_mem_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1         # Solicita 1 GPU
#SBATCH --time=00:05:00

module purge
module load nvhpc/24.9     # Ajusta según tu sistema

make >> make.out || exit 1   # Asegúrate de tener una regla para compilar vecadd_cuda

./partis_oacc_uni_mem 100000 0
