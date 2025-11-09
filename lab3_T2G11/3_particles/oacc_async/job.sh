#!/bin/bash

#SBATCH --job-name=ex3.4
#SBATCH -p gpu               # Usa la partición GPU
#SBATCH --output=out_partis_oacc_async_%j.out
#SBATCH --error=out_partis_oacc_async_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1         # Solicita 1 GPU
#SBATCH --time=00:05:00

module purge
module load nvhpc/24.9     # Ajusta según tu sistema

make >> make.out || exit 1  

nsys profile ./partis_oacc_async 10000 0
