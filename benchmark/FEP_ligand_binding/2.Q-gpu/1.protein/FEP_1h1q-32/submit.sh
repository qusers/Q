#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 10
#SBATCH --time 0-12:00:00
#SBATCH --account SNIC2019-2-1
module load CUDA/10.1.243
module load iccifort/2019.5.281
module load impi/2018.5.288
module load PyCUDA/2019.1.2-Python-3.7.4
python run.py