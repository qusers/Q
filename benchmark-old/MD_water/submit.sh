#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#              d-hh:mm:ss
#SBATCH --time=0-10:00:00
#SBATCH -A  SNIC2019-2-1
#SBATCH --gres=gpu:k80:1
module load CUDA/10.1.243
module load iccifort/2019.5.281
module load impi/2018.5.288
module load PyCUDA/2019.1.2-Python-3.7.4

# run the tests
python run_test.py
