#!/bin/bash
#SBATCH --job-name=q-standard-mpi-gate
#SBATCH --account=ugsei19097
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --chdir=/projects/prjs2157/astra-charge-change-perturbation
#SBATCH --output=/projects/prjs2157/astra-charge-change-perturbation/logs/standard-mpi-%j.out
set -eo pipefail
module purge
module load 2023 OpenMPI/4.1.5-GCC-12.3.0
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPATH=/projects/prjs2157/astra-charge-change-perturbation/releases/no-softcore-787e76eb-20260908/source/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/home3/dvidal/bash/envs/qligfep_new/bin/python "$1/standard_mpi_gate.py" "$1"
