#!/bin/bash
#SBATCH --job-name=q-charge-target-smoke
#SBATCH --account=ugsei19097
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --chdir=/projects/prjs2157/astra-charge-change-perturbation
#SBATCH --output=/projects/prjs2157/astra-charge-change-perturbation/logs/smoke-%j.out

set -euo pipefail
release=$1
export PYTHONPATH="$release/source/src"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export TMPDIR="$release/tmp"
mkdir "$TMPDIR"
/gpfs/home3/dvidal/bash/envs/qligfep_new/bin/python \
  "$release/source/experiments/analytical-charge-validation/hpc_smoke.py" "$release"
