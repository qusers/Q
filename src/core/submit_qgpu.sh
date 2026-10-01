#!/bin/bash
#SBATCH --job-name=qgpu-test
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_MI250X_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=qgpu-%j.out
#SBATCH --error=qgpu-%j.err

set -euo pipefail
module use /appl/local/containers/test-modules/
ml rocm/7.2.4

QDYN="$HOME/code/Q/src/core/.build/qgpu/hip/gfx90a/double/qdyn"

cd "$SLURM_SUBMIT_DIR"

INPUT="eq5.inp"

test -x "$QDYN" || {
    echo "Can't not find the executable file: $QDYN" >&2
    exit 1
}

test -f "$INPUT" || {
    echo "Can't not find the input: $INPUT" >&2
    exit 1
}

echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Directory: $PWD"

srun --ntasks=1 "$QDYN" --gpu "$INPUT"