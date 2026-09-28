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

module load YOUR_ROCM_MODULE

QDYN="$HOME/code/Q/src/core/.build/qgpu/hip/gfx90a/double/qdyn"

cd "$SLURM_SUBMIT_DIR"

INPUT="eq5.inp"

test -x "$QDYN" || {
    echo "找不到可执行文件：$QDYN" >&2
    exit 1
}

test -f "$INPUT" || {
    echo "找不到输入文件：$INPUT" >&2
    exit 1
}

echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Directory: $PWD"

srun --ntasks=1 "$QDYN" --gpu "$INPUT"