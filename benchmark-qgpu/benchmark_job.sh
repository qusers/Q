#!/bin/bash
#SBATCH --job-name=qgpu_benchmark
#SBATCH --output=logs/qgpu_benchmark_%j.log
#SBATCH --time=01:00:00
#SBATCH --gpus-per-node=a100:1
#SBATCH --mem=48G

mkdir -p logs

module load Miniconda3

# Proper conda initialization for batch jobs
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home6/p323093/.conda/envs/Q

cd /home6/p323093/code/Q

python ./benchmark-qgpu/main.py \
  --input /home6/p323093/code/Q/test/cdk2/eq4.inp \
  --bin /home6/p323093/code/Q/bin/qdyn \
  --replicates 1 2 4 6 8 10
