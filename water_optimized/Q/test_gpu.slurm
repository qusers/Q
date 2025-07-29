#!/bin/bash
#SBATCH --job-name=big_water_test_gpu
#SBATCH --mail-user="goodstudyqaq@gmail.com"
#SBATCH --mail-type="ALL"
#SBATCH --time=02:00:00
#SBATCH --partition=gpu-long
#SBATCH --output=%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --gres=gpu:a100:1


# making sure we start with a clean module environment
module purge
echo "## Starting GPU test on $HOSTNAME"
echo "## Loading module"
module load ALICE/default
module load slurm
module load CUDA/12.3.2
module load GCC/11.3.0


TEST_DIR=$(pwd)
echo "## Current directory $TEST_DIR"
echo "## Available CUDA devices: $CUDA_VISIBLE_DEVICES"
echo "## Checking status of CUDA device with nvidia-smi"
nvidia-smi

# Record start time
START_TIME=$(date +%s)

echo "## Running test"
$HOME/code/Q/bin/qdyn --gpu $HOME/code/Q/test/water_30

# Record end time
END_TIME=$(date +%s)
echo "## Start time: $(START_TIME)"
echo "## End time: $(END_TIME)"

# Calculate and display elapsed time in seconds
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "## Elapsed time: ${ELAPSED_TIME} seconds"

echo "## Test finished. Goodbye"
