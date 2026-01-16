#!/bin/bash
#SBATCH --job-name=fsched
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH  --mem=16gb
source ~/.bashrc
micromamba activate qligfep_new
export OMP_NUM_THREADS=1  # Each Python job gets a single thread

#this is an analyzer file
#it automatically executes the analysis singlegroup script in parallel for all folders

# Inner loop: replicate 0–9,
	for rep in {0..9}; do
    		echo "analysing $rep"
    		python analysis_singlegroup.py $rep $prefix &
	done


# Wait for all background jobs to complete
wait
echo "All jobs finished."


wait
echo "All jobs finished."
