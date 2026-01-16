#!/bin/bash
#SBATCH --job-name=fsched
#SBATCH --time=2-13:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=rome
#SBATCH  --mem=64gb
source ~/.bashrc

export OMP_NUM_THREADS=1  # Each Python job gets a single thread
mkdir logs


#This is the submission file for the water side of the FEP
#Adjust the cpus per task to roughly match the number of your edges * number of replicates
#in this case, 22 edges * 10 replicates = 220 cpus.
#so we just set the cpus per task to maximum allowable quantity of 128


# Outer loop: iterate over subfolders in ./water
for folder in ./water/*/; do
    folder_name=$(basename "$folder")  # Extract folder name
    echo "Processing folder: $folder_name"

    # Middle loop: water / protein
    for system in water; do
        # Inner loop: replicate 0–9
        for rep in {0..9}; do
            echo "Launching: $folder_name - $system replicate $rep"
            python folder_execution_neq.py  "$folder_name" "$system" "$rep" > "logs/${folder_name}_${system}_${rep}.log" &
	    
        done
    done
done

# Wait for all background jobs to complete
wait
echo "All jobs finished."


wait
echo "All jobs finished."
