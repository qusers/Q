#!/bin/bash
#SBATCH --job-name=fsched
#SBATCH --time=2-13:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH  --mem=1gb
source ~/.bashrc

#to save space, the tarer file is used separately to tarball the produced neq folders

tar -cvf neq_protein.tgz neq_protein &
tar -cvf neq_water.tgz neq_water &


wait
echo "All tar finished."

rm neq_protein -r &
rm neq_water -r &
wait
echo "all rm finished"

