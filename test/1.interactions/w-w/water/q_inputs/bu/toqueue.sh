#!/bin/bash
#sbatch -n 8  -t 6:00:00 -J jobname script.sh 
module load openmpi-x86_64
sbatch -n 16 -t 6:00:00 -J fepEQ runFEP.sh run
