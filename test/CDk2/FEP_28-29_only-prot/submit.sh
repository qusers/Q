#!/bin/sh

#SBATCH --gres=gpu:v100:1
#SBATCH --nodes 1
#SBATCH --time  0-12:00:00  # d-hh:mm:ss
#SBATCH --account SNIC2019-2-1
time python ../../../bin/qdyn.py --verbose -t dualtop.top -m eq1.inp -d TEST > log
