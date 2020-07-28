#!/bin/bash -l 
#module load gcc/5.3.0
#module load openmpi-1.8-x86_64
module load openmpi-x86_64
export qdyn=/home/apps/q-5.06/qdynp


#export qdyn=/home/jespers/software/qsource/bin/qdyn5p

time mpirun -np 16 $qdyn eq1.inp > eq1.log
time mpirun -np 16 $qdyn eq2.inp > eq2.log
time mpirun -np 16 $qdyn eq3.inp > eq3.log
time mpirun -np 16 $qdyn eq4.inp > eq4.log
time mpirun -np 16 $qdyn eq5.inp > eq5.log
#time mpirun -np 16 $qdyn eq6.inp > eq6.log
#time mpirun -np 16 $qdyn eq7.inp > eq7.log
#time mpirun -np 16 $qdyn eq8.inp > eq8.log
#time mpirun -np 16 $qdyn eq9.inp > eq9.log
#time mpirun -np 16 $qdyn eq10.inp > eq10.log
#time mpirun -np 16 $qdyn eq11.inp > eq11.log
