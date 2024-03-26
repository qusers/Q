#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH -A ACCOUNT 
#              d-hh:mm:ss
#SBATCH --time=TIME
#SBATCH -J JOBNAME
#SBATCH -o slurm.%N.%j.out # STDOUT

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
MODULES

## define qdynp location
QDYN
fepfiles=(FEPS)
temperature=298
run=10
finalMDrestart=md_0000_1000.re
seed="$[1 + $[RANDOM % 32767]]"
starttime=$(date)

workdir=$(pwd)
inputfiles=$workdir/inputfiles
length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length);do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir

rundir=$tempdir/$run
mkdir -p $rundir
cd $rundir

cp $inputfiles/eq*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .



sed -i s/SEED_VAR/$seed/ eq1.inp
if [ $index -lt 1 ]; then
sed -i s/'0.000 1.000'/'1.000 0.000'/ eq*inp
cp $inputfiles/md*_F.inp .
else
sed -i s/'1.000 0.000'/'0.000 1.000'/ eq*inp
cp $inputfiles/md*_R.inp .
fi  

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
#RUN_FILES
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ] && echo "##### WARNING - QFEP TIMED OUT #####"
else
#EQ_FILES
#RUN_FILES_R
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ] && echo "##### WARNING - QFEP TIMED OUT #####"
done
#CLEANUP

endtime=$(date)
echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
echo "#    Starting time: $starttime"
echo "#    Ending time: $endtime"
echo "#    Random seed: $seed"
echo "#    Replicate Number: $run"
echo "#    Working Directory: $workdir"