#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=128  # more than enough for 25A sphere size FEP
#SBATCH -A ACCOUNT 
#              d-hh:mm:ss
#SBATCH --time=TIME
#SBATCH -J JOBNAME
#SBATCH -o slurm.%N.%j.out # STDOUT

## Load modules for qdynp
MODULES

## define qdynp location
QDYN
fepfiles=(FEPS)
temperature=298
run=10
finalMDrestart=md_0000_1000.re
seed="$[1 + $[RANDOM % 32767]]"
starttime=$(date +%s)
starttime_readable=$(date)
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

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/$seed/ eq1.inp # change the random seed to custom
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
fi
#RUN_FILES
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ]
done
#CLEANUP

endtime=$(date +%s)
endtime_readable=$(date)
# Calculate runtime
runtime=$((endtime - starttime))

# Convert runtime to hours:minutes:seconds
hours=$(($runtime / 3600))
minutes=$(($runtime % 3600 / 60))
seconds=$(($runtime % 60))

echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
echo "#    Slurm tasks: $SLURM_NTASKS"
echo "#    Starttime: $starttime_readable"
echo "#    Endtime: $endtime_readable"
echo "#    Runtime: ${hours}h:${minutes}m:${seconds}s"
echo "#    Random seed: $seed"
echo "#    Replicate Number: $run"
echo "#    Working Directory: $workdir"