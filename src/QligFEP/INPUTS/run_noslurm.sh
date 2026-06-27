#!/bin/bash
#
# No-SLURM variant of run.sh. The SLURM job array that would normally fan out the
# individual temperature/replicate jobs is replaced by a sequential loop, so the
# multi-seed/replicate structure is preserved without a scheduler.

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
ntasks=NTASKS  # number of MPI ranks (replaces $SLURM_NTASKS)
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEPS)

# Debug prints
echo "Number of temperatures: ${#temperatures[@]}"
echo "Number of runs: $runs"

# Validate inputs
if [ -z "$runs" ] || [ "$runs" -eq 0 ]; then
    echo "Error: 'runs' variable is not set or is zero"
    exit 1
fi

if [ ${#temperatures[@]} -eq 0 ]; then
    echo "Error: No temperatures specified"
    exit 1
fi

if [ ${#seeds[@]} -eq 0 ]; then
    echo "Error: No seeds specified"
    exit 1
fi

## Load modules for qdynp
MODULES

## define qdynp location
QDYN

# Total number of (temperature x replicate) jobs the SLURM array would have run
NUMTEMPS=${#temperatures[@]}
TOTAL_JOBS=$((NUMTEMPS * runs))

starttime=$(date +%s)
starttime_readable=$(date)

# Loop over every temperature/replicate combination sequentially (the work a SLURM
# array task would otherwise pick up via $SLURM_ARRAY_TASK_ID).
for TID in $(seq 0 $((TOTAL_JOBS - 1))); do

# Calculate which temperature and replicate this iteration corresponds to
temp_idx=$((TID / runs))
run_num=$((TID % runs + 1))

# Validate indices
if [ $temp_idx -ge ${#temperatures[@]} ]; then
    echo "Error: Temperature index ($temp_idx) out of bounds"
    exit 1
fi

if [ $((run_num - 1)) -ge ${#seeds[@]} ]; then
    echo "Error: Run number ($run_num) exceeds number of seeds"
    exit 1
fi

# Get the actual values for this job
temperature=${temperatures[$temp_idx]}
seed=${seeds[$run_num-1]}

length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length); do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir || exit
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir || exit

rundir=$tempdir/$run_num
mkdir -p $rundir
cd $rundir || exit

echo "Running job in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed"
echo

echo -e "\n=== CPU Model Information ==="
lscpu | grep -E "Model name|Architecture|CPU op|Thread|Core|Socket|NUMA|CPU(s)"
echo

echo -e "\n=== Available CPU List ==="
cpu_list=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2}')
echo "CPU list: $cpu_list"
echo

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i "s/SEED_VAR/$seed/" eq1.inp # change the random seed to custom
else
    lastfep=FEP$index
    cp $workdir/$lastfep/$temperature/$run_num/$finalMDrestart $rundir/eq5.re
fi
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
fi
#RUN_FILES
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ]
done
#CLEANUP

done

endtime=$(date +%s)
endtime_readable=$(date)
# Calculate runtime
runtime=$((endtime - starttime))

# Convert runtime to hours:minutes:seconds
hours=$(($runtime / 3600))
minutes=$(($runtime % 3600 / 60))
seconds=$(($runtime % 60))

echo "#    Starttime: $starttime_readable"
echo "#    Endtime: $endtime_readable"
echo "#    Runtime: ${hours}h:${minutes}m:${seconds}s"
echo "#    Working Directory: $workdir"
