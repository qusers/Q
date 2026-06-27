#!/bin/bash
#
# No-SLURM variant of run.sh. The SLURM job array that would normally fan out the
# individual temperature/replicate jobs is replaced by a sequential loop, so the
# multi-seed/replicate structure is preserved without a scheduler.

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
# Number of MPI ranks per replicate (replaces $SLURM_NTASKS). Set via the
# Q_NCORES environment variable; otherwise default to 8 (Q's sweet spot) or fewer
# if the machine has fewer physical cores. Count physical cores (not hyperthreads):
# mpirun --bind-to core needs ntasks <= physical cores or it aborts with
# "not enough slots".
if [ -n "$Q_NCORES" ]; then
    ntasks=$Q_NCORES
else
    phys_cores=$(lscpu -p=Core 2>/dev/null | grep -v '^#' | sort -u | wc -l)
    if ! [ "$phys_cores" -ge 1 ] 2>/dev/null; then
        phys_cores=$(nproc 2>/dev/null || echo 1)
    fi
    ntasks=$(( phys_cores >= 8 ? 8 : phys_cores ))
fi
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

# Which job(s) to run. By default run every temperature/replicate combination
# sequentially (the work a SLURM array would fan out via $SLURM_ARRAY_TASK_ID).
# Set Q_TID=<index> to run a single combination, so a host-side orchestrator can
# launch one container per replicate and run several in parallel.
if [ -n "$Q_TID" ]; then
    if [ "$Q_TID" -lt 0 ] || [ "$Q_TID" -ge "$TOTAL_JOBS" ]; then
        echo "Error: Q_TID ($Q_TID) out of range 0..$((TOTAL_JOBS - 1))"
        exit 1
    fi
    TID_LIST="$Q_TID"
else
    TID_LIST=$(seq 0 $((TOTAL_JOBS - 1)))
fi

starttime=$(date +%s)
starttime_readable=$(date)

for TID in $TID_LIST; do

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
