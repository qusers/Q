#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=512  # more than enough for 25A sphere size FEP
#SBATCH -A ACCOUNT
#              d-hh:mm:ss
#SBATCH --time=TIME
#SBATCH -J JOBNAME
#SBATCH --array=1-TOTAL_JOBS
#SBATCH -o slurm.run%a.%N.%j.out
#
# QmapFEP + SLURM variant.
#
# One array task is one replicate, and nothing else. Where run.sh treats the
# array index as a (temperature, replicate) pair and launches qdyn under MPI,
# this script:
#
#   - runs a single temperature, so the array index IS the replicate number;
#   - runs qdyn serially, because the target is a Raspberry Pi cluster where MPI
#     is not the cost-effective way to use the hardware. Concurrency comes from
#     SLURM running several array tasks at once (sbatch --array=1-N%C), not from
#     ranks within one task;
#   - writes .finished / .failed sentinel files, so an external orchestrator can
#     tell a successful run from a crashed one without parsing logs. SLURM alone
#     cannot: a run whose energies went to NaN still exits 0 and is recorded as
#     COMPLETED.
#
# Everything else — the directory layout, the seeds array, the restart chaining —
# matches run.sh, so results land where the existing SLURM tooling expects them.

set -eE
set -o pipefail

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEPS)

# Debug prints
echo "Number of runs: $runs"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

# Validate inputs
if [ -z "$runs" ] || [ "$runs" -eq 0 ]; then
    echo "Error: 'runs' variable is not set or is zero"
    exit 1
fi

if [ ${#seeds[@]} -eq 0 ]; then
    echo "Error: No seeds specified"
    exit 1
fi

# A single temperature is a precondition, not a coincidence: the array index maps
# straight onto the replicate number, which is only correct while there is one
# temperature to run at.
if [ ${#temperatures[@]} -ne 1 ]; then
    echo "Error: this variant expects exactly one temperature, got ${#temperatures[@]}"
    exit 1
fi
temperature=${temperatures[0]}

# One array task per replicate. SLURM arrays are 1-based and so are the replicate
# directories, so no offset arithmetic is needed.
run_num=$SLURM_ARRAY_TASK_ID

if [ "$run_num" -lt 1 ] || [ "$run_num" -gt "$runs" ]; then
    echo "Error: array task ID ($run_num) outside 1..$runs"
    exit 1
fi

# Indexed from the generated seeds array, so replicates are reproducible and
# distinct. Drawing a fresh $RANDOM here instead would make a rerun of one
# replicate a different experiment.
seed=${seeds[$run_num-1]}

## Load modules for qdynp
MODULES

## define qdyn location
QDYN
# The scheduler may point at a different build than the cluster default above
# (SLURM --export=Q_QDYN=...,Q_QFEP=...); fall back to the default when unset.
qdyn=${Q_QDYN:-$qdyn}
qfep_bin=${Q_QFEP:-QFEP}

# check_nan fails the step when Q's energies stopped being finite.
#
# qdyn exits 0 on a NaN'd trajectory, so without this the run would continue
# through every remaining window and be reported as a success. Returning non-zero
# trips the ERR trap and writes the .failed sentinel instead.
check_nan() {
    local logfile=$1
    if [ ! -f "$logfile" ]; then
        echo "ERROR: expected log file $logfile does not exist"
        return 1
    fi
    if grep -qE '(^|[[:space:]])[-+]?([Nn][Aa][Nn]|[Ii][Nn][Ff]([Ii][Nn][Ii][Tt][Yy])?)([[:space:]]|$)' "$logfile"; then
        echo "ERROR: non-finite energy detected in $logfile"
        return 1
    fi
    return 0
}

starttime=$(date +%s)
starttime_readable=$(date)

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

# Registered once the replicate directory exists, so the sentinel lands beside
# the run it describes. The path is absolute so a later cd cannot misplace it.
trap 'touch "$rundir/.failed"' ERR

echo "Running job in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed"
echo

echo "=== Process Binding Information ==="
echo "slurm tasks per node: $SLURM_TASKS_PER_NODE"
echo "slurm cpus per node: $SLURM_JOB_CPUS_PER_NODE"
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
    # Chain from the previous FEP's final restart. run.sh references a variable
    # named finalMDrestart here, which is never assigned anywhere - the one that
    # actually holds this value is restartfile.
    lastfep=FEP$index
    cp $workdir/$lastfep/$temperature/$run_num/$restartfile $rundir/eq5.re
fi
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
fi
#RUN_FILES
timeout 3m $qfep_bin < qfep.inp > qfep.out || [ $? -eq 124 ]
done
#CLEANUP

# Written last, after every FEP and its qfep analysis have completed without
# tripping the trap. Its presence is what marks the replicate as succeeded.
touch "$rundir/.finished"

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
echo "#    Replicate Number: $run_num"
echo "#    Working Directory: $workdir"
