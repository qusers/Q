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

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEPS)

# Debug prints
echo "Number of temperatures: ${#temperatures[@]}"
echo "Number of runs: $runs"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

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

# Calculate indices based on array task ID
NUMTEMPS=${#temperatures[@]}
TID=$((SLURM_ARRAY_TASK_ID - 1))  # Convert to 0-based indexing

# Calculate which temperature and run this job corresponds to
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

## Load modules for qdynp
MODULES

## define qdynp location
QDYN

# Fail on the first error, so a failed mpirun can't exit 0 and be recorded COMPLETED. Set after
# `module`, which returns non-zero on some stacks just for unload warnings.
set -eo pipefail

starttime=$(date +%s)
starttime_readable=$(date)

# Emitted from an EXIT trap so the footer is written even when set -e aborts the job midway.
# Keep in sync with run_neq.sh; parsed by QligFEP.IO.read_slurm_diagnostics.
write_footer() {
    rc=$?
    endtime=$(date +%s)
    runtime=$((endtime - starttime))
    echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
    echo "#    Slurm tasks: $SLURM_NTASKS"
    echo "#    Starttime: $starttime_readable"
    echo "#    Endtime: $(date)"
    echo "#    Runtime: $((runtime / 3600))h:$((runtime % 3600 / 60))m:$((runtime % 60))s"
    echo "#    Random seed: ${seed:-unknown}"
    echo "#    Replicate Number: ${run_num:-unknown}"
    echo "#    Working Directory: ${workdir:-unknown}"
    echo "#    Exit status: $rc"
}
trap write_footer EXIT

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

echo "=== Process Binding Information ==="
echo "slurm tasks per node: $SLURM_TASKS_PER_NODE"
echo "slurm cpus per node: $SLURM_JOB_CPUS_PER_NODE"
echo

echo -e "\n=== CPU Model Information ==="
lscpu | grep -E "Model name|Architecture|CPU op|Thread|Core|Socket|NUMA|CPU(s)" || true
echo

# Diagnostics only -- must never be fatal under set -e, hence the `|| true`.
echo -e "\n=== Available CPU List ==="
cpu_list=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}') || true
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

# qdyn's serial abort still exits 0, so a completed run is confirmed by "terminated normally" in
# its log, not by the exit code. Failing here also stops FEP N+1 from starting on a bad eq5.re.
incomplete=""
logs_to_check=(md*.log)  # only the first FEP equilibrates; later ones inherit eq5.re
if [ $index -lt 1 ]; then
    logs_to_check+=(eq*.log)
fi
for log in "${logs_to_check[@]}"; do
    grep -q "terminated normally" "$log" 2>/dev/null || incomplete="$incomplete $log"
done
if [ -n "$incomplete" ]; then
    echo "ERROR: FEP$((index+1)) replicate $run_num (T=$temperature, seed=$seed) is incomplete."
    echo "ERROR: logs without a normal qdyn termination:$incomplete"
    exit 1
fi
done
#CLEANUP
