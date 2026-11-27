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
# One array task runs one replicate, or several side by side. Where run.sh treats
# the array index as a (temperature, replicate) pair and launches qdyn under MPI,
# this script:
#
#   - runs a single temperature. With Q_REPLICATES_PER_TASK=k (default 1), task t
#     runs replicates (t-1)k+1 .. tk, so with k=1 the array index IS the
#     replicate number. Submit ceil(replicates/k) tasks with --array: the
#     #SBATCH --array line above counts one task per replicate;
#   - runs qdyn serially, one process per replicate. Concurrency comes from SLURM
#     running several array tasks at once (sbatch --array=1-N%C) and, with k > 1,
#     from the replicates of a task running side by side - on a GPU engine such
#     as QGPU they share the task's GPU, under MPS when it can be started;
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

# Which replicates this task runs. SLURM arrays are 1-based and so are the
# replicate directories.
per_task=${Q_REPLICATES_PER_TASK:-1}
if ! [[ "$per_task" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Q_REPLICATES_PER_TASK must be a positive whole number, got '$per_task'"
    exit 1
fi
first=$(( (SLURM_ARRAY_TASK_ID - 1) * per_task + 1 ))
last=$(( SLURM_ARRAY_TASK_ID * per_task ))
if [ "$last" -gt "$runs" ]; then
    last=$runs
fi
if [ "$first" -lt 1 ] || [ "$first" -gt "$runs" ]; then
    echo "Error: array task $SLURM_ARRAY_TASK_ID has no replicates: $per_task per task, $runs replicates"
    exit 1
fi
echo "Replicates in this task: $first..$last"

## Load modules for qdynp
MODULES

## define qdyn location
QDYN
# The scheduler may point at a different build than the cluster default above
# (SLURM --export=Q_QDYN=...,Q_QFEP=...); fall back to the default when unset.
qdyn=${Q_QDYN:-$qdyn}
qfep_bin=${Q_QFEP:-QFEP}
# Arguments for qdyn from the scheduler, such as --gpu for the QGPU engine, which
# reads the same inputs and is otherwise called the same way. Left unquoted where
# it is used, so several arguments split and none expands to nothing.
qdyn_args=${Q_QDYN_ARGS:-}

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

# Marking this replicate failed, wherever the failure came from.
#
# $rundir is only set once the loop reaches its replicate directory. The guard is
# what makes this safe to register before that point: a signal arriving during
# setup must not touch /.failed.
mark_failed() {
    if [ -n "${rundir:-}" ]; then
        touch "$rundir/.failed"
    fi
}

# run_replicate runs one replicate through every FEP, in a subshell of its own:
# its traps, its sentinels and its footer are its own, whether it runs alone or
# beside the other replicates of its task.
#
# SLURM sends SIGTERM KillWait seconds before SIGKILL when a job exceeds its time
# limit or is scancelled, and that is the only chance to record the outcome: the
# ERR trap does not fire for a signal. It matters because the sentinels are the
# whole status contract - on a cluster without accounting storage, squeue forgets
# the job within MinJobAge, and after that nothing distinguishes a killed
# replicate from one that was never submitted. SIGKILL cannot be caught, so a
# node failure or the OOM killer still leaves no sentinel behind.
run_replicate() {
    run_num=$1
    # Indexed from the generated seeds array, so replicates are reproducible and
    # distinct. Drawing a fresh $RANDOM here instead would make a rerun of one
    # replicate a different experiment.
    seed=${seeds[$run_num-1]}
    rundir=""
    trap 'echo "Terminated by signal - marking replicate $run_num failed"; mark_failed; exit 143' TERM INT

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
    trap 'mark_failed' ERR

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
    # The same start and end as seconds since the epoch, which read the same in any
    # locale or timezone; the analysis takes the run's elapsed time from them.
    echo "#    Start epoch: $starttime"
    echo "#    End epoch: $endtime"
    echo "#    Random seed: $seed"
    echo "#    Replicate Number: $run_num"
    echo "#    Working Directory: $workdir"
}

# Several replicates on one GPU run best under MPS, which lets their kernels run
# at the same time instead of taking turns. The daemon is the task's own, with
# its pipe and log directories under the job, so it cannot meet another job's.
# Without a GPU, with one replicate, or where the cluster already provides MPS,
# nothing is started; when it cannot start, the replicates still share the GPU,
# only less efficiently. Q_MPS=off turns it off.
mps_started=""
start_mps() {
    if [ "${Q_MPS:-auto}" = "off" ] || [ "$per_task" -eq 1 ]; then
        return 0
    fi
    if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then
        echo "No GPU in this task: not starting MPS"
        return 0
    fi
    if [ -n "${CUDA_MPS_PIPE_DIRECTORY:-}" ]; then
        echo "MPS is provided by the cluster at $CUDA_MPS_PIPE_DIRECTORY"
        return 0
    fi
    if ! command -v nvidia-cuda-mps-control >/dev/null 2>&1; then
        echo "WARNING: nvidia-cuda-mps-control not found: the replicates share the GPU without MPS"
        return 0
    fi
    mps_dir="${TMPDIR:-/tmp}/qmapfep-mps-${SLURM_JOB_ID}"
    export CUDA_MPS_PIPE_DIRECTORY="$mps_dir/pipe"
    export CUDA_MPS_LOG_DIRECTORY="$mps_dir/log"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
    if nvidia-cuda-mps-control -d; then
        mps_started=1
        echo "Started MPS for GPU(s) $CUDA_VISIBLE_DEVICES in $mps_dir"
    else
        echo "WARNING: could not start MPS: the replicates share the GPU without it"
        unset CUDA_MPS_PIPE_DIRECTORY CUDA_MPS_LOG_DIRECTORY
    fi
}

stop_mps() {
    if [ -n "$mps_started" ]; then
        echo quit | nvidia-cuda-mps-control || true
        rm -rf "$mps_dir"
    fi
}

if [ "$per_task" -eq 1 ]; then
    # One replicate: it runs in the foreground and writes to this task's log, as
    # it always has.
    ( run_replicate "$first" )
    exit $?
fi

start_mps
trap stop_mps EXIT

# Several replicates: side by side, each writing its output and footer to a log
# of its own, which the analysis reads per replicate like a task's log.
pids=()
for run_num in $(seq "$first" "$last"); do
    ( run_replicate "$run_num" ) > "$workdir/slurm.run${run_num}.${SLURM_JOB_ID}.replicate.out" 2>&1 &
    pids+=($!)
    echo "Replicate $run_num started (pid $!)"
done

# A time limit or scancel reaches this shell; pass it on, so every replicate
# writes its .failed before SIGKILL.
trap 'echo "Terminated by signal - stopping replicates $first..$last"; kill -TERM "${pids[@]}" 2>/dev/null; wait; exit 143' TERM INT

failed=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        failed=$((failed + 1))
    fi
done
if [ "$failed" -gt 0 ]; then
    echo "$failed of $((last - first + 1)) replicates failed"
    exit 1
fi
echo "Replicates $first..$last finished"
