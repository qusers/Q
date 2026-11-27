#! /bin/bash
#
# simple-vm: QmapFEP on one machine reached over SSH, with no scheduler.
#
# The run script is the qmapfep-slurm one, unchanged: replicates, sentinels, NaN
# checks, several replicates per task sharing a GPU under MPS. What a scheduler
# would do around it, this script does for one leg, plainly: it runs the leg's
# array tasks one after another, giving each the environment SLURM would.
#
#   SLURM_ARRAY_TASK_ID  the task, 1..ceil(replicates / Q_REPLICATES_PER_TASK)
#   SLURM_JOB_ID         Q_JOB_ID when the orchestrator sets one, else the time
#   CUDA_VISIBLE_DEVICES Q_GPUS (such as 0) when set, for the run script's MPS
#                        start; a GPU the job is not given is a GPU it ignores
#
# Each task writes to slurm.run<task>.<host>.<job>.out in the leg directory, the
# name SLURM would give it, so collection and the analysis read it as they read
# a SLURM task's output. A signal ends this script and, through the process
# group, the task it is running, whose replicates then write their .failed.

workdir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE

if [ ! -f "$submitfile" ]; then
    echo "Error: no run script at $submitfile"
    exit 1
fi
cd "$workdir" || exit 1

# The replicate count is the length of the seeds array the run script was
# generated with.
seeds_line=$(grep -m1 '^seeds=(' "$submitfile")
seeds_line=${seeds_line#seeds=(}
seeds_line=${seeds_line%)}
replicates=$(wc -w <<< "$seeds_line")
if [ "$replicates" -lt 1 ]; then
    echo "Error: no seeds in $submitfile"
    exit 1
fi

per_task=${Q_REPLICATES_PER_TASK:-1}
if ! [[ "$per_task" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Q_REPLICATES_PER_TASK must be a positive whole number, got '$per_task'"
    exit 1
fi
tasks=$(( (replicates + per_task - 1) / per_task ))
job_id=${Q_JOB_ID:-$(date +%s)}
host=$(hostname -s)
if [ -n "${Q_GPUS:-}" ]; then
    export CUDA_VISIBLE_DEVICES=$Q_GPUS
fi

echo "Running $replicates replicates of $workdir as $tasks task(s), $per_task per task"
failed=0
for task in $(seq 1 "$tasks"); do
    if ! SLURM_ARRAY_TASK_ID=$task SLURM_JOB_ID=$job_id SLURM_NTASKS=$per_task \
        bash "$submitfile" > "$workdir/slurm.run${task}.${host}.${job_id}.out" 2>&1; then
        echo "Task $task had failed replicates"
        failed=$((failed + 1))
    fi
done

if [ "$failed" -gt 0 ]; then
    echo "$failed of $tasks task(s) had failed replicates"
    exit 1
fi
echo "All $tasks task(s) finished"
