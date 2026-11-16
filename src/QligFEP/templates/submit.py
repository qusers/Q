"""Job script pieces that depend on the execution profile.

The run scripts in ``INPUTS/`` are scheduler independent: they read the task index and
core count from ``QLIGFEP_TASK_ID`` / ``QLIGFEP_NTASKS`` and fall back to the Slurm
array variables. What differs between execution types is only how those tasks are
started, which is what ``FEP_submit.sh`` does.
"""

from __future__ import annotations

from ..config.schema import Profile

SLURM_SUBMIT = """\
#!/bin/bash
# Submit this FEP as a Slurm job array (profile: {profile}).
#   ./FEP_submit.sh          submit every task
#   ./FEP_submit.sh 1 3      submit only array tasks 1 and 3
cd "$(cd -P "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)" || exit 1
submitfile=inputfiles/{runfile}

if [ $# -eq 0 ]; then # if no arguments provided, submit the whole job array
    echo "No arguments provided - submitting full job array (all replicates)"
    sbatch $submitfile
else
    array_indices=$(echo "$@" | tr ' ' ',') # convert space-separated arguments to comma-separated
    echo "Submitting job array indexes: $array_indices"
    sbatch --array=$array_indices $submitfile
fi
"""

LOCAL_SUBMIT = """\
#!/bin/bash
# Run this FEP on the local machine (profile: {profile}), without a scheduler.
#   ./FEP_submit.sh          run every task
#   ./FEP_submit.sh 1 3      run only tasks 1 and 3
# Each task writes slurm.run<task>.local.0.out, the same log the Slurm jobs produce,
# so the analysis tools read local and cluster runs alike.
cd "$(cd -P "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)" || exit 1
runfile=inputfiles/{runfile}
total_tasks={total_tasks}
max_parallel=${{QLIGFEP_MAX_PARALLEL:-{max_parallel}}}
export QLIGFEP_NTASKS=${{QLIGFEP_NTASKS:-{ntasks}}}

if [ $# -eq 0 ]; then
    tasks=$(seq 1 $total_tasks)
else
    tasks="$*"
fi

echo "Running tasks: $(echo $tasks | tr '\\n' ' ')"
echo "$max_parallel at a time, $QLIGFEP_NTASKS core(s) each"
for task in $tasks; do
    while [ "$(jobs -rp | wc -l)" -ge "$max_parallel" ]; do
        sleep 1
    done
    echo "[$(date '+%H:%M:%S')] starting task $task"
    QLIGFEP_TASK_ID=$task bash "$runfile" > "slurm.run${{task}}.local.0.out" 2>&1 &
done
wait

failed=0
for task in $tasks; do
    if grep -q "^#    Exit status: 0" "slurm.run${{task}}.local.0.out"; then
        echo "task $task finished"
    else
        echo "task $task FAILED (see slurm.run${{task}}.local.0.out)"
        failed=1
    fi
done
exit $failed
"""


def render_submit_script(profile: Profile, total_tasks: int, ntasks: int | None = None) -> str:
    """Render ``FEP_submit.sh`` for the profile's execution type.

    Args:
        profile: Resolved execution profile.
        total_tasks: Number of tasks (replicates times temperatures).
        ntasks: Cores per task when it differs from the profile's ``ntasks``.
    """
    if profile.is_local:
        return LOCAL_SUBMIT.format(
            profile=profile.name,
            runfile=profile.runfile_name,
            total_tasks=total_tasks,
            max_parallel=profile.resources.max_parallel_replicates,
            ntasks=ntasks if ntasks is not None else profile.resources.ntasks,
        )
    return SLURM_SUBMIT.format(profile=profile.name, runfile=profile.runfile_name)


def engine_block(profile: Profile, method: str) -> tuple[str, str]:
    """Shell assignments for the Q binaries and the variables to check before running.

    Returns:
        ``(assignments, variables)``; ``variables`` is a space-separated list of quoted
        shell variables for the run script's executable check.
    """
    if method == "nonequilibrium":
        lines = [f"qdynp={profile.binary('qdynp')}", f"qdyn={profile.binary('qdyn')}"]
        return "\n".join(lines) + "\n", '"$qdynp" "$qdyn"'
    lines = [f"qdyn={profile.binary(profile.engine)}", f"qfep={profile.binary('qfep')}"]
    return "\n".join(lines) + "\n", '"$qdyn" "$qfep"'
