"""Minimal launch scripts for the platform-independent production runner."""

from __future__ import annotations

import shlex
from dataclasses import dataclass


@dataclass(frozen=True)
class ProductionRunScriptConfig:
    qdyn: str
    qfep: str
    ntasks: int
    temperatures: tuple[str, ...]
    seeds: tuple[int, ...]
    fep_file: str
    job_name: str
    use_mpi: bool = True
    slurm: bool = True
    nodes: int = 1
    time: str = "0-12:00:00"
    modules: str = ""
    account: str | None = None
    partition: str | None = None


def _runner_invocation(config: ProductionRunScriptConfig, launcher: str) -> str:
    arguments = [
        '"$runner"',
        '--input-dir "$inputfiles"',
        '--work-dir "$rundir"',
        '--temperature "$temperature"',
        '--seed "$seed"',
        f"--fep-file {shlex.quote(config.fep_file)}",
        '--qdyn "$qdyn"',
        '--qfep "$qfep"',
    ]
    if launcher:
        arguments.append(f'--launcher "{launcher}"')
    return " \\\n    ".join(arguments)


def render_production_run_script(config: ProductionRunScriptConfig) -> str:
    """Render a Slurm-array or sequential local launcher.

    The generated script only maps scheduler/local replicate parameters to the
    Python runner.  Scientific stage orchestration stays in tested Python code.
    """
    if not config.temperatures:
        raise ValueError("at least one temperature is required")
    if not config.seeds:
        raise ValueError("at least one seed is required")

    temperatures = " ".join(shlex.quote(value) for value in config.temperatures)
    seeds = " ".join(str(value) for value in config.seeds)
    common = f'''runner=${{QLIGFEP_RUNNER:-qligfep-run}}
qdyn={shlex.quote(config.qdyn)}
qfep={shlex.quote(config.qfep)}
temperatures=({temperatures})
seeds=({seeds})
inputfiles="$(cd -P "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
workdir="$(dirname "$inputfiles")"
'''

    if config.slurm:
        total_jobs = len(config.temperatures) * len(config.seeds)
        header = [
            "#!/bin/bash",
            f"#SBATCH --nodes={config.nodes}",
            f"#SBATCH --ntasks-per-node={config.ntasks}",
            "#SBATCH --mem-per-cpu=512",
            f"#SBATCH --time={config.time}",
            f"#SBATCH -J {config.job_name}",
            f"#SBATCH --array=1-{total_jobs}",
            "#SBATCH -o slurm.run%a.%N.%j.out",
        ]
        if config.account:
            header.insert(4, f"#SBATCH -A {config.account}")
        if config.partition:
            header.insert(4, f"#SBATCH -p {config.partition}")
        launcher = "mpirun -n $SLURM_NTASKS --map-by core --bind-to core" if config.use_mpi else ""
        invocation = _runner_invocation(config, launcher)
        return (
            "\n".join(header)
            + "\n\n"
            + config.modules.rstrip()
            + "\n\nset -eo pipefail\n\n"
            + common
            + '''
runs=${#seeds[@]}
task_index=$((SLURM_ARRAY_TASK_ID - 1))
temp_index=$((task_index / runs))
run_num=$((task_index % runs + 1))
temperature=${temperatures[$temp_index]}
seed=${seeds[$((run_num - 1))]}
rundir="$workdir/FEP1/$temperature/$run_num"
mkdir -p "$rundir"

starttime=$(date +%s)
starttime_readable=$(date)
write_footer() {
    rc=$?
    runtime=$(($(date +%s) - starttime))
    echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
    echo "#    Slurm tasks: $SLURM_NTASKS"
    echo "#    Starttime: $starttime_readable"
    echo "#    Endtime: $(date)"
    echo "#    Runtime: $((runtime / 3600))h:$((runtime % 3600 / 60))m:$((runtime % 60))s"
    echo "#    Random seed: $seed"
    echo "#    Replicate Number: $run_num"
    echo "#    Working Directory: $workdir"
    echo "#    Exit status: $rc"
}
trap write_footer EXIT

'''
            + invocation
            + "\n"
        )

    launcher = f"mpirun -n {config.ntasks}" if config.use_mpi else ""
    invocation = _runner_invocation(config, launcher)
    return (
        "#!/bin/bash\nset -eo pipefail\n\n"
        + common
        + '''
for temperature in "${temperatures[@]}"; do
    for run_index in "${!seeds[@]}"; do
        run_num=$((run_index + 1))
        seed=${seeds[$run_index]}
        rundir="$workdir/FEP1/$temperature/$run_num"
        mkdir -p "$rundir"
'''
        + "        "
        + invocation.replace("\n", "\n        ")
        + '''
    done
done
'''
    )
