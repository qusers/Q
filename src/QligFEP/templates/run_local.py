"""Local FEP run script template.

Generates bash scripts for running FEP simulations locally (without SLURM),
with proper error handling: stops on crash, writes FAILED.log, skips cleanup.
"""

from dataclasses import dataclass
from textwrap import dedent


@dataclass
class LocalRunConfig:
    """Configuration for generating a local FEP run script.

    Attributes:
        qdyn_path: Absolute path to qdyn (serial) or qdynp (MPI) binary.
        qfep_path: Absolute path to qfep binary.
        use_mpi: If True, run qdyn via mpirun (LOCALP). If False, run directly (LOCAL).
        ntasks: Number of MPI tasks (only used when use_mpi=True).
        temperatures: List of temperature strings, e.g. ["298"].
        seeds: Random seeds, one per replicate.
        eq_files: Equilibration stage basenames, e.g. ["eq1", "eq2", ..., "eq5"].
        md_files: Production MD basenames in execution order.
        fep_files: FEP file names, e.g. ["FEP1.fep"].
        final_md_restart: Restart file from final MD for multi-FEP continuation.
        cleanup_patterns: File extensions to remove on success, e.g. ["dcd"]. None to skip.
    """

    qdyn_path: str
    qfep_path: str
    use_mpi: bool
    ntasks: int
    temperatures: list[str]
    seeds: list[int]
    eq_files: list[str]
    md_files: list[str]
    fep_files: list[str]
    final_md_restart: str
    cleanup_patterns: list[str] | None


def _build_step_block(stages: list[str], qdyn_cmd: str) -> str:
    """Build bash if/then blocks for each MD stage with error checking."""
    blocks = []
    for stage in stages:
        blocks.append(
            f"if [ $failed -eq 0 ]; then\n"
            f"    {qdyn_cmd} {stage}.inp > {stage}.log\n"
            f'    if [ $? -ne 0 ]; then echo "FAILED at {stage}" > FAILED.log; failed=1; fi\n'
            f"fi"
        )
    return "\n".join(blocks)


def render_local_run(config: LocalRunConfig) -> str:
    """Render a bash script for local FEP execution.

    The generated script:
    - Loops sequentially over all (temperature, replicate) combinations
    - Checks exit codes after each MD step
    - Writes FAILED.log and stops the MD pipeline on any failure
    - Only runs qfep analysis on success
    - Only runs cleanup on success
    - Always logs runtime information
    """
    qdyn_cmd = (
        f"mpirun -n {config.ntasks} {config.qdyn_path}"
        if config.use_mpi
        else config.qdyn_path
    )

    eq_block = _build_step_block(config.eq_files, qdyn_cmd)
    md_block = _build_step_block(config.md_files, qdyn_cmd)

    cleanup_block = ""
    if config.cleanup_patterns:
        patterns = " ".join(f"*{ext}" for ext in config.cleanup_patterns)
        cleanup_block = (
            f"if [ $failed -eq 0 ]; then\n"
            f"    rm -r {patterns}\n"
            f"fi\n"
            f"\n"
        )

    # dedent strips the Python-level indentation from the template, then
    # .replace() injects multi-line blocks without breaking the indentation.
    # (f-string multi-line substitutions break dedent because injected lines
    # have 0 indent, making the common prefix 0.)
    template = dedent(f"""\
    #!/bin/bash
    set -uo pipefail

    temperatures=({" ".join(config.temperatures)})
    seeds=({" ".join(str(s) for s in config.seeds)})
    fepfiles=({" ".join(config.fep_files)})
    qdyn="{config.qdyn_path}"
    qfep="{config.qfep_path}"
    finalMDrestart={config.final_md_restart}

    workdir="$(cd -P "$(dirname "${{BASH_SOURCE[0]}}")/.." && pwd)"
    inputfiles=$workdir/inputfiles

    global_starttime=$(date +%s)

    length=${{#fepfiles[@]}}
    length=$((length-1))

    for temperature in ${{temperatures[@]}}; do

    for run_idx in $(seq 0 $((${{#seeds[@]}}-1))); do
    seed=${{seeds[$run_idx]}}
    run_num=$((run_idx+1))

    starttime=$(date +%s)
    starttime_readable=$(date)

    for index in $(seq 0 $length); do
    fepfile=${{fepfiles[$index]}}
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

    cp $inputfiles/md*.inp .
    cp $inputfiles/*.top .
    cp $inputfiles/qfep.inp .
    cp $inputfiles/$fepfile .

    if [ $index -lt 1 ]; then
    cp $inputfiles/eq*.inp .
    sed -i "" "s/SEED_VAR/$seed/" eq1.inp
    else
        lastfep=FEP$index
        cp $workdir/$lastfep/$temperature/$run_num/$finalMDrestart $rundir/eq5.re
    fi
    sed -i "" "s/T_VAR/$temperature/" *.inp
    sed -i "" "s/FEP_VAR/$fepfile/" *.inp

    failed=0

    if [ $index -lt 1 ]; then
    __EQ_BLOCK__
    fi

    __MD_BLOCK__

    if [ $failed -eq 0 ]; then
        $qfep < qfep.inp > qfep.out
    fi

    __CLEANUP_BLOCK__\
    done

    endtime=$(date +%s)
    endtime_readable=$(date)
    runtime=$((endtime - starttime))
    hours=$((runtime / 3600))
    minutes=$((runtime % 3600 / 60))
    seconds=$((runtime % 60))

    echo "#    Replicate $run_num complete"
    echo "#    Starttime: $starttime_readable"
    echo "#    Endtime: $endtime_readable"
    echo "#    Runtime: ${{hours}}h:${{minutes}}m:${{seconds}}s"
    echo "#    Random seed: $seed"
    echo "#    Replicate Number: $run_num"
    echo "#    Working Directory: $workdir"
    echo

    done
    done

    global_endtime=$(date +%s)
    global_runtime=$((global_endtime - global_starttime))
    global_hours=$((global_runtime / 3600))
    global_minutes=$((global_runtime % 3600 / 60))
    global_seconds=$((global_runtime % 60))
    echo "#    Total Runtime: ${{global_hours}}h:${{global_minutes}}m:${{global_seconds}}s"
    """)

    return (
        template
        .replace("__EQ_BLOCK__", eq_block)
        .replace("__MD_BLOCK__", md_block)
        .replace("__CLEANUP_BLOCK__", cleanup_block)
    )
