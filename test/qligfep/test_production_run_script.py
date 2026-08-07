"""Tests for thin platform adapters generated for the Python runner."""

from QligFEP.qligfep import QligFEP
from QligFEP.templates.production_run import (
    ProductionRunScriptConfig,
    render_production_run_script,
)


def _config(**overrides):
    values = {
        "qdyn": "/opt/Q/qdynp",
        "qfep": "/opt/Q/qfep",
        "ntasks": 16,
        "temperatures": ("298", "310"),
        "seeds": (1234, 5678),
        "fep_file": "FEP1.fep",
        "job_name": "p_lig1_lig2",
        "modules": "module purge\nmodule load OpenMPI",
        "account": "project1",
        "partition": "compute",
    }
    values.update(overrides)
    return ProductionRunScriptConfig(**values)


def test_slurm_script_only_adapts_array_parameters_to_runner():
    script = render_production_run_script(_config())

    assert "#SBATCH --array=1-4" in script
    assert "#SBATCH --ntasks-per-node=16" in script
    assert "#SBATCH -A project1" in script
    assert "#SBATCH -p compute" in script
    assert "module load OpenMPI" in script
    assert 'runner=${QLIGFEP_RUNNER:-qligfep-run}' in script
    assert '--input-dir "$inputfiles"' in script
    assert '--work-dir "$rundir"' in script
    assert '--launcher "mpirun -n $SLURM_NTASKS --map-by core --bind-to core"' in script
    assert 'rundir="$workdir/FEP1/$temperature/$run_num"' in script
    assert "TMPDIR" not in script
    assert "cp " not in script
    assert "\nrm " not in script


def test_serial_local_script_loops_without_scheduler_or_mpi():
    script = render_production_run_script(
        _config(slurm=False, use_mpi=False, ntasks=1, modules="", account=None, partition=None)
    )

    assert "#SBATCH" not in script
    assert 'for temperature in "${temperatures[@]}"' in script
    assert 'for run_index in "${!seeds[@]}"' in script
    assert "--launcher" not in script
    assert '--qdyn "$qdyn"' in script
    assert "TMPDIR" not in script


def test_qligfep_dispatches_production_runfile_generation(tmp_path):
    run = object.__new__(QligFEP)
    run.production = True
    run.cluster = "LOCAL"
    run.system = "water"
    run.lig1 = "lig1"
    run.lig2 = "lig2"
    run.temperature = "298"
    run.seeds = [1234, 5678]

    run.write_runfile(str(tmp_path), file_list=[])

    script = (tmp_path / "runLOCAL.sh").read_text()
    assert "qligfep-run" in script
    assert "temperatures=(298)" in script
    assert "seeds=(1234 5678)" in script
    assert "cp " not in script


def test_parallel_local_script_uses_configured_process_count():
    script = render_production_run_script(
        _config(slurm=False, use_mpi=True, ntasks=8, modules="", account=None, partition=None)
    )

    assert '--launcher "mpirun -n 8"' in script
    assert "SLURM" not in script
