"""Tests for non-equilibrium FEP input generation and log parsing.

These exercise the file-generation logic of QligFEP in NEQ mode without needing the compiled
Q binaries: a QligFEP instance is built directly and the writers are pointed at a temp dir.
"""

import pytest

from QligFEP.analyze_neq import collect_works, read_final_work
from QligFEP.qligfep import QligFEP


def make_run(**overrides):
    """Build a QligFEP instance with only the attributes the NEQ writers need."""
    run = QligFEP.__new__(QligFEP)
    run.timestep = "2fs"
    run.atomoffset = 0
    run.sphereradius = "25"
    run.system = "water"
    run.ABS = False
    run.dr_force = 0.5
    run.neq = True
    run.neq_reps = 3
    run.neq_steps = 20000
    run.neq_eq_steps = 500
    run.neq_relax_steps = 2500
    run.neq_L = 8.0
    run.neq_schedule = "sigmoidal"
    run.cluster = "SNELLIUS"
    run.replicates = "4"
    run.temperature = "298"
    run.seeds = [11, 22, 33, 44]
    run.lig1 = "lig1"
    run.lig2 = "lig2"
    run.to_clean = None
    for key, value in overrides.items():
        setattr(run, key, value)
    return run


@pytest.fixture
def neq_inputs(tmp_path):
    run = make_run()
    overlapping = [(1, 2), (3, 4)]
    files = run.write_MD_neq(str(tmp_path), lig_size1=10, lig_size2=12, overlapping_atoms=overlapping)
    return run, tmp_path, files


def test_writes_equilibration_endpoint_and_switch_files(neq_inputs):
    _, tmp_path, files = neq_inputs
    for name in [
        "eq1.inp",
        "eq2.inp",
        "eq3.inp",
        "eq4.inp",
        "eq5.inp",
        "eq6_0.inp",
        "eq6_1.inp",
        "neq_0.inp",
        "neq_1.inp",
    ]:
        assert (tmp_path / name).exists(), f"missing {name}"
        assert name in files


def test_writes_endpoint_relaxation_files(neq_inputs):
    _, tmp_path, files = neq_inputs
    for name in ["relax_0.inp", "relax_1.inp"]:
        assert (tmp_path / name).exists(), f"missing {name}"
        assert name in files


def test_relaxation_files_use_relax_steps_endpoint_lambda_no_scaling(neq_inputs):
    # The one-time endpoint relaxation is a plain (non-switching) equilibration at the
    # endpoint lambda, run for the longer neq_relax_steps (here 2500) rather than the tEQ
    # spacing of eq6 (here 500).
    _, tmp_path, _ = neq_inputs
    for state, lam in [("0", "0.000 1.000"), ("1", "1.000 0.000")]:
        text = (tmp_path / f"relax_{state}.inp").read_text()
        assert "steps                     2500" in text
        assert f"[lambdas]\n{lam}" in text
        assert "[lambda_scaling]" not in text
        for placeholder in ["RESTARTFILE", "FINALFILE", "T_VAR"]:
            assert placeholder in text


def test_runfile_relaxes_first_iteration_then_uses_spacing(tmp_path):
    run = make_run()
    run.write_MD_neq(str(tmp_path), 10, 12, [(1, 2)])
    run.write_neq_runfile(str(tmp_path), [])
    script = (tmp_path / "runSNELLIUS.sh").read_text()
    # the relaxation inputs must be staged into the per-replicate rundir
    assert "cp $inputfiles/relax_*.inp" in script
    # the first endpoint iteration runs the longer relaxation; later ones the tEQ spacing
    assert "relax_${s}.inp" in script
    assert "eq6_${s}.inp" in script


def test_no_windowed_md_files_in_neq_mode(neq_inputs):
    _, tmp_path, _ = neq_inputs
    # NEQ must not emit the ~100 fixed-lambda md_XXXX_YYYY.inp windows of equilibrium FEP.
    assert not list(tmp_path.glob("md_0*.inp"))
    assert not list(tmp_path.glob("md_1*.inp"))


def test_endpoint_lambdas(neq_inputs):
    _, tmp_path, _ = neq_inputs
    assert "[lambdas]\n0.000 1.000" in (tmp_path / "eq6_0.inp").read_text()
    assert "[lambdas]\n1.000 0.000" in (tmp_path / "eq6_1.inp").read_text()
    assert "[lambdas]\n0.000 1.000" in (tmp_path / "neq_0.inp").read_text()
    assert "[lambdas]\n1.000 0.000" in (tmp_path / "neq_1.inp").read_text()


def test_only_switch_files_carry_lambda_scaling(neq_inputs):
    _, tmp_path, _ = neq_inputs
    for switch in ["neq_0.inp", "neq_1.inp"]:
        text = (tmp_path / switch).read_text()
        assert "[lambda_scaling]" in text
        assert "scaling_parameter          sigmoidal" in text
        assert "L_sigmoid        8.0" in text
    for endpoint in ["eq6_0.inp", "eq6_1.inp"]:
        assert "[lambda_scaling]" not in (tmp_path / endpoint).read_text()


def test_step_counts_and_restraints(neq_inputs):
    _, tmp_path, _ = neq_inputs
    assert "steps                     20000" in (tmp_path / "neq_0.inp").read_text()
    assert "steps                     500" in (tmp_path / "eq6_0.inp").read_text()
    # overlapping atoms injected as distance restraints with the configured force constant
    assert "1 2 0.0 0.1 0.5 0" in (tmp_path / "neq_0.inp").read_text()


def test_per_replicate_placeholders_left_for_runscript(neq_inputs):
    _, tmp_path, _ = neq_inputs
    text = (tmp_path / "neq_0.inp").read_text()
    for placeholder in ["RESTARTFILE", "FINALFILE", "T_VAR", "FEP_VAR"]:
        assert placeholder in text


def test_restart_final_placeholders_survive_temperature_sed(neq_inputs):
    # The run script runs `sed s/T_VAR/<temp>/` over every .inp; the restart/final
    # placeholders must not contain "T_VAR" as a substring or they get corrupted
    # (regression: "RESTART_VAR" -> "RESTAR298").
    _, tmp_path, _ = neq_inputs
    for name in ["eq6_0.inp", "neq_0.inp"]:
        lines = (tmp_path / name).read_text().splitlines()
        restart_line = next(line for line in lines if line.startswith("restart"))
        final_line = next(line for line in lines if line.startswith("final"))
        assert "T_VAR" not in restart_line and "T_VAR" not in final_line
        assert restart_line.split()[1] == "RESTARTFILE"
        assert final_line.split()[1] == "FINALFILE"


def test_files_section_declares_trajectory(neq_inputs):
    # Q aborts with "Invalid data in input file" if [intervals] enables trajectory output
    # but [files] has no trajectory entry (regression).
    _, tmp_path, _ = neq_inputs
    for name in ["eq6_0.inp", "neq_0.inp"]:
        files_section = (tmp_path / name).read_text().split("[files]")[1].split("\n[")[0]
        assert any(line.startswith("trajectory") for line in files_section.splitlines())


def test_neq_runfile_parallelizes_switches(tmp_path):
    run = make_run()
    run.write_MD_neq(str(tmp_path), 10, 12, [(1, 2)])
    run.write_neq_runfile(str(tmp_path), [])
    script = (tmp_path / "runSNELLIUS.sh").read_text()
    # equilibration uses the parallel engine across all cores; switches use the serial engine
    assert "qdynp=" in script and "$qdynp " in script  # qdynp (parallel) for equilibration
    assert "qdyn=" in script and "$qdyn " in script  # qdyn (serial) for switches
    # switches are packed one-per-core via mpirun binding (Snellius bills the whole node)
    assert "mpirun" in script
    assert "--bind-to core" in script
    assert "#SBATCH --ntasks-per-node=16" in script  # use the billed cores
    assert "neq_reps=3" in script
    assert "#SBATCH --array=1-4" in script
    assert "qfep" not in script  # NEQ uses BAR, not the windowed qfep step


def test_read_final_work_parses_engine_log_format(tmp_path):
    # Mirrors the qdyn NEQ-mode output: work value at split index 6; completion line at the end.
    log = tmp_path / "neq_1_0.log"
    log.write_text(
        "Initializing dynamics\n"
        "At step 19990, work accumulated was   0.3100000000000000E+01 and dU and dlambda were   0.1E+00   0.1E-03\n"
        "At step 20000, work accumulated was   0.3210000000000000E+01 and dU and dlambda were   0.1E+00   0.1E-03\n"
        "qdyn version 6.0.1 terminated normally.\n"
    )
    assert read_final_work(str(log)) == pytest.approx(3.21)


def test_read_final_work_rejects_incomplete_log(tmp_path):
    log = tmp_path / "neq_1_0.log"
    log.write_text(
        "At step 10000, work accumulated was   0.5000000000000000E+01 and dU and dlambda were 0.1 0.1\n"
        "ABNORMAL TERMINATION of qdyn\n"
    )
    assert read_final_work(str(log)) is None


def test_collect_works_separates_forward_and_reverse(tmp_path):
    def write_log(name, work):
        (tmp_path / name).write_text(
            f"At step 20000, work accumulated was   {work:.16E} and dU and dlambda were 0.1 0.1\n"
            "qdyn version 6.0.1 terminated normally.\n"
        )

    write_log("neq_1_0.log", 2.0)  # forward
    write_log("neq_1_1.log", 2.5)  # forward
    write_log("neq_0_0.log", -1.0)  # reverse
    forward, reverse = collect_works(str(tmp_path))
    assert sorted(forward) == pytest.approx([2.0, 2.5])
    assert sorted(reverse) == pytest.approx([-1.0])
