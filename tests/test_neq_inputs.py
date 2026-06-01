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
    run.replacements = {}
    run.timestep = "2fs"
    run.set_timestep()
    run.atomoffset = 0
    run.sphereradius = "25"
    run.system = "water"
    run.ABS = False
    run.dr_force = 0.5
    run.neq = True
    run.neq_reps = 3
    run.neq_steps = 20000
    run.neq_eq_steps = 500
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
    for placeholder in ["RESTART_VAR", "FINAL_VAR", "T_VAR", "FEP_VAR"]:
        assert placeholder in text


def test_neq_runfile_uses_serial_qdyn_neq(tmp_path):
    run = make_run()
    run.write_MD_neq(str(tmp_path), 10, 12, [(1, 2)])
    run.write_neq_runfile(str(tmp_path), [])
    script = (tmp_path / "runSNELLIUS.sh").read_text()
    assert "$qdyn_neq" in script
    assert "qdyn_neq=" in script  # the QDYN_NEQ binary path was substituted in
    assert "mpirun" not in script  # qdyn_neq is serial
    assert "neq_reps=3" in script
    assert "#SBATCH --array=1-4" in script
    assert "qfep" not in script  # NEQ uses BAR, not the windowed qfep step


def test_read_final_work_parses_engine_log_format(tmp_path):
    # Mirrors the qdyn_neq output: work value at split index 6; completion line at the end.
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
