"""Short engine regressions; build qdyn/qdynp with src/q6/makefile first.

The existing benzene fixtures supply geometry/bonds. Their LJ tables are
converted to arithmetic Rmin/2, epsilon for a synthetic neutral annihilation.
These are numerical regressions, not a solvation free-energy benchmark.
"""

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from QligFEP.analyze_fhe import read_fhe

ROOT = Path(__file__).resolve().parents[2]


def arithmetic_topology(source):
    text = source.read_text()
    text = re.sub(r"1 = vdW combination rule", "2 = vdW combination rule", text)
    pattern = r"sqrt \(Aii\) ([^:]+):\n(.*?)sqrt \(Bii\) \1:\n(.*?)(?=sqrt \(Aii\)|\s+0 = No. of type-2)"

    def convert(match):
        a = np.fromstring(match[2], sep=" ")
        b = np.fromstring(match[3], sep=" ")
        radius, epsilon = np.zeros_like(a), np.zeros_like(b)
        real = (a > 0) & (b > 0)
        radius[real] = (2 * (a[real] / b[real]) ** 2) ** (1 / 6) / 2
        epsilon[real] = b[real] ** 4 / (4 * a[real] ** 2)
        return (f"Rmin/2 {match[1]}:\n" + " ".join(map(str, radius)) + "\n"
                + f"epsilon {match[1]}:\n" + " ".join(map(str, epsilon)) + "\n")

    text, count = re.subn(pattern, convert, text, flags=re.S)
    assert count == 3
    return text


@pytest.fixture
def engine_case(tmp_path):
    def setup(solvent=True, source=0.5, mode="validate", constraints="lincs settle"):
        top = ROOT / "test/data/topology" / ("benzene-water.top" if solvent else "benzene-vacuum.top")
        (tmp_path / "test.top").write_text(arithmetic_topology(top))
        (tmp_path / "test.fep").write_text("""[FEP]
states 2
softcore_method gapsys_v2
softcore_use_max_potential off
[atoms]
1 1
2 4
3 7
4 10
[atom_types]
C 1.9 0.1 0 0 1.9 0.05 12.011
H 0 0 0 0 0 0 1.008
DUM 0 0 0 0 0 0 1.008
[change_atoms]
1 C DUM
2 C DUM
3 H DUM
4 H DUM
[change_charges]
1 -0.115 0
2 -0.115 0
3 0.115 0
4 0.115 0
[softcore]
1 1 0
2 1 0
3 1 0
4 1 0
""")
        (tmp_path / "test.inp").write_text(f"""[MD]
steps 6
stepsize 0.01
temperature 1
initial_temperature 1
random_seed 1
bath_coupling 1
constrain_solvent on
constrain_hydrogens on
constraint_algorithm {constraints}
lrf on
foreign_hamiltonian_mode {mode}
[cut-offs]
solute_solute 10
solute_solvent 10
solvent_solvent 10
q_atom 99
lrf 99
[sphere]
shell_radius 18
shell_force 10
[intervals]
non_bond 2
output 2
energy 2
[files]
topology test.top
fep test.fep
final test.re
energy test.en
[lambdas]
{1-source:.12f} {source:.12f}
[foreign_lambdas]
0
0.0001
0.5
0.9999
1
""")
        return tmp_path
    return setup


def run_engine(directory, mpi=False):
    binary = ROOT / "src/q6" / ("qdynp" if mpi else "qdyn")
    if not binary.exists():
        pytest.skip(f"build {binary.name} first")
    command = [str(binary), "test.inp"]
    if mpi:
        if not shutil.which("mpirun"):
            pytest.skip("mpirun unavailable")
        command = ["mpirun", "-n", "2"] + command
    result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=60)
    (directory / ("mpi.log" if mpi else "serial.log")).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout[-5000:] + result.stderr
    assert "terminated normally" in result.stdout, result.stdout[-5000:]
    return read_fhe(directory / "test.en.fhe")


@pytest.mark.parametrize("solvent,source", [(False, 0.5), (True, 0), (True, 0.0001), (True, 0.5), (True, 0.9999), (True, 1)])
def test_correction_matches_complete_hamiltonian(engine_case, solvent, source):
    samples = run_engine(engine_case(solvent=solvent, source=source))
    assert samples.values.shape == (2, 5)
    assert np.isfinite(samples.values).all()
    np.testing.assert_allclose(samples.work_to(source), 0, atol=1e-7)


def test_validation_does_not_change_dynamics(engine_case):
    directory = engine_case(mode="correction")
    correction = run_engine(directory)
    restart = (directory / "test.re").read_bytes()
    engine_case(mode="validate")
    validated = run_engine(directory)
    assert restart == (directory / "test.re").read_bytes()
    np.testing.assert_array_equal(correction.values, validated.values)


@pytest.mark.skipif(os.environ.get("Q_TEST_MPI") != "1", reason="set Q_TEST_MPI=1 to run MPI regression")
def test_mpi_matches_serial(engine_case):
    directory = engine_case()
    serial = run_engine(directory)
    mpi = run_engine(directory, mpi=True)
    np.testing.assert_allclose(mpi.values, serial.values, atol=1e-7, rtol=1e-7)


@pytest.mark.parametrize("entry", ["nan", "1.1", "0.5"])
def test_invalid_foreign_targets_fail(engine_case, entry):
    directory = engine_case()
    path = directory / "test.inp"
    path.write_text(path.read_text() + entry + "\n")
    binary = ROOT / "src/q6/qdyn"
    if not binary.exists():
        pytest.skip("build qdyn first")
    result = subprocess.run([str(binary), "test.inp"], cwd=directory, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert "[foreign_lambdas] entry" in result.stdout


def test_qfep_rejects_linear_reconstruction_of_foreign_work(engine_case):
    directory = engine_case(solvent=False)
    run_engine(directory)
    binary = ROOT / "src/q6/qfep"
    if not binary.exists():
        pytest.skip("build qfep first")
    inputs = "1\n2 0\n0.592 0\n1\n1\n1\n0\n0\n1 0\ntest.en\n"
    result = subprocess.run([str(binary)], input=inputs, cwd=directory,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert "Nonlinear foreign work found" in result.stderr
