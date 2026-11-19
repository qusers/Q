"""Engine regressions for charge-changing nonlinear soft-core FEP.

Build qdyn/qdynp with src/q6/makefile first. The benzene case from
test_foreign_work is turned into a synthetic 0 -> -1 e perturbation.
"""

import os
import re
import subprocess

import numpy as np
import pytest

from test_foreign_work import ROOT, engine_case, run_engine  # noqa: F401

BORN_SOLVENT = {
    "polarization": "on",
    "charge_correction": "off",
    "polarization_adaptation": "off",
    "polarization_zero_offsets": "on",
    "perstate_born_correction": "on",
    "born_dielectric": "80.0",
}


def charge_changing_case(engine_case, solvent_settings, steps=6, restart=None, **kwargs):
    directory = engine_case(**kwargs)
    fep = directory / "test.fep"
    fep.write_text(fep.read_text().replace("1 -0.115 0\n2 -0.115 0\n", "1 -0.115 -0.5\n2 -0.115 -0.5\n"))
    section = "[solvent]\n" + "".join(f"{key} {value}\n" for key, value in solvent_settings.items())
    text = (directory / "test.inp").read_text().replace("[intervals]", section + "[intervals]")
    text = text.replace("steps 6\n", f"steps {steps}\n")
    if restart:
        text = text.replace("initial_temperature 1\nrandom_seed 1\n", "")
        text = text.replace("final test.re", f"restart {restart}\nfinal test.re")
    (directory / "test.inp").write_text(text)
    return directory


def run_failing(directory):
    binary = ROOT / "src/q6/qdyn"
    if not binary.exists():
        pytest.skip("build qdyn first")
    result = subprocess.run([str(binary), "test.inp"], cwd=directory, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 or "terminated normally" not in result.stdout
    return result.stdout + result.stderr


def initial_offsets(log):
    return np.array(re.search(r"Initial polarization offsets \(rad\)\s*=(.*)", log)[1].split(), dtype=float)


@pytest.mark.parametrize("source", [0, 0.0001, 0.5, 0.9999, 1])
def test_born_term_is_carried_by_foreign_work(engine_case, source):
    directory = charge_changing_case(engine_case, BORN_SOLVENT, source=source)
    samples = run_engine(directory)
    assert np.isfinite(samples.values).all()
    np.testing.assert_allclose(samples.work_to(source), 0, atol=1e-7)

    log = (directory / "serial.log").read_text()
    radius = float(re.search(r"Born solvent radius\s*=\s*([\d.]+)", log)[1])
    coefficient = float(re.search(r"Born coefficient\s*=\s*([\d.]+)", log)[1])
    np.testing.assert_allclose(coefficient, 332.0637 * (1 - 1 / 80) / (2 * radius), rtol=1e-4)
    states = re.findall(r"Boundary state\s+\d+ Q-region charge =\s*(\S+) total charge =\s*(\S+)", log)
    region = [float(q) for q, _ in states]
    total = [float(q) for _, q in states]
    np.testing.assert_allclose(region[1] - region[0], -1, atol=1e-4)
    difference = float(re.search(r"Born state difference E\(2\)-E\(1\)\s*=\s*(\S+)", log)[1])
    np.testing.assert_allclose(difference, -coefficient * (total[1] ** 2 - total[0] ** 2), atol=1e-3)


def born_difference(directory):
    log = (directory / "serial.log").read_text()
    return float(re.search(r"Born state difference E\(2\)-E\(1\)\s*=\s*(\S+)", log)[1])


def test_born_term_shifts_endpoint_work(engine_case):
    """The Born term is a per-state constant, so it changes work but not dynamics."""
    directory = charge_changing_case(engine_case, BORN_SOLVENT, source=0)
    high = run_engine(directory).work_to(1)
    high_difference = born_difference(directory)
    charge_changing_case(engine_case, {**BORN_SOLVENT, "born_dielectric": "2.0"}, source=0)
    low = run_engine(directory).work_to(1)
    np.testing.assert_allclose(high - low, high_difference - born_difference(directory), atol=2e-3)


@pytest.mark.parametrize("settings,message", [
    ({**BORN_SOLVENT, "perstate_born_correction": "off"}, "requires perstate_born_correction on"),
    ({**BORN_SOLVENT, "charge_correction": "on"}, "requires charge_correction off"),
    ({**BORN_SOLVENT, "born_dielectric": "1.0"}, "born_dielectric must be > 1"),
    ({**BORN_SOLVENT, "polarization_adaptation": "on"}, "requires polarization_adaptation off"),
    ({**BORN_SOLVENT, "polarization": "off"}, "requires polarization on"),
    ({**BORN_SOLVENT, "polarization_zero_offsets": "maybe"}, "polarization_zero_offsets must be on or off"),
    ({**BORN_SOLVENT, "polarization_adaptation": "maybe"}, "polarization_adaptation must be on or off"),
])
def test_unsupported_settings_fail(engine_case, settings, message):
    assert message in run_failing(charge_changing_case(engine_case, settings))


def test_zero_offsets_override_restart_and_stay_frozen(engine_case):
    adaptive = {**BORN_SOLVENT, "polarization_adaptation": "on", "polarization_zero_offsets": "off"}
    directory = charge_changing_case(engine_case, adaptive, steps=200, mode="correction")
    run_engine(directory)
    (directory / "adapted.re").write_bytes((directory / "test.re").read_bytes())

    charge_changing_case(engine_case, adaptive, restart="adapted.re", mode="correction")
    run_engine(directory)
    inherited = initial_offsets((directory / "serial.log").read_text())
    assert np.any(inherited != 0)

    charge_changing_case(engine_case, BORN_SOLVENT, steps=200, restart="adapted.re", mode="correction")
    run_engine(directory)
    log = (directory / "serial.log").read_text()
    assert "ZERO (restart offsets overridden)" in log
    assert "OFF (frozen)" in log
    assert np.all(initial_offsets(log) == 0)
    blocks = re.findall(r"shell\s+<n>\s+<theta>\s+theta_0\s+theta_corr\n((?:\s*\d+\s[^\n]*\n)+)", log)
    assert blocks
    assert all(float(row.split()[-1]) == 0 for block in blocks for row in block.splitlines())


@pytest.mark.parametrize("minimizer", ["fire", "steepest_descent"])
def test_minimization_does_not_adapt_offsets(engine_case, minimizer):
    adaptive = {**BORN_SOLVENT, "polarization_adaptation": "on", "polarization_zero_offsets": "off"}
    directory = charge_changing_case(engine_case, adaptive, mode="correction")
    inp = directory / "test.inp"
    inp.write_text(inp.read_text().replace(
        "lrf on\n", f"lrf on\nminimize on\nminimize_method {minimizer}\nmax_minimize_steps 250\n", 1))
    run_engine(directory)
    log = (directory / "serial.log").read_text()
    assert "Energy minimization enabled" in log
    assert "Water polarization restraint data" not in log

    charge_changing_case(engine_case, adaptive, restart="test.re", mode="correction")
    run_engine(directory)
    assert np.all(initial_offsets((directory / "serial.log").read_text()) == 0)


@pytest.mark.skipif(os.environ.get("Q_TEST_MPI") != "1", reason="set Q_TEST_MPI=1 to run MPI regression")
def test_mpi_matches_serial(engine_case):
    directory = charge_changing_case(engine_case, BORN_SOLVENT)
    serial = run_engine(directory)
    mpi = run_engine(directory, mpi=True)
    np.testing.assert_allclose(mpi.values, serial.values, atol=1e-7, rtol=1e-7)
