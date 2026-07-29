"""Compile and run the pure spherical-boundary helper regression test."""

from pathlib import Path
import shutil
import subprocess

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_SOURCE = PROJECT_ROOT / "src" / "q6" / "boundary_corrections.f90"
DRIVER_SOURCE = Path(__file__).with_suffix(".f90")


def test_boundary_correction_helpers(tmp_path):
    compiler = shutil.which("gfortran-11") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("No gfortran compiler available")

    executable = tmp_path / "test_boundary_corrections"
    subprocess.run(
        [
            compiler,
            "-std=legacy",
            "-ffree-line-length-none",
            str(MODULE_SOURCE),
            str(DRIVER_SOURCE),
            "-o",
            str(executable),
        ],
        cwd=tmp_path,
        check=True,
    )
    result = subprocess.run([str(executable)], capture_output=True, text=True, check=False)

    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASS: spherical-boundary correction helpers" in result.stdout
