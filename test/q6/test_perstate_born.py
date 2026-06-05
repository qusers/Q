"""Pytest wrapper for the Fortran-level per-state Born self-energy regression test.

The Fortran driver lives at test/q6/test_perstate_born.f90 and is linked against
the qatom module's pure helpers (born_coefficient, born_self_energy). It asserts
the continuum Born coefficient C = k_e(1-1/eps)/(2R), the self-energy -C*Q^2, the
per-state difference -C*(Q2^2-Q1^2) with the correct sign across both signs of
in-sphere charge, and exact cancellation on neutral edges.
"""
from pathlib import Path
import subprocess
import shutil
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
Q6_SRC = PROJECT_ROOT / "src" / "q6"
TEST_DRIVER_SRC = Path(__file__).with_suffix(".f90")
TEST_BIN = Path(__file__).with_suffix("")  # binary alongside .py


def _find_fortran() -> str:
    for cand in ("gfortran-11", "gfortran"):
        if shutil.which(cand):
            return cand
    pytest.skip("No gfortran available for per-state Born regression test")


def _build():
    """Build the Fortran test driver against the existing src/q6 module objects.

    Requires src/q6 to have been compiled (qatom.mod / .o etc. present either in
    src/q6/ directly or under src/q6/obj/).
    """
    fc = _find_fortran()
    needed_objs = [
        "qatom.o", "misc.o", "nrgy.o", "prmfile.o", "sizes.o",
        "index.o", "topo.o", "mask.o", "mpiglob.o", "trj.o",
    ]
    obj_dir = None
    for cand in (Q6_SRC, Q6_SRC / "obj"):
        if all((cand / o).exists() for o in needed_objs):
            obj_dir = cand
            break
    if obj_dir is None:
        pytest.skip("src/q6 not built; objects not found in src/q6/ or src/q6/obj/")

    obj = TEST_DRIVER_SRC.with_suffix(".o")
    flags = ["-O2", "-ffast-math", "-funroll-loops", "-std=legacy", "-w",
             "-fPIC", "-cpp", f"-I{obj_dir}"]
    subprocess.check_call(
        [fc, *flags, "-c", str(TEST_DRIVER_SRC), "-o", str(obj)],
        cwd=obj_dir,
    )
    subprocess.check_call(
        [fc, *flags, str(obj), *needed_objs, "-o", str(TEST_BIN)],
        cwd=obj_dir,
    )


def test_perstate_born_helpers():
    """C = k_e(1-1/eps)/(2R); E_born = -C*Q^2; ddG_Born = -C*(Q2^2-Q1^2) with
    correct sign for both charge signs; neutral edges cancel exactly."""
    _build()
    result = subprocess.run([str(TEST_BIN)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, (
        f"per-state Born regression test failed:\nstdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "PASS: all per-state Born" in result.stdout, result.stdout
