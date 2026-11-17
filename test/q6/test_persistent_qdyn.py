"""qdyn --serve: the generated reset routines are up to date, and the protocol works.

That a persistent process gives the same bytes as one process per stage is tested through
FortranEngine (src/core/python/test/test_fortran_engine.py, both modes) and the run plan
(test/qligfep/test_run_plan.py).
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
Q6 = ROOT / "src" / "q6"
QDYN = Q6 / "bin" / "q6" / "qdyn"


def test_reset_routines_cover_every_module_variable():
    """A module variable added without regenerating would leak from one stage into the next."""
    result = subprocess.run([sys.executable, str(Q6 / "tools" / "reset_state.py"), "--check"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.skipif(not QDYN.is_file(), reason=f"{QDYN} is not built")
def test_serve_protocol(tmp_path):
    """Commands, answers, errors, and one topology per process."""
    shutil.copy(ROOT / "test" / "data" / "topology" / "benzene-water.top", tmp_path / "top.top")
    shutil.copy(ROOT / "test" / "data" / "inputs" / "FEP_benzene.fep", tmp_path / "fep.fep")
    stage = ("[MD]\nsteps 10\nstepsize 0.5\ntemperature 1\nbath_coupling 1\nrandom_seed 7\n"
             "initial_temperature 1\n[files]\ntopology {top}\nfep fep.fep\n[lambdas]\n1.0 0.0\n")
    shutil.copy(tmp_path / "top.top", tmp_path / "other.top")
    commands = "\n".join([
        "save nothing.re",                          # no state yet
        "bogus",                                    # unknown command
        "run", stage.format(top="top.top"), "@end",
        "save first.re",
        "load first.re",
        "load missing.re",
        "run", stage.format(top="other.top"), "@end",  # another topology: refused
    ]) + "\n"
    result = subprocess.run([str(QDYN), "--serve", "."], input=commands, cwd=tmp_path,
                            capture_output=True, text=True)
    answers = [line for line in result.stdout.splitlines() if line.startswith("@qdyn")]
    assert answers == [
        "@qdyn ready",
        "@qdyn error no state to save",
        "@qdyn error unknown command: bogus",
        "@qdyn done",
        "@qdyn done",
        "@qdyn done",
        "@qdyn error cannot read missing.re",
    ]
    assert "runs one topology" in result.stdout and result.returncode != 0
    assert (tmp_path / "first.re").is_file() and (tmp_path / "recovery.re").is_file()
    assert (tmp_path / "recovery.re").read_bytes() == (tmp_path / "first.re").read_bytes()
