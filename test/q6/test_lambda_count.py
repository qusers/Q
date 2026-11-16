"""qdyn must refuse an input whose number of lambdas differs from the FEP file's states.

It used to run: the energies of the states without a lambda were weighted with
whatever was in memory behind the lambda array.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
QDYN = Path(os.environ.get("QDYN", PROJECT_ROOT / "src" / "q6" / "bin" / "q6" / "qdyn"))
TOPOLOGY = PROJECT_ROOT / "test" / "data" / "topology" / "benzene-water.top"
FEP = PROJECT_ROOT / "test" / "data" / "inputs" / "FEP_benzene.fep"  # two states

INPUT_TEMPLATE = """\
[MD]
steps                     10
stepsize                  1.0
temperature               300
random_seed               1
initial_temperature       300

[files]
topology                  top.top
fep                       fep.fep
final                     final.re

[lambdas]
{lambdas}
"""


def _run(tmp_path: Path, lambdas: str) -> subprocess.CompletedProcess:
    shutil.copy(TOPOLOGY, tmp_path / "top.top")
    shutil.copy(FEP, tmp_path / "fep.fep")
    (tmp_path / "md.inp").write_text(INPUT_TEMPLATE.format(lambdas=lambdas))
    return subprocess.run([str(QDYN), "md.inp"], cwd=tmp_path, text=True, capture_output=True, timeout=300)


@pytest.mark.skipif(not QDYN.is_file(), reason=f"{QDYN} is not built")
@pytest.mark.parametrize("lambdas", ["1.0", "0.3 0.3 0.4"])
def test_wrong_number_of_lambdas_is_an_error(tmp_path, lambdas):
    result = _run(tmp_path, lambdas)
    assert result.returncode != 0
    assert f"{len(lambdas.split())} lambda value(s) in the input file, but the FEP file has 2 state(s)" in result.stdout


@pytest.mark.skipif(not QDYN.is_file(), reason=f"{QDYN} is not built")
def test_matching_lambdas_run(tmp_path):
    result = _run(tmp_path, "0.5 0.5")
    assert result.returncode == 0
    assert "terminated normally" in result.stdout
