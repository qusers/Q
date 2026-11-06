"""Regression test for issue #131: heating of restrained atoms outside the sphere.

Excluded atoms are integrated under their positional restraints but do not
contribute to the free-atom temperature. Before the fix they were rescaled with
the factor computed from the free atoms only, which pumps kinetic energy into
the excluded region. ``excluded_thermostat separate`` (the default) couples them
to their own bath at the target temperature.
"""

import os
import re
import shutil
import subprocess
from pathlib import Path
from statistics import mean

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
QDYN = Path(os.environ.get("QDYN", PROJECT_ROOT / "src" / "q6" / "qdyn"))
TOPOLOGY = PROJECT_ROOT / "test" / "data" / "topology" / "ala_wat.top"
TARGET_K = 298.0

INPUT_TEMPLATE = """\
[MD]
steps                     {steps}
stepsize                  1.0
temperature               {target}
bath_coupling             10
random_seed               971
initial_temperature       {target}
separate_scaling          on
{thermostat_line}
shake_solvent             on
shake_hydrogens           off

[cut_offs]
solute_solute             10
solute_solvent            10
solvent_solvent           10
q_atom                    99
lrf                       99

[intervals]
non_bond                  25
output                    1000
temperature               10
energy                    0
trajectory                0

[files]
topology                  ala_wat.top
final                     final.re
"""


def _run(tmp_path: Path, mode: str | None, steps: int = 8000) -> str:
    workdir = tmp_path / (mode or "default")
    workdir.mkdir()
    shutil.copy(TOPOLOGY, workdir / TOPOLOGY.name)
    line = f"excluded_thermostat       {mode}" if mode else ""
    (workdir / "md.inp").write_text(
        INPUT_TEMPLATE.format(steps=steps, target=TARGET_K, thermostat_line=line)
    )
    result = subprocess.run(
        [str(QDYN), "md.inp"],
        cwd=workdir,
        text=True,
        capture_output=True,
        timeout=300,
    )
    log = result.stdout + result.stderr
    # Older engines may exit 0 on failure; require the normal-termination marker.
    assert "terminated normally" in log, log[-2000:]
    return log


def _late_excluded_temperature(log: str) -> float:
    values = [float(v) for v in re.findall(r"T_excl=\s*([-\d.]+)", log)]
    assert len(values) > 20, "excluded-atom temperature was not reported"
    return mean(values[len(values) // 2 :])


pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not QDYN.is_file(), reason="qdyn must be compiled in src/q6 (or set QDYN)"),
]


def test_default_couples_excluded_atoms_to_their_own_bath(tmp_path):
    log = _run(tmp_path, None)
    assert "Excluded atoms coupled to their own heat bath." in log
    assert abs(_late_excluded_temperature(log) - TARGET_K) < 25.0


@pytest.mark.slow
def test_legacy_scaling_still_overheats_excluded_atoms(tmp_path):
    """Keeps the reproducer: the old coupling drives excluded atoms far above target."""
    legacy = _late_excluded_temperature(_run(tmp_path, "legacy", steps=20000))
    separate = _late_excluded_temperature(_run(tmp_path, "separate", steps=20000))
    assert legacy > 1.5 * TARGET_K
    assert legacy > separate + 150.0


def test_frozen_excluded_atoms_have_no_kinetic_energy(tmp_path):
    log = _run(tmp_path, "frozen", steps=2000)
    assert "Excluded atoms are frozen at their current positions." in log
    values = [float(v) for v in re.findall(r"T_excl=\s*([-\d.]+)", log)]
    assert values and max(values) == 0.0


def test_invalid_excluded_thermostat_is_rejected(tmp_path):
    workdir = tmp_path / "bogus"
    workdir.mkdir()
    shutil.copy(TOPOLOGY, workdir / TOPOLOGY.name)
    (workdir / "md.inp").write_text(
        INPUT_TEMPLATE.format(
            steps=10, target=TARGET_K, thermostat_line="excluded_thermostat       bogus"
        )
    )
    result = subprocess.run(
        [str(QDYN), "md.inp"], cwd=workdir, text=True, capture_output=True, timeout=60
    )
    assert "terminated normally" not in result.stdout
    assert "excluded_thermostat must be separate, none, frozen or legacy" in result.stdout
