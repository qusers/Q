"""Native regression: Born constants follow topology units without changing dynamics."""
from pathlib import Path
import struct
import subprocess

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
QDYN = PROJECT_ROOT / "src/q6/qdyn"
DATA = PROJECT_ROOT / "test/data"


def _real_to_dummy_fep() -> str:
    return """[FEP]
states 2
[atoms]
1 13

[change_charges]
1 1.0 0.0

[atom_types]
SOD 1036.9613 6.8099 0.0 0.0 733.2424 4.8153 22.9898
DUM 0.0 0.0 0 0 0.0 0.0 1.0080

[softcore]
1 0 0

[change_atoms]
1 SOD DUM
"""


def _md_input(topology: Path, fep: Path, final: Path, lambda1: float) -> str:
    return f"""[MD]
steps 1
stepsize 0.001
temperature 1
bath_coupling 1
random_seed 112
initial_temperature 1
shake_solvent off
shake_hydrogens off
shake_solute off
lrf off

[cut-offs]
solute_solvent 99.0
solute_solute 99.0
solvent_solvent 99.0
q_atom 99.0
lrf 99.0

[sphere]
shell_force 10.0
shell_radius 20

[solvent]
radial_force 60.0
polarization on
polarization_force 20.0
charge_correction off

[intervals]
output 1
non_bond 1
energy 1

[files]
topology {topology}
fep {fep}
final {final}
energy audit.en

[lambdas]
{lambda1:.6f} {1.0 - lambda1:.6f}
"""


def _restart_record(handle) -> bytes:
    size_bytes = handle.read(4)
    if not size_bytes:
        raise EOFError
    size = struct.unpack("=i", size_bytes)[0]
    payload = handle.read(size)
    assert struct.unpack("=i", handle.read(4))[0] == size
    return payload


def _read_fortran_records(path: Path) -> list[bytes]:
    records = []
    with path.open("rb") as handle:
        while True:
            try:
                records.append(_restart_record(handle))
            except EOFError:
                return records


def _energy_state_records(path: Path) -> list[dict[int, tuple[float, ...]]]:
    records = _read_fortran_records(path)
    steps = []
    for offset in range(0, len(records), 3):
        state_records = records[offset : offset + 2]
        assert len(state_records) == 2
        states = {}
        for payload in state_records:
            assert len(payload) == 4 + 15 * 8
            state = struct.unpack("=i", payload[:4])[0]
            states[state] = struct.unpack("=15d", payload[4:])
        assert states.keys() == {1, 2}
        steps.append(states)
    return steps


def _last_born_block(stdout: str) -> dict[int, float]:
    states: dict[int, float] = {}
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) == 8 and fields[0] == "State" and fields[5] == "E_Born":
            states[int(fields[1])] = float(fields[7])
    assert states.keys() == {1, 2}
    return states


@pytest.mark.parametrize("ke,override,lambda1,q_charge", [
    (332.0716, None, .5, 1),
    (250.0, None, .0001, -1),
    (332.0716, 3.5, .9999, 1),
])
def test_born_changes_only_serialized_state_constants(tmp_path, ke, override, lambda1, q_charge):
    """Compare Born on/off on the same coordinates, including a coefficient override."""
    if not QDYN.is_file():
        pytest.skip("Build src/q6/qdyn to run the native Born regression")
    newest_source = max(path.stat().st_mtime for path in (PROJECT_ROOT / "src/q6").glob("*.f90"))
    assert QDYN.stat().st_mtime >= newest_source, "Rebuild stale src/q6/qdyn"

    # All other topology fields stay identical between these two runs.
    template = (DATA / "topology/Na-benzene-water.top").read_text()
    assert "332.0716 = Electrostatic" in template
    topology = tmp_path / "system.top"
    topology.write_text(template.replace("332.0716 = Electrostatic", f"{ke:.4f} = Electrostatic"))
    records, restarts = {}, {}
    born_states = None
    for enabled in (False, True):
        label = "born" if enabled else "baseline"
        run = tmp_path / label
        run.mkdir()
        # Q atom 1 leaves included non-Q solute charge +1.115 e.
        fep = _real_to_dummy_fep().replace("1 13", "1 1")
        fep = fep.replace("1 1.0 0.0", f"1 {q_charge:.1f} 0.0")
        (run / "audit.fep").write_text(fep)
        # Q's legacy filename fields are only 80 characters wide.
        inp = _md_input(Path("../system.top"), Path("audit.fep"), Path("audit.re"), lambda1)
        inp = inp.replace("steps 1\n", "steps 10\n")
        if enabled:
            extra = "perstate_born_correction on\nborn_dielectric 80.0"
            if override is not None:
                extra += f"\nborn_coefficient {override}"
            inp = inp.replace("charge_correction off", "charge_correction off\n"+extra)
        (run / "audit.inp").write_text(inp)
        result = subprocess.run([str(QDYN), "audit.inp"], cwd=run,
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        records[label] = _energy_state_records(run / "audit.en")
        restarts[label] = (run / "audit.re").read_bytes()
        if enabled:
            born_states = _last_born_block(result.stdout)

    assert restarts["baseline"] == restarts["born"]
    assert len(records["baseline"]) == len(records["born"]) == 9
    # Independent fixture values; no use of qdyn's accumulator or printed coefficient.
    coefficient = override if override is not None else ke*(1-1/80)/(2*20.1)
    expected = {1: -coefficient*(1.115+q_charge)**2, 2: -coefficient*1.115**2}
    shifts = {1: [], 2: []}
    for old, new in zip(records["baseline"], records["born"]):
        for state in (1, 2):
            shifts[state].append(new[state][1]-old[state][1])
            assert new[state][0] == old[state][0]
            assert new[state][2:] == pytest.approx(old[state][2:], rel=0, abs=1e-12)
    for state in (1, 2):
        assert max(shifts[state])-min(shifts[state]) < 1e-12
        # Atomic charges are serialized in single precision in this topology.
        assert shifts[state][0] == pytest.approx(expected[state], rel=0, abs=1e-6)
        assert shifts[state][0] == pytest.approx(born_states[state], rel=0, abs=.000051)
    assert shifts[2][0]-shifts[1][0] == pytest.approx(expected[2]-expected[1], rel=0, abs=2e-6)
