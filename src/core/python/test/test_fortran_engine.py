"""qgpu.FortranEngine against Fortran Q run by hand on .inp files.

Needs the Fortran binaries (cd src/q6 && make all COMP=gcc && make mpi COMP=gcc), no GPU:

    PYTHONPATH=src/core/python/src pytest src/core/python/test/test_fortran_engine.py

QDYN / QDYNP select other binaries, QGPU_TEST_MPI the launcher (default "mpirun -np 2").
"""

import os
from pathlib import Path
import shlex
import shutil
import subprocess

import numpy as np
import pytest

from qgpu import Config, DistRestraint, FortranEngine, InputError, SeqRestraint
from qgpu.fortran import read_restart, write_restart

ROOT = Path(__file__).resolve().parents[4]
BIN = ROOT / "src" / "q6" / "bin" / "q6"
QDYN = Path(os.environ.get("QDYN", BIN / "qdyn"))
QDYNP = Path(os.environ.get("QDYNP", BIN / "qdynp"))
MPI = shlex.split(os.environ.get("QGPU_TEST_MPI", "mpirun -np 2"))
TOPOLOGY = ROOT / "test" / "data" / "topology" / "benzene-water.top"
FEP = ROOT / "test" / "data" / "inputs" / "FEP_benzene.fep"

serial = pytest.param((QDYN, ()), id="serial", marks=pytest.mark.skipif(
    not QDYN.is_file(), reason=f"{QDYN} is not built"))
mpi = pytest.param((QDYNP, MPI), id="mpi", marks=pytest.mark.skipif(
    not QDYNP.is_file() or shutil.which(MPI[0]) is None, reason=f"{QDYNP} or {MPI[0]} is missing"))


@pytest.fixture(params=[serial, mpi])
def launch(request):
    """(qdyn executable, launcher prefix)"""
    return request.param


def base_config(**changes):
    values = dict(
        topology=TOPOLOGY, fep=FEP, seed=7, initial_temperature=1, lrf=True,
        steps=100, stepsize=0.5, temperature=1, bath_coupling=1,
        constrain_solvent=True, constrain_hydrogens=False, lambdas=(1.0, 0.0), output=50,
        seq_restraints=(SeqRestraint(1, 12, 5.0),), dist_restraints=(DistRestraint(1, 4, 2.0, 3.0, 1.5),),
    )
    return Config(**{**values, **changes})


# The same three stages as the API calls in test_matches_inp_files, written by hand.
STAGES = {
    "eq1": """[MD]
steps 100
stepsize 0.5
temperature 1
bath_coupling 1
random_seed 7
initial_temperature 1
constrain_solvent on
constrain_hydrogens off
lrf on
[intervals]
output 50
[files]
topology top.top
fep fep.fep
final eq1.re
[lambdas]
1.0 0.0
[sequence_restraints]
1 12 5.0 0 0
[distance_restraints]
1 4 2.0 3.0 1.5 0
""",
    "eq2": """[MD]
steps 150
stepsize 1.0
temperature 300
bath_coupling 10
constrain_solvent on
constrain_hydrogens on
lrf on
[intervals]
output 50
[files]
topology top.top
fep fep.fep
restart eq1.re
final eq2.re
[lambdas]
1.0 0.0
[distance_restraints]
1 4 2.0 3.0 0.5 0
""",
    "md": """[MD]
steps 100
stepsize 1.0
temperature 300
bath_coupling 10
constrain_solvent on
constrain_hydrogens on
lrf on
[intervals]
output 50
energy 10
trajectory 50
[files]
topology top.top
fep fep.fep
restart eq2.re
final md.re
energy md.en
trajectory md.dcd
[trajectory_atoms]
not excluded
[lambdas]
0.5 0.5
[distance_restraints]
1 4 2.0 3.0 0.5 0
""",
}


def run_by_hand(folder, launch):
    qdyn, prefix = launch
    folder.mkdir()
    shutil.copy(TOPOLOGY, folder / "top.top")
    shutil.copy(FEP, folder / "fep.fep")
    for name, text in STAGES.items():
        (folder / f"{name}.inp").write_text(text)
        with open(folder / f"{name}.log", "w") as log:
            subprocess.run([*prefix, str(qdyn), f"{name}.inp"], cwd=folder, stdout=log, check=True)
    return folder


def test_matches_inp_files(tmp_path, launch):
    """Three stages through the API give the same restart and energy files, byte for byte."""
    qdyn, prefix = launch
    reference = run_by_hand(tmp_path / "by_hand", launch)

    engine = FortranEngine(base_config(topology=reference / "top.top"), qdyn=qdyn, mpi=prefix,
                           workdir=tmp_path / "api")
    engine.run()
    engine.update_config(steps=150, stepsize=1.0, temperature=300, bath_coupling=10, constrain_hydrogens=True,
                         seq_restraints=(), dist_restraints=(DistRestraint(1, 4, 2.0, 3.0, 0.5),))
    engine.run()
    engine.update_config(steps=100, lambdas=(0.5, 0.5), energy=10, energy_file=tmp_path / "out" / "md.en",
                         final_file=tmp_path / "out" / "md.re",
                         trajectory=50, trajectory_file=tmp_path / "out" / "md.dcd", trajectory_atoms="not excluded")
    engine.run()

    assert engine.n_atoms == 3420
    assert (tmp_path / "out" / "md.re").read_bytes() == (reference / "md.re").read_bytes()
    assert (tmp_path / "out" / "md.en").read_bytes() == (reference / "md.en").read_bytes()
    assert (tmp_path / "out" / "md.dcd").read_bytes() == (reference / "md.dcd").read_bytes()


# What the C++ parser fills in for the settings base_config() leaves out (inp_parser.cpp).
CPP_DEFAULTS = dict(
    shell_radius=0.85, exclude_bonded=False, shell_force=10.0, solute_solute=10.0, solvent_solvent=10.0,
    solute_solvent=10.0, q_atom=99.0, lrf_cutoff=99.0, radial_force=60.0, polarisation=True,
    polarisation_force=20.0, charge_correction=True, non_bond=10, hmr=False, hmr_target_mass=3.024,
    solute_constraint_algorithm="lincs", solvent_constraint_algorithm="settle", thermostat="berendsen",
    constrain_solute=False, separate_scaling=False, energy=0, trajectory=0,
)


def test_explicit_defaults_match_fortran_defaults(tmp_path):
    """With the CUDA modules installed, Config holds the C++ parser's defaults and FortranEngine
    writes them all out. That must be the same simulation as leaving them to Fortran Q."""
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    finals = []
    for name, extra in (("unset", {}), ("explicit", CPP_DEFAULTS)):
        final = tmp_path / f"{name}.re"
        engine = FortranEngine(base_config(final_file=final, **extra), qdyn=QDYN, workdir=tmp_path / name,
                                   keep_files=True)
        engine.run()
        engine.update_config(stepsize=1.0, temperature=300, bath_coupling=10, constrain_hydrogens=True)
        engine.run()
        finals.append(final.read_bytes())
    assert "polarization_force" in (tmp_path / "explicit" / "run_0001.inp").read_text()
    assert finals[0] == finals[1]


def test_state_branches(tmp_path, launch):
    """set_state() returns to a saved state: the same window run twice from it is identical."""
    qdyn, prefix = launch
    engine = FortranEngine(base_config(), qdyn=qdyn, mpi=prefix, workdir=tmp_path / "work")
    engine.run()
    midpoint = engine.get_state()
    assert midpoint["coords"].shape == midpoint["velocities"].shape == (engine.n_atoms, 3)

    def window(name, lambdas):
        engine.update_config(lambdas=lambdas, energy=10, energy_file=tmp_path / name)
        engine.run()
        return (tmp_path / name).read_bytes()

    first = window("a.en", (0.5, 0.5))
    moved = engine.get_state()
    assert not np.array_equal(moved["coords"], midpoint["coords"])
    assert window("b.en", (0.5, 0.5)) != first       # continues from the moved state
    engine.set_state(midpoint)
    assert window("c.en", (0.5, 0.5)) == first       # back at the midpoint


def test_output_goes_to_stdout(tmp_path, capsys, launch):
    """qdyn's output is printed, as the CUDA engine's is, and is also in the log of the run."""
    qdyn, prefix = launch
    engine = FortranEngine(base_config(steps=10), qdyn=qdyn, mpi=prefix, workdir=tmp_path, keep_files=True)
    engine.run()
    first = capsys.readouterr().out
    assert "terminated normally" in first and first == engine.log_file.read_text()
    engine.update_config(temperature=0)
    with pytest.raises(InputError):
        engine.run()
    assert "No dynamics at zero temperature" in capsys.readouterr().out


def test_leaves_no_files_behind(tmp_path):
    """Only what the Config names stays on disk; keep_files keeps the engine's own files too."""
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    work = tmp_path / "work"
    with FortranEngine(base_config(steps=10), qdyn=QDYN, workdir=work) as engine:
        engine.run()
        state = engine.get_state()
        engine.update_config(energy=5, energy_file=tmp_path / "out" / "a.en", final_file=tmp_path / "out" / "a.re")
        engine.run()
        engine.set_state(state)
        engine.run()
        assert sorted(p.name for p in work.iterdir()) == ["benzene-water.top", "fep.fep", "run_0004.re"]
    assert not work.exists()
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["a.en", "a.re"]

    kept = tmp_path / "kept"
    with FortranEngine(base_config(steps=10), qdyn=QDYN, workdir=kept, keep_files=True) as engine:
        engine.run()
    assert {"run_0001.inp", "run_0001.log", "run_0001.re"} <= {p.name for p in kept.iterdir()}


def test_restart_round_trip(tmp_path):
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    engine = FortranEngine(base_config(steps=10, final_file=tmp_path / "final.re"), qdyn=QDYN)
    engine.run()
    write_restart(tmp_path / "copy.re", read_restart(tmp_path / "final.re"))
    assert (tmp_path / "copy.re").read_bytes() == (tmp_path / "final.re").read_bytes()


def test_from_inp(tmp_path):
    """Config.from_inp() reads a hand-written file (old shake_* keys included) and runs it."""
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    # the file names its outputs next to itself: keep the test from writing into the repository
    config = Config.from_inp(ROOT / "test" / "q6" / "test1" / "eq1.inp", trajectory=0, trajectory_file="",
                             final_file=tmp_path / "eq1.re")
    assert (config.steps, config.stepsize, config.seed, config.lrf) == (100, 0.1, 1, True)
    assert config.constrain_solvent is True and config.shell_radius == 16
    assert config.seq_restraints == (SeqRestraint(1, 46, 0.5, False, 1),)
    assert Path(config.topology).name == "lig_w.top"
    assert Path(Config.from_inp(ROOT / "test" / "q6" / "test1" / "eq1.inp").final_file).name == "eq1.re"
    engine = FortranEngine(config, qdyn=QDYN, workdir=tmp_path / "work")
    engine.run()
    assert (tmp_path / "eq1.re").is_file() and engine.log_file is None


def test_errors(tmp_path):
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    with pytest.raises(InputError, match="thermostat"):
        FortranEngine(base_config(thermostat="langevin"), qdyn=QDYN)
    with pytest.raises(InputError, match="seed"):
        FortranEngine(base_config(initial_temperature=None), qdyn=QDYN)
    with pytest.raises(InputError, match="not found"):
        FortranEngine(base_config(), qdyn=tmp_path / "no_such_qdyn")

    engine = FortranEngine(base_config(steps=10), qdyn=QDYN, workdir=tmp_path)
    engine.run()
    with pytest.raises(InputError, match="Expected 2 lambdas"):
        engine.update_config(lambdas=(1.0,))
    before = engine.get_state()
    engine.update_config(temperature=0)              # qdyn refuses this
    with pytest.raises(InputError, match="zero temperature"):
        engine.run()
    assert np.array_equal(engine.get_state()["coords"], before["coords"])
    engine.update_config(temperature=1)              # still usable, at the state before the failure
    engine.run()


def test_extra_keywords(tmp_path, capsys):
    """Keywords without a Config field reach qdyn through Config.extra and can change between runs."""
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    with pytest.raises(InputError, match="own Config field"):
        base_config(extra={"MD": {"steps": 5}})
    with pytest.raises(InputError, match="own Config field"):
        base_config(extra={"solvent": {"polarization": True}})

    config = base_config(steps=10, extra={"MD": {"excluded_thermostat": "bogus"}})
    assert config.extra == {"md": {"excluded_thermostat": "bogus"}}
    engine = FortranEngine(config, qdyn=QDYN, workdir=tmp_path, keep_files=True)
    with pytest.raises(InputError):
        engine.run()                                  # qdyn read the keyword and refused its value
    assert "excluded_thermostat must be" in capsys.readouterr().out
    engine.update_config(extra={"md": {"excluded_thermostat": "none"}, "PBC": {}})
    engine.run()
    assert "excluded_thermostat" in (tmp_path / "run_0001.inp").read_text()
    engine.update_config(extra={})
    engine.run()
    assert "excluded_thermostat" not in (tmp_path / "run_0002.inp").read_text()


def test_from_sections_and_lossless_inp(tmp_path):
    """Config.from_sections() equals Config.from_inp() of the same content; unknown keys are kept."""
    from qgpu.fortran import read_inp

    (tmp_path / "a.inp").write_text(
        "[MD]\nsteps 10\nstepsize 1.0\ntemperature 300\nminimize on\nshake_solvent on\n"
        "constraint_algorithm lincs settle\n[PBC]\npressure 1.0\n[files]\ntopology t.top\nfinal a.re\n"
        "[lambdas]\n0.5 0.5\n[trajectory_atoms]\nnot excluded\n[sequence_restraints]\n1 12 5.0 0 1\n")
    values = read_inp(tmp_path / "a.inp")
    assert values["extra"] == {"md": {"minimize": "on"}, "pbc": {"pressure": "1.0"}}
    sections = {
        "MD": {"steps": 10, "stepsize": 1.0, "temperature": 300, "minimize": True, "constrain_solvent": True,
               "constraint_algorithm": "lincs settle"},
        "PBC": {"pressure": 1.0}, "files": {"topology": "t.top", "final": "a.re"},
        "lambdas": [0.5, 0.5], "trajectory_atoms": "not excluded", "sequence_restraints": [[1, 12, 5.0, 0, 1]],
    }
    from qgpu import _native
    if not _native.is_any_available():
        assert Config.from_sections(sections, folder=tmp_path) == Config.from_inp(tmp_path / "a.inp")
    config = Config.from_sections(sections, folder=tmp_path, steps=20)
    assert config.steps == 20 and config.extra["md"] == {"minimize": "on"}
    assert config.seq_restraints == (SeqRestraint(1, 12, 5.0, False, 1),)
    assert config.final_file == str(tmp_path.resolve() / "a.re")


def test_energies_in_memory(tmp_path, launch):
    """get_energies() holds what energy_file would hold; without a file none is left anywhere."""
    from qgpu.fortran import read_energy_file

    qdyn, prefix = launch
    results = []
    for name, energy_file in (("file", tmp_path / "md.en"), ("memory", "")):
        work = tmp_path / name
        engine = FortranEngine(base_config(steps=40), qdyn=qdyn, mpi=prefix, workdir=work)
        with pytest.raises(InputError, match="No energies"):
            engine.get_energies()
        engine.run()
        engine.update_config(lambdas=(0.5, 0.5), energy=10, energy_file=energy_file)
        engine.run()
        results.append(engine.get_energies())
        assert sorted(p.name for p in work.iterdir()) == ["benzene-water.top", "fep.fep", "run_0002.re"]
        engine.update_config(energy=0, energy_file="")
        engine.run()
        with pytest.raises(InputError, match="No energies"):
            engine.get_energies()
        engine.close()
    from_file, in_memory = results
    assert from_file == in_memory == read_energy_file(tmp_path / "md.en")
    assert in_memory.to_bytes() == (tmp_path / "md.en").read_bytes()
    assert (in_memory.n_frames, in_memory.n_states, in_memory.n_offdiag) == (3, 2, 0)  # steps 10, 20, 30: qdyn counts steps from 0
    rows = in_memory.q_energies
    assert list(rows["frame"]) == [1, 1, 2, 2, 3, 3] and list(rows["state"]) == [1, 2] * 3
    assert set(rows["lambda"]) == {0.5}


def test_energy_file_codec(tmp_path):
    """Frames with off-diagonal elements, an empty file, and what is refused."""
    import struct

    from qgpu.fortran import ENERGY_TERMS, parse_energies

    def record(payload):
        return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))

    data = b""
    for frame in range(3):
        for state in (1, 2, 3):
            data += record(struct.pack("<i", state) + struct.pack(f"<{len(ENERGY_TERMS)}d",
                                                                    *(frame + state + k / 7 for k in range(15))))
        data += record(b"".join(struct.pack("<iidd", i, i + 1, frame * 0.1, -i / 3) for i in (1, 2)))
    energies = parse_energies(data)
    assert (energies.n_frames, energies.n_states, energies.n_offdiag) == (3, 3, 2)
    assert energies.to_bytes() == data
    assert list(energies.offdiag["j"]) == [2, 3] * 3 and energies.offdiag["Hij"][-1] == 0.2
    assert energies.q_energies["q_bond"][4] == 1 + 2 + 2 / 7
    assert parse_energies(b"").to_bytes() == b"" and parse_energies(b"").n_frames == 0
    with pytest.raises(InputError, match="Truncated"):
        parse_energies(data[:-5])
    with pytest.raises(InputError, match="Not a Q energy file"):
        parse_energies(record(b"x" * 10))
    with pytest.raises(InputError, match="losslessly"):
        parse_energies(data[:-4] + struct.pack("<i", 47))


def test_cuda_engine_refuses_energies_in_memory():
    """Refused before the CUDA modules are loaded, so this needs no GPU."""
    from qgpu import Engine

    with pytest.raises(InputError, match="energy > 0 needs energy_file"):
        Engine(base_config(energy=10))
    with pytest.raises(InputError, match="cannot return energies"):
        Engine.__new__(Engine).get_energies()
    assert base_config(energy=10).energy_file == ""  # Config itself accepts it
