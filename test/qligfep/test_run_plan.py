"""The run plan: one file per perturbation leg instead of one .inp file per stage.

The plan must describe exactly the simulation the .inp files describe. That is tested
at three levels: the plan against the rendered inputs, the engine settings of every
stage against the ones read from the .inp files, and the output of a short run of both
routes with Fortran Q, byte for byte.

The runs need the Fortran binaries (cd src/q6 && make all COMP=gcc && make mpi COMP=gcc).
QDYN / QDYNP select other binaries, QGPU_TEST_MPI the launcher (default "mpirun -np 2").
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import textwrap
from dataclasses import replace
from pathlib import Path

import pytest

from QligFEP.CLI.config_args import resolve_run_settings
from QligFEP.config import ConfigError, StageOverride, load_config
from QligFEP.energies import export, read_energies
from QligFEP.energies import main as energies_main
from QligFEP.qligfep import QligFEP
from QligFEP.run_plan import (
    _STAGE_HEADER,
    PlanError,
    _import_qgpu,
    build_plan,
    expand,
    load_plan,
    main,
    parse_inp,
    window_name,
    window_order,
    write_plan,
)

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "src" / "q6" / "bin" / "q6"
QDYN = Path(os.environ.get("QDYN", BIN / "qdyn"))
QDYNP = Path(os.environ.get("QDYNP", BIN / "qdynp"))
QFEP = Path(os.environ.get("QFEP", BIN / "qfep"))
MPI = shlex.split(os.environ.get("QGPU_TEST_MPI", "mpirun -np 2"))
TOPOLOGY = ROOT / "test" / "data" / "topology" / "benzene-water.top"
FEP = ROOT / "test" / "data" / "inputs" / "FEP_benzene.fep"
SEEDS = [11, 12]

serial = pytest.param((QDYN, []), id="serial", marks=pytest.mark.skipif(
    not QDYN.is_file(), reason=f"{QDYN} is not built"))  # fmt: skip
mpi = pytest.param((QDYNP, MPI), id="mpi", marks=pytest.mark.skipif(
    not QDYNP.is_file() or shutil.which(MPI[0]) is None, reason=f"{QDYNP} or {MPI[0]} is missing"))  # fmt: skip


@pytest.fixture(autouse=True)
def isolated_config(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("QLIGFEP_CONFIG", raising=False)
    monkeypatch.delenv("QLIGFEP_Q_BIN", raising=False)
    monkeypatch.chdir(tmp_path)


def short(steps: int, **inp) -> StageOverride:
    return StageOverride(steps=steps, inp=inp)


#: Stages short enough to run in seconds on the benzene test system.
SHORT_STAGES = {
    **{f"eq{i}": short(20) for i in range(1, 6)},
    "production": short(20, intervals={"energy": 5, "trajectory": 10}),
}


def setup(run_plan: bool, lambdas=None, start="0.5", system="water", stages=None, production=False) -> QligFEP:
    """A QligFEP object that can write the MD inputs of the benzene test system.

    Only the attributes the input writers use are set: building the topology is not
    what is tested here.
    """
    run = QligFEP.__new__(QligFEP)
    run.timestep, run.sphereradius, run.minimize = "2fs", "20", False
    run.system, run.atomoffset, run.dr_force, run.start = system, 0, 0.5, start
    run.production, run.counter_water_atoms = production, []
    run.stage_overrides = SHORT_STAGES if stages is None else stages
    run.temperature, run.seeds = "298", SEEDS
    run.run_plan, run._plan_stages = run_plan, []
    run.custom_lambdas = lambdas
    return run


def write_inputs(run: QligFEP, folder: Path, windows=4, sampling="sigmoidal") -> list[str]:
    """inputfiles/ of the benzene system as QligFEP writes it; returns the lambda list."""
    folder.mkdir(parents=True)
    shutil.copy(TOPOLOGY, folder / "dualtop.top")
    shutil.copy(FEP, folder / "FEP1.fep")
    lambdas = run.get_lambdas(windows, sampling)
    run.write_md_files(lambdas, str(folder), 6, 6, [[1, 4], [2, 5]], wall_restraints="7 12 15.0 1.0 0 0 0")
    return lambdas


def stage_inputs(folder: Path) -> dict[str, str]:
    return {path.stem: path.read_text() for path in folder.glob("*.inp")}


def fill(text: str, temperature=298, seed=SEEDS[0]) -> str:
    """What the run script's sed lines do to an .inp file."""
    return text.replace("T_VAR", str(temperature)).replace("SEED_VAR", str(seed)).replace("FEP_VAR", "FEP1.fep")


# ---------------------------------------------------------------------------------
# The plan against the .inp files
# ---------------------------------------------------------------------------------

CASES = {
    "default": {},
    "start-at-endpoint": {"start": "1"},
    "protein-leg": {"system": "protein"},
    "no-trajectories": {"production": True},
    "explicit-lambdas": {"lambdas": [0.0, 0.25, 0.5, 0.9, 0.97, 1.0]},
    "overrides": {
        "stages": {
            **SHORT_STAGES,
            "eq1": short(20, MD={"minimize": True, "max_minimize_steps": 5}),
            "eq5": short(20, MD={"excluded_thermostat": "none"}, sphere={"shell_force": 5.0}),
            "production": short(20, **{"cut-offs": {"solute_solvent": 9}, "intervals": {"energy": 5}}),
        }
    },
}


@pytest.mark.parametrize("case", CASES)
def test_plan_describes_the_same_stages_as_the_inp_files(tmp_path, case):
    """Every stage of the plan has the engine settings and the restart of its .inp file."""
    qgpu = _import_qgpu()
    options = CASES[case]
    write_inputs(setup(False, **options), tmp_path / "inp")
    write_inputs(setup(True, **options), tmp_path / "plan")

    assert not list((tmp_path / "plan").glob("*.inp"))
    plan = load_plan(tmp_path / "plan" / "run_plan.yaml")
    assert plan["replicates"] == {"temperatures": [298], "seeds": SEEDS}
    inputs = stage_inputs(tmp_path / "inp")
    stages = expand(plan, 310, SEEDS[1], inputfiles=tmp_path / "inp", rundir=tmp_path / "inp")
    assert sorted(stage.name for stage in stages) == sorted(inputs)

    done = set()
    for stage in stages:
        (tmp_path / "inp" / "filled.inp").write_text(fill(inputs[stage.name], 310, SEEDS[1]))
        expected = qgpu.Config.from_inp(tmp_path / "inp" / "filled.inp")
        restart = Path(expected.restart).stem if expected.restart else None
        assert stage.parent == restart and (restart is None or restart in done), stage.name
        assert qgpu.Config.from_sections(stage.sections, folder=tmp_path) == replace(expected, restart=""), stage.name
        done.add(stage.name)


def test_plan_is_compact_and_states_per_replicate_values(tmp_path):
    write_inputs(setup(True), tmp_path / "plan")
    plan = load_plan(tmp_path / "plan" / "run_plan.yaml")
    assert plan["system"] == {"topology": "dualtop.top", "fep": "FEP1.fep"}
    assert plan["settings"]["cut-offs"]["q_atom"] == 99 and "temperature" not in plan["settings"]["MD"]
    first, last = plan["equilibration"]["stages"][0], plan["equilibration"]["stages"][-1]
    assert first["MD"]["random_seed"] == "replicate" and first["MD"]["temperature"] == 1
    assert last["MD"]["temperature"] == "replicate" and "random_seed" not in last["MD"]
    assert plan["restraints"]["distance"] == {"range": [0.0, 0.1], "state": 0, "pairs": [[1, 4], [2, 5]]}
    assert (first["distance_force"], last["distance_force"]) == (1.5, 0.5)
    assert plan["windows"]["start"] == 0.5
    assert plan["windows"]["lambdas"] == [[1.0, 0.0], [0.958, 0.042], [0.5, 0.5], [0.042, 0.958], [0.0, 1.0]]
    assert plan["windows"]["settings"]["intervals"]["energy"] == 5


def test_window_order():
    ladder = [[1.0, 0.0], [0.9, 0.1], [0.5, 0.5], [0.25, 0.75], [0.0, 1.0]]
    # the start window, then away from it on one side, then on the other
    assert window_order(ladder, 0.5) == [(2, None), (1, 2), (0, 1), (3, 2), (4, 3)]
    assert window_order(ladder, 1.0) == [(0, None), (1, 0), (2, 1), (3, 2), (4, 3)]
    with pytest.raises(PlanError, match="not one of the lambda windows"):
        window_order(ladder, 0.3)
    assert window_name([0.5, 0.5]) == "md_0500_0500" and window_name([1, 0.067]) == "md_1000_0067"


def test_explicit_lambdas_need_not_be_symmetric(tmp_path):
    """State 2 runs at 1 - lambda, in the inputs and in qfep.inp alike."""
    run = setup(False, lambdas=[0.0, 0.5, 0.9, 1.0])
    lambdas = write_inputs(run, tmp_path / "inp")
    assert lambdas == ["1.000", "0.900", "0.500", "0.000"]
    assert run.lambda_pairs(lambdas)[1] == ("0.900", "0.100")
    inputs = stage_inputs(tmp_path / "inp")
    assert {name for name in inputs if name.startswith("md_")} == {
        "md_1000_0000", "md_0900_0100", "md_0500_0500", "md_0000_1000",
    }  # fmt: skip
    restarts = {name: parse_inp(text)["files"]["restart"] for name, text in inputs.items() if name.startswith("md_")}
    assert restarts == {
        "md_0500_0500": "eq5.re",
        "md_0900_0100": "md_0500_0500.re",
        "md_1000_0000": "md_0900_0100.re",
        "md_0000_1000": "md_0500_0500.re",
    }


def test_start_in_the_middle_needs_a_middle_window(tmp_path):
    with pytest.raises(ValueError, match="needs a window at lambda 0.500"):
        write_inputs(setup(False), tmp_path / "inp", windows=5, sampling="linear")


def test_what_a_plan_cannot_express_is_refused(tmp_path):
    write_inputs(setup(False), tmp_path / "inp")
    inputs = stage_inputs(tmp_path / "inp")
    stages = sorted(inputs.items())  # eq1..eq5, then the windows
    with pytest.raises(PlanError, match="differ only in their lambdas"):
        build_plan([(n, t.replace("steps                     20", "steps 30") if n == "md_0000_1000" else t)
                    for n, t in stages], [298], SEEDS, "FEP1.fep")  # fmt: skip
    with pytest.raises(PlanError, match="final file must be eq1.re"):
        build_plan([(n, t.replace("eq1.re", "other.re") if n == "eq1" else t) for n, t in stages], [298], SEEDS, "f")
    with pytest.raises(PlanError, match="needs equilibration stages and lambda windows"):
        build_plan([(n, t) for n, t in stages if n.startswith("eq")], [298], SEEDS, "f")


# ---------------------------------------------------------------------------------
# Lambda windows in the protocol
# ---------------------------------------------------------------------------------


def protocol_file(tmp_path: Path, text: str) -> None:
    (tmp_path / "qligfep.yaml").write_text(textwrap.dedent(text))


def cli_args(**given) -> argparse.Namespace:
    names = [
        "cluster", "protocol", "config", "q_bin", "sphereradius", "start", "temperature", "replicates",
        "sampling", "windows", "timestep", "neq", "neq_reps", "neq_steps", "neq_eq_steps",
        "neq_relax_steps", "neq_L", "neq_schedule",
    ]  # fmt: skip
    return argparse.Namespace(**{**dict.fromkeys(names), **given})


class TestProtocolLambdas:
    def test_explicit_list(self, tmp_path):
        protocol_file(tmp_path, "protocols: {p: {equilibrium: {lambdas: [0, 0.5, 0.9004, 1]}}}")
        protocol = load_config().get_protocol("p")
        assert protocol.lambdas == (1.0, 0.9, 0.5, 0.0)  # from 1 down, 3 decimals
        assert protocol.as_dict()["equilibrium"]["lambdas"] == [1.0, 0.9, 0.5, 0.0]
        assert load_config().get_protocol().lambdas is None

    @pytest.mark.parametrize(
        "equilibrium, message",
        [
            ("{lambdas: [0, 0.5, 1], windows: 10}", "either the list or windows/sampling"),
            ("{lambdas: [0, 0.5, 1], sampling: linear}", "either the list or windows/sampling"),
            ("{lambdas: [0, 0.5, 1.2]}", "not a number between 0 and 1"),
            ("{lambdas: [0, 0.5, 0.5001, 1]}", "more than once"),
            ("{lambdas: [0.1, 0.5, 1]}", "both end points"),
            ("{lambdas: [0, 0.4, 1]}", "start window"),
            ("{lambdas: 5}", "expected a list"),
        ],
    )
    def test_invalid_lists(self, tmp_path, equilibrium, message):
        protocol_file(tmp_path, f"protocols: {{p: {{equilibrium: {equilibrium}}}}}")
        with pytest.raises(ConfigError, match=message):
            load_config().get_protocol("p")

    def test_start_at_the_endpoint_needs_no_middle_window(self, tmp_path):
        protocol_file(tmp_path, "protocols: {p: {equilibrium: {lambdas: [0, 0.4, 1], start: 1}}}")
        assert load_config().get_protocol("p").lambdas == (1.0, 0.4, 0.0)

    def test_a_protocol_that_sets_windows_drops_an_inherited_list(self, tmp_path):
        protocol_file(
            tmp_path,
            """
            protocols:
              listed: {equilibrium: {lambdas: [0, 0.5, 1]}}
              sampled: {extends: listed, equilibrium: {windows: 20}}
            """,
        )
        protocol = load_config().get_protocol("sampled")
        assert protocol.lambdas is None and protocol.windows == 20

    def test_command_line(self, tmp_path):
        protocol_file(tmp_path, "protocols: {p: {equilibrium: {lambdas: [0, 0.5, 0.9, 1]}}}")
        args = cli_args(cluster="HABROK", protocol="p")
        resolve_run_settings(args)
        assert args.lambdas == [1.0, 0.9, 0.5, 0.0] and args.windows == "3"
        args = cli_args(cluster="HABROK", protocol="p", windows="10")  # the flag wins
        resolve_run_settings(args)
        assert args.lambdas is None and args.windows == "10"


# ---------------------------------------------------------------------------------
# Run script
# ---------------------------------------------------------------------------------


def test_run_script_calls_the_runner_once(tmp_path):
    run = setup(True)
    run.profile = load_config().get_profile("HABROK")
    run.cluster, run.replicates, run.lig1, run.lig2, run.to_clean = "HABROK", "2", "a", "b", ["dcd"]
    write_inputs(run, tmp_path / "plan")
    run.write_runfile(str(tmp_path / "plan"), [])
    script = (tmp_path / "plan" / "runHABROK.sh").read_text()
    assert script.count("-m QligFEP.run_plan $inputfiles/run_plan.yaml") == 1
    assert '--replicate $run_num --temperature $temperature --qdyn $qdyn --mpi "mpirun -n $ncores' in script
    # the job may not inherit the setup environment: the script names the checkout itself
    assert f'export PYTHONPATH="{ROOT / "src"}${{PYTHONPATH:+:$PYTHONPATH}}"' in script
    assert "eq1.inp" not in script and "sedi" not in script and "T_VAR" not in script
    # what the script keeps: seeds, footer, binary check, qfep, the completion check, cleanup
    for kept in ("seeds=(11 12)", "trap write_footer EXIT", "ERROR: Q binary missing", "*dcd"):
        assert kept in script, kept
    # qfep and the completion check moved into the runner; qfep.inp is not copied
    assert "--qfep $qfep --qfep-input $inputfiles/qfep.inp" in script
    for gone in ("qfep.inp .", "< qfep.inp", "grep -q", "md*.log", "eq1.minimization"):
        assert gone not in script, gone
    assert subprocess.run(["bash", "-n", tmp_path / "plan" / "runHABROK.sh"]).returncode == 0


def test_plan_is_the_default_for_windowed_fep_only():
    """No flag gives the run plan; a non-equilibrium setup keeps .inp files without one."""
    profile = load_config().get_profile("HABROK")
    for neq, expected in ((False, True), (True, False)):
        run = QligFEP.__new__(QligFEP)
        QligFEP.__init__(run, "a", "b", "AMBER14sb", "water", profile=profile, neq=neq)
        assert run.run_plan is expected
    run = QligFEP.__new__(QligFEP)
    QligFEP.__init__(run, "a", "b", "AMBER14sb", "water", profile=profile, run_plan=False)
    assert run.run_plan is False


def test_command_line_route_flags():
    from QligFEP.CLI.parser_base import parse_arguments
    from QligFEP.CLI.setupFEP import create_call

    base = ["qligfep", "-l1", "a", "-l2", "b", "-FF", "AMBER14sb", "-s", "water", "-c", "LOCAL"]
    for extra, inp_files in (([], False), (["--inp-files"], True), (["--run-plan"], False)):
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr("sys.argv", [*base, *extra])
            assert parse_arguments("QligFEP").inp_files is inp_files
    kwargs = dict(
        lig1="a", lig2="b", FF="AMBER14sb", system="water", cluster="LOCAL", replicates=1, sampling="sigmoidal",
        sphereradius=25, start="0.5", windows=10, temperature=298, timestep="2fs", rest="overlap", dr_force=0.5,
        log="info",
    )  # fmt: skip
    assert "--inp-files" not in create_call(**kwargs) and "--run-plan" not in create_call(**kwargs)
    assert "--inp-files" in create_call(**kwargs, inp_files=True)


# ---------------------------------------------------------------------------------
# Both routes, run
# ---------------------------------------------------------------------------------


def run_inp_route(inputfiles: Path, rundir: Path, order: list[str], launch, temperature, seed) -> None:
    """What run.sh does for one replicate: copy, fill in the placeholders, run stage by stage."""
    qdyn, prefix = launch
    rundir.mkdir(parents=True)
    for name in ("dualtop.top", "FEP1.fep"):
        shutil.copy(inputfiles / name, rundir / name)
    for name in order:
        (rundir / f"{name}.inp").write_text(fill((inputfiles / f"{name}.inp").read_text(), temperature, seed))
        with open(rundir / f"{name}.log", "w") as log:
            subprocess.run([*prefix, str(qdyn), f"{name}.inp"], cwd=rundir, stdout=log, check=True)
        assert "terminated normally" in (rundir / f"{name}.log").read_text(), name


def run_plan_route(plan_file: Path, rundir: Path, launch, *options) -> int:
    return main([str(plan_file), "--replicate", "2", "--temperature", "298", "--rundir", str(rundir),
                 "--qdyn", str(launch[0]), "--mpi", shlex.join(launch[1]), *options])  # fmt: skip


def stage_headers(log: Path) -> list[str]:
    return [line.split()[3] for line in log.read_text().splitlines() if line.startswith(_STAGE_HEADER)]


@pytest.mark.parametrize("launch", [serial, mpi])
@pytest.mark.parametrize("case", ["default", "start-at-endpoint", "explicit-lambdas", "overrides", "no-trajectories"])
def test_plan_run_equals_inp_run(tmp_path, launch, case):
    """Energies and trajectories of a replicate are the same bytes on both routes; with
    --keep-stage-files also the restart file and energy file of every stage."""
    options = CASES[case]
    write_inputs(setup(False, **options), tmp_path / "inp")
    write_inputs(setup(True, **options), tmp_path / "plan")
    plan_file = tmp_path / "plan" / "run_plan.yaml"
    order = [stage.name for stage in expand(load_plan(plan_file), 298, SEEDS[1])]
    windows = [name for name in order if name.startswith("md_")]

    reference = tmp_path / "reference"
    run_inp_route(tmp_path / "inp", reference, order, launch, 298, SEEDS[1])
    energies = sorted(p.name for p in reference.glob("*.en"))
    trajectories = sorted(p.name for p in reference.glob("*.dcd"))
    assert energies == sorted(f"{name}.en" for name in windows)
    assert bool(trajectories) == (case != "no-trajectories")

    # lean: one log, one energy file, the trajectories; nothing else
    rundir = tmp_path / "FEP1" / "298" / "2"
    assert run_plan_route(plan_file, rundir, launch) == 0
    assert sorted(p.name for p in rundir.iterdir()) == sorted([*trajectories, "energies.parquet", "replicate.log"])
    assert stage_headers(rundir / "replicate.log") == order
    assert (rundir / "replicate.log").read_text().count("terminated normally") == len(order)
    assert [w.name for w in read_energies(rundir / "energies.parquet")] == windows
    export(rundir / "energies.parquet", tmp_path / "exported")
    assert sorted(p.name for p in (tmp_path / "exported").iterdir()) == energies
    for name in energies:
        assert (tmp_path / "exported" / name).read_bytes() == (reference / name).read_bytes(), name
    for name in trajectories:
        assert (rundir / name).read_bytes() == (reference / name).read_bytes(), name

    # everything, as the .inp route writes it
    kept = tmp_path / "kept"
    assert run_plan_route(plan_file, kept, launch, "--keep-stage-files") == 0
    restarts = [f"{name}.re" for name in order]
    logs = [f"{name}.log" for name in order]
    assert sorted(p.name for p in kept.iterdir()) == sorted(
        [*energies, *trajectories, *restarts, *logs, "energies.parquet", "replicate.log"]
    )
    for name in [*energies, *trajectories, *restarts]:
        assert (kept / name).read_bytes() == (reference / name).read_bytes(), name
    assert (kept / "energies.parquet").read_bytes() == (rundir / "energies.parquet").read_bytes()


def write_qfep_input(run: QligFEP, folder: Path, lambdas: list[str]) -> Path:
    run.write_dir = str(folder.parent)
    run.write_qfep(len(lambdas) - 1, lambdas)
    return folder / "qfep.inp"


def test_qfep_runs_on_the_energy_file(tmp_path):
    """--qfep gives the qfep.out and energies.csv of qfep run on the .en files, and leaves no .en."""
    if not (QDYN.is_file() and QFEP.is_file()):
        pytest.skip(f"{QDYN} or {QFEP} is not built")
    # qfep skips as many points per window as there are lambda steps (4): record 7
    run = setup(True, stages={**SHORT_STAGES, "production": short(40, intervals={"energy": 5, "trajectory": 10})})
    inputfiles = tmp_path / "setup" / "inputfiles"
    lambdas = write_inputs(run, inputfiles)
    qfep_input = write_qfep_input(run, inputfiles, lambdas)
    rundir = tmp_path / "run"
    assert run_plan_route(inputfiles / "run_plan.yaml", rundir, (QDYN, []), "--qfep", str(QFEP)) == 0
    assert sorted(p.name for p in rundir.iterdir()) == [
        "energies.csv", "energies.parquet", "eq1.dcd", "eq2.dcd", "eq3.dcd", "eq4.dcd", "eq5.dcd",
        *sorted(f"{window_name(pair)}.dcd" for pair in load_plan(inputfiles / "run_plan.yaml")["windows"]["lambdas"]),
        "qfep.out", "replicate.log",
    ]  # fmt: skip

    by_hand = tmp_path / "by_hand"
    export(rundir / "energies.parquet", by_hand)
    with open(qfep_input) as stdin, open(by_hand / "qfep.out", "w") as stdout:
        subprocess.run([str(QFEP)], stdin=stdin, stdout=stdout, cwd=by_hand, check=True)
    for name in ("qfep.out", "energies.csv"):
        assert (rundir / name).read_bytes() == (by_hand / name).read_bytes(), name
    assert "Part 3" in (rundir / "qfep.out").read_text()

    assert energies_main(["qfep", str(rundir / "energies.parquet"), "--qfep", str(QFEP), "--input",
                          str(qfep_input), "--out", str(tmp_path / "again")]) == 0  # fmt: skip
    assert (tmp_path / "again" / "qfep.out").read_bytes() == (rundir / "qfep.out").read_bytes()


def test_failed_stage_keeps_what_helps(tmp_path, capsys):
    """The log up to the failure, the state the failed stage started from, the finished windows."""
    if not QDYN.is_file():
        pytest.skip(f"{QDYN} is not built")
    stages = {**SHORT_STAGES, "eq3": short(20, MD={"excluded_thermostat": "bogus"})}
    write_inputs(setup(True, stages=stages), tmp_path / "plan")
    arguments = [str(tmp_path / "plan" / "run_plan.yaml"), "--temperature", "298", "--rundir", str(tmp_path / "run"),
                 "--qdyn", str(QDYN)]  # fmt: skip
    assert main([*arguments, "--replicate", "1"]) == 1
    error = capsys.readouterr().err
    assert "Stage eq3 failed" in error and "excluded_thermostat" in error
    assert f"state at the end of eq2: {tmp_path / 'run' / 'eq2.re'}" in error
    assert sorted(p.name for p in (tmp_path / "run").iterdir()) == ["eq1.dcd", "eq2.dcd", "eq2.re", "replicate.log"]
    assert stage_headers(tmp_path / "run" / "replicate.log") == ["eq1", "eq2", "eq3"]
    assert main([*arguments, "--replicate", "3"]) == 2
    assert "Replicate 3 is not in the plan" in capsys.readouterr().err

    # a failure among the windows: the energy file holds the windows finished before it
    stages = {**SHORT_STAGES, "production": short(20, intervals={"energy": 5, "trajectory": 10})}
    write_inputs(setup(True, stages=stages), tmp_path / "plan2")
    plan = load_plan(tmp_path / "plan2" / "run_plan.yaml")
    plan["windows"]["lambdas"][0] = [1.0, 0.0, 0.5]  # three lambdas for two states: qdyn refuses it
    write_plan(tmp_path / "plan2" / "run_plan.yaml", plan)
    rundir = tmp_path / "run2"
    assert main([str(tmp_path / "plan2" / "run_plan.yaml"), "--replicate", "1", "--temperature", "298",
                 "--rundir", str(rundir), "--qdyn", str(QDYN)]) == 1  # fmt: skip
    error = capsys.readouterr().err
    done = [w.name for w in read_energies(rundir / "energies.parquet")]
    assert done == ["md_0500_0500", "md_0958_0042"] and "energies of 2 completed window(s)" in error
    assert (rundir / "md_0958_0042.re").is_file()


def test_analysis_counts_lambdas_from_the_plan(tmp_path):
    """qligfep_analyze finds the number of lambda steps without md*.inp files."""
    from QligFEP.analyze_FEP import FepReader

    write_inputs(setup(True), tmp_path / "plan" / "inputfiles")
    write_inputs(setup(False), tmp_path / "inp" / "inputfiles")
    assert FepReader._input_n_lambdas(tmp_path / "plan", 1) == FepReader._input_n_lambdas(tmp_path / "inp", 1) == 4
    assert FepReader._input_n_lambdas(tmp_path, 1) is None


def test_energy_file_round_trip(tmp_path):
    """Windows with different frame counts, off-diagonal elements and an empty window come back
    unchanged; the file is replaced in one step."""
    import struct

    qgpu = _import_qgpu()
    from qgpu.fortran import ENERGY_TERMS, parse_energies

    from QligFEP.energies import EnergyStoreError, Window, write_energies

    def en_file(frames, offdiag, seed):
        def record(payload):
            return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))

        data = b""
        for frame in range(frames):
            for state in (1, 2):
                values = [seed + frame / 3 + state + k / 11 for k in range(len(ENERGY_TERMS))]
                data += record(struct.pack("<i", state) + struct.pack(f"<{len(values)}d", *values))
            data += record(b"".join(struct.pack("<iidd", 1, 2, seed * frame, 0.25) for _ in range(offdiag)))
        return data

    raw = {"md_0500_0500": en_file(5, 0, 1.0), "md_1000_0000": en_file(3, 1, 2.0), "md_0000_1000": b""}
    windows = [Window(name, (0.5, 0.5), parse_energies(data)) for name, data in raw.items()]
    path = write_energies(tmp_path / "energies.parquet", windows)
    assert [p.name for p in tmp_path.iterdir()] == ["energies.parquet"]
    assert [w.name for w in read_energies(path)] == list(raw)
    export(path, tmp_path / "out")
    for name, data in raw.items():
        assert (tmp_path / "out" / f"{name}.en").read_bytes() == data, name
    assert read_energies(path, ["md_1000_0000"])[0].energies == windows[1].energies

    import pandas as pd

    table = pd.read_parquet(path)
    assert list(table.columns[:6]) == ["window", "lambda1", "lambda2", "frame", "state", "lambda"]
    assert len(table) == 2 * (5 + 3) and table["offdiag"].map(lambda cells: cells is not None and len(cells) > 0).sum() == 3
    with pytest.raises(EnergyStoreError, match="no window md_0100_0900"):
        read_energies(path, ["md_0100_0900"])
    assert qgpu  # the engine package is what reads and writes the records
