"""Run a windowed FEP replicate from one run plan instead of one .inp file per stage.

At setup, ``qligfep --run-plan`` writes ``inputfiles/run_plan.yaml``: the settings every
stage shares, what each equilibration stage changes, the lambda windows, and what is
filled in per replicate (temperature, random seed). The run script then calls

    qligfep_run inputfiles/run_plan.yaml --replicate 1 --temperature 298 \\
        --qdyn /path/qdynp --mpi "mpirun -n 8"

which runs all stages through the qgpu engine API with one engine, and so one persistent
``qdyn --serve`` process, for the replicate. It
leaves one log (``replicate.log``), one energy file (``energies.parquet``, see
:mod:`QligFEP.energies`) and the trajectories, and with ``--qfep`` runs qfep on the
energies in a temporary folder.

The plan is built from the text of the .inp files QligFEP would otherwise write
(:func:`build_plan`), so both routes describe the same simulation, and
:func:`expand` turns it back into the content of each stage.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import re
import shlex
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, fields, replace
from pathlib import Path

import yaml

PLAN_VERSION = 1
#: Value of a setting that the runner fills in for each replicate.
REPLICATE = "replicate"
_PLACEHOLDERS = {"T_VAR": REPLICATE, "SEED_VAR": REPLICATE}
#: Per-replicate settings: [MD] key -> what the runner puts there.
_REPLICATE_KEYS = {"temperature": "temperature", "random_seed": "seed"}
_ROW_SECTIONS = ("sequence_restraints", "distance_restraints", "wall_restraints", "angle_restraints", "atom_restraints")
_WINDOW_PREFIX = "md_"

_HEADER = """\
# QligFEP run plan: every stage of one perturbation leg, in place of one .inp file per stage.
# Run a replicate with:  qligfep_run run_plan.yaml --replicate N --temperature T --qdyn ... --mpi "..."
#
#   settings       qdyn input keywords shared by every stage, as section -> key
#   equilibration  stages run in order at fixed lambdas; each lists what it does not share
#   windows        lambda windows, each with windows.settings. The window at `start` runs first,
#                  from the end of equilibration; every other window continues from its
#                  neighbour toward `start`.
#   replicates     `replicate` as a value below means: filled in per replicate
#                  (temperature, and the random seed of the first stage)
#   restraints     distance restraint pairs; each stage gives their force constant (distance_force)
# A replicate folder receives replicate.log (every stage's output), energies.parquet (the
# energies of every window) and <name>.dcd for each stage whose trajectory interval is > 0.
# Windows are named md_<lambda1>_<lambda2>.
"""


class PlanError(ValueError):
    """A run plan cannot be built, read or run as written."""


# ----------------------------------------------------------------------------------
# Building the plan at setup
# ----------------------------------------------------------------------------------


def _scalar(text: str):
    """A value as written in an .inp file -> YAML scalar."""
    if text in _PLACEHOLDERS:
        return _PLACEHOLDERS[text]
    if text.lower() in ("on", "off"):
        return text.lower() == "on"
    if re.fullmatch(r"[-+]?\d+", text):
        return int(text)
    try:
        return float(text)
    except ValueError:
        return text


def parse_inp(text: str) -> dict:
    """The content of a qdyn .inp file as data.

    Key-value sections become ``{key: value}``, restraint sections lists of rows,
    ``lambdas`` a list and ``trajectory_atoms`` a string. Empty sections are left out.
    """
    sections: dict = {}
    rows = None
    for line in text.splitlines():
        line = re.split(r"[!#*]", line)[0].strip()
        if line.startswith("["):
            rows = sections.setdefault(line.strip("[]").strip(), [])
        elif line and rows is not None:
            rows.append(line.split())
    parsed: dict = {}
    for name, rows in sections.items():
        if not rows:
            continue
        lower = name.lower()
        if lower == "lambdas":
            parsed["lambdas"] = [float(x) for row in rows for x in row]
        elif lower == "trajectory_atoms":
            parsed["trajectory_atoms"] = " ".join(rows[0])
        elif lower in _ROW_SECTIONS:
            parsed[lower] = [[_scalar(x) for x in row] for row in rows]
        else:
            parsed[name] = {row[0]: _scalar(" ".join(row[1:])) for row in rows}
    return parsed


def _interval(stage: dict, key: str) -> int:
    return int((stage.get("intervals") or {}).get(key) or 0)


def _take_files(name: str, stage: dict) -> tuple[str, str | None]:
    """Remove ``files`` from a parsed stage after checking it follows the naming rules.

    Returns (topology, name of the stage it restarts from).
    """
    files = stage.pop("files", {})
    if files.get("final") != f"{name}.re":
        raise PlanError(f"Stage {name}: the final file must be {name}.re, not {files.get('final')}.")
    for key, suffix in (("trajectory", "dcd"), ("energy", "en")):
        expected = f"{name}.{suffix}" if _interval(stage, key) > 0 else None
        if files.get(key) != expected and not (expected is None and key == "trajectory"):
            raise PlanError(
                f"Stage {name}: {key} interval {_interval(stage, key)} does not go with "
                f"{key} file {files.get(key)!r} (expected {expected!r})."
            )
    restart = files.get("restart")
    if restart is not None and not str(restart).endswith(".re"):
        raise PlanError(f"Stage {name}: unexpected restart file {restart}.")
    return files.get("topology"), (str(restart)[:-3] if restart else None)


def _take_distance(name: str, stage: dict) -> tuple[list, float | None]:
    """Replace the distance restraint rows of a stage by their force constant."""
    rows = stage.pop("distance_restraints", [])
    if not rows:
        return [], None
    forces = {row[4] for row in rows}
    if any(len(row) != 6 for row in rows) or len(forces) != 1:
        raise PlanError(f"Stage {name}: distance restraints must share one force constant.")
    stage["distance_force"] = forces.pop()
    return [[row[0], row[1], row[2], row[3], row[5]] for row in rows], stage["distance_force"]


def _items(stage: dict) -> dict:
    """Flat view of a stage: (section, key) -> value; whole-value entries have key None."""
    flat = {}
    for section, content in stage.items():
        if isinstance(content, dict):
            flat.update({(section, key): value for key, value in content.items()})
        else:
            flat[(section, None)] = content
    return flat


def _nested(flat: dict) -> dict:
    stage: dict = {}
    for (section, key), value in flat.items():
        if key is None:
            stage[section] = value
        else:
            stage.setdefault(section, {})[key] = value
    return stage


def window_name(lambdas) -> str:
    """md_1000_0000 for (1.0, 0.0): the name QligFEP, qfep.inp and the analysis use."""
    return _WINDOW_PREFIX + "_".join(f"{value:.3f}".replace(".", "") for value in lambdas)


def window_order(lambdas: list, start: float) -> list[tuple[int, int | None]]:
    """Run order of the windows as (index, index of the window it continues from).

    The start window comes first (parent None: the end of equilibration), then the
    windows above it in the list walking away from it, then the ones below it.
    """
    firsts = [pair[0] for pair in lambdas]
    if start not in firsts:
        raise PlanError(f"windows.start {start} is not one of the lambda windows.")
    middle = firsts.index(start)
    order: list[tuple[int, int | None]] = [(middle, None)]
    order += [(i, i + 1) for i in range(middle - 1, -1, -1)]
    order += [(i, i - 1) for i in range(middle + 1, len(lambdas))]
    return order


def build_plan(stages: list[tuple[str, str]], temperatures: list, seeds: list, fep: str) -> dict:
    """The run plan for stages given as (name, .inp text), as QligFEP renders them.

    Raises:
        PlanError: If the stages do not form what a plan can express: equilibration
            stages in one chain, then windows that differ only in their lambdas.
    """
    parsed, parents, topologies, pairs = {}, {}, set(), None
    for name, text in stages:
        stage = parse_inp(text)
        topology, parents[name] = _take_files(name, stage)
        topologies.add(topology)
        stage_pairs, _ = _take_distance(name, stage)
        if pairs is None or not pairs:
            pairs = stage_pairs
        if stage_pairs and stage_pairs != pairs:
            raise PlanError(f"Stage {name}: its distance restraints differ from the other stages'.")
        parsed[name] = stage
    if len(topologies) != 1 or None in topologies:
        raise PlanError("All stages must use the same topology.")

    equilibration = [name for name in parsed if not name.startswith(_WINDOW_PREFIX)]
    windows = [name for name in parsed if name.startswith(_WINDOW_PREFIX)]
    if not equilibration or not windows:
        raise PlanError("A run plan needs equilibration stages and lambda windows.")
    for previous, name in zip([None, *equilibration], equilibration):
        if parents[name] != previous:
            raise PlanError(f"Stage {name} must continue from {previous}, not {parents[name]}.")

    flat = {name: _items(stage) for name, stage in parsed.items()}
    lambdas = {name: items.pop(("lambdas", None)) for name, items in flat.items()}
    first = next(iter(flat.values()))
    shared = {key: value for key, value in first.items() if all(key in f and f[key] == value for f in flat.values())}
    own = {name: {k: v for k, v in items.items() if k not in shared} for name, items in flat.items()}

    if len({tuple(lambdas[name]) for name in equilibration}) != 1:
        raise PlanError("All equilibration stages must run at the same lambdas.")
    if any(own[name] != own[windows[0]] for name in windows):
        raise PlanError("Lambda windows may differ only in their lambdas.")
    for name in windows:
        if name != window_name(lambdas[name]):
            raise PlanError(f"Window {name} is not named after its lambdas {lambdas[name]}.")

    # ladder order; the start window is the one that continues from equilibration
    ladder = sorted(windows, key=lambda name: -lambdas[name][0])
    starts = [name for name in windows if parents[name] == equilibration[-1]]
    if len(starts) != 1:
        raise PlanError("Exactly one lambda window must continue from the last equilibration stage.")
    start = lambdas[starts[0]][0]
    ladder_lambdas = [lambdas[name] for name in ladder]
    for index, parent in window_order(ladder_lambdas, start)[1:]:
        if parents[ladder[index]] != ladder[parent]:
            raise PlanError(
                f"Window {ladder[index]} continues from {parents[ladder[index]]}, "
                f"not from its neighbour {ladder[parent]}."
            )

    plan: dict = {
        "version": PLAN_VERSION,
        "system": {"topology": topologies.pop(), "fep": fep},
        "replicates": {"temperatures": [_number(t) for t in temperatures], "seeds": [int(s) for s in seeds]},
        "settings": _nested(shared),
    }
    if pairs:
        geometry = {tuple(row[2:]) for row in pairs}
        if len(geometry) != 1:
            raise PlanError("Distance restraints must share their distance range and state.")
        d1, d2, state = geometry.pop()
        plan["restraints"] = {"distance": {"range": [d1, d2], "state": state, "pairs": _Flow(_Flow(r[:2]) for r in pairs)}}
    plan["equilibration"] = {
        "lambdas": lambdas[equilibration[0]],
        "stages": [{"name": name, **_nested(own[name])} for name in equilibration],
    }
    plan["windows"] = {
        "start": start,
        "settings": _nested(own[windows[0]]),
        "lambdas": _Flow(_Flow(pair) for pair in ladder_lambdas),
    }
    return plan


def _number(value):
    value = float(value)
    return int(value) if value.is_integer() else value


class _Flow(list):
    """A list written on one (wrapped) line."""


class _Dumper(yaml.SafeDumper):
    pass


_Dumper.add_representer(_Flow, lambda dumper, data: dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True))


def write_plan(path, plan: dict) -> None:
    text = yaml.dump(plan, Dumper=_Dumper, sort_keys=False, default_flow_style=None, width=100)
    Path(path).write_text(_HEADER + text)


# ----------------------------------------------------------------------------------
# Reading the plan back: the content of every stage of one replicate
# ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class PlannedStage:
    name: str
    #: The stage's input as ``qgpu.Config.from_sections`` reads it (no restart file:
    #: the runner carries the state from ``parent``).
    sections: dict
    #: Stage this one continues from; None for the first.
    parent: str | None


def load_plan(path) -> dict:
    plan = yaml.safe_load(Path(path).read_text())
    if not isinstance(plan, dict) or plan.get("version") != PLAN_VERSION:
        raise PlanError(f"{path}: not a version {PLAN_VERSION} QligFEP run plan.")
    for key in ("system", "replicates", "settings", "equilibration", "windows"):
        if key not in plan:
            raise PlanError(f"{path}: missing '{key}'.")
    return plan


def _merge(base: dict, override: dict) -> dict:
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = {**merged[key], **value}
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _stage_sections(plan: dict, name: str, own: dict, lambdas, temperature, seed, inputfiles: Path, rundir: Path):
    sections = _merge(plan["settings"], own)
    md = next((content for section, content in sections.items() if section.lower() == "md"), {})
    values = {"temperature": temperature, "seed": seed}
    for section, content in sections.items():
        if not isinstance(content, dict):
            continue
        for key, value in content.items():
            if value == REPLICATE:
                if content is not md or key not in _REPLICATE_KEYS:
                    raise PlanError(f"Stage {name}: [{section}] {key} cannot be set per replicate.")
                content[key] = values[_REPLICATE_KEYS[key]]

    force = sections.pop("distance_force", None)
    distance = (plan.get("restraints") or {}).get("distance")
    if distance and force is not None:
        d1, d2 = distance["range"]
        sections["distance_restraints"] = [[a, b, d1, d2, force, distance.get("state", 0)] for a, b in distance["pairs"]]

    files = {"topology": str(inputfiles / plan["system"]["topology"]), "final": str(rundir / f"{name}.re")}
    if plan["system"].get("fep"):
        files["fep"] = str(inputfiles / plan["system"]["fep"])
    for key, suffix in (("trajectory", "dcd"), ("energy", "en")):
        if _interval(sections, key) > 0:
            files[key] = str(rundir / f"{name}.{suffix}")
    sections["files"] = files
    sections["lambdas"] = [float(x) for x in lambdas]
    return sections


def expand(plan: dict, temperature: float, seed: int, inputfiles=".", rundir=".") -> list[PlannedStage]:
    """Every stage of one replicate, in run order.

    Args:
        plan: A loaded run plan.
        temperature: Target temperature of this replicate (K).
        seed: Random seed of this replicate.
        inputfiles: Folder holding the topology and FEP file the plan names.
        rundir: Folder the outputs are written to.
    """
    inputfiles, rundir = Path(inputfiles).resolve(), Path(rundir).resolve()
    temperature = _number(temperature)
    stages: list[PlannedStage] = []
    equilibration = plan["equilibration"]
    for own in equilibration["stages"]:
        own = dict(own)
        name = own.pop("name")
        sections = _stage_sections(plan, name, own, equilibration["lambdas"], temperature, seed, inputfiles, rundir)
        stages.append(PlannedStage(name, sections, stages[-1].name if stages else None))
    if not stages:
        raise PlanError("The plan has no equilibration stage.")

    windows = plan["windows"]
    last_equilibration = stages[-1].name
    names = [window_name(pair) for pair in windows["lambdas"]]
    if len(set(names)) != len(names):
        raise PlanError("Two lambda windows have the same name (lambdas are told apart to 3 decimals).")
    for index, parent in window_order(windows["lambdas"], windows["start"]):
        sections = _stage_sections(
            plan, names[index], windows["settings"], windows["lambdas"][index], temperature, seed, inputfiles, rundir
        )
        stages.append(PlannedStage(names[index], sections, last_equilibration if parent is None else names[parent]))
    return stages


# ----------------------------------------------------------------------------------
# Running a replicate
# ----------------------------------------------------------------------------------

#: Fixed engine settings that only describe how the first stage starts.
_START_ONLY = {"restart", "seed", "initial_temperature"}


def _import_qgpu():
    """The engine API. It lives in this repository (src/core/python) and need not be installed."""
    try:
        import qgpu
    except ImportError:
        from QligFEP import SRC

        bundled = Path(SRC) / "core" / "python" / "src"
        if not (bundled / "qgpu").is_dir():
            raise PlanError(
                "The qgpu package (the Q engine API) is not installed and was not found in "
                f"{bundled}. Install it from src/core/python: QGPU_PRECISIONS=none pip install ."
            ) from None
        sys.path.insert(0, str(bundled))
        import qgpu
    return qgpu


def _open_engine(qgpu, kind: str, config, qdyn, mpi):
    """Start the engine of the chosen kind. The GPU engine plugs in here."""
    if kind == "fortran":
        return qgpu.FortranEngine(config, qdyn=qdyn, mpi=mpi)
    raise PlanError(f"Engine '{kind}' is not supported by the plan runner yet (only 'fortran').")


def run_replicate(
    plan_path,
    replicate: int,
    temperature: float,
    qdyn="qdyn",
    mpi=(),
    rundir=".",
    engine_kind: str = "fortran",
    keep_stage_files: bool = False,
    qfep=None,
    qfep_input=None,
    report=lambda message: print(message, flush=True),
) -> None:
    """Run every stage of one replicate of a run plan, with one engine.

    What stays in ``rundir``: ``replicate.log`` (qdyn's output of every stage, in run
    order), ``energies.parquet`` (the energies of every window, see :mod:`QligFEP.energies`),
    the trajectories the plan asks for, and ``qfep.out`` and ``energies.csv`` when ``qfep``
    is given. ``keep_stage_files`` also writes what each stage wrote before: ``<stage>.log``,
    ``<stage>.re`` and the windows' ``<stage>.en``.

    The engine is only replaced when a stage changes a setting that is fixed for an
    engine's lifetime (cut-offs, shell radius, ...); the new one starts from the
    previous stage's state, through a restart file in a temporary folder.

    Raises:
        PlanError: The plan or the replicate number is invalid.
        qgpu.QgpuError: A stage failed. ``replicate.log`` holds the output up to the failure,
            ``<parent>.re`` the state the failed stage started from, and ``energies.parquet``
            the windows completed so far.
    """
    plan_path = Path(plan_path).resolve()
    plan = load_plan(plan_path)
    seeds = plan["replicates"]["seeds"]
    if not 1 <= replicate <= len(seeds):
        raise PlanError(f"Replicate {replicate} is not in the plan: it has {len(seeds)} seed(s).")
    rundir = Path(rundir).resolve()
    rundir.mkdir(parents=True, exist_ok=True)
    stages = expand(plan, temperature, seeds[replicate - 1], plan_path.parent, rundir)

    qgpu = _import_qgpu()
    from qgpu.config import FIXED
    from qgpu.fortran import write_restart

    from .energies import FILE_NAME, Window, run_qfep, write_energies

    fixed = [name for name in FIXED if name not in _START_ONLY]
    changeable = [f.name for f in fields(qgpu.Config) if f.name not in FIXED]
    # stages that are continued from later, after another stage has run in between
    branch_points = {s.parent for previous, s in zip(stages, stages[1:]) if s.parent != previous.name}
    energy_file = rundir / FILE_NAME

    engine, states, previous, windows, stage, start = None, {}, None, [], None, None
    report(f"Run plan {plan_path}: replicate {replicate}, T={temperature:g}, seed={seeds[replicate - 1]}, {len(stages)} stages")
    with tempfile.TemporaryDirectory(prefix="qligfep_run_") as scratch, open(rundir / LOG_NAME, "w") as log:
        try:
            for number, stage in enumerate(stages, 1):
                config = qgpu.Config.from_sections(stage.sections, folder=rundir)
                if not keep_stage_files:
                    config = replace(config, final_file="", energy_file="")
                if engine is not None and any(getattr(config, n) != getattr(engine.config, n) for n in fixed):
                    # a new engine, from the state the old one holds (or the branch point's)
                    start = states[stage.parent] if stage.parent != previous else engine.get_state()
                    engine.close()
                    engine = None
                    restart = Path(scratch) / f"{stage.parent}.re"
                    write_restart(restart, start)
                    config = replace(config, restart=str(restart))
                    previous = stage.parent
                if engine is None:
                    engine = _open_engine(qgpu, engine_kind, config, qdyn, mpi)
                    start = None
                else:
                    if stage.parent != previous:
                        engine.set_state(states[stage.parent])
                    engine.update_config(**{name: getattr(config, name) for name in changeable})
                started = time.time()
                log.write(f"{_STAGE_HEADER} {stage.name} ({number}/{len(stages)})\n")
                log.flush()
                output = _StageOutput(log, rundir / f"{stage.name}.log" if keep_stage_files else None)
                with output, contextlib.redirect_stdout(output):
                    engine.run()
                if not output.terminated_normally:
                    raise qgpu.QgpuError(f"qdyn did not report a normal termination for stage {stage.name}.")
                if (config.energy or 0) > 0:
                    windows.append(Window(stage.name, tuple(config.lambdas or ()), engine.get_energies()))
                report(f"{stage.name}: {time.time() - started:.1f} s")
                if stage.name in branch_points:
                    states[stage.name] = engine.get_state()
                previous = stage.name
        except BaseException as error:  # a failed stage, or the job being stopped
            kept = [f"log: {rundir / LOG_NAME}"]
            if stage is not None and stage.parent is not None:
                try:  # the engine is at the state the failed stage started from
                    state = start if engine is None else engine.get_state()
                    if state is not None:
                        write_restart(rundir / f"{stage.parent}.re", state)
                        kept.append(f"state at the end of {stage.parent}: {rundir / stage.parent}.re")
                except Exception:  # noqa: BLE001 - keep reporting the original failure
                    pass
            if windows:
                write_energies(energy_file, windows)
                kept.append(f"energies of {len(windows)} completed window(s): {energy_file}")
            where = f"Stage {stage.name} failed" if stage is not None else "The replicate failed"
            if isinstance(error, qgpu.QgpuError):
                raise type(error)(f"{where}: {error} Kept: {'; '.join(kept)}.") from error
            report(f"{where} ({type(error).__name__}). Kept: {'; '.join(kept)}.")
            raise
        finally:
            if engine is not None:
                engine.close()
    write_energies(energy_file, windows)
    report(f"Energies of {len(windows)} windows: {energy_file}")
    if qfep is not None:
        run_qfep(energy_file, qfep, qfep_input or plan_path.parent / "qfep.inp", rundir)
        report(f"qfep: {rundir / 'qfep.out'}")


#: The log of a replicate: every stage's qdyn output, each under a header line.
LOG_NAME = "replicate.log"
_STAGE_HEADER = "=== qligfep_run stage"


class _StageOutput:
    """Where qdyn's output of one stage goes: the replicate log, and the stage's own log if kept."""

    def __init__(self, log, stage_log=None):
        self._log, self._path, self._stage = log, stage_log, None
        self.terminated_normally = False

    def __enter__(self):
        if self._path is not None:
            self._stage = open(self._path, "w")
        return self

    def __exit__(self, *error):
        self._log.flush()
        if self._stage is not None:
            self._stage.close()

    def write(self, text):
        if "terminated normally" in text:
            self.terminated_normally = True
        self._log.write(text)
        if self._stage is not None:
            self._stage.write(text)
        return len(text)

    def flush(self):
        self._log.flush()


class _Stopped(SystemExit):
    """The job is being stopped (SIGTERM, as Slurm sends at the time limit)."""


def _stop(signum, frame):
    raise _Stopped(128 + signum)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="qligfep_run",
        description="Run one replicate of a QligFEP run plan (inputfiles/run_plan.yaml): all "
        "equilibration stages and lambda windows, with one engine.",
    )
    parser.add_argument("plan", help="Path to run_plan.yaml. The topology and FEP file are read from its folder.")
    parser.add_argument("--replicate", type=int, required=True, help="Replicate number, from 1; selects the seed.")
    parser.add_argument("--temperature", type=float, required=True, help="Target temperature (K).")
    parser.add_argument("--qdyn", default="qdyn", help="The qdyn or qdynp executable. Defaults to qdyn on PATH.")
    parser.add_argument(
        "--mpi", default="", help='MPI launcher prefix for qdynp, e.g. "mpirun -n 8". Empty runs qdyn directly.'
    )
    parser.add_argument("--rundir", default=".", help="Folder for the outputs. Defaults to the current folder.")
    parser.add_argument(
        "--engine", default="fortran", choices=["fortran"], help="Engine running the plan. Only Fortran Q for now."
    )
    parser.add_argument(
        "--qfep", help="Run this qfep on the energies at the end; writes qfep.out and energies.csv to --rundir."
    )
    parser.add_argument(
        "--qfep-input", help="The qfep input. Defaults to qfep.inp next to the plan."
    )
    parser.add_argument(
        "--keep-stage-files", action="store_true",
        help="Also write each stage's .log and .re and each window's .en into --rundir, as the .inp route does "
        "(for debugging). By default only replicate.log, energies.parquet and trajectories are written.",
    )
    args = parser.parse_args(argv)
    signal.signal(signal.SIGTERM, _stop)  # stopped jobs still save the windows finished so far
    try:
        run_replicate(
            args.plan, args.replicate, args.temperature, qdyn=args.qdyn, mpi=shlex.split(args.mpi),
            rundir=args.rundir, engine_kind=args.engine, keep_stage_files=args.keep_stage_files,
            qfep=args.qfep, qfep_input=args.qfep_input,
        )  # fmt: skip
    except PlanError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    except (RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as error:
        # qgpu.QgpuError: a stage failed; EnergyStoreError; qfep failed
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    finally:
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
    return 0


def main_exe() -> None:
    sys.exit(main())


if __name__ == "__main__":
    main_exe()
