"""Platform-agnostic production runner for windowed Q FEP calculations.

The runner executes generated Q input definitions without copying every input,
topology, and FEP file into each replicate directory.  It discovers the restart
dependency graph, renders one transient input at a time, and assigns stages to a
bounded set of rolling checkpoint lanes.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

_SECTION_RE = re.compile(r"^\s*\[([^]]+)]\s*$")
_EQ_NUMBER_RE = re.compile(r"^eq(\d+)$")


class RunnerConfigurationError(ValueError):
    """Raised when generated inputs do not describe a runnable stage graph."""


class StageExecutionError(RuntimeError):
    """Raised when qdyn does not complete a stage normally."""


@dataclass(frozen=True)
class Stage:
    """One Q dynamics invocation and its logical file dependencies."""

    name: str
    source: Path
    kind: str
    restart: str | None
    final: str
    energy: str | None
    lambdas: tuple[float, ...]


@dataclass(frozen=True)
class PlannedStage:
    """A stage assigned to a reusable physical checkpoint lane."""

    stage: Stage
    lane: str

    @property
    def checkpoint(self) -> str:
        return f"checkpoint.{self.lane}.re"

    @property
    def next_checkpoint(self) -> str:
        return f"checkpoint.{self.lane}.next.re"


@dataclass(frozen=True)
class RunPlan:
    """Validated execution order and dependency metadata."""

    stages: tuple[PlannedStage, ...]
    children: dict[str, tuple[str, ...]]
    terminal_restarts: tuple[str, ...]

    @property
    def energy_files(self) -> tuple[str, ...]:
        return tuple(item.stage.energy for item in self.stages if item.stage.energy is not None)


def _normalized_section(name: str) -> str:
    return name.strip().lower().replace("_", "-")


def read_input_sections(path: str | Path) -> dict[str, list[tuple[str, str]]]:
    """Read key/value and free-form lines from a generated Q input file."""
    sections: dict[str, list[tuple[str, str]]] = defaultdict(list)
    section = ""
    with Path(path).open(encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.split("!", 1)[0].strip()
            if not line:
                continue
            match = _SECTION_RE.match(line)
            if match:
                section = _normalized_section(match.group(1))
                continue
            fields = line.split(maxsplit=1)
            key = fields[0].lower()
            value = fields[1].strip() if len(fields) == 2 else ""
            sections[section].append((key, value))
    return dict(sections)


def _section_value(sections: dict[str, list[tuple[str, str]]], section: str, key: str) -> str | None:
    for candidate, value in sections.get(section, []):
        if candidate == key:
            return value
    return None


def _read_stage(path: Path, kind: str) -> Stage:
    sections = read_input_sections(path)
    final = _section_value(sections, "files", "final")
    if final is None:
        raise RunnerConfigurationError(f"{path}: [files] final is required")

    lambda_fields = sections.get("lambdas", [])
    lambdas: tuple[float, ...] = ()
    if lambda_fields:
        first, remainder = lambda_fields[0]
        try:
            lambdas = tuple(float(value) for value in (first + " " + remainder).split())
        except ValueError as exc:
            raise RunnerConfigurationError(f"{path}: invalid [lambdas] values") from exc

    return Stage(
        name=path.stem,
        source=path,
        kind=kind,
        restart=_section_value(sections, "files", "restart"),
        final=final,
        energy=_section_value(sections, "files", "energy"),
        lambdas=lambdas,
    )


def _equilibration_sort_key(path: Path) -> int:
    match = _EQ_NUMBER_RE.match(path.stem)
    if match is None:
        raise RunnerConfigurationError(f"unsupported equilibration input name: {path.name}")
    return int(match.group(1))


def _validate_unique_outputs(stages: Iterable[Stage]) -> dict[str, Stage]:
    producers: dict[str, Stage] = {}
    for stage in stages:
        if stage.final in producers:
            raise RunnerConfigurationError(
                f"{stage.source}: final restart {stage.final!r} is also produced by "
                f"{producers[stage.final].source}"
            )
        producers[stage.final] = stage
    return producers


def _order_production_stages(
    production: list[Stage],
    available_restarts: set[str],
) -> tuple[list[Stage], dict[str, list[Stage]]]:
    """Return a deterministic depth-first order for the production restart DAG."""
    by_parent: dict[str, list[Stage]] = defaultdict(list)
    roots: list[Stage] = []
    production_outputs = {stage.final for stage in production}
    for stage in production:
        if stage.restart in production_outputs:
            assert stage.restart is not None
            by_parent[stage.restart].append(stage)
        elif stage.restart is None or stage.restart in available_restarts:
            roots.append(stage)
        else:
            raise RunnerConfigurationError(
                f"{stage.source}: restart {stage.restart!r} is not produced by any stage"
            )

    for children in by_parent.values():
        children.sort(key=lambda item: item.name)
    roots.sort(key=lambda item: item.name)
    if len(roots) != 1:
        raise RunnerConfigurationError(
            f"expected one production root, found {len(roots)}: "
            + ", ".join(stage.name for stage in roots)
        )

    ordered: list[Stage] = []
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(stage: Stage) -> None:
        if stage.final in visiting:
            raise RunnerConfigurationError(f"restart dependency cycle at {stage.source}")
        if stage.final in visited:
            return
        visiting.add(stage.final)
        ordered.append(stage)
        for child in by_parent.get(stage.final, []):
            visit(child)
        visiting.remove(stage.final)
        visited.add(stage.final)

    visit(roots[0])
    if len(visited) != len(production):
        missing = sorted(stage.name for stage in production if stage.final not in visited)
        raise RunnerConfigurationError(f"unreachable production stages: {', '.join(missing)}")
    return ordered, by_parent


def _assign_checkpoint_lanes(
    equilibration: list[Stage],
    production: list[Stage],
    production_children: dict[str, list[Stage]],
) -> dict[str, str]:
    lanes = {stage.final: "equilibration" for stage in equilibration}
    if not production:
        return lanes

    branch_number = 0
    for stage in production:
        if stage.restart not in {item.final for item in production}:
            lane = "production"
        else:
            parent_lane = lanes[stage.restart]
            siblings = production_children.get(stage.restart, [])
            if len(siblings) <= 1:
                lane = parent_lane
            else:
                branch_number += 1
                lane = f"branch-{branch_number}"
        lanes[stage.final] = lane
    return lanes


def build_run_plan(input_dir: str | Path) -> RunPlan:
    """Discover generated inputs and construct a bounded-checkpoint run plan."""
    root = Path(input_dir).resolve()
    equilibration_paths = sorted(root.glob("eq[0-9]*.inp"), key=_equilibration_sort_key)
    production_paths = sorted(root.glob("md_*.inp"))
    if not equilibration_paths:
        raise RunnerConfigurationError(f"{root}: no eq*.inp files found")
    if not production_paths:
        raise RunnerConfigurationError(f"{root}: no md_*.inp files found")

    equilibration = [_read_stage(path, "equilibration") for path in equilibration_paths]
    production = [_read_stage(path, "production") for path in production_paths]
    producers = _validate_unique_outputs([*equilibration, *production])

    previous_final: str | None = None
    for index, stage in enumerate(equilibration):
        if index == 0:
            if stage.restart is not None:
                raise RunnerConfigurationError(f"{stage.source}: first equilibration stage has a restart")
        elif stage.restart != previous_final:
            raise RunnerConfigurationError(
                f"{stage.source}: expected restart {previous_final!r}, found {stage.restart!r}"
            )
        previous_final = stage.final

    production_order, production_children = _order_production_stages(
        production,
        {stage.final for stage in equilibration},
    )
    lanes = _assign_checkpoint_lanes(equilibration, production_order, production_children)
    planned = tuple(
        PlannedStage(stage, lanes[stage.final]) for stage in [*equilibration, *production_order]
    )

    children: dict[str, list[str]] = defaultdict(list)
    for stage in [*equilibration, *production]:
        if stage.restart in producers:
            assert stage.restart is not None
            children[stage.restart].append(stage.final)
    normalized_children = {
        final: tuple(sorted(values)) for final, values in children.items()
    }
    terminal_restarts = tuple(
        stage.final for stage in production if not normalized_children.get(stage.final)
    )
    return RunPlan(planned, normalized_children, terminal_restarts)


def _runtime_path(value: str, input_dir: Path) -> str:
    path = Path(value)
    return str(path if path.is_absolute() else (input_dir / path).resolve())


def render_runtime_input(
    planned: PlannedStage,
    input_dir: str | Path,
    temperature: str | float,
    seed: int,
    fep_file: str,
    restart_file: str | None,
    trajectory: bool = False,
) -> str:
    """Render one transient Q input with runtime paths and checkpoint names."""
    input_root = Path(input_dir).resolve()
    replacements = {
        "T_VAR": str(temperature),
        "SEED_VAR": str(seed),
        "FEP_VAR": fep_file,
    }
    lines = planned.stage.source.read_text(encoding="utf-8").splitlines()
    output: list[str] = []
    section = ""
    saw_restart = False

    for raw_line in lines:
        line = raw_line
        for token, replacement in replacements.items():
            line = line.replace(token, replacement)
        match = _SECTION_RE.match(line.strip())
        if match:
            section = _normalized_section(match.group(1))
            output.append(line)
            continue

        content = line.split("!", 1)[0].strip()
        fields = content.split(maxsplit=1)
        key = fields[0].lower() if fields else ""
        value = fields[1].strip() if len(fields) == 2 else ""

        if section == "intervals" and key == "trajectory" and not trajectory:
            output.append("trajectory                0")
            continue
        if section == "files":
            if key == "topology":
                output.append(f"topology                  {_runtime_path(value, input_root)}")
                continue
            if key == "fep":
                output.append(f"fep                       {_runtime_path(fep_file, input_root)}")
                continue
            if key == "trajectory" and not trajectory:
                continue
            if key == "restart":
                saw_restart = True
                if restart_file is not None:
                    output.append(f"restart                   {restart_file}")
                continue
            if key == "final":
                if restart_file is not None and not saw_restart:
                    output.append(f"restart                   {restart_file}")
                    saw_restart = True
                output.append(f"final                     {planned.next_checkpoint}")
                continue
        output.append(line)

    if any(token in "\n".join(output) for token in replacements):
        raise RunnerConfigurationError(f"{planned.stage.source}: unresolved runtime placeholder")
    return f"! QligFEP stage: {planned.stage.name}\n" + "\n".join(output) + "\n"


class ProductionRunner:
    """Execute a run plan directly in a replicate directory.

    Progress is recorded after every successful stage.  A restart dependency is
    promoted from ``checkpoint.<lane>.next.re`` only after qdyn returns and
    prints its normal-termination marker.
    """

    state_filename = "run-state.json"
    log_filename = "run.log"
    current_input_filename = "current.inp"

    def __init__(
        self,
        input_dir: str | Path,
        work_dir: str | Path,
        temperature: str | float,
        seed: int,
        fep_file: str,
        qdyn: str,
        launcher: str = "",
        trajectories: bool = False,
    ) -> None:
        self.input_dir = Path(input_dir).resolve()
        self.work_dir = Path(work_dir).resolve()
        self.temperature = str(temperature)
        self.seed = int(seed)
        self.fep_file = fep_file
        self.qdyn = qdyn
        self.launcher = launcher
        self.trajectories = trajectories
        self.plan = build_run_plan(self.input_dir)

    @property
    def state_path(self) -> Path:
        return self.work_dir / self.state_filename

    @property
    def log_path(self) -> Path:
        return self.work_dir / self.log_filename

    @property
    def current_input_path(self) -> Path:
        return self.work_dir / self.current_input_filename

    def _parameters(self) -> dict:
        return {
            "input_dir": str(self.input_dir),
            "temperature": self.temperature,
            "seed": self.seed,
            "fep_file": self.fep_file,
            "qdyn": self.qdyn,
            "launcher": self.launcher,
            "trajectories": self.trajectories,
            "stages": [item.stage.name for item in self.plan.stages],
        }

    def _new_state(self) -> dict:
        return {
            "schema_version": 1,
            "status": "pending",
            "parameters": self._parameters(),
            "completed_stages": [],
            "checkpoints": {},
            "terminal_restarts": list(self.plan.terminal_restarts),
            "error": None,
        }

    def _load_state(self) -> dict:
        if not self.state_path.exists():
            return self._new_state()
        try:
            state = json.loads(self.state_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RunnerConfigurationError(f"cannot read runner state {self.state_path}") from exc
        if state.get("schema_version") != 1:
            raise RunnerConfigurationError(f"unsupported runner state in {self.state_path}")
        if state.get("parameters") != self._parameters():
            raise RunnerConfigurationError(
                f"{self.state_path}: run parameters differ from the existing resumable run"
            )
        expected = [item.stage.name for item in self.plan.stages]
        completed = state.get("completed_stages", [])
        if completed != expected[: len(completed)]:
            raise RunnerConfigurationError(
                f"{self.state_path}: completed stages are not a valid execution prefix"
            )
        return state

    def _write_state(self, state: dict) -> None:
        partial = self.state_path.with_name(self.state_path.name + ".partial")
        with partial.open("w", encoding="utf-8") as stream:
            json.dump(state, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(partial, self.state_path)

    def _run_qdyn(self, stage: PlannedStage, log_stream) -> None:
        command = launcher_command(self.launcher, self.qdyn, self.current_input_filename)
        log_stream.write(f"\n===== START {stage.stage.name} =====\n")
        log_stream.write(f"command: {shlex.join(command)}\n")
        log_stream.flush()
        try:
            process = subprocess.Popen(  # noqa: S603 - command is an explicit argument list
                command,
                cwd=self.work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as exc:
            raise StageExecutionError(f"could not start {stage.stage.name}: {exc}") from exc

        terminated_normally = False
        assert process.stdout is not None
        for line in process.stdout:
            log_stream.write(line)
            if "terminated normally" in line.lower():
                terminated_normally = True
        return_code = process.wait()
        log_stream.write(
            f"===== END {stage.stage.name}: returncode={return_code}, "
            f"normal_termination={terminated_normally} =====\n"
        )
        log_stream.flush()
        if return_code != 0 or not terminated_normally:
            raise StageExecutionError(
                f"{stage.stage.name} failed: return code {return_code}, "
                f"normal termination marker={terminated_normally}"
            )

    def _parent_checkpoint(self, stage: Stage, state: dict) -> str | None:
        if stage.restart is None:
            return None
        checkpoint = state["checkpoints"].get(stage.restart)
        if checkpoint is None:
            raise RunnerConfigurationError(
                f"cannot run {stage.name}: no live checkpoint for {stage.restart}"
            )
        if not (self.work_dir / checkpoint).is_file():
            raise RunnerConfigurationError(
                f"cannot run {stage.name}: checkpoint {checkpoint} is missing"
            )
        return checkpoint

    def _parent_has_pending_children(self, parent: str, state: dict) -> bool:
        completed_names = set(state["completed_stages"])
        final_by_name = {item.stage.name: item.stage.final for item in self.plan.stages}
        completed_finals = {final_by_name[name] for name in completed_names}
        return any(child not in completed_finals for child in self.plan.children.get(parent, ()))

    def _promote_checkpoint(self, planned: PlannedStage, state: dict) -> None:
        next_path = self.work_dir / planned.next_checkpoint
        checkpoint_path = self.work_dir / planned.checkpoint
        if not next_path.is_file():
            raise StageExecutionError(
                f"{planned.stage.name} terminated normally but did not write {planned.next_checkpoint}"
            )
        os.replace(next_path, checkpoint_path)
        state["checkpoints"][planned.stage.final] = planned.checkpoint

        parent = planned.stage.restart
        if parent is None or self._parent_has_pending_children(parent, state):
            return
        old_checkpoint = state["checkpoints"].pop(parent, None)
        if old_checkpoint and old_checkpoint != planned.checkpoint:
            old_path = self.work_dir / old_checkpoint
            if old_checkpoint not in state["checkpoints"].values() and old_path.exists():
                old_path.unlink()

    def _publish_terminal_restarts(self, state: dict) -> None:
        for logical_restart in self.plan.terminal_restarts:
            checkpoint = state["checkpoints"].get(logical_restart)
            if checkpoint is None:
                raise RunnerConfigurationError(
                    f"terminal checkpoint for {logical_restart} is not available"
                )
            source = self.work_dir / checkpoint
            destination = self.work_dir / logical_restart
            if source != destination:
                os.replace(source, destination)
            state["checkpoints"][logical_restart] = logical_restart

    def run_dynamics(self) -> dict:
        """Run or resume all qdyn stages, returning the persisted run state."""
        self.work_dir.mkdir(parents=True, exist_ok=True)
        state = self._load_state()
        if state["status"] in {"dynamics_complete", "complete"}:
            return state

        state["status"] = "running"
        state["error"] = None
        self._write_state(state)
        completed_count = len(state["completed_stages"])

        with self.log_path.open("a", encoding="utf-8") as log_stream:
            for planned in self.plan.stages[completed_count:]:
                try:
                    parent_checkpoint = self._parent_checkpoint(planned.stage, state)
                    next_path = self.work_dir / planned.next_checkpoint
                    if next_path.exists():
                        next_path.unlink()
                    rendered = render_runtime_input(
                        planned,
                        self.input_dir,
                        self.temperature,
                        self.seed,
                        self.fep_file,
                        parent_checkpoint,
                        trajectory=self.trajectories,
                    )
                    self.current_input_path.write_text(rendered, encoding="utf-8")
                    self._run_qdyn(planned, log_stream)
                    state["completed_stages"].append(planned.stage.name)
                    self._promote_checkpoint(planned, state)
                    self._write_state(state)
                except Exception as exc:
                    state["status"] = "failed"
                    state["error"] = str(exc)
                    self._write_state(state)
                    raise

        self._publish_terminal_restarts(state)
        state["status"] = "dynamics_complete"
        state["error"] = None
        self._write_state(state)
        return state


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse runner CLI arguments. Execution is added by the runner integration layer."""
    parser = argparse.ArgumentParser(prog="qligfep-run")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--work-dir", default=".")
    parser.add_argument("--temperature", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--fep-file", default="FEP1.fep")
    parser.add_argument("--qdyn", required=True)
    parser.add_argument("--qfep", required=True)
    parser.add_argument(
        "--launcher",
        default="",
        help="Optional platform launcher prefix, for example 'mpirun -n 8'",
    )
    parser.add_argument("--trajectories", action="store_true")
    return parser.parse_args(argv)


def launcher_command(launcher: str, executable: str, *arguments: str) -> list[str]:
    """Build a subprocess argument list without shell interpolation."""
    return [*shlex.split(launcher), executable, *arguments]
