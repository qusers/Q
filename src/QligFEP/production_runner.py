"""Platform-agnostic production runner for windowed Q FEP calculations.

The runner executes generated Q input definitions without copying every input,
topology, and FEP file into each replicate directory.  It discovers the restart
dependency graph, renders one transient input at a time, and assigns stages to a
bounded set of rolling checkpoint lanes.
"""

from __future__ import annotations

import argparse
import re
import shlex
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

_SECTION_RE = re.compile(r"^\s*\[([^]]+)]\s*$")
_EQ_NUMBER_RE = re.compile(r"^eq(\d+)$")


class RunnerConfigurationError(ValueError):
    """Raised when generated inputs do not describe a runnable stage graph."""


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
    return "\n".join(output) + "\n"


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
