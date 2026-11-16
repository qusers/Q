"""Validated configuration objects: execution profiles and simulation protocols.

A *profile* says where and how a job runs (scheduler, engine, modules, resources).
A *protocol* says what is simulated (method, lambda schedule, stage lengths, qdyn
input overrides). Keeping them apart lets the same science run on any machine.
"""

from __future__ import annotations

import difflib
import os
import re
from dataclasses import dataclass, field
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any

EXECUTIONS = ("slurm-array", "local")
ENGINES = ("qdynp", "qdyn")
METHODS = ("equilibrium", "nonequilibrium")
TIMESTEPS = ("1fs", "2fs")
SAMPLINGS = ("linear", "sigmoidal", "exponential", "reverse_exponential")
NEQ_SCHEDULES = ("sigmoidal", "linear")
EQUILIBRATION_STAGES = ("eq1", "eq2", "eq3", "eq4", "eq5")
STAGES = {
    "equilibrium": EQUILIBRATION_STAGES + ("production",),
    "nonequilibrium": EQUILIBRATION_STAGES + ("relax", "spacing", "switch"),
}
ALL_STAGES = tuple(dict.fromkeys(STAGES["equilibrium"] + STAGES["nonequilibrium"]))

TOP_KEYS = ("version", "default_profile", "default_protocol", "profiles", "protocols")
PROFILE_KEYS = (
    "extends",
    "description",
    "execution",
    "engine",
    "q_bin",
    "modules",
    "setup",
    "mpi_command",
    "resources",
    "slurm",
    "systems",
)
RESOURCE_KEYS = ("nodes", "ntasks", "time", "mem_per_cpu", "max_parallel_replicates", "neq_cores")
SLURM_KEYS = ("account", "partition", "export", "extra_directives")
SYSTEMS = ("water", "protein", "vacuum")
PROTOCOL_KEYS = (
    "extends",
    "description",
    "method",
    "timestep",
    "temperature",
    "replicates",
    "sphere_radius",
    "equilibrium",
    "nonequilibrium",
    "stages",
)
EQUILIBRIUM_KEYS = ("windows", "sampling", "start")
NONEQUILIBRIUM_KEYS = ("switch_pairs", "schedule", "steepness")
STAGE_KEYS = ("time", "steps", "inp")

#: Values every profile starts from before its own keys are applied.
PROFILE_BASE = {
    "description": "",
    "execution": "slurm-array",
    "engine": "qdynp",
    "q_bin": None,
    "modules": [],
    "setup": "",
    "mpi_command": None,
    "resources": {
        "nodes": 1,
        "ntasks": 8,
        "time": "0-12:00:00",
        "mem_per_cpu": "auto",
        "max_parallel_replicates": 1,
        "neq_cores": "auto",
    },
    "slurm": {"account": None, "partition": None, "export": "ALL", "extra_directives": []},
    "systems": {},
}

#: Memory per core (MB) requested when a profile leaves ``mem_per_cpu: auto``.
DEFAULT_MEM_PER_CPU = {"equilibrium": 512, "nonequilibrium": 1700}

BUILD_INFO_NAME = "build-info.json"
_TIME_UNITS_FS = {"fs": 1.0, "ps": 1e3, "ns": 1e6}


class ConfigError(ValueError):
    """A configuration file is invalid or asks for something that cannot be done."""


def check_keys(mapping: Any, allowed: tuple[str, ...], where: str) -> None:
    """Reject unknown keys, suggesting the closest valid one."""
    if not isinstance(mapping, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(mapping).__name__} ({mapping!r})")
    for key in mapping:
        if key not in allowed:
            close = difflib.get_close_matches(str(key), allowed, n=1)
            hint = f" Did you mean '{close[0]}'?" if close else ""
            raise ConfigError(f"{where}: unknown key '{key}'.{hint} Valid keys: {', '.join(allowed)}.")


def _choice(value: Any, choices: tuple[str, ...], where: str) -> str:
    if value not in choices:
        raise ConfigError(f"{where}: {value!r} is not one of {', '.join(choices)}.")
    return value


def _positive_int(value: Any, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigError(f"{where}: expected a positive integer, got {value!r}.")
    return value


def parse_duration_fs(value: Any, where: str) -> float:
    """Parse a duration such as ``100ps`` or ``2.5 ns`` into femtoseconds."""
    match = re.fullmatch(r"\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*(fs|ps|ns)\s*", str(value))
    if not match:
        raise ConfigError(f"{where}: {value!r} is not a duration. Use a number with a unit: fs, ps or ns.")
    return float(match.group(1)) * _TIME_UNITS_FS[match.group(2)]


def deep_merge(base: dict, override: dict) -> dict:
    """Return ``base`` with ``override`` applied; mappings merge, everything else replaces."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


# ----------------------------------------------------------------------------------
# Profiles
# ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class Resources:
    nodes: int
    ntasks: int
    time: str
    mem_per_cpu: int | str | None  # MB, "auto" (method default) or None (omit the directive)
    max_parallel_replicates: int
    neq_cores: int | str  # "auto" or an explicit core count for non-equilibrium jobs


@dataclass(frozen=True)
class SlurmOptions:
    account: str | None
    partition: str | None
    export: str | None
    extra_directives: tuple[str, ...]


@dataclass(frozen=True)
class Profile:
    """A resolved execution profile."""

    name: str
    description: str
    execution: str
    engine: str
    q_bin: Path
    q_bin_source: str
    modules: tuple[str, ...]
    setup: str
    mpi_command: str | None
    resources: Resources
    slurm: SlurmOptions
    sources: tuple[str, ...] = ()

    @property
    def is_local(self) -> bool:
        return self.execution == "local"

    @property
    def use_mpi(self) -> bool:
        return self.engine == "qdynp"

    def binary(self, name: str) -> Path:
        return self.q_bin / name

    @property
    def runfile_name(self) -> str:
        return f"run{self.name}.sh"

    def environment_block(self) -> str:
        """Shell lines preparing the job environment: free-form setup, then module loads."""
        lines = [line for line in self.setup.splitlines() if line.strip()]
        lines += [f"module load {module}" for module in self.modules]
        return "\n".join(lines) + "\n"

    def mpi_binding(self) -> str:
        """Core-binding flags for mpirun; only meaningful inside a scheduler cpuset."""
        return "" if self.is_local else "--map-by core --bind-to core"

    def launch_prefix(self) -> str:
        """Command prefix that runs the profile's engine; ``$ncores`` is a shell variable."""
        if not self.use_mpi:
            return "time"
        if self.mpi_command:
            return f"time {self.mpi_command}"
        binding = self.mpi_binding()
        return "time mpirun -n $ncores" + (f" {binding}" if binding else "")

    def mem_per_cpu(self, method: str) -> int | None:
        if self.resources.mem_per_cpu == "auto":
            return DEFAULT_MEM_PER_CPU[method]
        return self.resources.mem_per_cpu

    def neq_cores(self, n_switches: int) -> int:
        """Cores for a non-equilibrium job that fires ``n_switches`` concurrent switches.

        ``auto`` grows the allocation so every switch runs in one wave and never drops
        below ``ntasks``. An explicit ``neq_cores`` is used as is; with fewer cores than
        switches the run script fires them in successive waves.
        """
        if self.resources.neq_cores == "auto":
            return max(n_switches, self.resources.ntasks)
        return int(self.resources.neq_cores)

    def sbatch_header(self, method: str, jobname: str, array_size: int, ntasks: int | None = None) -> str:
        """``#SBATCH`` directives for this profile; a comment for local execution."""
        if self.is_local:
            return "# Local execution profile: no scheduler directives. Launch with FEP_submit.sh.\n"
        res, slurm = self.resources, self.slurm
        lines = [
            f"#SBATCH --nodes={res.nodes}",
            f"#SBATCH --ntasks-per-node={ntasks if ntasks is not None else res.ntasks}",
        ]
        mem = self.mem_per_cpu(method)
        if mem is not None:
            lines.append(f"#SBATCH --mem-per-cpu={mem}")
        if slurm.account:
            lines.append(f"#SBATCH -A {slurm.account}")
        if slurm.partition:
            lines.append(f"#SBATCH -p {slurm.partition}")
        if slurm.export:
            lines.append(f"#SBATCH --export={slurm.export}")
        lines += [
            "#              d-hh:mm:ss",
            f"#SBATCH --time={res.time}",
            f"#SBATCH -J {jobname}",
            f"#SBATCH --array=1-{array_size}",
            "#SBATCH -o slurm.run%a.%N.%j.out",
        ]
        lines += [f"#SBATCH {directive}" for directive in slurm.extra_directives]
        return "\n".join(lines) + "\n"

    def as_dict(self) -> dict:
        """Plain-data view, used for the provenance snapshot and ``qligfep_config show``."""
        return {
            "name": self.name,
            "description": self.description,
            "execution": self.execution,
            "engine": self.engine,
            "q_bin": str(self.q_bin),
            "q_bin_source": self.q_bin_source,
            "modules": list(self.modules),
            "setup": self.setup,
            "mpi_command": self.mpi_command,
            "resources": dict(vars(self.resources)),
            "slurm": {**vars(self.slurm), "extra_directives": list(self.slurm.extra_directives)},
            "sources": list(self.sources),
        }

    def as_legacy_dict(self) -> dict:
        """The historical ``CLUSTER_DICT`` entry, for code not yet ported to profiles."""
        qdynp, qdyn = str(self.binary("qdynp")), str(self.binary("qdyn"))
        legacy = {
            "NODES": str(self.resources.nodes),
            "NTASKS": str(self.resources.ntasks),
            "TIME": self.resources.time,
            "MODULES": self.environment_block(),
            "QPREP": str(self.binary("qprep")),
            "QFEP": str(self.binary("qfep")),
            "QCALC": str(self.binary("qcalc")),
            "EXECUTION": self.execution,
        }
        if self.is_local:
            legacy.update({"QDYN": qdynp if self.use_mpi else qdyn, "USE_MPI": self.use_mpi})
        else:
            legacy.update({"QDYN": f"qdyn={qdynp}", "QDYNP": f"qdynp={qdynp}", "QDYN_NEQ": f"qdyn={qdyn}"})
        if self.slurm.account:
            legacy["ACCOUNT"] = self.slurm.account
        return legacy


def resolve_q_bin(configured: Any, where: str, override: str | os.PathLike | None = None) -> tuple[Path, str]:
    """Pick the directory holding the Q binaries and say where the choice came from.

    Order: explicit ``override`` (CLI), ``$QLIGFEP_Q_BIN``, the profile's ``q_bin``, then
    the binaries compiled inside the imported package. The package default is what keeps
    several clones independent: each checkout runs its own engine unless told otherwise.
    """
    from QligFEP import SRC

    if override:
        return Path(override).expanduser().resolve(), "command line"
    env = os.environ.get("QLIGFEP_Q_BIN")
    if env:
        return Path(env).expanduser().resolve(), "$QLIGFEP_Q_BIN"
    if configured:
        path = Path(os.path.expandvars(str(configured))).expanduser()
        if not path.is_absolute():
            raise ConfigError(f"{where}.q_bin: '{configured}' must be an absolute path.")
        return path, "profile"
    return (SRC / "q6" / "bin" / "q6").absolute(), "package"


def read_build_info(q_bin: Path) -> dict | None:
    """Read the toolchain manifest written next to the binaries at build time, if any."""
    import json

    path = Path(q_bin) / BUILD_INFO_NAME
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise ConfigError(f"Could not read {path}: {error}") from error


def build_profile(
    name: str, raw: dict, sources: tuple[str, ...] = (), q_bin_override: str | os.PathLike | None = None
) -> Profile:
    """Validate a fully merged profile mapping and turn it into a :class:`Profile`."""
    where = f"profiles.{name}"
    check_keys(raw, PROFILE_KEYS, where)
    execution = _choice(raw["execution"], EXECUTIONS, f"{where}.execution")
    engine = _choice(raw["engine"], ENGINES, f"{where}.engine")

    resources = raw["resources"]
    check_keys(resources, RESOURCE_KEYS, f"{where}.resources")
    ntasks = resources["ntasks"]
    if ntasks == "auto":
        if execution != "local":
            raise ConfigError(
                f"{where}.resources.ntasks: 'auto' is only valid for local execution; "
                "a Slurm job needs an explicit core count."
            )
        ntasks = cpu_count() if engine == "qdynp" else 1
    ntasks = _positive_int(ntasks, f"{where}.resources.ntasks")
    if engine == "qdyn" and ntasks != 1 and execution == "local":
        raise ConfigError(f"{where}: the serial engine 'qdyn' runs on one core; set ntasks to 1 or use qdynp.")
    mem = resources["mem_per_cpu"]
    if mem not in ("auto", None):
        mem = _positive_int(mem, f"{where}.resources.mem_per_cpu")
    neq_cores = resources["neq_cores"]
    if neq_cores != "auto":
        neq_cores = _positive_int(neq_cores, f"{where}.resources.neq_cores")
    if not isinstance(resources["time"], str) or not resources["time"].strip():
        raise ConfigError(f"{where}.resources.time: expected a quoted Slurm time such as \"0-08:00:00\".")

    slurm = raw["slurm"]
    check_keys(slurm, SLURM_KEYS, f"{where}.slurm")
    extra = slurm["extra_directives"] or []
    if not isinstance(extra, list) or not all(isinstance(item, str) for item in extra):
        raise ConfigError(f"{where}.slurm.extra_directives: expected a list of strings such as \"--qos=short\".")
    extra = [re.sub(r"^\s*#SBATCH\s+", "", item).strip() for item in extra]

    q_bin, q_bin_source = resolve_q_bin(raw["q_bin"], where, q_bin_override)

    modules = raw["modules"]
    if modules == "from-build":
        info = read_build_info(q_bin)
        if info is None:
            raise ConfigError(
                f"{where}.modules: 'from-build' needs {q_bin / BUILD_INFO_NAME}, which does not exist. "
                "Run `qligfep_config build-info` in the environment the binaries were compiled in, "
                "or list the modules explicitly."
            )
        modules = info.get("modules", [])
    if modules is None:
        modules = []
    if not isinstance(modules, list) or not all(isinstance(item, (str, int, float)) for item in modules):
        raise ConfigError(f"{where}.modules: expected a list of module names or 'from-build'.")

    setup = raw["setup"] or ""
    if not isinstance(setup, str):
        raise ConfigError(f"{where}.setup: expected a string of shell lines.")
    mpi_command = raw["mpi_command"]
    if mpi_command is not None and not isinstance(mpi_command, str):
        raise ConfigError(f"{where}.mpi_command: expected a string such as \"srun\".")

    return Profile(
        name=name,
        description=str(raw.get("description") or ""),
        execution=execution,
        engine=engine,
        q_bin=q_bin,
        q_bin_source=q_bin_source,
        modules=tuple(str(item) for item in modules),
        setup=setup,
        mpi_command=mpi_command,
        resources=Resources(
            nodes=_positive_int(resources["nodes"], f"{where}.resources.nodes"),
            ntasks=ntasks,
            time=resources["time"].strip(),
            mem_per_cpu=mem,
            max_parallel_replicates=_positive_int(
                resources["max_parallel_replicates"], f"{where}.resources.max_parallel_replicates"
            ),
            neq_cores=neq_cores,
        ),
        slurm=SlurmOptions(
            account=str(slurm["account"]) if slurm["account"] else None,
            partition=str(slurm["partition"]) if slurm["partition"] else None,
            export=str(slurm["export"]) if slurm["export"] else None,
            extra_directives=tuple(extra),
        ),
        sources=sources,
    )


# ----------------------------------------------------------------------------------
# Protocols
# ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class StageOverride:
    """User changes to one simulation stage."""

    time_fs: float | None = None
    steps: int | None = None
    inp: dict = field(default_factory=dict)

    def n_steps(self, stepsize_fs: float, where: str) -> int | None:
        """Length of the stage in MD steps, or None when the protocol leaves it alone."""
        if self.steps is not None:
            return self.steps
        if self.time_fs is None:
            return None
        steps = self.time_fs / float(stepsize_fs)
        if abs(steps - round(steps)) > 1e-6 or round(steps) < 1:
            raise ConfigError(
                f"{where}: a duration of {self.time_fs:g} fs is not a whole number of "
                f"{stepsize_fs:g} fs steps."
            )
        return int(round(steps))

    def as_dict(self) -> dict:
        data: dict = {}
        if self.time_fs is not None:
            data["time"] = f"{self.time_fs:g}fs"
        if self.steps is not None:
            data["steps"] = self.steps
        if self.inp:
            data["inp"] = self.inp
        return data


@dataclass(frozen=True)
class Protocol:
    """A resolved simulation protocol."""

    name: str
    description: str
    method: str
    timestep: str
    temperature: str  # comma-separated when several temperatures are requested
    replicates: int
    sphere_radius: int
    windows: int
    sampling: str
    start: str
    switch_pairs: int
    schedule: str
    steepness: float
    stages: dict
    sources: tuple[str, ...] = ()

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "method": self.method,
            "timestep": self.timestep,
            "temperature": self.temperature,
            "replicates": self.replicates,
            "sphere_radius": self.sphere_radius,
            "equilibrium": {"windows": self.windows, "sampling": self.sampling, "start": self.start},
            "nonequilibrium": {
                "switch_pairs": self.switch_pairs,
                "schedule": self.schedule,
                "steepness": self.steepness,
            },
            "stages": {name: stage.as_dict() for name, stage in self.stages.items()},
            "sources": list(self.sources),
        }


def _build_stage(name: str, raw: Any, where: str) -> StageOverride:
    check_keys(raw, STAGE_KEYS, where)
    if "time" in raw and "steps" in raw:
        raise ConfigError(f"{where}: give either 'time' or 'steps', not both.")
    inp = raw.get("inp") or {}
    if not isinstance(inp, dict):
        raise ConfigError(f"{where}.inp: expected a mapping of section -> key -> value.")
    for section, keys in inp.items():
        if not isinstance(keys, dict):
            raise ConfigError(f"{where}.inp.{section}: expected a mapping of key -> value.")
        for key, value in keys.items():
            if isinstance(value, (dict, list)):
                raise ConfigError(f"{where}.inp.{section}.{key}: expected a single value.")
        if "steps" in keys and section.lower() == "md" and ("time" in raw or "steps" in raw):
            raise ConfigError(f"{where}: the stage length is set twice (inp.MD.steps and time/steps).")
    return StageOverride(
        time_fs=parse_duration_fs(raw["time"], f"{where}.time") if "time" in raw else None,
        steps=_positive_int(raw["steps"], f"{where}.steps") if "steps" in raw else None,
        inp={str(section): dict(keys) for section, keys in inp.items()},
    )


def validate_protocol_definition(name: str, raw: Any, source: str) -> None:
    """Check the keys of one protocol as written in one file (before inheritance)."""
    where = f"{source}: protocols.{name}"
    check_keys(raw, PROTOCOL_KEYS, where)
    if "equilibrium" in raw:
        check_keys(raw["equilibrium"], EQUILIBRIUM_KEYS, f"{where}.equilibrium")
    if "nonequilibrium" in raw:
        check_keys(raw["nonequilibrium"], NONEQUILIBRIUM_KEYS, f"{where}.nonequilibrium")
    stages = raw.get("stages") or {}
    check_keys(stages, ALL_STAGES, f"{where}.stages")


def build_protocol(name: str, raw: dict, own: dict, sources: tuple[str, ...] = ()) -> Protocol:
    """Validate a merged protocol mapping and turn it into a :class:`Protocol`.

    Args:
        name: Protocol name.
        raw: The protocol after inheritance and layering.
        own: The keys written under this protocol's own name (no inherited ones); used to
            reject settings that belong to the method the protocol does not use.
        sources: Files that contributed.
    """
    where = f"protocols.{name}"
    check_keys(raw, PROTOCOL_KEYS, where)
    method = _choice(raw["method"], METHODS, f"{where}.method")
    other = "nonequilibrium" if method == "equilibrium" else "equilibrium"
    if own.get(other):
        raise ConfigError(
            f"{where}: method is '{method}', but the protocol sets '{other}' options "
            f"({', '.join(own[other])}). They would be ignored; remove them or change the method."
        )
    for stage in own.get("stages") or {}:
        if stage not in STAGES[method]:
            raise ConfigError(
                f"{where}.stages.{stage}: not a stage of the '{method}' method. "
                f"Its stages are: {', '.join(STAGES[method])}."
            )

    equilibrium, nonequilibrium = raw["equilibrium"], raw["nonequilibrium"]
    check_keys(equilibrium, EQUILIBRIUM_KEYS, f"{where}.equilibrium")
    check_keys(nonequilibrium, NONEQUILIBRIUM_KEYS, f"{where}.nonequilibrium")

    start = str(equilibrium["start"])
    start = {"1.0": "1", "1": "1", "0.5": "0.5"}.get(start)
    if start is None:
        raise ConfigError(f"{where}.equilibrium.start: {equilibrium['start']!r} is not one of 0.5, 1.")

    temperature = raw["temperature"]
    temperatures = temperature if isinstance(temperature, list) else str(temperature).split(",")
    try:
        temperature = ",".join(f"{float(str(t).strip()):g}" for t in temperatures)
    except ValueError:
        raise ConfigError(f"{where}.temperature: expected a number or a list of numbers.") from None

    try:
        steepness = float(nonequilibrium["steepness"])
    except (TypeError, ValueError):
        raise ConfigError(f"{where}.nonequilibrium.steepness: expected a number.") from None

    stages_raw = raw.get("stages") or {}
    check_keys(stages_raw, ALL_STAGES, f"{where}.stages")
    stages = {
        stage: _build_stage(stage, stage_raw or {}, f"{where}.stages.{stage}")
        for stage, stage_raw in stages_raw.items()
    }

    return Protocol(
        name=name,
        description=str(raw.get("description") or ""),
        method=method,
        timestep=_choice(str(raw["timestep"]), TIMESTEPS, f"{where}.timestep"),
        temperature=temperature,
        replicates=_positive_int(raw["replicates"], f"{where}.replicates"),
        sphere_radius=_positive_int(raw["sphere_radius"], f"{where}.sphere_radius"),
        windows=_positive_int(equilibrium["windows"], f"{where}.equilibrium.windows"),
        sampling=_choice(equilibrium["sampling"], SAMPLINGS, f"{where}.equilibrium.sampling"),
        start=start,
        switch_pairs=_positive_int(nonequilibrium["switch_pairs"], f"{where}.nonequilibrium.switch_pairs"),
        schedule=_choice(nonequilibrium["schedule"], NEQ_SCHEDULES, f"{where}.nonequilibrium.schedule"),
        steepness=steepness,
        stages=stages,
        sources=sources,
    )
