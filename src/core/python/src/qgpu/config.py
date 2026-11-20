"""qgpu.Config: all settings of one Engine, with the same names and meaning as the .inp keys.

    config = Config.from_inp("eq1.inp")   # or Config(topology=..., steps=..., stepsize=..., temperature=...)
    config = replace(config, steps=1000)  # or engine.update_config(steps=1000)
"""

from dataclasses import dataclass, fields
from pathlib import Path

from . import _native
from .errors import InputError
from .restraints import AngRestraint, DistRestraint, PosRestraint, SeqRestraint, WallRestraint

# Config field -> restraint class
_RESTRAINTS = {
    "seq_restraints": SeqRestraint,
    "dist_restraints": DistRestraint,
    "ang_restraints": AngRestraint,
    "wall_restraints": WallRestraint,
    "pos_restraints": PosRestraint,
}

# Settings used to build the Engine. They cannot change afterwards: create a new Engine.
FIXED = frozenset({
    "topology", "fep", "restart", "seed", "initial_temperature",
    "shell_radius", "exclude_bonded",
    "solute_solute", "solvent_solvent", "solute_solvent", "q_atom", "lrf", "lrf_cutoff",
    "radial_force", "polarisation", "polarisation_force", "charge_correction",
    "non_bond", "hmr", "hmr_target_mass",
    "solute_constraint_algorithm", "solvent_constraint_algorithm",
})


def _onoff(value):
    return "on" if value else "off"


def _abs(path):
    return str(Path(path).expanduser().resolve()) if path else ""


# Config field -> (.inp section, key, type). The C++ parser is the only place with
# defaults: unset fields are filled in by it, and the FIXED ones are passed to the engine this way.
_INP_KEYS = {
    "seed": ("md", "random_seed", int),
    "initial_temperature": ("md", "initial_temperature", float),
    "lrf": ("md", "lrf", bool),
    "hmr": ("md", "hmr", bool),
    "hmr_target_mass": ("md", "hmr_target_mass", float),
    "steps": ("md", "steps", int),
    "stepsize": ("md", "stepsize", float),
    "temperature": ("md", "temperature", float),
    "thermostat": ("md", "thermostat", str),
    "bath_coupling": ("md", "bath_coupling", float),
    "constrain_solvent": ("md", "constrain_solvent", bool),
    "constrain_solute": ("md", "constrain_solute", bool),
    "constrain_hydrogens": ("md", "constrain_hydrogens", bool),
    "separate_scaling": ("md", "separate_scaling", bool),
    "solute_solute": ("cut-offs", "solute_solute", float),
    "solvent_solvent": ("cut-offs", "solvent_solvent", float),
    "solute_solvent": ("cut-offs", "solute_solvent", float),
    "q_atom": ("cut-offs", "q_atom", float),
    "lrf_cutoff": ("cut-offs", "lrf", float),
    "shell_radius": ("sphere", "shell_radius", float),
    "exclude_bonded": ("sphere", "exclude_bonded", bool),
    "shell_force": ("sphere", "shell_force", float),
    "radial_force": ("solvent", "radial_force", float),
    "polarisation": ("solvent", "polarisation", bool),
    "polarisation_force": ("solvent", "polarisation_force", float),
    "charge_correction": ("solvent", "charge_correction", bool),
    "non_bond": ("intervals", "non_bond", int),
    "output": ("intervals", "output", int),
    "energy": ("intervals", "energy", int),
    "trajectory": ("intervals", "trajectory", int),
}
# One .inp key: constraint_algorithm <solute> <solvent>
_ALGORITHMS = ("solute_constraint_algorithm", "solvent_constraint_algorithm")
# Fields that the parser fills in when they are None
_RESOLVED = (*_INP_KEYS, *_ALGORITHMS, "lambdas")
# initial_temperature has no default: None means "not given" (as in Fortran Q)
_DEFAULTED = tuple(n for n in _RESOLVED if n != "initial_temperature")


def _inp_value(value, kind):
    if isinstance(value, str):
        return value  # as in an .inp file, e.g. "on"
    if kind is bool:
        return _onoff(value)
    return repr(kind(value))


@dataclass(frozen=True, kw_only=True)
class Config:
    """All settings of one engine, with the same names and meaning as the .inp keys.

    Leave a setting out (None) to get the .inp default: the C++ parser fills it in,
    so Config(topology=..., steps=...) behaves like an .inp file with only those keys.
    After creation every field holds the value the engine uses.

    Settings in FIXED are used when the Engine is created; all others can be changed
    between runs with Engine.update_config().
    """

    # ---- fixed: used when the Engine is created
    topology: str
    fep: str = ""
    restart: str = ""                   # "" = fresh start from the topology coordinates
    seed: int | None = None             # for the velocities from initial_temperature
    # Given: new Maxwell velocities at the start, also with a restart (as in Fortran Q).
    # None: keep the restart velocities; a fresh start (no restart) then fails.
    initial_temperature: float | None = None
    shell_radius: float | None = None
    exclude_bonded: bool | None = None
    solute_solute: float | None = None
    solvent_solvent: float | None = None
    solute_solvent: float | None = None
    q_atom: float | None = None
    lrf: bool | None = None
    lrf_cutoff: float | None = None
    radial_force: float | None = None
    polarisation: bool | None = None
    polarisation_force: float | None = None
    charge_correction: bool | None = None
    non_bond: int | None = None
    hmr: bool | None = None
    hmr_target_mass: float | None = None
    # shake | lincs; set both or neither
    solute_constraint_algorithm: str | None = None
    solvent_constraint_algorithm: str | None = None  # shake | lincs | settle

    # ---- changeable with Engine.update_config()
    steps: int | None = None
    stepsize: float | None = None
    temperature: float | None = None
    thermostat: str | None = None
    bath_coupling: float | None = None
    constrain_solvent: bool | None = None
    constrain_solute: bool | None = None
    constrain_hydrogens: bool | None = None
    separate_scaling: bool | None = None
    shell_force: float | None = None
    lambdas: tuple | None = None
    seq_restraints: tuple = ()          # of qgpu.SeqRestraint
    dist_restraints: tuple = ()         # of qgpu.DistRestraint
    ang_restraints: tuple = ()          # of qgpu.AngRestraint
    wall_restraints: tuple = ()         # of qgpu.WallRestraint
    pos_restraints: tuple = ()          # of qgpu.PosRestraint
    output: int | None = None           # log interval
    energy: int | None = None           # energy file interval; > 0 needs energy_file
    # trajectory interval; > 0 needs trajectory_file
    trajectory: int | None = None
    energy_file: str = ""               # "" = none; files are overwritten by every run()
    trajectory_file: str = ""
    final_file: str = ""                # restart file written at the end
    trajectory_atoms: str = ""          # "" = default mask

    def __post_init__(self):
        if any(getattr(self, name) is None for name in _DEFAULTED):
            self._resolve()

        object.__setattr__(self, "lambdas", tuple(float(x)
                           for x in self.lambdas))
        for name, cls in _RESTRAINTS.items():
            items = tuple(getattr(self, name))
            for r in items:
                if not isinstance(r, cls):
                    raise TypeError(
                        f"{name} must contain qgpu.{cls.__name__}, got {type(r).__name__}")
            object.__setattr__(self, name, items)

        if self.steps < 0:
            raise InputError("steps must be >= 0")
        if self.stepsize <= 0:
            raise InputError("stepsize must be > 0")
        if self.energy > 0 and not self.energy_file:
            raise InputError("energy > 0 needs energy_file (or set energy=0)")
        if self.trajectory > 0 and not self.trajectory_file:
            raise InputError(
                "trajectory > 0 needs trajectory_file (or set trajectory=0)")

    def _resolve(self):
        """Fill in every setting from the C++ parser, given the ones that are set."""
        solute, solvent = (getattr(self, name) for name in _ALGORITHMS)
        if (solute is None) != (solvent is None):
            raise InputError("Set solute_constraint_algorithm and solvent_constraint_algorithm "
                             "together (one .inp key: constraint_algorithm <solute> <solvent>).")
        lists = {}
        if self.lambdas is not None:
            lists["lambdas"] = [[_inp_value(x, float) for x in self.lambdas]]
        values = _native.load_any().resolve_settings(self._sections(), lists)
        for name in _RESOLVED:
            object.__setattr__(self, name, values[name])

    @classmethod
    def from_inp(cls, path, **overrides):
        """Read a Q .inp file; `overrides` replace values from the file.

        Relative paths in the file are resolved against the file's folder.
        """
        values = _native.load_any().read_inp(str(path))
        for name, rcls in _RESTRAINTS.items():
            values[name] = tuple(rcls(*row) for row in values[name])
        return cls(**{**values, **overrides})

    def _sections(self):
        """The settings that are set, as .inp sections."""
        names = [n for n in _INP_KEYS if getattr(self, n) is not None]
        sections = {}
        for name in names:
            section, key, kind = _INP_KEYS[name]
            sections.setdefault(section, {})[key] = _inp_value(
                getattr(self, name), kind)
        if self.solute_constraint_algorithm is not None:
            sections.setdefault("md", {})["constraint_algorithm"] = \
                f"{self.solute_constraint_algorithm} {self.solvent_constraint_algorithm}"
        files = {"fep": _abs(self.fep), "restart": _abs(self.restart)}
        files = {key: path for key, path in files.items() if path}
        if files:
            # fep: one lambda by default; restart: start state
            sections["files"] = files
        return sections


_FIELDS = frozenset(f.name for f in fields(Config))
assert FIXED <= _FIELDS and set(_RESOLVED) <= _FIELDS
