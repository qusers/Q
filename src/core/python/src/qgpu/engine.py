"""In-process QGPU engine: load the system once, run many times.

    config = Config.from_inp("eq1.inp")   # or Config(topology=..., steps=..., stepsize=..., temperature=...)
    engine = Engine(config)
    engine.run()
    engine.update_config(lambdas=(0.5, 0.5), energy_file="fep_050.en")
    engine.run()                          # continues from the current coordinates/velocities
"""

from collections.abc import Sequence
from dataclasses import replace
import inspect
from typing import Any

from . import _native
from .config import _FIELDS, FIXED, Config, _abs
from .errors import InputError
from .restraints import AngRestraint, DistRestraint, PosRestraint, SeqRestraint, WallRestraint


class _Keep:
    """Default of every update_config() parameter: keep the current value."""

    def __repr__(self):
        return "<keep>"


_KEEP: Any = _Keep()


class Engine:
    def __init__(self, config: Config, device=0, precision="double"):
        """precision: "double" or "spfp" (single-precision forces, fixed-point accumulation)."""
        if not isinstance(config, Config):
            raise TypeError(
                f"config must be a qgpu.Config, got {type(config).__name__}")
        self._mod = _native.load(precision)
        self._precision = precision
        lists = {"lambdas": [[str(x) for x in config.lambdas]]}
        self._native = self._mod.Engine(_abs(config.topology), _abs(config.fep),
                                        config._sections(), lists, device)
        self._base = self._native.base_stage()
        # fresh start: velocities generated here
        self._native.configure(self._to_native(config))
        self._config = config

    @property
    def config(self):
        """The settings in use. Change them with update_config()."""
        return self._config

    @property
    def precision(self):
        return self._precision

    @property
    def n_atoms(self):
        return self._native.n_atoms

    def update_config(
        self, *,
        # [MD]
        steps: int = _KEEP,
        stepsize: float = _KEEP,
        temperature: float = _KEEP,
        thermostat: str = _KEEP,
        bath_coupling: float = _KEEP,
        constrain_solvent: bool = _KEEP,
        constrain_solute: bool = _KEEP,
        constrain_hydrogens: bool = _KEEP,
        separate_scaling: bool = _KEEP,
        # [sphere]
        shell_force: float = _KEEP,
        # [lambdas]
        lambdas: Sequence[float] = _KEEP,
        # restraints
        seq_restraints: Sequence[SeqRestraint] = _KEEP,
        dist_restraints: Sequence[DistRestraint] = _KEEP,
        ang_restraints: Sequence[AngRestraint] = _KEEP,
        wall_restraints: Sequence[WallRestraint] = _KEEP,
        pos_restraints: Sequence[PosRestraint] = _KEEP,
        # [intervals]
        output: int = _KEEP,
        energy: int = _KEEP,
        trajectory: int = _KEEP,
        # [files]
        energy_file: str = _KEEP,
        trajectory_file: str = _KEEP,
        final_file: str = _KEEP,
        trajectory_atoms: str = _KEEP,
    ):
        """Change settings for the next run(), e.g. update_config(lambdas=(0.5, 0.5)).

        Settings that are not passed keep their current value. They have the same
        meaning as in Config. Settings used to create the Engine (shell_radius, lrf,
        restart, ...) are not parameters here: create a new Engine to change them.
        Everything is checked before the engine is touched.
        """
        changes = {name: value for name, value in locals().items()
                   if name != "self" and value is not _KEEP}
        config = replace(self._config, **changes)  # validates
        self._native.configure(self._to_native(config))
        self._config = config

    def run(self):
        """Run config.steps steps from the current state."""
        self._native.run()

    def get_state(self):
        return self._native.get_state()

    def set_state(self, state):
        self._native.set_state(
            state["coords"], state["velocities"], list(state["theta_corr"]))

    def _to_native(self, config: Config):
        if len(config.lambdas) != self._native.n_states:
            raise InputError(
                f"Expected {self._native.n_states} lambdas, got {len(config.lambdas)}")

        # copy: FIXED fields stay as the engine loaded them
        s = self._mod.StageSettings(self._base)
        md = s.md
        md.steps = config.steps
        md.stepsize = config.stepsize
        md.temperature = config.temperature
        md.thermostat = config.thermostat
        md.bath_coupling = config.bath_coupling
        md.constraint_solvent = md.shake_solvent = config.constrain_solvent
        md.constraint_solute = md.shake_solute = config.constrain_solute
        md.constraint_hydrogens = md.shake_hydrogens = config.constrain_hydrogens
        md.separate_scaling = config.separate_scaling
        md.shell_force = config.shell_force
        md.output = config.output
        md.energy = config.energy
        md.trajectory = config.trajectory

        s.lambdas = list(config.lambdas)
        # Convert to this engine's module; assign whole lists (vectors are copied).
        s.restrseqs = [r._native(self._mod) for r in config.seq_restraints]
        s.restrdists = [r._native(self._mod) for r in config.dist_restraints]
        s.restrangs = [r._native(self._mod) for r in config.ang_restraints]
        s.restrwalls = [r._native(self._mod) for r in config.wall_restraints]
        s.restrspos = [r._native(self._mod) for r in config.pos_restraints]

        s.outputs.energy_file = _abs(config.energy_file)
        s.outputs.final_file = _abs(config.final_file)
        s.outputs.trajectory_file = _abs(config.trajectory_file)
        s.outputs.trajectory_atoms = config.trajectory_atoms
        return s


# update_config() lists the changeable Config fields by hand, for autocompletion.
# Fail at import if the two drift apart.
_UPDATABLE = set(inspect.signature(Engine.update_config).parameters) - {"self"}
if _UPDATABLE != _FIELDS - FIXED:
    raise RuntimeError(
        "Engine.update_config() parameters do not match the changeable Config fields: "
        f"missing {sorted(_FIELDS - FIXED - _UPDATABLE)}, extra {sorted(_UPDATABLE - (_FIELDS - FIXED))}")
