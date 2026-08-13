"""Equilibration stage configurations for FEP simulations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .md_template import MDParameters


@dataclass
class EquilibrationConfig:
    """Data class for a single equilibration stage. These parameters are used to
    override the defaults in MDParameters (see method get_equilibration_configs).

    To render a full configuration:

    >>> from QligFEP.templates import get_equilibration_configs, render_md_input
    >>> configs = get_equilibration_configs("2fs", shell_radius=25)
    >>> eq5 = configs[4]
    >>> content = render_md_input(
    ...     params=eq5.params,
    ...     lambda1="0.500",
    ...     lambda2="0.500",
    ...     trajectory_file="eq5.dcd",
    ...     final_file="eq5.re",
    ...     restart_file="eq4.re",
    ... )

    Attributes:
        name: Stage name (eq1, eq2, etc.)
        params: MD simulation parameters
        sequence_restraint_force: Force constant for sequence restraints (eq1-eq4)
        distance_restraint_force: Force constant for distance restraints. If None,
            use the production dr_force value (eq5 uses production force).
        use_water_restraint: If True, use WATER_RESTRAINT placeholder instead of
            atom-range sequence restraints (eq5 only)
    """

    name: str
    params: MDParameters
    sequence_restraint_force: float
    distance_restraint_force: float | None
    use_water_restraint: bool = False


# ======================================================================
# Equilibration parameter sets used to override defaults in MDParameters
# ======================================================================

# eq1 is identical for both timesteps (fixed small timestep for initial equilibration)
# Minimization is enabled for eq1 to relax initial geometries before MD
EQ1_PARAMS = dict(
    steps=5000,  # with minimization, I expect we can use fewer steps here
    stepsize=0.2,
    temperature=1,
    bath_coupling=0.2,
    shake_hydrogens=False,
    interval_output=5,
    minimize=True,
)

# 2fs timestep variants
EQ2_2FS_PARAMS = dict(
    steps=5000,
    stepsize=2.0,
    temperature=50,
    bath_coupling=2.0,
    shake_hydrogens=True,
    interval_output=5,
)

EQ3_2FS_PARAMS = dict(
    steps=5000,
    stepsize=2.0,
    temperature=150,
    bath_coupling=2.0,
    shake_hydrogens=True,
    interval_output=5,
)

EQ4_2FS_PARAMS = dict(
    steps=5000,
    stepsize=2.0,
    temperature=275,
    bath_coupling=2.0,
    shake_hydrogens=True,
    interval_output=5,
)

EQ5_2FS_PARAMS = dict(
    steps=50000,
    stepsize=2.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    shake_hydrogens=True,
    interval_output=25,
)

# 1fs timestep variants
EQ2_1FS_PARAMS = dict(
    steps=10000,
    stepsize=1.0,
    temperature=50,
    bath_coupling=1.0,
    shake_hydrogens=False,
    interval_output=5,
)

EQ3_1FS_PARAMS = dict(
    steps=10000,
    stepsize=1.0,
    temperature=150,
    bath_coupling=1.0,
    shake_hydrogens=False,
    interval_output=5,
)

EQ4_1FS_PARAMS = dict(
    steps=10000,
    stepsize=1.0,
    temperature=275,
    bath_coupling=1.0,
    shake_hydrogens=False,
    interval_output=5,
)

EQ5_1FS_PARAMS = dict(
    steps=100000,
    stepsize=1.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    shake_hydrogens=False,
    interval_output=25,
)

# ======================================================================================
# Config tuples: (name, params_dict, seq_restraint_force, dr_force, use_water_restraint)
# ======================================================================================

_CONFIGS_2FS = [
    ("eq1", EQ1_PARAMS, 10.0, 1.5, False),
    ("eq2", EQ2_2FS_PARAMS, 10.0, 1.5, False),
    ("eq3", EQ3_2FS_PARAMS, 5.0, 1.5, False),
    ("eq4", EQ4_2FS_PARAMS, 2.5, 1.5, False),
    ("eq5", EQ5_2FS_PARAMS, 0.0, None, True),  # dr_force=None means use production dr_force
]

_CONFIGS_1FS = [
    ("eq1", EQ1_PARAMS, 10.0, 1.5, False),
    ("eq2", EQ2_1FS_PARAMS, 10.0, 1.5, False),
    ("eq3", EQ3_1FS_PARAMS, 5.0, 1.5, False),
    ("eq4", EQ4_1FS_PARAMS, 2.5, 1.5, False),
    ("eq5", EQ5_1FS_PARAMS, 0.0, None, True),  # dr_force=None means use production dr_force
]


# ======================================================================================
# QresFEP equilibration
#
# Amino-acid mutations keep their own protocol, as published in Jespers et al. (2019)
# and Koenekoop et al. (2025). It differs from the ligand protocol above in three ways
# that matter: eq1-eq4 run at a fixed 1 fs regardless of the production timestep, eq5 is
# two orders of magnitude longer (the hybrid side chain needs to settle in the pocket
# before sampling), and the heat bath is coupled to the whole system rather than to
# solute and solvent separately.
# ======================================================================================

_RESFEP_SHARED = dict(
    lrf=True,
    # Coupled to one bath, Q derives a single Berendsen scaling factor from the whole
    # system's temperature and applies it to both, so it enforces only the total. The
    # spherical boundary does continuous work on the outer water shell, and removing
    # that excess globally over-cools the solute: measured on barnase, the solute sat
    # near 250 K against a solvent near 318 K, with the total on 298 K as asked. No
    # length of equilibration closes the gap, because the thermostat cannot see it.
    # The published protocol left this off and its 5 ns eq5 is plausibly compensating.
    separate_scaling=True,
)

RESFEP_EQ1_PARAMS = dict(
    steps=10000, stepsize=0.1, temperature=1, bath_coupling=0.2,
    shake_hydrogens=False, interval_output=5, **_RESFEP_SHARED,
)

RESFEP_EQ2_PARAMS = dict(
    steps=10000, stepsize=1.0, temperature=50, bath_coupling=1.0,
    shake_hydrogens=False, interval_output=5, **_RESFEP_SHARED,
)

RESFEP_EQ3_PARAMS = dict(
    steps=10000, stepsize=1.0, temperature=150, bath_coupling=1.0,
    shake_hydrogens=False, interval_output=5, **_RESFEP_SHARED,
)

RESFEP_EQ4_PARAMS = dict(
    steps=10000, stepsize=1.0, temperature=275, bath_coupling=1.0,
    shake_hydrogens=False, interval_output=5, **_RESFEP_SHARED,
)

#: eq5 is the only stage that follows the production timestep: (steps, stepsize,
#: shake_hydrogens) per timestep. The default 2 fs protocol runs 1,250,000 steps,
#: or 2.5 ns: half the 2,500,000-step equilibration used by the published runs.
_RESFEP_EQ5_BY_TIMESTEP = {
    "1fs": (500000, 1.0, False),
    "2fs": (1250000, 2.0, True),
}


def get_resfep_equilibration_configs(
    timestep: Literal["1fs", "2fs"],
    shell_radius: int,
    eq5_steps: int | None = None,
    separate_scaling: bool = True,
) -> list[EquilibrationConfig]:
    """Return the QresFEP equilibration configurations.

    Args:
        timestep: Production timestep ("1fs" or "2fs"); only eq5 follows it.
        shell_radius: Inner radius of the restrained boundary shell.
        eq5_steps: Length of the final equilibration, overriding the 2 fs protocol's
            own 2.5 ns. Charge-changing mutations have been run with longer
            ones, and reproducing such a set means matching it.
        separate_scaling: Drive solute and solvent to the heat bath separately.
            Set false only to reproduce the published legacy protocol.

    Returns:
        List of EquilibrationConfig for eq1 through eq5.
    """
    steps, stepsize, shake_hydrogens = _RESFEP_EQ5_BY_TIMESTEP[timestep]
    if eq5_steps is not None:
        steps = int(eq5_steps)
    eq5_params = dict(
        steps=steps,
        stepsize=stepsize,
        temperature="T_VAR",
        bath_coupling=10.0,
        shake_hydrogens=shake_hydrogens,
        interval_output=5,
        **_RESFEP_SHARED,
    )

    # (name, params, sequence restraint force). The hybrid residue is held by distance
    # restraints throughout, so the whole-solute sequence restraint is released at eq5,
    # where only the reference tripeptide keeps an anchor.
    raw_configs = [
        ("eq1", RESFEP_EQ1_PARAMS, 10.0),
        ("eq2", RESFEP_EQ2_PARAMS, 10.0),
        ("eq3", RESFEP_EQ3_PARAMS, 5.0),
        ("eq4", RESFEP_EQ4_PARAMS, 5.0),
        ("eq5", eq5_params, 0.0),
    ]

    # Override on a copy: the module-level dictionaries are shared constants and
    # a legacy-compatible call must not affect later corrected-protocol setups.
    return [
        EquilibrationConfig(
            name=name,
            params=MDParameters(
                **{**params, "separate_scaling": separate_scaling},
                shell_radius=shell_radius,
            ),
            sequence_restraint_force=seq_force,
            distance_restraint_force=None,
            use_water_restraint=(name == "eq5"),
        )
        for name, params, seq_force in raw_configs
    ]


def get_equilibration_configs(
    timestep: Literal["1fs", "2fs"],
    shell_radius: int,
) -> list[EquilibrationConfig]:
    """Return equilibration configurations for the given timestep.

    The equilibration protocol consists of 5 stages:
    - eq1: Initial equilibration at 1K with small timestep (0.1fs)
    - eq2: Heating to 50K
    - eq3: Heating to 150K
    - eq4: Heating to 275K
    - eq5: Final equilibration at target temperature (T_VAR)

    Args:
        timestep: Simulation timestep ("1fs" or "2fs")
        shell_radius: Spherical boundary radius

    Returns:
        List of EquilibrationConfig for eq1 through eq5
    """
    raw_configs = _CONFIGS_2FS if timestep == "2fs" else _CONFIGS_1FS

    return [
        EquilibrationConfig(
            name=name,
            params=MDParameters(**params, shell_radius=shell_radius),
            sequence_restraint_force=seq_force,
            distance_restraint_force=dr_force,
            use_water_restraint=use_water,
        )
        for name, params, seq_force, dr_force, use_water in raw_configs
    ]
