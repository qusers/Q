"""Production MD configurations for FEP simulations."""

from dataclasses import dataclass
from typing import Literal

from .md_template import MDParameters


@dataclass
class ProductionConfig:
    """Configuration for production MD runs.

    Attributes:
        params: MD simulation parameters
        distance_restraint_force: Force constant for distance restraints
    """

    params: MDParameters
    distance_restraint_force: float


# ============================================================
# Explicit parameter dictionaries - can be imported and inspected
# ============================================================

PRODUCTION_2FS_PARAMS = dict(
    steps=5000,
    stepsize=2.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    constrain_hydrogens=True,
    interval_output=25,
    interval_energy=10,
    interval_trajectory=100,
    interval_non_bond=25,
)

PRODUCTION_1FS_PARAMS = dict(
    steps=10000,
    stepsize=1.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    constrain_hydrogens=False,
    interval_output=25,
    interval_energy=10,
    interval_trajectory=100,
    interval_non_bond=25,
)


def get_production_config(
    timestep: Literal["1fs", "2fs"],
    shell_radius: int,
    distance_restraint_force: float = 0.5,
    **overrides,
) -> ProductionConfig:
    """Return production MD configuration for the given timestep.

    Production MD files include energy output and use WATER_RESTRAINT
    for sequence restraints (empty for protein systems).

    Args:
        timestep: Simulation timestep ("1fs" or "2fs")
        shell_radius: Spherical boundary radius
        distance_restraint_force: Force constant for distance restraints
        **overrides: Override any base parameter by name (e.g., steps=10000,
            bath_coupling=5.0, interval_energy=5)

    Returns:
        ProductionConfig for md_*.inp files

    Examples:
        # Standard production run
        config = get_production_config("2fs", shell_radius=25)

        # Extended sampling for difficult targets
        config = get_production_config("2fs", shell_radius=25, steps=10000)

        # Fine-grained energy output for analysis
        config = get_production_config("2fs", shell_radius=25, interval_energy=5)
    """
    base_params = PRODUCTION_2FS_PARAMS if timestep == "2fs" else PRODUCTION_1FS_PARAMS
    merged_params = {**base_params, **overrides, "shell_radius": shell_radius}

    return ProductionConfig(
        params=MDParameters(**merged_params),
        distance_restraint_force=distance_restraint_force,
    )
