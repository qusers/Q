"""Non-equilibrium endpoint MD configurations for fast-switching FEP.

These drive the relax_*, eq6_* and neq_* input files: MD held at an endpoint
lambda, with the neq_* switching runs additionally carrying a [lambda_scaling]
section (added by the caller).
"""

from typing import Literal

from .md_template import MDParameters

NEQ_ENDPOINT_4FS_PARAMS = dict(
    stepsize=4.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    constrain_hydrogens=True,
    constrain_solute=True,
    interval_output=5,
    interval_trajectory=100000000,
    interval_non_bond=12,
)

NEQ_ENDPOINT_2FS_PARAMS = dict(
    stepsize=2.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    constrain_hydrogens=True,
    constrain_solute=True,
    interval_output=10,
    interval_trajectory=100000000,
)


NEQ_ENDPOINT_1FS_PARAMS = dict(
    stepsize=1.0,
    temperature="T_VAR",
    bath_coupling=10.0,
    constrain_hydrogens=False,
    constrain_solute=False,
    interval_output=10,
    interval_trajectory=100000000,
)


def get_neq_endpoint_config(
    timestep: Literal["1fs", "2fs", "4fs"],
    shell_radius: int,
    steps: int,
) -> MDParameters:
    """Return the MD parameters for a non-equilibrium endpoint input file.

    Shared by the relax_*, eq6_* and neq_* files, which differ only in step
    count (and, for neq_*, an added [lambda_scaling] section handled by the
    caller). Temperature stays as the "T_VAR" placeholder the run script fills.

    Args:
        timestep: Simulation timestep ("1fs", "2fs", or HMR-only "4fs")
        shell_radius: Spherical boundary radius
        steps: Number of MD steps for this endpoint file

    Returns:
        MDParameters for a relax/eq6/neq .inp file
    """
    params_by_timestep = {
        "1fs": NEQ_ENDPOINT_1FS_PARAMS,
        "2fs": NEQ_ENDPOINT_2FS_PARAMS,
        "4fs": NEQ_ENDPOINT_4FS_PARAMS,
    }
    try:
        base_params = params_by_timestep[timestep]
    except KeyError as exc:
        raise ValueError(f"unsupported timestep: {timestep}") from exc
    return MDParameters(**base_params, steps=steps, shell_radius=shell_radius)
