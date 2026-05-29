"""Base MD template and render function for Q input files."""

from dataclasses import dataclass


@dataclass
class MDParameters:
    """Parameters for MD input files.

    This dataclass captures all the configurable parameters for Q MD simulations.
    The temperature field accepts "T_VAR" for runtime substitution by the job
    submission script.
    """

    # Simulation core
    steps: int
    stepsize: float
    temperature: int | str  # Can be int or "T_VAR" placeholder
    bath_coupling: float

    # Shake settings
    shake_solvent: bool = True  # Should be always on - we use rigid water models
    shake_hydrogens: bool = True
    shake_solute: bool = False

    # Other MD settings
    lrf: bool = True
    separate_scaling: bool = True

    # Cut-offs (usually constant)
    cutoff_solute_solvent: int = 10
    cutoff_solute_solute: int = 10
    cutoff_solvent_solvent: int = 10
    cutoff_q_atom: int = 99
    cutoff_lrf: int = 99

    # Sphere
    shell_force: float = 10.0
    shell_radius: int = 25

    # Solvent
    radial_force: float = 60.0
    polarisation: bool = True
    polarisation_force: float = 20.0

    # Intervals
    interval_output: int = 5
    interval_trajectory: int = 100
    interval_non_bond: int = 25
    interval_energy: int | None = None  # Only for production MD

    # Files
    topology: str = "dualtop.top"
    fep_file: str = "FEP_VAR"

    # Energy minimization (steepest descent before MD)
    minimize: bool = False
    max_minimize_steps: int = 1000
    minimize_tolerance: float = 0.1  # kcal/mol/Å
    minimize_step_size: float = 0.001  # Å


def onoff(val: bool) -> str:
    """Convert boolean to 'on'/'off' string."""
    return "on" if val else "off"


def render_md_input(
    params: MDParameters,
    lambda1: str,
    lambda2: str,
    trajectory_file: str,
    final_file: str,
    distance_restraints: str = "",
    sequence_restraints: str = "",
    wall_restraints: str = "",
    restart_file: str | None = None,
    energy_file: str | None = None,
    correction_file: str | None = None,
    correction_exclude_last: int = 0,
    is_eq1: bool = False,
) -> str:
    """Render an MD input file from parameters.

    Args:
        params: MD simulation parameters
        lambda1: First lambda value (e.g., "0.500" or "FLOAT_LAMBDA1")
        lambda2: Second lambda value (e.g., "0.500" or "FLOAT_LAMBDA2")
        trajectory_file: Trajectory output filename (e.g., "eq1.dcd")
        final_file: Final restart filename (e.g., "eq1.re")
        distance_restraints: Pre-formatted distance restraints section
        sequence_restraints: Pre-formatted sequence restraints section
        wall_restraints: Pre-formatted wall restraints section (for counter-ions)
        restart_file: Restart input filename (None for eq1)
        energy_file: Energy output filename (None for equilibration)
        correction_file: Charge-correction observable log filename (None to disable
            the optional [correction] section); sampled at the energy interval.
        correction_exclude_last: Number of trailing Q-atoms — the co-alchemical
            counter-water, appended last in the FEP [atoms] list — to drop from
            the phi_cog centroid (they sit far from the ligand and bias it). 0
            (default) disables it, preserving byte-identical output.
        is_eq1: True for eq1.inp which has random_seed and initial_temperature

    Returns:
        Complete .inp file content as string
    """
    # Pre-build optional sections (empty string = omitted from output).
    # The template literal below is written flush-left so the rendered output
    # carries no leading whitespace regardless of how multi-line restraint
    # sections are passed in.
    equilibration_start = (
        "random_seed               SEED_VAR\ninitial_temperature       1\n" if is_eq1 else ""
    )

    minimization_settings = (
        (
            f"minimize                  on\n"
            f"max_minimize_steps        {params.max_minimize_steps}\n"
            f"minimize_tolerance        {params.minimize_tolerance}\n"
            f"minimize_step_size        {params.minimize_step_size}\n"
        )
        if params.minimize
        else ""
    )

    energy_interval = (
        f"energy                    {params.interval_energy}\n"
        if params.interval_energy is not None
        else ""
    )

    restart_file_name = f"restart                   {restart_file}\n" if restart_file else ""
    energy_file_name = f"energy                    {energy_file}\n" if energy_file else ""

    # Optional [correction] block; "\n" alone reproduces the blank line that
    # otherwise separates [files] from [trajectory_atoms], so output is
    # byte-identical when correction logging is disabled.
    if correction_file is not None:
        correction_interval = (
            params.interval_energy if params.interval_energy is not None else params.interval_output
        )
        exclude_line = (
            f"exclude_last              {correction_exclude_last}\n"
            if correction_exclude_last > 0
            else ""
        )
        correction_block = (
            f"\n[correction]\n"
            f"interval                  {correction_interval}\n"
            f"file                      {correction_file}\n"
            f"kernel                    1\n"
            f"{exclude_line}"
            f"\n"
        )
    else:
        correction_block = "\n"

    return f"""\
[MD]
steps                     {params.steps}
stepsize                  {params.stepsize}
temperature               {params.temperature}
bath_coupling             {params.bath_coupling}
{equilibration_start}\
shake_solvent             {onoff(params.shake_solvent)}
shake_hydrogens           {onoff(params.shake_hydrogens)}
shake_solute              {onoff(params.shake_solute)}
lrf                       {onoff(params.lrf)}
separate_scaling          {onoff(params.separate_scaling)}
{minimization_settings}\

[cut-offs]
solute_solvent            {params.cutoff_solute_solvent}
solute_solute             {params.cutoff_solute_solute}
solvent_solvent           {params.cutoff_solvent_solvent}
q_atom                    {params.cutoff_q_atom}
lrf                       {params.cutoff_lrf}

[sphere]
shell_force               {params.shell_force}
shell_radius              {params.shell_radius}

[solvent]
radial_force              {params.radial_force}
polarisation              {onoff(params.polarisation)}
polarisation_force        {params.polarisation_force}

[intervals]
output                    {params.interval_output}
{energy_interval}\
trajectory                {params.interval_trajectory}
non_bond                  {params.interval_non_bond}

[files]
topology                  {params.topology}
trajectory                {trajectory_file}
{restart_file_name}\
{energy_file_name}\
final                     {final_file}
fep                       {params.fep_file}
{correction_block}[trajectory_atoms]
not excluded

[lambdas]
{lambda1} {lambda2}

[sequence_restraints]
{sequence_restraints}

[distance_restraints]
{distance_restraints}

[wall_restraints]
{wall_restraints}
"""
