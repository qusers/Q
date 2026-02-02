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


def _bool_to_onoff(val: bool) -> str:
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
    restart_file: str | None = None,
    energy_file: str | None = None,
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
        restart_file: Restart input filename (None for eq1)
        energy_file: Energy output filename (None for equilibration)
        is_eq1: True for eq1.inp which has random_seed and initial_temperature

    Returns:
        Complete .inp file content as string
    """
    lines = []

    # [MD] section
    lines.append("[MD]")
    lines.append(f"steps                     {params.steps}")
    lines.append(f"stepsize                  {params.stepsize}")
    lines.append(f"temperature               {params.temperature}")
    lines.append(f"bath_coupling             {params.bath_coupling}")

    # eq1 has random_seed and initial_temperature (set by job submission script)
    if is_eq1:
        lines.append("random_seed               SEED_VAR")
        lines.append("initial_temperature       1")

    lines.append(f"shake_solvent             {_bool_to_onoff(params.shake_solvent)}")
    lines.append(f"shake_hydrogens           {_bool_to_onoff(params.shake_hydrogens)}")
    lines.append(f"shake_solute              {_bool_to_onoff(params.shake_solute)}")
    lines.append(f"lrf                       {_bool_to_onoff(params.lrf)}")
    lines.append(f"separate_scaling          {_bool_to_onoff(params.separate_scaling)}")
    lines.append("")

    # [cut-offs] section
    lines.append("[cut-offs]")
    lines.append(f"solute_solvent            {params.cutoff_solute_solvent}")
    lines.append(f"solute_solute             {params.cutoff_solute_solute}")
    lines.append(f"solvent_solvent           {params.cutoff_solvent_solvent}")
    lines.append(f"q_atom                    {params.cutoff_q_atom}")
    lines.append(f"lrf                       {params.cutoff_lrf}")
    lines.append("")

    # [sphere] section
    lines.append("[sphere]")
    lines.append(f"shell_force               {params.shell_force}")
    lines.append(f"shell_radius              {params.shell_radius}")
    lines.append("")

    # [solvent] section
    lines.append("[solvent]")
    lines.append(f"radial_force              {params.radial_force}")
    lines.append(f"polarisation              {_bool_to_onoff(params.polarisation)}")
    lines.append(f"polarisation_force        {params.polarisation_force}")
    lines.append("")

    # [intervals] section
    lines.append("[intervals]")
    lines.append(f"output                    {params.interval_output}")
    if params.interval_energy is not None:
        lines.append(f"energy                    {params.interval_energy}")
    lines.append(f"trajectory                {params.interval_trajectory}")
    lines.append(f"non_bond                  {params.interval_non_bond}")
    lines.append("")

    # [files] section
    lines.append("[files]")
    lines.append(f"topology                  {params.topology}")
    lines.append(f"trajectory                {trajectory_file}")
    if restart_file:
        lines.append(f"restart                   {restart_file}")
    if energy_file:
        lines.append(f"energy                    {energy_file}")
    lines.append(f"final                     {final_file}")
    lines.append(f"fep                       {params.fep_file}")
    lines.append("")

    # [trajectory_atoms] section
    lines.append("[trajectory_atoms]")
    lines.append("not excluded")
    lines.append("")

    # [lambdas] section
    lines.append("[lambdas]")
    lines.append(f"{lambda1} {lambda2}")
    lines.append("")

    # [sequence_restraints] section
    lines.append("[sequence_restraints]")
    lines.append(sequence_restraints)
    lines.append("")

    # [distance_restraints] section
    lines.append("[distance_restraints]")
    lines.append(distance_restraints)
    lines.append("")

    return "\n".join(lines)
