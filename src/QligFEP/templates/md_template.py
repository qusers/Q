"""Base MD template and render function for Q input files."""

from dataclasses import dataclass


@dataclass
class MDParameters:
    """Parameters for a Q molecular-dynamics input (.inp) file.

    Each field maps to a key the Q engine (src/q6/md.f90) reads from a named
    [section]. The sections are:

      - [MD]       core dynamics: integration, thermostat, SHAKE, minimiser.
      - [cut-offs] hard spherical non-bonded cutoffs (A); q_atom ~ infinite.
      - [sphere]   solute boundary-shell restraint (SCAAS droplet).
      - [solvent]  water boundary restraints + opt-in per-state Born correction.
      - [intervals]/[files] output cadence and filenames.

    The temperature field accepts "T_VAR" for runtime substitution by the job
    submission script.
    """

    # Simulation core
    steps: int
    stepsize: float  # fs; engine converts to internal units (dt = 0.020462 * stepsize)
    temperature: int | str  # K target; or "T_VAR" placeholder
    bath_coupling: float  # fs; Berendsen bath relaxation time (tau_T), must be >= stepsize

    # Shake settings
    shake_solvent: bool = True  # Should be always on - we use rigid water models
    shake_hydrogens: bool = True
    shake_solute: bool = False

    # Other MD settings
    lrf: bool = True  # Local Reaction Field Taylor expansion beyond the cutoff
    separate_scaling: bool = True  # couple solute and solvent to the heat bath separately

    # Cut-offs in Angstrom (usually constant)
    cutoff_solute_solvent: int = 10
    cutoff_solute_solute: int = 10
    cutoff_solvent_solvent: int = 10
    cutoff_q_atom: int = 99  # ~infinite: Q-atoms see all environment (FEP consistency)
    cutoff_lrf: int = 99  # outer LRF radius; engine requires >= the cutoffs above

    # Sphere
    shell_force: float = 10.0  # kcal/mol/A^2; harmonic restraint of boundary-shell solute atoms
    shell_radius: int = 25  # A; inner radius of the restrained shell (rexcl_i)

    # Solvent
    radial_force: float = 60.0  # kcal/mol/A^2; radial restraint keeping surface waters in
    polarization: bool = True  # SCAAS surface-water dipole-orientation restraint
    polarization_force: float = 20.0  # kcal/mol/rad^2; force constant for the orientation restraint

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
        f"energy                    {params.interval_energy}\n" if params.interval_energy is not None else ""
    )

    restart_file_name = f"restart                   {restart_file}\n" if restart_file else ""
    energy_file_name = f"energy                    {energy_file}\n" if energy_file else ""

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
polarization              {onoff(params.polarization)}
polarization_force        {params.polarization_force}

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

[trajectory_atoms]
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
