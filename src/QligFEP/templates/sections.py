"""Formatters for restraint sections in MD input files."""


def format_distance_restraints(
    atom_pairs: list[tuple[int, int]],
    force: float,
) -> str:
    """Format the [distance_restraints] section content.

    Args:
        atom_pairs: List of (ligand1_atom, ligand2_atom) index pairs
        force: Force constant for the restraints

    Returns:
        Formatted string for distance_restraints section
    """
    lines = []
    for atom1, atom2 in atom_pairs:
        lines.append(f"{atom1} {atom2} 0.0 0.1 {force:.1f} 0")
    return "\n".join(lines)


def format_sequence_restraint(
    atom_start: int,
    atom_end: int,
    force: float,
) -> str:
    """Format a sequence restraint line for equilibration files (eq1-eq4).

    Args:
        atom_start: First atom index (ATOM_START_LIG1)
        atom_end: Last atom index (ATOM_END)
        force: Force constant

    Returns:
        Formatted sequence restraint line
    """
    # Format matches original template: ATOM_START_LIG1 ATOM_END 10.0 0  0
    # where atoms are left-padded to 6 chars, force is right-padded to 4 chars
    return f"{atom_start:<6} {atom_end:<6} {force:>4.1f} 0  0"


def format_water_restraint(
    atom_start: int,
    atom_end: int,
    force: float = 1.0,
) -> str:
    """Format a water system restraint line (eq5 and production MD).

    Used only in water systems to restrain ligands in place.
    For protein systems, this should return an empty string.

    Args:
        atom_start: First ligand atom index
        atom_end: Last ligand atom index (both ligands)
        force: Force constant (default 1.0)

    Returns:
        Formatted water restraint line
    """
    return f"{atom_start:<7}{atom_end:<7} {force:.1f} 0 1   "
