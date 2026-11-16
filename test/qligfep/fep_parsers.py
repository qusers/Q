"""Parser functions for FEP setup output files.

Extracts structural data from Q/FEP output files for validation and testing.
"""

import re


def parse_distance_restraints(md_content: str) -> list[tuple[int, int, float]]:
    """Extract (lig1_atom, lig2_atom, force) from [distance_restraints] section."""
    restraints = []
    in_section = False

    for line in md_content.splitlines():
        line = line.strip()

        if line.startswith("[distance_restraints]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            parts = line.split()
            if len(parts) >= 5:
                lig1_atom = int(parts[0])
                lig2_atom = int(parts[1])
                force = float(parts[4])
                restraints.append((lig1_atom, lig2_atom, force))

    return restraints


def parse_sequence_restraints(md_content: str) -> bool:
    """Check if [sequence_restraints] section has content."""
    in_section = False

    for line in md_content.splitlines():
        line = line.strip()

        if line.startswith("[sequence_restraints]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            return True

    return False


def parse_topology_header(top_content: str) -> tuple[int, int]:
    """Extract (total_atoms, solute_atoms) from topology header.

    The header line format is:
    '    7603    4753 = Total no. of atoms, no. of solute atoms.'
    """
    for line in top_content.splitlines():
        match = re.match(r"\s*(\d+)\s+(\d+)\s*=\s*Total no\. of atoms", line)
        if match:
            return int(match.group(1)), int(match.group(2))

    raise ValueError("Could not parse topology header")


def parse_fep_atom_count(fep_content: str) -> int:
    """Count entries in [atoms] section of FEP file."""
    count = 0
    in_section = False

    for line in fep_content.splitlines():
        line = line.strip()

        if line.startswith("[atoms]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            count += 1

    return count


def get_lig_start_indices(pdb_content: str) -> tuple[int, int]:
    """Find first atom index for LIG and LID residues in combined PDB."""
    lig_start = None
    lid_start = None

    for line in pdb_content.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        atom_num = int(line[6:11].strip())
        resname = line[17:20].strip()

        if resname == "LIG" and lig_start is None:
            lig_start = atom_num
        elif resname == "LID" and lid_start is None:
            lid_start = atom_num

        if lig_start is not None and lid_start is not None:
            break

    if lig_start is None or lid_start is None:
        raise ValueError("Could not find LIG and LID residues in PDB")

    return lig_start, lid_start


def normalize_restraints(
    restraints: list[tuple[int, int, float]], lig1_start: int, lig2_start: int
) -> tuple[list[int], list[int]]:
    """Convert absolute indices to relative (0-based from each ligand start).

    Returns (lig1_relative_indices, lig2_relative_indices).
    """
    lig1_indices = [r[0] - lig1_start for r in restraints]
    lig2_indices = [r[1] - lig2_start for r in restraints]
    return lig1_indices, lig2_indices
