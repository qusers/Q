"""Module for atom mapping utilities, to be used to set FEP restraints. Main function is `process_rings_separately`"""

from rdkit import Chem


def get_surrounding_idxs(atom, mol):
    return [a.GetIdx() for a in mol.GetAtomWithIdx(atom).GetNeighbors()]


def identify_and_enumerate_rings(mol):
    """Identify distinct ring structures and enumerate them."""
    ring_info = mol.GetRingInfo()
    rings = {}  # Maps ring index to a set of atom indices in the ring
    for idx, ring in enumerate(ring_info.AtomRings()):
        rings[idx] = set(ring)
    return rings


def map_rings_between_molecules(
    ringsA, ringsB, atom_mapping
):  # TODO: test when you don't have respective ring
    """Map rings between molecules based on atom mapping."""
    ring_mapping = {}  # Maps ring indices in A to ring indices in B
    for idxA, atomsA in ringsA.items():
        for idxB, atomsB in ringsB.items():
            if all(atom_mapping.get(a) in atomsB for a in atomsA):
                ring_mapping[idxA] = idxB
                break
    return ring_mapping


def process_rings_separately(molA, molB, atom_mapping):
    """Process each ring separately based on ring equivalency and their substituents.
    The returned dictionary is the main output for the restraint setting process. It contains
    the mapping between rings, the atoms in each ring, and the substituents for each ring
    in both molecules.

    Args:
        molA: RDKit molecule object for molecule A.
        molB: RDKit molecule object for molecule B.
        atom_mapping: Kartograf mapping between atoms in molecule A and molecule B.

    Returns:
        dict: Dictionary containing the ring mapping, the atoms in each ring, and the substituents
            for each ring in both molecules.
    Example:
        >>> {'Ring mapping': {0: 0, 1: 1, 2: 2, 3: 3},
        >>> 'Ring 0': {'ringAtomsA': {set of AtomIndexes in ring 0 of molecule A},
        >>> 'ringAtomsB': {set of AtomIndexes in ring 0 of molecule B},
        >>> 'substituentsA': {AtominRing0 for molecule A: [substituent atom indexes]
        >>>    ...},
        >>> 'substituentsB': {AtominRing0 for molecule B: [substituent atom indexes]
        >>>    ...},
        >>> 'Ring 1': ...,
        >>> }
    """
    ringsA = identify_and_enumerate_rings(molA)
    ringsB = identify_and_enumerate_rings(molB)
    rings = {}
    ring_mapping = map_rings_between_molecules(ringsA, ringsB, atom_mapping)
    rings.update({"Ring mapping": ring_mapping})

    for ring_idx, (idxA, idxB) in enumerate(ring_mapping.items()):
        atomsA = ringsA[idxA]
        atomsB = ringsB[idxB]
        substituentsA = get_substituents(molA, atomsA)
        substituentsB = get_substituents(molB, atomsB)
        rings.update(
            {
                f"Ring {ring_idx}": {
                    "ringAtomsA": atomsA,
                    "ringAtomsB": atomsB,
                    "substituentsA": substituentsA,  # here we have ringAtom : [substituent atoms]
                    "substituentsB": substituentsB,
                }
            }
        )
    return rings


def remove_from_mapping(mapping, to_remove):
    """Remove specified atoms from the mapping."""
    for atom in to_remove:
        if atom in mapping:
            del mapping[atom]


def is_atom_in_other_ring(atom, current_ring_atoms, mol):
    """Check if the atom is part of a ring that is not the current one being processed."""
    # Check all rings the atom is a part of, if any is fully outside current_ring_atoms, it's another ring
    for ring in mol.GetRingInfo().AtomRings():
        if atom.GetIdx() in ring and not all(atom_idx in current_ring_atoms for atom_idx in ring):
            return True
    return False


def walk_bonds(atom, found_atoms, visited_atoms, mol, current_ring_atoms):
    if atom.GetIdx() in visited_atoms:
        return found_atoms, visited_atoms  # Return immediately if atom has been visited

    visited_atoms.add(atom.GetIdx())  # Mark this atom as visited
    found_atoms.add(atom.GetIdx())  # Add this atom to the found set

    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in visited_atoms and neighbor.GetIdx() not in current_ring_atoms:
            # Recursively walk bonds from this neighbor
            found_atoms, visited_atoms = walk_bonds(
                neighbor, found_atoms, visited_atoms, mol, current_ring_atoms
            )

    return found_atoms, visited_atoms


def get_bond_symbol(bond):
    """Return the symbol for the bond type."""
    bond_type = bond.GetBondType()
    if bond_type == Chem.BondType.SINGLE:
        return "-"
    elif bond_type == Chem.BondType.DOUBLE:
        return "="
    elif bond_type == Chem.BondType.TRIPLE:
        return "#"
    elif bond_type == Chem.BondType.AROMATIC:
        return ":"
    return ""  # Default or unrecognized bond type


def get_substituents(mol, ring_atoms):
    substituents = {}
    visited_atoms = set(ring_atoms)  # Start with ring atoms as visited to avoid re-traversing the ring

    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)

        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited_atoms:
                # Explore substituents starting from this neighbor
                new_found, _ = walk_bonds(neighbor, set(), visited_atoms.copy(), mol, ring_atoms)
                if new_found:  # If any new atoms were found, add them to the substituents list
                    substituents[atom_idx] = list(
                        new_found  # Convert set to list if needed for downstream processing
                    )
        # after completing the for loop, add itself if no substituents were found
        if atom_idx not in substituents:
            substituents[atom_idx] = [atom_idx]

    return substituents
