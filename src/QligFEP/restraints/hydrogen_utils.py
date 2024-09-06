"""Utility functions for setting hydrogen indexes at the end of the RDKit.Mol object."""

from rdkit.Chem import rdmolops


def are_hydrogens_at_end(mol):
    """Created to certify that heavy atom indexes will be conserved upon removal of hydrogens."""
    atoms = mol.GetAtoms()
    non_h_found = False
    for atom in reversed(atoms):
        if atom.GetAtomicNum() == 1:  # If hydrogen
            if non_h_found:
                return False
        else:
            non_h_found = True
    return True


def reindex_hydrogens_to_end(mol):
    """Modify the indexes of the atoms in the molecule to have H's at the end.

    Args:
        mol: rdkit molecule object

    Returns:
        rdkit molecule object with hydrogens at the end of the atom list.
    """
    non_h_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    h_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1]
    new_order = non_h_atoms + h_atoms
    mol_reordered = rdmolops.RenumberAtoms(mol, new_order)
    return mol_reordered
