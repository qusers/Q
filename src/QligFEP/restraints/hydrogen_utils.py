"""Utility functions for setting hydrogen indexes at the end of the RDKit.Mol object."""

from rdkit.Chem import rdmolops


def are_hydrogens_at_end(mol):
    """Created to certify that heavy atom indexes will be conserved upon removal of hydrogens."""
    atoms = mol.GetAtoms()
    h_block_ended = False
    for atom in reversed(atoms):
        # logger.trace(f"Atom {atom.GetIdx()} atomicNumber {atom.GetAtomicNum()}")
        if atom.GetAtomicNum() == 1:
            if h_block_ended:
                return False
        else:
            h_block_ended = True
    return True  # All checks passed


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
    for prop_name in mol.GetPropNames():  # assert properties are conserved
        prop_value = mol.GetProp(prop_name)
        mol_reordered.SetProp(prop_name, prop_value)
    return mol_reordered
