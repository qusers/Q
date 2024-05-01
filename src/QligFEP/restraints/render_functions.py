"""Module with functions to extract submolecules from atom mapping & to render restraints"""

from rdkit import Chem
from rdkit.Chem import RWMol, rdDepictor


def compute2Dcoords(mol):
    """Because I can never get this import right on the first try..."""
    rdDepictor.Compute2DCoords(mol)
    return mol


def extract_sub_molecule(atom_indices, input_molecule):
    """Helper function to extract a sub molecule from an input molecule based on a set of
    atom indices.

    Args:
        atom_indices: atom indices to use for the extracted submolecule.
        input_molecule: molecule from which to substract the submolecule.

    Returns:
        rdkit.Chem.Mol: the extracted submolecule.
    """
    if not atom_indices:
        raise ValueError("No atom indices provided.")

    # Handle the case of a single atom
    if len(atom_indices) == 1:
        return Chem.MolFragmentToSmiles(input_molecule, atomsToUse=atom_indices, isomericSmiles=True)

    editable_molecule = RWMol(input_molecule)  # Create the sub-molecule
    atoms_to_remove = set(range(editable_molecule.GetNumAtoms())) - set(atom_indices)
    for atom_idx in sorted(atoms_to_remove, reverse=True):
        editable_molecule.RemoveAtom(atom_idx)
    return editable_molecule.GetMol()


def remove_hydrogens(mol):
    """Forcefully removes hydrogens from a molecule by creating a new one from scratch
    and removing any hydrogen atoms. This is useful for creating a pattern that only
    matches heavy atoms for the substructure matching.

    Args:
        mol: RDKit molecule object to remove hydrogens from.

    Returns:
        rdchem.Mol: the molecule with hydrogens removed.
    """
    emol = Chem.EditableMol(Chem.Mol())
    index_map = {}  # Map from old atom indices to new atom indices

    # Add non-hydrogen atoms to the new molecule and build the index map
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:  # Exclude hydrogens
            new_idx = emol.AddAtom(atom)
            index_map[atom.GetIdx()] = new_idx

    for bond in mol.GetBonds():  # Add bonds between non-H atoms
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() != 1 and end_atom.GetAtomicNum() != 1:
            begin_idx = index_map[begin_atom.GetIdx()]
            end_idx = index_map[end_atom.GetIdx()]
            emol.AddBond(begin_idx, end_idx, bond.GetBondType())

    new_mol = emol.GetMol()  # get final molecule
    return new_mol
