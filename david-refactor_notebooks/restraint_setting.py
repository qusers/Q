from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from copy import deepcopy
from rdkit.Chem import RWMol
from kartograf.atom_aligner import align_mol_shape
from kartograf.atom_mapping_scorer import (
    # MappingRMSDScorer,
    # MappingShapeOverlapScorer,
    MappingVolumeRatioScorer,
)
from kartograf import KartografAtomMapper, SmallMoleculeComponent
from QligFEP.chemIO import MoleculeIO

root_path = Path(__file__).parent

molio = MoleculeIO(root_path / "BACE_ligands.sdf")
molio.setup_mols_and_names()

# there's an error currently because the atoms aren't correctly aligned (1A between rings)
mol1 = molio.molecules[molio.lig_names.index("CAT-17c")]
mol2 = molio.molecules[molio.lig_names.index("CAT-17d")]

img = Draw.MolsToGridImage([mol1.to_rdkit(), mol2.to_rdkit()])

rdmols = [mol1.to_rdkit(), mol2.to_rdkit()]
molA, molB = [SmallMoleculeComponent.from_rdkit(m) for m in rdmols]

# Align the mols first - this might not needed, depends on input.
a_molB = align_mol_shape(molB, ref_mol=molA)

# Build Kartograf Atom Mapper
mapper = KartografAtomMapper(atom_map_hydrogens=True)

# Get Mapping
kartograf_mapping = next(mapper.suggest_mappings(molA, a_molB))

# Score Mapping
rmsd_scorer = MappingVolumeRatioScorer()
score = rmsd_scorer(mapping=kartograf_mapping)
print(f"RMSD Score: {score}")

kartograf_mapping


atom_mapping = deepcopy(kartograf_mapping.to_dict()["componentA_to_componentB"])

def get_surrounding_idxs(atom, mol):
    return [a.GetIdx() for a in mol.GetAtomWithIdx(atom).GetNeighbors()]

def identify_and_enumerate_rings(mol):
    """Identify distinct ring structures and enumerate them."""
    ring_info = mol.GetRingInfo()
    rings = {}  # Maps ring index to a set of atom indices in the ring
    for idx, ring in enumerate(ring_info.AtomRings()):
        rings[idx] = set(ring)
    return rings


def map_rings_between_molecules(ringsA, ringsB, atom_mapping): # TODO: test when you don't have respective ring
    """Map rings between molecules based on atom mapping."""
    ring_mapping = {}  # Maps ring indices in A to ring indices in B
    for idxA, atomsA in ringsA.items():
        for idxB, atomsB in ringsB.items():
            if all(atom_mapping.get(a) in atomsB for a in atomsA):
                ring_mapping[idxA] = idxB
                break
    return ring_mapping


def process_rings_separately(molA, molB, atom_mapping):
    """Process each ring separately based on ring equivalency."""
    ringsA = identify_and_enumerate_rings(molA)
    ringsB = identify_and_enumerate_rings(molB)
    rings = {}
    ring_mapping = map_rings_between_molecules(ringsA, ringsB, atom_mapping)
    rings.update({"Ring mapping": ring_mapping})

    print("Ring mapping: ", ring_mapping)
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
        if atom.GetIdx() in ring and not all(
            atom_idx in current_ring_atoms for atom_idx in ring
        ):
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
            found_atoms, visited_atoms = walk_bonds(neighbor, found_atoms, visited_atoms, mol, current_ring_atoms)
    
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
                    substituents[atom_idx] = list(new_found)  # Convert set to list if needed for downstream processing

    return substituents

def extract_sub_molecule(atom_indices, input_molecule):
    """Helper function to extract a sub molecule from an input molecule based on a set of
    atom indices.

    Args:
        atom_indices: atom indices to use for the extracted submolecule.
        input_molecule: molecule from which to substract the submolecule.

    Returns:
        rdkit.Chem.Mol: the extracted submolecule.
    """    
    # Validate that the atom indices list is not empty
    if not atom_indices:
        raise ValueError("No atom indices provided.")

    # Handle the case of a single atom
    if len(atom_indices) == 1:
        return Chem.MolFragmentToSmiles(
            input_molecule, atomsToUse=atom_indices, isomericSmiles=True
        )

    # Create the sub-molecule
    editable_molecule = RWMol(input_molecule)
    atoms_to_remove = set(range(editable_molecule.GetNumAtoms())) - set(atom_indices)
    for atom_idx in sorted(atoms_to_remove, reverse=True):
        editable_molecule.RemoveAtom(atom_idx)
    return editable_molecule.GetMol()

def submolecule_to_img(sub):
    Chem.rdDepictor.Compute2DCoords(sub)
    return Draw.MolToImage(sub, size=(300, 300))

ring_restraint_struct = process_rings_separately(
    molA.to_rdkit(), molB.to_rdkit(), atom_mapping
)
# print("Updated mapping:", atom_mapping)


for ringA, ringB in ring_restraint_struct["Ring mapping"].items():
    subsA = ring_restraint_struct[f"Ring {ringA}"]["substituentsA"]
    subsB = ring_restraint_struct[f"Ring {ringA}"]["substituentsB"]
    print(f"For Ring {ringA}:")
    ringA = extract_sub_molecule(list(subsA.keys()), molA.to_rdkit())
    ringB = extract_sub_molecule(list(subsB.keys()), molB.to_rdkit())
    img = submolecule_to_img([ringA, ringB])
    # display(img) # for the jupyter notebook
    print("for the substituents:")
    for key in subsA.keys():
        keyB = atom_mapping[key]
        print(f"from : Atom {key} -> Atom {keyB}")
        atomsA = [key] + subsA[key]
        atomsB = [keyB] + subsB[keyB]
        subA = extract_sub_molecule(atomsA, molA.to_rdkit())
        subB = extract_sub_molecule(atomsB, molB.to_rdkit())
        img = submolecule_to_img([subA, subB])
        # display(img)