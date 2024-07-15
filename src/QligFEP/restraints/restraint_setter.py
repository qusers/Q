from copy import deepcopy
from itertools import zip_longest
from pathlib import Path

import numpy as np
from kartograf import KartografAtomMapper, SmallMoleculeComponent
from kartograf.atom_aligner import align_mol_shape
from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer
from rdkit import Chem

from ..logger import logger
from .atom_mapping import AtomMapperHelper


class RestraintSetter:
    """
    TODO: currently, this class correctly sets the restraints for the ring structures,
    but it takes into account:
        1) the ring structure itself
        2) the substitutions
    In the future, it should also account for: # see FIXME below...
        3) the fact that the ring, despite being the same, has different decorations and
        should therefore be left unrestrained
    """

    def __init__(self, molA: str, molB: str) -> None:
        self.ligname1 = Path(molA).stem
        self.ligname2 = Path(molB).stem
        self._load_molecules(molA, molB)
        self._align_and_map_molecules()

    def _path_to_mol(self, _path):
        if Path(_path).suffix == ".pdb":
            logger.error(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost. PDBs can be used along with a valid smiles string with RDKit using "
                "the constructor from `openff.toolkit` Molecule.from_pdb_and_smiles(file_path, smiles) "
                "but this is not yet implemented in this class."
            )
            raise ValueError("PDBs are not supported yet.")
            return Chem.MolFromPDBFile(_path)
        elif Path(_path).suffix == ".sdf":
            return Chem.SDMolSupplier(str(_path), removeHs=False)[0]
        return None

    def _load_molecules(self, molA, molB):
        rdmolA = self._path_to_mol(molA)
        rdmolB = self._path_to_mol(molB)
        if any([rdmolA is None, rdmolB is None]):
            raise ValueError("Molecule could not be loaded!")
        self.molA = SmallMoleculeComponent.from_rdkit(rdmolA)
        self.molB = SmallMoleculeComponent.from_rdkit(rdmolB)

    def _align_and_map_molecules(self):
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        a_molB = align_mol_shape(self.molB, ref_mol=self.molA)
        self.kartograf_mapping = next(mapper.suggest_mappings(self.molA, a_molB))
        # Score Mapping
        rmsd_scorer = MappingVolumeRatioScorer()
        score = rmsd_scorer(mapping=self.kartograf_mapping)
        logger.info(f"Volume ratio score between ligands: {score}")
        self.atom_mapping = deepcopy(self.kartograf_mapping.to_dict()["componentA_to_componentB"])

    @staticmethod
    def get_ring_and_ring_immediates(mol: Chem.Mol, connect_struct: dict):
        """Function to get the ring atoms and their immediate neighbors as an array of atomic
        numbers. The list is ordered as [ring_atom1, bound_atom1, ring_atom2, bound_atom2, ...].

        Args:
            mol: rdkit molecule object
            connect_struct: dictionary containing the ring atoms (keys) and their bound atoms (values)

        Returns:
            np.array: array of atomic numbers for the ring atoms and their immediate neighbors.
        """
        ringAtom_subsAtom = []  # list containing the ring atoms & bound atoms
        for ringAtom, subsAtom in connect_struct.items():
            ring_atomicnum = mol.GetAtomWithIdx(ringAtom).GetAtomicNum()
            if subsAtom[0] != ringAtom:
                subs_atomicnum = mol.GetAtomWithIdx(subsAtom[0]).GetAtomicNum()
                ringAtom_subsAtom.extend([ring_atomicnum, subs_atomicnum])
            else:
                # no 0 atomic number so we use for missing atoms
                ringAtom_subsAtom.extend([ring_atomicnum, 0])
        return np.array(ringAtom_subsAtom)

    @staticmethod
    def are_atoms_equivalent(atom_a, atom_b):
        return atom_a.GetAtomicNum() == atom_b.GetAtomicNum()

    def are_substituents_equivalent(self, subsA, subsB, atom_mapping, mol_a, mol_b):
        logger.trace(f"Substituents A: {subsA}")
        logger.trace(f"Substituents B: {subsB}")
        is_same = set()
        for ring_atomA, subs_atomsA in subsA.items():
            all_atomsA = [ring_atomA] + subs_atomsA
            ring_atomB = atom_mapping.get(ring_atomA)
            if ring_atomB is None or ring_atomB not in subsB:
                continue
            all_atomsB = [ring_atomB] + subsB[ring_atomB]
            for atomA, atomB in zip_longest(all_atomsA, all_atomsB):
                if atomA is None or atomB is None:
                    break
                atomA_obj = mol_a.GetAtomWithIdx(atomA)
                atomB_obj = mol_b.GetAtomWithIdx(atomB)

                # The iteration should stop when it finds a ring structure
                if any([atomA != ring_atomA, atomB != ring_atomB]) and any(
                    [atomA in self.atom_mapper.ringIdxsA, atomB in self.atom_mapper.ringIdxsB]
                ):
                    if self.are_atoms_equivalent(atomA_obj, atomB_obj):
                        is_same.add(atomA)
                    break

                if self.are_atoms_equivalent(atomA_obj, atomB_obj):
                    is_same.add(atomA)  # Add original ring atom if all checks pass
                else:
                    break
        return is_same

    def is_ring_equivalent(self, ring_data, mol_a, mol_b, atom_mapping) -> bool:
        # Enhanced ring equivalence check
        ring_atoms_a_indices = ring_data["ringAtomsA"]
        ring_atoms_b_indices = [
            atom_mapping.get(a) for a in ring_atoms_a_indices if atom_mapping.get(a) is not None
        ]

        # Check if mapped indices fully match the ringAtomsB set
        if set(ring_atoms_b_indices) != ring_data["ringAtomsB"]:
            return False

        # Further check the atomic number and connectivity
        for a, b in zip(ring_atoms_a_indices, ring_atoms_b_indices):
            atom_a = mol_a.GetAtomWithIdx(a)
            atom_b = mol_b.GetAtomWithIdx(b)
            if not self.are_atoms_equivalent(atom_a, atom_b):
                return False

        return True

    def compare_molecule_rings(
        self, data: dict, atom_mapping: dict, mol_a: Chem.Mol, mol_b: Chem.Mol, strict: bool = True
    ) -> dict:
        """Compares the ring structures based on the data output from `process_rings_separately`
        and returns a dictionary with the final atoms to be restrained.

        Args:
            data: dictionary output from `process_rings_separately`.
            atom_mapping: atom mapping dictionary output from Kartograf.
            mol_a: molecule A in RDKit format.
            mol_b: molecule B in RDKit format.
            strict: if true, it will purge rings with differing immediate-surrounds
                (e.g. a methil decoration v.s. a Cl atom). Defaults to True.

        Returns:
            a dictionary with the final atoms to be restrained.
        """
        matching_atoms = {}
        purge_list = []  # make this for the ring atoms that aren't equal including neighbors
        data = deepcopy(data)
        data.pop("Mapped Rings")
        connections = {k: v for k, v in data.items() if "->" in k}

        for ring_key, ring_data in data.items():
            if "->" in ring_key:
                continue
            if "Ring" not in ring_key:
                continue
            if self.is_ring_equivalent(ring_data, mol_a, mol_b, atom_mapping):
                matching_atoms[ring_key] = set(ring_data["ringAtomsA"])
                if strict:
                    molA_ring_and_subs = self.get_ring_and_ring_immediates(mol_a, ring_data["substituentsA"])
                    molB_ring_and_subs = self.get_ring_and_ring_immediates(mol_b, ring_data["substituentsB"])
                    is_same = np.equal(molA_ring_and_subs, molB_ring_and_subs).all()
                    if not is_same:
                        ring_connections = {k: v for k, v in connections.items() if k.startswith(ring_key)}
                        conserve_atoms = []  # we want to conserve the atoms that connect to other rings
                        for conserve in ring_connections.values():
                            conserve_atoms.extend(conserve)
                        to_purge = np.unique(
                            [item for k, vals in ring_data["substituentsA"].items() for item in [k] + vals]
                        )
                        # purge everything else that is not connected to another ring
                        purge_list.extend(np.setdiff1d(to_purge, conserve_atoms).tolist())
            else:
                purge_list.extend(ring_data["ringAtomsA"])

        # Process and compare substituents only for equivalent rings
        for ring_key in matching_atoms:
            ring_data = data[ring_key]
            subsA = ring_data["substituentsA"]
            subsB = ring_data["substituentsB"]

            matched_subs = self.are_substituents_equivalent(subsA, subsB, atom_mapping, mol_a, mol_b)
            matching_atoms[ring_key].update(matched_subs)

        logger.trace(f"PURGE: {purge_list}")
        all_atomA_idxs = set(
            idx for ring_atoms in matching_atoms.values() for idx in ring_atoms if idx not in purge_list
        )
        restraints = {k: atom_mapping[k] for k in all_atomA_idxs if k in atom_mapping}
        return restraints

    def set_restraints(self) -> dict:
        self.atom_mapper = AtomMapperHelper()
        ringStruc_compareDict = self.atom_mapper.process_rings_separately(
            Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit()), self.atom_mapping
        )
        restraints = self.compare_molecule_rings(
            ringStruc_compareDict, self.atom_mapping, self.molA.to_rdkit(), self.molB.to_rdkit()
        )
        logger.debug(f"Atoms to restrain: {restraints}")
        self.restraints = restraints
        return restraints
