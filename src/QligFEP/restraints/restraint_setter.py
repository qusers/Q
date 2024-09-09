"""Contain the `RestraintSetter` class - map atoms to be restrained during hybrid topology FEP calculations"""

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
from .hydrogen_utils import are_hydrogens_at_end, reindex_hydrogens_to_end


class RestraintSetter:
    """
    Class for mapping atom indexes between molecules A and B for force restraints during
    hybrid topology FEP simulations. The class uses Kartograf to define an initial atom
    mapping and compares both molecules based on their ring structures and respective decorations.

    The main method for this class is `set_restraints`, which returns a dictionary with the atoms
    to be restrained.

    Example:
        >>> from QligFEP.restraints.restraint_setter import RestraintSetter
        >>> rsetter = RestraintSetter('molecule_A.sdf', 'molecule_B.sdf')
        >>> restraints = rsetter.set_restraints(strict=False, ignore_substituent_atom_type=True)

    Attributes:
        self.molA: Molecule A (arg[0]) in RDKit format
        self.molB: Molecule B (arg[1]) in RDKit format
        self.kartograf_mapping: kartograf mapping object
        self.atom_mapping: kartagraf's `componentA_to_componentB` mapping dictionary
    """

    def __init__(self, molA: str, molB: str) -> None:
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
            mol = Chem.SDMolSupplier(str(_path), removeHs=False)[0]
            if not are_hydrogens_at_end(mol):
                raise ValueError(
                    "Hydrogens are not in the end of the atom list. Please reindex them "
                    "before proceeding and assert that your .lib, .prm and .pdb files match "
                    "the atom indexes of the molecule you're setting restraints for."
                )
            return mol
        return None

    def _load_molecules(self, molA, molB):
        rdmolA = self._path_to_mol(molA)
        rdmolB = self._path_to_mol(molB)
        if any([rdmolA is None, rdmolB is None]):
            raise ValueError("Molecule could not be loaded!")
        self.molA = SmallMoleculeComponent.from_rdkit(rdmolA)
        self.molB = SmallMoleculeComponent.from_rdkit(rdmolB)

    def _align_and_map_molecules(self):
        # TODO: now we run this alignment but the aligned molecule isn't saved
        mapper = KartografAtomMapper(atom_map_hydrogens=False, map_exact_ring_matches_only=True)
        a_molB = align_mol_shape(self.molB, ref_mol=self.molA)
        self.kartograf_mapping = next(mapper.suggest_mappings(self.molA, a_molB))
        # Score Mapping
        rmsd_scorer = MappingVolumeRatioScorer()
        score = rmsd_scorer(mapping=self.kartograf_mapping)
        logger.debug(f"Volume ratio score between ligands: {score}")
        self.atom_mapping = deepcopy(self.kartograf_mapping.to_dict()["componentA_to_componentB"])

    @staticmethod
    def get_ring_and_ringbound(
        mol: Chem.Mol,
        connect_struct: dict,
        ring_compare_method: str = "element",
        ignore_surround_elements: bool = False,
    ):
        """Function to get the ring atoms and their immediate neighbors as an array of atomic
        numbers. The list is ordered as [ring_atom1, bound_atom1, ring_atom2, bound_atom2, ...].

        Args:
            mol: rdkit molecule object
            connect_struct: dictionary containing the ring atoms (keys) and their bound atoms (values)
            ring_compare_method: method to compare the ring atoms. Choices are `element`, `hybridization`,
                and `aromaticity`.
            ignore_surround_elements: if true, it will ignore the atomic number of the substituents

        Returns:
            np.array: array of atomic numbers for the ring atoms and their immediate neighbors.
        """
        ringAtom_subsAtom = []  # list containing the ring atoms & bound atoms
        for ringAtom, subsAtom in connect_struct.items():

            if ring_compare_method == "hybridization":
                ringAtomInfo = mol.GetAtomWithIdx(ringAtom).GetHybridization().name
            elif ring_compare_method == "element":
                ringAtomInfo = mol.GetAtomWithIdx(ringAtom).GetAtomicNum()
            elif ring_compare_method == "aromaticity":
                ringAtomInfo = mol.GetAtomWithIdx(ringAtom).GetIsAromatic()

            if subsAtom[0] != ringAtom:
                subs_atomicnum = mol.GetAtomWithIdx(subsAtom[0]).GetAtomicNum()
                ringAtom_subsAtom.extend(
                    [ringAtomInfo, (subs_atomicnum if not ignore_surround_elements else 1)]
                )
            else:
                # 0 is missing atom (can't do nan on array equality)
                ringAtom_subsAtom.extend([ringAtomInfo, 0])
        return np.array(ringAtom_subsAtom)

    @staticmethod
    def are_atoms_equivalent(atom_a, atom_b, compare_method: str = "element") -> bool:
        """check if two atoms are equivalent based on the `compare_method`.

        Args:
            atom_a: rdkit Atom object A.
            atom_b: rdkit Atom object B.
            compare_method: method to compare the atoms A & B. Defaults to 'element'.

        Returns:
            bool: ture if the atoms are equivalent, false otherwise.
        """
        logger.trace(
            f"Comparing atoms {atom_a.GetSymbol()}:{atom_a.GetIdx()} and {atom_b.GetSymbol()}:{atom_b.GetIdx()}"
        )
        if compare_method == "hybridization":
            result = atom_a.GetHybridization().name == atom_b.GetHybridization().name
        elif compare_method == "element":
            result = atom_a.GetAtomicNum() == atom_b.GetAtomicNum()
        elif compare_method == "aromaticity":
            result = atom_a.GetIsAromatic() == atom_b.GetIsAromatic()
        logger.trace(f"Result: {result}")
        return result

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

    def is_ring_equivalent(self, ring_data, mol_a, mol_b, atom_mapping, compare_method="element") -> bool:
        """takes as input the mapped ring data from `AtomMapperHelper.process_rings_separately`
        compares the ring atoms between molecule A and molecule B.

        Args:
            ring_data: dictionary containing the ring data as output from `AtomMapperHelper.process_rings_separately`.
            mol_a: molecule A in RDKit format.
            mol_b: molecule B in RDKit format.
            atom_mapping: dictionary with the atom index mapping between molecules A and B.
            compare_method: how to compare the atoms within the ring. Defaults to "element".

        Returns:
            bool: True if the ring atoms are equivalent, False otherwise.
        """
        ring_atoms_a_indices = ring_data["ringAtomsA"]
        ring_atoms_b_indices = [
            atom_mapping.get(a) for a in ring_atoms_a_indices if atom_mapping.get(a) is not None
        ]

        # Check if mapped indices fully match the ringAtomsB set
        if set(ring_atoms_b_indices) != ring_data["ringAtomsB"]:
            return False

        # Further check the atomic number and connectivity
        for a, b in zip(ring_atoms_a_indices, ring_atoms_b_indices):
            logger.trace(f"Comparison method: {compare_method}")
            atom_a = mol_a.GetAtomWithIdx(a)
            atom_b = mol_b.GetAtomWithIdx(b)
            if not self.are_atoms_equivalent(atom_a, atom_b, compare_method=compare_method):
                return False

        return True

    def compare_molecule_rings(
        self,
        data: dict,
        atom_mapping: dict,
        mol_a: Chem.Mol,
        mol_b: Chem.Mol,
        ring_compare_method: str = "element",
        strict_surround: bool = False,
        ignore_surround_atom_type: bool = True,
    ) -> dict:
        """Compares the ring structures based on the data output from `AtomMapperHelper.process_rings_separately`
        and returns a dictionary with the final atoms to be restrained.

        Args:
            data: dictionary output from `process_rings_separately`.
            atom_mapping: atom mapping dictionary output from Kartograf.
            mol_a: molecule A in RDKit format.
            mol_b: molecule B in RDKit format.
            ring_compare_method: method to compare the ring atoms. Choices are `element`, `hybridization`,
                and `aromaticity`.
            strict_surround: apply a strict method of `compare_molecule_rings()`. If True, it will
                purge rings and their not ring-linker substituents based on both ring structure
                & immediate neighbors. Defaults to False.
            ignore_substituent_atom_type: if true, it will ignore the atomic number of the substituents
                when comparing them for the `strict_surround` method. Defaults to False.

        Returns:
            a dictionary with the final atoms to be restrained.
        """
        matching_atoms = {}
        purge_list = []  # make this for the ring atoms that aren't equal including neighbors
        data = deepcopy(data)
        data.pop("Mapped Rings")
        connections = {k: v for k, v in data.items() if "->" in k}

        conserve_atoms = []  # we want to conserve the atoms that connect to other rings
        for ring_key, ring_data in data.items():
            if "->" in ring_key:
                continue
            if "Ring" not in ring_key:
                continue
            if self.is_ring_equivalent(ring_data, mol_a, mol_b, atom_mapping, ring_compare_method):
                matching_atoms[ring_key] = set(ring_data["ringAtomsA"])
                if strict_surround:
                    molA_ring_and_subs = self.get_ring_and_ringbound(
                        mol_a, ring_data["substituentsA"], ring_compare_method, ignore_surround_atom_type
                    )
                    compareB = {}
                    for key in ring_data["substituentsA"]:
                        compareB[atom_mapping[key]] = ring_data["substituentsB"][atom_mapping[key]]
                    molB_ring_and_subs = self.get_ring_and_ringbound(
                        mol_b, compareB, ring_compare_method, ignore_surround_atom_type
                    )
                    is_same = np.equal(molA_ring_and_subs, molB_ring_and_subs).all()
                    # check for len(v) > 2 to avoid conserving the restraints for the ring itself (ring-ring bond)
                    # logger.trace(data)
                    ring_connections = {k: v for k, v in connections.items() if k.startswith(ring_key)}
                    for conserve in ring_connections.values():
                        conserve_atoms.extend(conserve)
                        logger.debug(f"conserving atoms: {conserve_atoms}")
                    if not is_same:
                        to_purge = np.unique(
                            [item for k, vals in ring_data["substituentsA"].items() for item in [k] + vals]
                        )
                        purge_list.extend(to_purge)  # purge all rings & substituents another ring
            else:
                logger.trace(f"`{ring_key}` is not equivalent, purging atoms: {ring_data['ringAtomsA']}")
                purge_list.extend(ring_data["ringAtomsA"])
        purge_list = np.setdiff1d(purge_list, conserve_atoms).tolist()  # conserve connection to rings

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

    def set_restraints(
        self,
        ring_compare_method: str = "element",
        strict_surround: bool = False,
        ignore_surround_atom_type: bool = False,
    ) -> dict:
        """Main method for the class `RestraintSetter`. Run all the methods necessary to compare two ligands
        to be transmuted, find the atoms to be restrained, and return a dictionary with the atoms to be restrained.

        Args:
            ring_compare_method: method to compare the ring atoms. Choices are `element`, `hybridization`,
                and `aromaticity`.
            strict_surround: apply a strict method of `compare_molecule_rings()`. If True, it will
                purge rings and their not ring-linker substituents based on both ring structure
                & immediate neighbors. Defaults to False.
            ignore_substituent_atom_type: if true, it will ignore the atomic number of the substituents
                when comparing them for the `strict_surround` method. Defaults to False.

        Returns:
            restraints: dictionary containing the atoms to be restrained.
        """
        self.atom_mapper = AtomMapperHelper()
        ringStruc_compareDict = self.atom_mapper.process_rings_separately(
            Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit()), self.atom_mapping
        )
        restraints = self.compare_molecule_rings(
            ringStruc_compareDict,
            self.atom_mapping,
            self.molA.to_rdkit(),
            self.molB.to_rdkit(),
            ring_compare_method=ring_compare_method,
            strict_surround=strict_surround,
            ignore_surround_atom_type=ignore_surround_atom_type,
        )
        logger.debug(f"Atoms to restrain: {restraints}")
        self.restraints = restraints
        return restraints
