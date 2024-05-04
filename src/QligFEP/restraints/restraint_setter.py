from copy import deepcopy
from itertools import zip_longest
from pathlib import Path

from kartograf import KartografAtomMapper, SmallMoleculeComponent
from kartograf.atom_aligner import align_mol_shape
from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer
from rdkit import Chem

from ..logger import logger
from .atom_mapping import process_rings_separately


class RestraintSetter:
    """
    TODO: currently, this class correctly sets the restraints for the ring structures,
    but it takes into account:
        1) the ring structure itself
        2) the substitutions
    In the future, it should also account for:
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
        mapper = KartografAtomMapper(atom_map_hydrogens=True)
        a_molB = align_mol_shape(self.molB, ref_mol=self.molA)
        self.kartograf_mapping = next(mapper.suggest_mappings(self.molA, a_molB))
        # Score Mapping
        rmsd_scorer = MappingVolumeRatioScorer()
        score = rmsd_scorer(mapping=self.kartograf_mapping)
        logger.info(f"Volume ratio score between ligands: {score}")
        self.atom_mapping = deepcopy(self.kartograf_mapping.to_dict()["componentA_to_componentB"])

    @staticmethod
    def are_atoms_equivalent(atom_a, atom_b):
        return atom_a.GetAtomicNum() == atom_b.GetAtomicNum()

    def are_substituents_equivalent(self, subsA, subsB, atom_mapping, mol_a, mol_b):
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
                if self.are_atoms_equivalent(atomA_obj, atomB_obj):
                    is_same.add(atomA)  # Add original ring atom if all checks pass
                else:
                    break
        return is_same

    def is_ring_equivalent(self, ring_data, mol_a, mol_b, atom_mapping):
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

    def compare_molecule_rings(self, data, atom_mapping, mol_a, mol_b):
        matching_atoms = {}
        purge_list = []
        data = deepcopy(data)
        data.pop("Ring mapping", None)  # Remove safely with default

        for ring_key, ring_data in data.items():
            if "Ring" not in ring_key:
                continue
            if self.is_ring_equivalent(ring_data, mol_a, mol_b, atom_mapping):
                matching_atoms[ring_key] = set(ring_data["ringAtomsA"])
            else:
                purge_list.extend(ring_data["ringAtomsA"])

        # Process and compare substituents only for equivalent rings
        for ring_key in matching_atoms:
            ring_data = data[ring_key]
            subsA = ring_data["substituentsA"]
            subsB = ring_data["substituentsB"]

            matched_subs = self.are_substituents_equivalent(subsA, subsB, atom_mapping, mol_a, mol_b)
            matching_atoms[ring_key].update(matched_subs)

        all_atomA_idxs = set(
            idx for ring_atoms in matching_atoms.values() for idx in ring_atoms if idx not in purge_list
        )
        restraints = {k: atom_mapping[k] for k in all_atomA_idxs if k in atom_mapping}
        return restraints

    def set_restraints(self):
        ringStruc_compareDict = process_rings_separately(
            Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit()), self.atom_mapping
        )
        restraints = self.compare_molecule_rings(
            ringStruc_compareDict, self.atom_mapping, self.molA.to_rdkit(), self.molB.to_rdkit()
        )

        # submol = extract_sub_molecule(np.unique(all_atomA_idxs).tolist(), self.molA.to_rdkit())
        # submol = remove_hydrogens(submol)
        # self.submol = submol

        # mols = [Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit())]
        # matches = [mol.GetSubstructMatches(submol)[0] for mol in mols]
        # restraints = {self.atom_mapping[k]: v for k, v in self.atom_mapping.items() if k in matches[0]}
        self.restraints = restraints
        return restraints
