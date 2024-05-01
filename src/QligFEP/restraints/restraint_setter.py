from copy import deepcopy
from importlib.util import find_spec
from itertools import zip_longest
from pathlib import Path

import numpy as np
from kartograf import KartografAtomMapper, SmallMoleculeComponent
from kartograf.atom_aligner import align_mol_shape
from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer
from rdkit import Chem

from ..logger import logger
from .atom_mapping import process_rings_separately
from .render_functions import extract_sub_molecule, remove_hydrogens


class RestraintSetter:
    def __init__(self, molA: str, molB: str) -> None:
        self.ligname1 = Path(molA).stem
        self.ligname2 = Path(molB).stem
        self._load_molecules(molA, molB)

    def _path_to_mol(self, _path):
        if Path(_path).suffix == ".pdb":
            return Chem.MolFromPDBFile(_path)
        elif Path(_path).suffix == ".sdf":
            return Chem.SDMolSupplier(_path)[0]
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
        self.kartograph_mapping = next(mapper.suggest_mappings(self.molA, a_molB))
        # Score Mapping
        rmsd_scorer = MappingVolumeRatioScorer()
        score = rmsd_scorer(mapping=self.kartograf_mapping)
        print(f"RMSD Score: {score}")
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
                continue  # Skip if no corresponding atom or no substituents list
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

    def compare_molecule_rings(self, data, atom_mapping, mol_a, mol_b):
        matching_atoms = {}

        for ring_key, ring_data in data.items():
            if ring_key == "Ring mapping" or "Ring" not in ring_key:
                continue
            matching_atoms[ring_key] = set()
            ring_atoms_a_indices = ring_data["ringAtomsA"]
            # get the ring atoms for mol_b based on the atom mapping
            ring_atoms_b_indices = [atom_mapping.get(a) for a in ring_atoms_a_indices]
            # check if both rings are the same
            eq_array = []
            for atom_idx_a, atom_idx_b in zip(ring_atoms_a_indices, ring_atoms_b_indices):
                atom_a = mol_a.GetAtomWithIdx(atom_idx_a)
                atom_b = mol_b.GetAtomWithIdx(atom_idx_b)
                eq_array.append(self.are_atoms_equivalent(atom_a, atom_b))
            if all(eq_array):
                matching_atoms[ring_key] = set(ring_atoms_a_indices)
            else:
                continue

            subsA = ring_data["substituentsA"]
            subsB = ring_data["substituentsB"]

            # Compare substituents
            matched_subs = self.are_substituents_equivalent(subsA, subsB, atom_mapping, mol_a, mol_b)
            matching_atoms[ring_key] = matching_atoms[ring_key].union(matched_subs)

        return matching_atoms

    def get_restraints(self):
        ringStruc_compareDict = process_rings_separately(self.molA, self.molB, self.atom_mapping)
        rest_per_ring_atomIdxs = self.compare_molecule_rings(
            ringStruc_compareDict, self.atom_mapping, self.molA, self.molB
        )

        all_atom_numbers = []
        for _, value in rest_per_ring_atomIdxs.items():
            all_atom_numbers.extend(value)

        submol = extract_sub_molecule(np.unique(all_atom_numbers).tolist(), self.molA.to_rdkit())
        submol = remove_hydrogens(submol)
        self.submol = submol

        mols = [Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit())]
        matches = [mol.GetSubstructMatches(submol)[0] for mol in mols]
        restraints = {self.atom_mapping[k]: v for k, v in self.atom_mapping.items() if k in matches[0]}
        self.restraints = restraints
        return restraints

    def render_restraints(self, restraints):
        if find_spec("chem-filters"):
            from chemFilters.img_render import MolGridPlotter
        else:
            logger.info(
                "chem-filters not found, install it from: "
                "https://github.com/David-Araripe/chemFilters/blob/master/pyproject.toml"
            )
        mols = [Chem.RemoveHs(self.molA.to_rdkit()), Chem.RemoveHs(self.molB.to_rdkit())]
        matches = [mol.GetSubstructMatches(self.submol)[0] for mol in mols]
        plotter = MolGridPlotter(from_smi=False, add_atom_indices=True, size=(400, 400))
        imgs = []
        imgs.append(plotter.render_mol(mols[0], highlightAtoms=matches[0]))
        imgs.append(plotter.render_mol(mols[1], highlightAtoms=matches[1]))
        img = plotter._images_to_grid(imgs, n_cols=2)
        return img
