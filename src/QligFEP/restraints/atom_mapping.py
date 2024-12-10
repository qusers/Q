"""Atom mapping utilities to be used to set FEP restraints. Main function is `process_rings_separately`"""

import numpy as np
from rdkit import Chem

from ..logger import logger


class AtomMapperHelper:
    """Helper class for atom mapping and ring processing, used by RestraintSetter class

    This class contains methods to identify and enumerate rings in a molecule, map rings
    between molecules based on their equivalency, and identify the respective decorations
    present in each of the molecule rings. The main method `process_rings_separately`
    returns a dictionary containing the mapped rings, the atoms in each ring, and the
    substituents for each ring in both molecules.
    """

    def __init__(self) -> None:
        self.mapped_ringIdxsA = None
        self.mapped_ringIdxsB = None
        self.all_ringIdxsA = None
        self.all_ringIdxsB = None

    def get_surrounding_idxs(self, atom: Chem.Atom, mol: Chem.Mol) -> list[int]:
        return [a.GetIdx() for a in mol.GetAtomWithIdx(atom).GetNeighbors()]

    def identify_and_enumerate_rings(self, mol: Chem.Mol) -> dict:
        """Identify ring structures and enumerate them into a dictionary"""

        def merge_sets(sets: list[set]) -> list[set]:
            """merge sets within a list based on shared numbers (rings that share atoms)."""
            changed = True
            while changed:
                changed = False
                new_sets = []
                for i, s in enumerate(sets):
                    merged = False
                    for j in range(i + 1, len(sets)):
                        if s & sets[j]:  # If there's an intersection
                            new_sets.append(s | sets[j])  # Merge the sets
                            sets[j] = set()  # Empty the merged set
                            merged = True
                            changed = True
                            break
                    if not merged and s:
                        new_sets.append(s)
                sets = [s for s in new_sets if s]  # Remove empty sets
            return sets

        ring_info = mol.GetRingInfo()
        rings = merge_sets([set(ring) for ring in ring_info.AtomRings()])
        rings = {idx: ring for idx, ring in enumerate(rings)}
        return rings

    def map_rings_between_molecules(self, ringsA: dict, ringsB: dict, atom_mapping: dict) -> dict:
        """Map rings between molecules based on atom mapping."""
        ring_mapping = {}  # Maps ring indices in A to ring indices in B
        for idxA, atomsA in ringsA.items():
            for idxB, atomsB in ringsB.items():
                if all(atom_mapping.get(a) in atomsB for a in atomsA):
                    ring_mapping[idxA] = idxB
                    break
        return ring_mapping

    def process_rings_separately(self, molA: Chem.Mol, molB: Chem.Mol, atom_mapping: dict) -> dict[dict]:
        """Process each ring separately based on ring equivalency and their substituents.
        The returned dictionary is the main output for the restraint setting process. It contains
        the mapping between rings, the atoms in each ring, and the substituents for each ring
        in both molecules.

        Args:
            molA: RDKit molecule object for molecule A.
            molB: RDKit molecule object for molecule B.
            atom_mapping: Kartograf mapping between atoms in molecule A and molecule B.

        Returns:
            dict: Dictionary containing the Mapped Rings, the atoms in each ring, and the substituents
                for each ring in both molecules.
        Example:
            >>> {'Mapped Rings': {0: 0, 1: 1, 2: 2, 3: 3},
            >>> 'Ring 0': {'ringAtomsA': {set of AtomIndexes in ring 0 of molecule A},
            >>> 'ringAtomsB': {set of AtomIndexes in ring 0 of molecule B},
            >>> 'substituentsA': {AtominRing0 for molecule A: [substituent atom indexes]
            >>>    ...},
            >>> 'substituentsB': {AtominRing0 for molecule B: [substituent atom indexes]
            >>>    ...},
            >>> 'Ring 1': ...,
            >>> }
        """
        ringsA = self.identify_and_enumerate_rings(molA)
        logger.debug(f"Rings in molecule A: {ringsA}")
        ringsB = self.identify_and_enumerate_rings(molB)
        logger.debug(f"Rings in molecule B: {ringsB}")
        self.all_ringIdxsA = [idx for set_values in ringsA.values() for idx in set_values]
        self.all_ringIdxsB = [idx for set_values in ringsB.values() for idx in set_values]
        self.mapped_ringIdxsA = []
        self.mapped_ringIdxsB = []

        for idxs in ringsA.values():  # only add rings that are fully mapped by kartograf
            if np.isin(list(idxs), list(atom_mapping.keys())).all():
                self.mapped_ringIdxsA.extend(idxs)

        for idxs in ringsB.values():
            if np.isin(list(idxs), list(atom_mapping.values())).all():
                self.mapped_ringIdxsB.extend(idxs)

        rings = {}
        ring_mapping = self.map_rings_between_molecules(ringsA, ringsB, atom_mapping)
        rings.update({"Mapped Rings": ring_mapping})

        for ring_idx, (idxA, idxB) in enumerate(ring_mapping.items()):
            atomsA = ringsA[idxA]
            atomsB = ringsB[idxB]
            substituentsA = self.get_substituents(molA, atomsA, a_or_b="a")
            substituentsB = self.get_substituents(molB, atomsB, a_or_b="b")
            rings.update(
                {
                    f"Ring {ring_idx}": {
                        "ringAtomsA": atomsA,  # set of atom indexes in ring
                        "ringAtomsB": atomsB,
                        "substituentsA": substituentsA,  #  ringAtom : [substituent atoms]
                        "substituentsB": substituentsB,
                    }
                }
            )
        ring_keys = [key for key in rings if key.startswith("Ring")]
        for key in ring_keys:
            for ringAtom, subsAtoms in rings[key]["substituentsA"].items():
                if ringAtom != subsAtoms[-1]:
                    neighbors = self.get_surrounding_idxs(subsAtoms[-1], molA)
                    for secondkey in ring_keys:
                        to_ring_atoms = rings[secondkey]["ringAtomsA"]
                        if secondkey == key:
                            continue
                        else:
                            if bool(set(neighbors) & set(to_ring_atoms)):  # add both link directions
                                rings.update({f"{key} -> {secondkey}": list(subsAtoms)})
                                rings.update({f"{secondkey} -> {key}": list(subsAtoms[::-1])})
                else:
                    continue
        return rings

    def is_atom_in_other_ring(self, atom, current_ring_atoms, mol):
        """Check if the atom is part of a ring that is not the current one being processed
        during  `get_substituents`"""
        for ring in mol.GetRingInfo().AtomRings():
            if atom.GetIdx() in ring and not all(atom_idx in current_ring_atoms for atom_idx in ring):
                return True
        return False

    def walk_bonds(
        self,
        atom: Chem.Atom,
        found_atoms: list,
        visited_atoms: set,
        mol: Chem.Mol,
        current_ring_atoms: list[int],
        a_or_b: str,
    ) -> tuple[list, set]:
        atom_idx = atom.GetIdx()
        rings_idxs = self.mapped_ringIdxsA if a_or_b == "a" else self.mapped_ringIdxsB
        if atom_idx in visited_atoms or atom_idx in rings_idxs:
            # Return immediately if atom has been visited or if it's part of a ring (e.g.: ring linked to another ring)
            return found_atoms, visited_atoms

        if atom_idx not in found_atoms:
            found_atoms.append(atom_idx)
        visited_atoms.add(atom_idx)  # Mark this atom as visited

        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in current_ring_atoms:
                continue

            if neighbor_idx in rings_idxs:
                continue

            if neighbor_idx not in visited_atoms:  # Recursively walk bonds from this neighbor
                found_atoms, visited_atoms = self.walk_bonds(
                    neighbor, found_atoms, visited_atoms, mol, current_ring_atoms, a_or_b
                )

        return found_atoms, visited_atoms

    def get_substituents(self, mol: Chem.Mol, ring_atoms: list[int], a_or_b: str) -> dict:
        substituents = {}
        visited_atoms = set(ring_atoms)  # Start with ring atoms as visited to avoid re-traversing the ring

        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            logger.trace(f"Processing atom {atom_idx} in {ring_atoms} for molecule {a_or_b}")

            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in ring_atoms:
                    continue

                # Explore substituents starting from this neighbor
                new_found, _ = self.walk_bonds(
                    neighbor, list(), visited_atoms.copy(), mol, ring_atoms, a_or_b
                )
                if new_found:  # If any new atoms were found, add them to the substituents list
                    substituents[atom_idx] = list(new_found)
                    # place the atom connecting to another ring in the end of the sequence
                    for aIdx in new_found[::-1]:
                        if self.is_atom_in_other_ring(mol.GetAtomWithIdx(aIdx), ring_atoms, mol):
                            substituents[atom_idx].remove(aIdx)
                            substituents[atom_idx].append(aIdx)
                            break
            if atom_idx not in substituents:  # add itself if no substituents were found
                substituents[atom_idx] = [atom_idx]

        return substituents
