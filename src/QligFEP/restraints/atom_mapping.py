"""Module for atom mapping utilities, to be used to set FEP restraints. Main function is `process_rings_separately`"""

from itertools import combinations

from rdkit import Chem

from ..logger import logger


class AtomMapperHelper:
    def __init__(self) -> None:
        self.ringIdxsA = None
        self.ringIdxsB = None

    def get_surrounding_idxs(self, atom, mol):
        return [a.GetIdx() for a in mol.GetAtomWithIdx(atom).GetNeighbors()]

    def identify_and_enumerate_rings(self, mol):
        """Identify distinct ring structures and enumerate them."""
        ring_info = mol.GetRingInfo()
        rings = []
        for ring in ring_info.AtomRings():
            rings.append(set(ring))

        for ring1, ring2 in combinations(rings, 2):  # check for ring overlap & merge
            if ring1.intersection(ring2):
                rings.remove(ring1)
                rings.remove(ring2)
                rings.append(ring1.union(ring2))

        rings = {idx: ring for idx, ring in enumerate(rings)}
        return rings

    def map_rings_between_molecules(
        self, ringsA, ringsB, atom_mapping
    ):  # TODO: test when you don't have respective ring
        """Map rings between molecules based on atom mapping."""
        ring_mapping = {}  # Maps ring indices in A to ring indices in B
        for idxA, atomsA in ringsA.items():
            for idxB, atomsB in ringsB.items():
                if all(atom_mapping.get(a) in atomsB for a in atomsA):
                    ring_mapping[idxA] = idxB
                    break
        return ring_mapping

    def process_rings_separately(self, molA: Chem.Mol, molB: Chem.Mol, atom_mapping):
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
        self.ringIdxsA = [item for sublist in molA.GetRingInfo().AtomRings() for item in sublist]
        ringsB = self.identify_and_enumerate_rings(molB)
        self.ringIdxsB = [item for sublist in molB.GetRingInfo().AtomRings() for item in sublist]
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
                        "ringAtomsA": atomsA,
                        "ringAtomsB": atomsB,
                        "substituentsA": substituentsA,  # here we have ringAtom : [substituent atoms]
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
                            if bool(set(neighbors) & set(to_ring_atoms)):
                                rings.update({f"{key} -> {secondkey}": list(subsAtoms)})
                else:
                    continue
        return rings

    def is_atom_in_other_ring(self, atom, current_ring_atoms, mol):
        """Check if the atom is part of a ring that is not the current one being processed."""
        # Check all rings the atom is a part of, if any is fully outside current_ring_atoms, it's another ring
        for ring in mol.GetRingInfo().AtomRings():
            if atom.GetIdx() in ring and not all(atom_idx in current_ring_atoms for atom_idx in ring):
                return True
        return False

    def walk_bonds(self, atom, found_atoms, visited_atoms, mol, current_ring_atoms, a_or_b):
        atom_idx = atom.GetIdx()
        if atom_idx in visited_atoms:
            return found_atoms, visited_atoms  # Return immediately if atom has been visited

        if atom_idx not in found_atoms:
            found_atoms.append(atom_idx)
        visited_atoms.add(atom_idx)  # Mark this atom as visited

        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in current_ring_atoms:
                continue

            rings_idxs = self.ringIdxsA if a_or_b == "a" else self.ringIdxsB
            if neighbor_idx in rings_idxs:
                continue

            if neighbor_idx not in visited_atoms:
                # Recursively walk bonds from this neighbor
                found_atoms, visited_atoms = self.walk_bonds(
                    neighbor, found_atoms, visited_atoms, mol, current_ring_atoms, a_or_b
                )

        return found_atoms, visited_atoms

    def get_substituents(self, mol, ring_atoms, a_or_b: str):
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
                    substituents[atom_idx] = list(
                        new_found  # Convert set to list if needed for downstream processing
                    )
            # after completing the for loop, add itself if no substituents were found
            if atom_idx not in substituents:
                substituents[atom_idx] = [atom_idx]

        return substituents
