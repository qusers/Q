"""Single-topology atom mapping for single-Hamiltonian FEP.

A single-topology hybrid represents the shared core of two congeneric ligands once
(those atoms morph A->B) and only the unique atoms grow/vanish. Building it needs an
atom map that, crucially, does NOT leave the substitution-point atoms split into an
overlapping vanishing+appearing pair: when one substituent replaces another at the
same position (e.g. a methyl C and an R-group C ~0.4 A apart), rdFMCS maps them as
separate unique atoms because their connectivity differs, and that overlap creates a
degenerate bonded geometry that blows up in dense solvent.

The mapper is therefore rdFMCS (connectivity) followed by an element-safe positional
post-process: any remaining unmapped atoms of the SAME element within `overlap_thresh`
of each other are paired, so the overlapping substitution-point atoms morph as one.
It deliberately never crosses elements (a relaxed positional mapper will map a methyl
H onto a ring C, which is chemically wrong).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdFMCS


@dataclass(frozen=True)
class SingleTopologyMap:
    """Result of mapping ligand A (endpoint, lambda_c=0) to ligand B (lambda_c=1).

    shared: list of (idxA, idxB) atom-index pairs that morph A->B (the core).
    vanish: idxA of atoms present only in A (real->DUM).
    appear: idxB of atoms present only in B (DUM->real).
    """
    shared: tuple[tuple[int, int], ...]
    vanish: tuple[int, ...]
    appear: tuple[int, ...]

    @property
    def a_to_b(self) -> dict[int, int]:
        return dict(self.shared)


def _conf_xyz(conf, i):
    p = conf.GetAtomPosition(i)
    return (p.x, p.y, p.z)


def map_single_topology(molA, molB, overlap_thresh: float = 0.6,
                        mcs_timeout: int = 30) -> SingleTopologyMap:
    """Map two RDKit mols (each with a 3D conformer) for single-topology FEP.

    rdFMCS connectivity map + element-safe positional post-process. Returns the
    shared / vanishing / appearing partition. `overlap_thresh` (Angstrom) is the
    distance under which two unmapped same-element atoms are treated as the same
    morphing atom (the substitution-point overlap); 0.6 A pairs the overlapping
    substituent atoms without touching normally-bonded neighbours (>~1 A).
    """
    if molA.GetNumConformers() == 0 or molB.GetNumConformers() == 0:
        raise ValueError("both molecules must have a 3D conformer for positional mapping")

    res = rdFMCS.FindMCS([molA, molB],
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         ringMatchesRingOnly=True, completeRingsOnly=True,
                         timeout=mcs_timeout)
    patt = Chem.MolFromSmarts(res.smartsString)
    mA = list(molA.GetSubstructMatch(patt))
    mB = list(molB.GetSubstructMatch(patt))
    shared = list(zip(mA, mB))

    # element-safe positional post-process for substitution-point overlaps
    cA, cB = molA.GetConformer(), molB.GetConformer()
    umA = [i for i in range(molA.GetNumAtoms()) if i not in mA]
    umB = [j for j in range(molB.GetNumAtoms()) if j not in mB]
    cand = sorted(
        (math.dist(_conf_xyz(cA, i), _conf_xyz(cB, j)), i, j)
        for i in umA for j in umB
        if molA.GetAtomWithIdx(i).GetAtomicNum() == molB.GetAtomWithIdx(j).GetAtomicNum()
    )
    usedA, usedB = set(), set()
    for d, i, j in cand:
        if d < overlap_thresh and i not in usedA and j not in usedB:
            shared.append((i, j))
            usedA.add(i)
            usedB.add(j)

    mappedA = {i for i, _ in shared}
    mappedB = {j for _, j in shared}
    vanish = tuple(i for i in range(molA.GetNumAtoms()) if i not in mappedA)
    appear = tuple(j for j in range(molB.GetNumAtoms()) if j not in mappedB)
    return SingleTopologyMap(shared=tuple(shared), vanish=vanish, appear=appear)
