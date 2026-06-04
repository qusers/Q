"""Tests for single-topology atom mapping (single-Hamiltonian FEP)."""

import math
from pathlib import Path

import pytest
from rdkit import Chem

from QligFEP.single_topology import map_single_topology

RESOURCES_DIR = Path(__file__).parent / "resources"


@pytest.fixture(scope="module")
def tyk2_mols():
    sdf = RESOURCES_DIR / "tyk2_ligands.sdf"
    return {m.GetProp("_Name"): m
            for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None}


def _xyz(conf, i):
    p = conf.GetAtomPosition(i)
    return (p.x, p.y, p.z)


def test_ejm31_ejm47_partition(tyk2_mols):
    """rdFMCS + positional post-process gives the validated 30/2/9 partition."""
    a, b = tyk2_mols["ejm_31"], tyk2_mols["ejm_47"]
    m = map_single_topology(a, b)
    assert len(m.shared) == 30
    assert len(m.vanish) == 2          # ejm_31 methyl H's lost to the R-group swap
    assert len(m.appear) == 9          # ejm_47 R-group atoms


def test_endpoints_reconstruct_each_ligand(tyk2_mols):
    """shared+vanish must be all of A, shared+appear all of B."""
    a, b = tyk2_mols["ejm_31"], tyk2_mols["ejm_47"]
    m = map_single_topology(a, b)
    assert len(m.shared) + len(m.vanish) == a.GetNumAtoms()   # 32 = ejm_31
    assert len(m.shared) + len(m.appear) == b.GetNumAtoms()   # 39 = ejm_47
    # every A index appears exactly once across shared+vanish, same for B
    a_idx = sorted([i for i, _ in m.shared] + list(m.vanish))
    b_idx = sorted([j for _, j in m.shared] + list(m.appear))
    assert a_idx == list(range(a.GetNumAtoms()))
    assert b_idx == list(range(b.GetNumAtoms()))


def test_no_residual_same_element_overlap(tyk2_mols):
    """The key property: the post-process leaves no vanishing/appearing same-element
    pair within the overlap threshold (the substitution-point overlap that crashes
    dense-solvent dynamics is mapped away)."""
    a, b = tyk2_mols["ejm_31"], tyk2_mols["ejm_47"]
    thresh = 0.6
    m = map_single_topology(a, b, overlap_thresh=thresh)
    ca, cb = a.GetConformer(), b.GetConformer()
    for i in m.vanish:
        for j in m.appear:
            if a.GetAtomWithIdx(i).GetAtomicNum() == b.GetAtomWithIdx(j).GetAtomicNum():
                assert math.dist(_xyz(ca, i), _xyz(cb, j)) >= thresh


def test_mapping_never_crosses_elements(tyk2_mols):
    """Shared pairs must be element-matched (no H->C, unlike relaxed kartograf)."""
    a, b = tyk2_mols["ejm_31"], tyk2_mols["ejm_47"]
    m = map_single_topology(a, b)
    for i, j in m.shared:
        assert a.GetAtomWithIdx(i).GetAtomicNum() == b.GetAtomWithIdx(j).GetAtomicNum()


def test_requires_conformers():
    a = Chem.MolFromSmiles("CCO")
    b = Chem.MolFromSmiles("CCC")
    with pytest.raises(ValueError):
        map_single_topology(a, b)
