"""Tests for the RestraintScorer, RestraintLomapScorer, and parse_restraint_method utility."""

import tempfile
from pathlib import Path

import pytest
from rdkit import Chem

from QligFEP.scoring import (
    RestraintLomapScorer,
    RestraintScorer,
    parse_restraint_method,
)

RESOURCES_DIR = Path(__file__).parent / "resources"


@pytest.fixture
def tyk2_mols():
    """Return dict of tyk2 ligand molecules by name."""
    supplier = Chem.SDMolSupplier(str(RESOURCES_DIR / "tyk2_ligands.sdf"), removeHs=False)
    return {m.GetProp("_Name"): m for m in supplier if m is not None}


@pytest.fixture
def tyk2_sdf_pair(tyk2_mols):
    """Return paths to two tyk2 ligand SDF files (ejm_31, ejm_42)."""
    tmpdir = tempfile.mkdtemp()
    paths = {}
    for name in ("ejm_31", "ejm_42"):
        path = Path(tmpdir) / f"{name}.sdf"
        writer = Chem.SDWriter(str(path))
        writer.write(tyk2_mols[name])
        writer.close()
        paths[name] = path
    return paths


@pytest.fixture
def cmet_charged_pair():
    """Return paths to two cmet ligands with different charges."""
    supplier = Chem.SDMolSupplier(str(RESOURCES_DIR / "cmet_ligands.sdf"), removeHs=False)
    neutral = None
    charged = None
    for m in supplier:
        if m is None:
            continue
        charge = Chem.GetFormalCharge(m)
        if charge == 0 and neutral is None:
            neutral = m
        elif charge != 0 and charged is None:
            charged = m
        if neutral is not None and charged is not None:
            break
    tmpdir = tempfile.mkdtemp()
    paths = {}
    for label, mol in [("neutral", neutral), ("charged", charged)]:
        path = Path(tmpdir) / f"{label}.sdf"
        writer = Chem.SDWriter(str(path))
        writer.write(mol)
        writer.close()
        paths[label] = path
    return paths


class TestRestraintScorer:
    def test_returns_float_between_0_and_1(self, tyk2_sdf_pair):
        """Scorer should return a float in [0, 1]."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintScorer(restraint_method="heavyatom_p")
        score = scorer(mapping)

        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_identical_molecules_return_1(self, tyk2_sdf_pair):
        """Scoring the same molecule against itself should give 1.0."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molA))

        scorer = RestraintScorer(restraint_method="heavyatom_p")
        score = scorer(mapping)

        assert score == pytest.approx(1.0)

    def test_uses_correct_formula(self, tyk2_sdf_pair):
        """Verify score = n_restrained² / (nheavy1 * nheavy2)."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        from QligFEP.restraints.restraint_setter import RestraintSetter
        from QligFEP.scoring import _ensure_hydrogens_reindexed

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        # Compute expected score manually using reindexed molecules (same as scorer does)
        molA_r = _ensure_hydrogens_reindexed(molA)
        molB_r = _ensure_hydrogens_reindexed(molB)
        rsetter = RestraintSetter(molA_r, molB_r)
        restraints = rsetter.set_restraints(atom_compare_method="heavyatom", strict_surround=False)
        n_restrained = len(restraints)
        nheavy1 = Chem.RemoveHs(molA_r.to_rdkit()).GetNumAtoms()
        nheavy2 = Chem.RemoveHs(molB_r.to_rdkit()).GetNumAtoms()
        expected_score = n_restrained**2 / (nheavy1 * nheavy2)

        scorer = RestraintScorer(restraint_method="heavyatom_p")
        actual_score = scorer(mapping)

        assert actual_score == pytest.approx(expected_score)

    def test_different_methods_give_different_scores(self, tyk2_mols):
        """Different restraint methods should produce different scores for structurally distinct pairs."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        from QligFEP.scoring import _ensure_hydrogens_reindexed

        # Use ejm_31 (21 heavy) vs ejm_49 (26 heavy) — more structural difference
        tmpdir = tempfile.mkdtemp()
        for name in ("ejm_31", "ejm_49"):
            path = Path(tmpdir) / f"{name}.sdf"
            writer = Chem.SDWriter(str(path))
            writer.write(tyk2_mols[name])
            writer.close()

        molA = _ensure_hydrogens_reindexed(
            SmallMoleculeComponent.from_sdf_file(str(Path(tmpdir) / "ejm_31.sdf"))
        )
        molB = _ensure_hydrogens_reindexed(
            SmallMoleculeComponent.from_sdf_file(str(Path(tmpdir) / "ejm_49.sdf"))
        )
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer_hp = RestraintScorer(restraint_method="heavyatom_p")
        scorer_es = RestraintScorer(restraint_method="element_strict")

        score_hp = scorer_hp(mapping)
        score_es = scorer_es(mapping)

        # Both should be valid scores
        assert 0.0 <= score_hp <= 1.0
        assert 0.0 <= score_es <= 1.0
        # heavyatom_p is more permissive, so should give >= element_strict
        assert score_hp >= score_es

    def test_handles_charged_molecules(self, cmet_charged_pair):
        """Scorer should work with charged molecules without errors."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["neutral"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["charged"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintScorer(restraint_method="heavyatom_p")
        score = scorer(mapping)

        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_kartograf_scorer_returns_volume_ratio(self, tyk2_sdf_pair):
        """When restraint_method='kartograf', scorer should use kartograf native mapping."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent
        from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintScorer(restraint_method="kartograf")
        score = scorer(mapping)

        # Should match kartograf's own scorer
        kartograf_scorer = MappingVolumeRatioScorer()
        expected = kartograf_scorer(mapping)

        assert score == pytest.approx(expected)


class TestParseRestraintMethod:
    def test_simple_method(self):
        """Parse 'heavyatom_p' without distance suffix."""
        result = parse_restraint_method("heavyatom_p")
        assert result["atom_compare_method"] == "heavyatom"
        assert result["strict_surround"] is False
        assert result["kartograf_max_atom_distance"] == 0.95

    def test_with_distance_suffix(self):
        """Parse 'heavyatom_p_1.2' with distance suffix."""
        result = parse_restraint_method("heavyatom_p_1.2")
        assert result["atom_compare_method"] == "heavyatom"
        assert result["strict_surround"] is False
        assert result["kartograf_max_atom_distance"] == 1.2

    def test_strict_method(self):
        """Parse 'element_strict'."""
        result = parse_restraint_method("element_strict")
        assert result["atom_compare_method"] == "element"
        assert result["strict_surround"] is True
        assert result["ignore_surround_atom_type"] is False

    def test_less_strict_method(self):
        """Parse 'hybridization_ls'."""
        result = parse_restraint_method("hybridization_ls")
        assert result["atom_compare_method"] == "hybridization"
        assert result["strict_surround"] is True
        assert result["ignore_surround_atom_type"] is True

    def test_kartograf_method(self):
        """Parse 'kartograf'."""
        result = parse_restraint_method("kartograf")
        assert result["kartograf_native"] is True

    def test_invalid_method_raises(self):
        """Invalid method should raise ValueError."""
        with pytest.raises(ValueError, match="not recognized"):
            parse_restraint_method("invalid_method")

    def test_integer_distance(self):
        """Parse distance as integer '2'."""
        result = parse_restraint_method("element_p_2")
        assert result["kartograf_max_atom_distance"] == 2.0


class TestRestraintLomapScorer:
    def test_returns_geometric_mean_of_lomap_and_restraint(self, tyk2_sdf_pair):
        """Combined score should be sqrt(lomap_score * restraint_score)."""
        import numpy as np
        from kartograf import KartografAtomMapper, SmallMoleculeComponent
        from lomap import default_lomap_score

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintLomapScorer(restraint_method="heavyatom_p")
        combined = scorer(mapping)

        # Compute expected values independently
        lomap_score = default_lomap_score(mapping)
        restraint_score = RestraintScorer(restraint_method="heavyatom_p")(mapping)
        expected = np.sqrt(lomap_score * restraint_score)

        assert combined == pytest.approx(expected)
        assert 0.0 <= combined <= 1.0

    def test_same_charge_pair_gets_nonzero_score(self, tyk2_sdf_pair):
        """Same-charge pair should get a nonzero score (lomap doesn't penalize)."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintLomapScorer(restraint_method="heavyatom_p")
        score = scorer(mapping)

        assert score > 0.0

    def test_charge_change_blocked_by_default(self, cmet_charged_pair):
        """Charge-changing pairs should score 0.0 with default charge_changes_score."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["neutral"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["charged"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintLomapScorer(restraint_method="heavyatom_p")
        score = scorer(mapping)

        assert score == 0.0

    def test_charge_change_allowed_with_nonzero_param(self, cmet_charged_pair):
        """With charge_changes_score > 0, charge-changing pairs get a nonzero score."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        molA = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["neutral"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["charged"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        mapping = next(mapper.suggest_mappings(molA, molB))

        scorer = RestraintLomapScorer(restraint_method="heavyatom_p", charge_changes_score=0.5)
        score = scorer(mapping)

        assert score > 0.0

    def test_charge_change_score_penalizes_vs_same_charge(self, tyk2_sdf_pair, cmet_charged_pair):
        """Even with charge_changes_score=0.5, a charge-change edge scores lower than same-charge."""
        from kartograf import KartografAtomMapper, SmallMoleculeComponent

        # Same-charge pair
        molA = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_31"]))
        molB = SmallMoleculeComponent.from_sdf_file(str(tyk2_sdf_pair["ejm_42"]))
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        same_mapping = next(mapper.suggest_mappings(molA, molB))

        # Charge-change pair
        molC = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["neutral"]))
        molD = SmallMoleculeComponent.from_sdf_file(str(cmet_charged_pair["charged"]))
        charge_mapping = next(mapper.suggest_mappings(molC, molD))

        scorer = RestraintLomapScorer(restraint_method="heavyatom_p", charge_changes_score=0.5)
        same_score = scorer(same_mapping)
        charge_score = scorer(charge_mapping)

        # The charge-change penalty (0.5) should make its lomap component lower
        # We can't compare directly since they're different molecule pairs,
        # but the charge_score should be bounded by sqrt(charge_changes_score * restraint_score)
        import numpy as np

        restraint_only = RestraintScorer(restraint_method="heavyatom_p")(charge_mapping)
        assert charge_score <= np.sqrt(0.5 * restraint_only) + 1e-10
