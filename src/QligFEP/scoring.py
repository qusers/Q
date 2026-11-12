"""Scoring functions for perturbation network edge quality based on atom mapping similarity."""

from functools import partial

import numpy as np
from gufe import AtomMapping
from kartograf import SmallMoleculeComponent
from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer
from lomap import default_lomap_score
from rdkit import Chem

from .CLI.utils import parse_restraint_method
from .logger import logger
from .restraints.hydrogen_utils import are_hydrogens_at_end, reindex_hydrogens_to_end
from .restraints.restraint_setter import RestraintSetter


def _ensure_hydrogens_reindexed(mol: SmallMoleculeComponent) -> SmallMoleculeComponent:
    """Reindex hydrogens to the end of the atom list if needed.

    RestraintSetter requires hydrogens at the end. Molecules coming from konnektor
    mappings may not satisfy this, so we reindex them here.
    """
    rdmol = mol.to_rdkit()
    if not are_hydrogens_at_end(rdmol):
        rdmol = reindex_hydrogens_to_end(rdmol)
        return SmallMoleculeComponent.from_rdkit(rdmol)
    return mol


class RestraintScorer:
    """Score atom mappings based on RestraintSetter atom pair count.

    Implements the konnektor scorer interface: `scorer(AtomMapping) -> float`.
    Score formula: n_restrained² / (nheavy_lig1 * nheavy_lig2)
    where n_restrained is the number of atom pairs from RestraintSetter.

    For the 'kartograf' method, delegates to kartograf's MappingVolumeRatioScorer.

    Args:
        restraint_method: restraint method string (e.g. 'heavyatom_p', 'element_strict_1.2', 'kartograf').
    """

    def __init__(self, restraint_method: str = "heavyatom_p"):
        self.parsed = parse_restraint_method(restraint_method)
        self._kartograf_native = self.parsed.pop("kartograf_native", False)
        self._kartograf_max_atom_distance = self.parsed.pop("kartograf_max_atom_distance", 0.95)
        if self._kartograf_native:
            self._kartograf_scorer = MappingVolumeRatioScorer()

    def __call__(self, mapping: AtomMapping) -> float:
        if self._kartograf_native:
            return self._kartograf_scorer(mapping)

        molA = _ensure_hydrogens_reindexed(mapping.componentA)
        molB = _ensure_hydrogens_reindexed(mapping.componentB)

        rsetter = RestraintSetter(molA, molB, kartograf_max_atom_distance=self._kartograf_max_atom_distance)
        restraints = rsetter.set_restraints(**self.parsed)
        n_restrained = len(restraints)

        nheavy_A = Chem.RemoveHs(molA.to_rdkit()).GetNumAtoms()
        nheavy_B = Chem.RemoveHs(molB.to_rdkit()).GetNumAtoms()

        if nheavy_A == 0 or nheavy_B == 0:
            logger.warning("One of the molecules has 0 heavy atoms, returning score 0.0")
            return 0.0

        score = n_restrained**2 / (nheavy_A * nheavy_B)
        logger.debug(
            f"RestraintScorer: n_restrained={n_restrained}, nheavy_A={nheavy_A}, nheavy_B={nheavy_B}, score={score:.4f}"
        )
        return score


class RestraintLomapScorer:
    """Combine LOMAP and restraint scores, then penalize charge changes.

    The unpenalized score is ``sqrt(lomap_score * restraint_score)``. LOMAP's
    charge-change rule is disabled while calculating that geometric mean so that
    ``charge_change_multiplier`` can be applied directly to the final score.

    Args:
        restraint_method: restraint method string for RestraintScorer.
        charge_change_multiplier: Multiplier applied to the final combined score
            when the ligands differ in net formal charge. ``0.5`` halves the final
            score, ``1.0`` applies no penalty, and ``0.0`` assigns a score of zero
            without excluding the edge. Default: 0.0.
    """

    def __init__(self, restraint_method: str = "heavyatom_p", charge_change_multiplier: float = 0.0):
        if not 0.0 <= charge_change_multiplier <= 1.0:
            raise ValueError("charge_change_multiplier must be between 0.0 and 1.0")
        self.charge_change_multiplier = charge_change_multiplier
        self._restraint_scorer = RestraintScorer(restraint_method=restraint_method)
        self._lomap_score = partial(default_lomap_score, charge_changes_score=1.0)

    @staticmethod
    def is_charge_change(mapping: AtomMapping) -> bool:
        """Return whether the mapped ligands have different net formal charges."""
        charge_a = Chem.GetFormalCharge(mapping.componentA.to_rdkit())
        charge_b = Chem.GetFormalCharge(mapping.componentB.to_rdkit())
        return charge_a != charge_b

    def base_score(self, mapping: AtomMapping) -> float:
        """Return the combined score before applying the charge-change multiplier."""
        lomap = self._lomap_score(mapping)
        if lomap == 0.0:
            return 0.0
        restraint = self._restraint_scorer(mapping)
        return float(np.sqrt(lomap * restraint))

    def __call__(self, mapping: AtomMapping) -> float:
        base_score = self.base_score(mapping)
        multiplier = self.charge_change_multiplier if self.is_charge_change(mapping) else 1.0
        score = base_score * multiplier
        logger.debug(
            f"RestraintLomapScorer: base={base_score:.4f}, multiplier={multiplier:.4f}, "
            f"combined={score:.4f}"
        )
        return score
