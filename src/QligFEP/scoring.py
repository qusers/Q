"""Scoring functions for perturbation network edge quality based on atom mapping similarity."""

import re

from gufe import AtomMapping
from kartograf import SmallMoleculeComponent
from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer
from rdkit import Chem

from .CLI.utils import get_avail_restraint_methods
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


def parse_restraint_method(restraint_method: str) -> dict:
    """Parse a restraint method string into RestraintSetter keyword arguments.

    Format: '{atom_compare}_{permissiveness}' with optional distance suffix '_{distance}'.
    Examples: 'heavyatom_p', 'element_strict_1.2', 'kartograf'

    Args:
        restraint_method: restraint method string to parse.

    Returns:
        Dictionary with keys appropriate for RestraintSetter.set_restraints() or
        {'kartograf_native': True} for the kartograf method.

    Raises:
        ValueError: if the method is not recognized.
    """
    if restraint_method == "kartograf":
        return {"kartograf_native": True, "kartograf_max_atom_distance": 0.95}

    # Check for optional atom max distance suffix
    distance_pattern = r"_(\d+\.?\d*)$"
    match = re.search(distance_pattern, restraint_method)
    if match:
        kartograf_max_atom_distance = float(match.group(1))
        restraint_method = re.sub(distance_pattern, "", restraint_method)
    else:
        kartograf_max_atom_distance = 0.95

    avail_methods = get_avail_restraint_methods()
    if restraint_method not in avail_methods:
        raise ValueError(f"Method {restraint_method} not recognized. Please use one of {avail_methods}")

    atom_compare_method, permissiveness_lvl = restraint_method.split("_")
    if permissiveness_lvl == "p":
        params = {"strict_surround": False, "ignore_surround_atom_type": False}
    elif permissiveness_lvl == "ls":
        params = {"strict_surround": True, "ignore_surround_atom_type": True}
    elif permissiveness_lvl == "strict":
        params = {"strict_surround": True, "ignore_surround_atom_type": False}

    return {
        "atom_compare_method": atom_compare_method,
        "kartograf_max_atom_distance": kartograf_max_atom_distance,
        **params,
    }


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
