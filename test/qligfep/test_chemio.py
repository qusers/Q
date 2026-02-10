"""Tests for MoleculeIO SDF loading and hydrogen handling."""

from pathlib import Path

import pytest
from loguru import logger

from QligFEP.chemIO import MoleculeIO

RESOURCES_DIR = Path(__file__).parent / "resources"

# Ground truth: H counts from tyk2_ligands.sdf (read with RDKit, removeHs=False)
EXPECTED_H_COUNTS = {
    "jmc_23": 12,
    "ejm_47": 15,
    "ejm_49": 13,
    "ejm_31": 11,
    "ejm_45": 15,
    "ejm_44": 17,
    "ejm_43": 15,
    "ejm_50": 11,
    "ejm_42": 13,
    "ejm_55": 11,
    "ejm_48": 17,
    "ejm_54": 14,
    "jmc_28": 15,
    "jmc_27": 12,
    "jmc_30": 12,
    "ejm_46": 13,
}


@pytest.fixture
def sdf_without_hydrogens():
    """Tyk2 ligands with all explicit hydrogens stripped."""
    return str(RESOURCES_DIR / "tyk2_ligands_no_hydrogens.sdf")


@pytest.fixture
def sdf_with_hydrogens():
    """Tyk2 ligands with explicit hydrogens as prepared by the user."""
    return str(RESOURCES_DIR / "tyk2_ligands.sdf")


@pytest.fixture
def capture_logs():
    """Capture loguru log messages into a list for assertion."""
    messages = []

    def sink(message):
        messages.append(message.record)

    handler_id = logger.add(sink, level="WARNING")
    yield messages
    logger.remove(handler_id)


class TestHydrogenWarning:
    """MoleculeIO should warn loudly when input molecules lack explicit hydrogens."""

    def test_warns_when_hydrogens_missing(self, sdf_without_hydrogens, capture_logs):
        """Loading an SDF without explicit H should produce a warning per molecule."""
        mio = MoleculeIO(sdf_without_hydrogens)

        hydrogen_warnings = [r for r in capture_logs if "hydrogen" in r["message"].lower()]
        assert len(hydrogen_warnings) > 0, (
            "Expected warnings about hydrogens being auto-added, but none were emitted"
        )

        # Each molecule without H should produce a warning
        for name in mio.lig_names:
            mol_warnings = [r for r in hydrogen_warnings if name in r["message"]]
            assert len(mol_warnings) > 0, (
                f"Expected a hydrogen warning for molecule '{name}', but none was found"
            )

    def test_no_warning_when_hydrogens_present(self, sdf_with_hydrogens, capture_logs):
        """Loading an SDF with explicit H should NOT produce hydrogen addition warnings."""
        MoleculeIO(sdf_with_hydrogens)

        hydrogen_addition_warnings = [
            r for r in capture_logs
            if "hydrogen" in r["message"].lower() and "added" in r["message"].lower()
        ]
        assert len(hydrogen_addition_warnings) == 0, (
            f"Did not expect hydrogen addition warnings for properly prepared molecules, "
            f"but got: {[r['message'] for r in hydrogen_addition_warnings]}"
        )

    def test_molecules_still_loaded_when_hydrogens_missing(self, sdf_without_hydrogens):
        """Even when H are missing, molecules should still be loaded (just with a warning)."""
        mio = MoleculeIO(sdf_without_hydrogens)
        assert len(mio.molecules) > 0, "Molecules should still be loaded despite missing H"
        # All molecules should have H atoms after loading (OpenFF adds them)
        for name, mol in mio:
            n_h = sum(1 for a in mol.atoms if a.atomic_number == 1)
            assert n_h > 0, f"Molecule {name} should have H atoms after OpenFF loading"


class TestHydrogenPreservation:
    """When input molecules have explicit H, MoleculeIO must preserve them exactly."""

    def test_preserves_exact_hydrogen_count(self, sdf_with_hydrogens):
        """Loaded molecules must have exactly the same H count as the input SDF."""
        mio = MoleculeIO(sdf_with_hydrogens)

        for name, mol in mio:
            expected_h = EXPECTED_H_COUNTS[name]
            actual_h = sum(1 for a in mol.atoms if a.atomic_number == 1)
            assert actual_h == expected_h, (
                f"Molecule '{name}': expected {expected_h} H atoms (from input), "
                f"but got {actual_h} after loading. Hydrogens must not be added or removed."
            )
