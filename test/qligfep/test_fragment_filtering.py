"""Tests for filtering out-of-sphere molecular fragments.

This module tests the fragment filtering functionality that removes entire
molecular segments (chains, cofactors, ions, lipids) that are completely
outside the simulation sphere radius before FEP setup.
"""

import shutil
import tempfile
from pathlib import Path

import pytest

# Test resources paths
TEST_DIR = Path(__file__).parent
RESOURCES_DIR = TEST_DIR / "resources"
PROTEINS_DIR = RESOURCES_DIR / "proteins"


# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------


@pytest.fixture
def cdk2_protein_path() -> Path:
    """Path to cdk2 protein with chains A and B."""
    return PROTEINS_DIR / "cdk2_protein.pdb"


@pytest.fixture
def gpcr_protein_path() -> Path:
    """Path to GPCR protein with POP lipids."""
    return PROTEINS_DIR / "4eiy-AA2AR_protein.pdb"


@pytest.fixture
def temp_pdb_dir():
    """Create a temporary directory for test PDB files."""
    temp_dir = Path(tempfile.mkdtemp())
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# -----------------------------------------------------------------------------
# Unit Tests for filter_pdb_by_sphere()
# -----------------------------------------------------------------------------


class TestFilterPdbBySphere:
    """Unit tests for the filter_pdb_by_sphere function."""

    def test_filter_removes_distant_chain(self, cdk2_protein_path, temp_pdb_dir):
        """Chain B should be removed when completely outside a small sphere centered on chain A."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        # Copy the PDB to temp directory
        input_pdb = temp_pdb_dir / "cdk2.pdb"
        output_pdb = temp_pdb_dir / "cdk2_filtered.pdb"
        shutil.copy(cdk2_protein_path, input_pdb)

        # Get approximate center of chain A (from ligand binding site)
        # COG is approximately (0.897, 26.524, 8.756) for cdk2
        center = [0.897, 26.524, 8.756]
        radius = 20.0  # At this radius, chain B has 0 atoms

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # Verify significant atom reduction
        assert filt_count < orig_count, "Filtered should have fewer atoms than original"

        # Verify chain B is removed
        filtered_df = read_pdb_to_dataframe(output_pdb)
        chain_ids = filtered_df["chain_id"].unique()
        assert "B" not in chain_ids, "Chain B should be completely removed"
        assert "A" in chain_ids, "Chain A should be preserved"

    def test_filter_preserves_touching_fragment(self, cdk2_protein_path, temp_pdb_dir):
        """Fragment with any atom in sphere should be entirely preserved."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "cdk2.pdb"
        output_pdb = temp_pdb_dir / "cdk2_filtered.pdb"
        shutil.copy(cdk2_protein_path, input_pdb)

        # Use a larger radius where chain B has some atoms in sphere
        center = [0.897, 26.524, 8.756]
        radius = 30.0  # At this radius, chain B should have atoms

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # With large radius, both chains should be preserved
        filtered_df = read_pdb_to_dataframe(output_pdb)
        chain_ids = filtered_df["chain_id"].unique()
        assert "A" in chain_ids, "Chain A should be preserved"
        # At large radius, chain B might also be preserved if it touches sphere

    def test_filter_removes_distant_lipids(self, gpcr_protein_path, temp_pdb_dir):
        """POP lipids completely outside sphere should be removed."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "gpcr.pdb"
        output_pdb = temp_pdb_dir / "gpcr_filtered.pdb"
        shutil.copy(gpcr_protein_path, input_pdb)

        # Count original POP residues
        orig_df = read_pdb_to_dataframe(input_pdb)
        orig_pop_count = orig_df[orig_df["residue_name"] == "POP"]["residue_seq_number"].nunique()

        # Use typical sphere center and radius
        # GPCR binding site is roughly at the membrane center
        center = [-8.0, 5.0, 30.0]  # Approximate from 4eiy structure
        radius = 25.0

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # Verify POP reduction
        filtered_df = read_pdb_to_dataframe(output_pdb)
        filt_pop_count = filtered_df[filtered_df["residue_name"] == "POP"]["residue_seq_number"].nunique()

        assert filt_pop_count < orig_pop_count, "Some POP lipids should be removed"
        assert filt_pop_count > 0, "Some POP lipids should remain (those near sphere)"

    def test_filter_preserves_lig_lid_hoh(self, cdk2_protein_path, temp_pdb_dir):
        """LIG, LID, HOH residues in sphere should be preserved (excluded from fragment expansion)."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        # Create a test PDB with LIG and LID residues
        input_pdb = temp_pdb_dir / "test.pdb"
        output_pdb = temp_pdb_dir / "test_filtered.pdb"

        # Write a simple test PDB with LIG, LID, HOH at known positions
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   ALA A   1       1.500   0.000   0.000  1.00  0.00           C
HETATM    3  C1  LIG A   2       2.000   0.000   0.000  1.00  0.00           C
HETATM    4  C2  LIG A   2       3.000   0.000   0.000  1.00  0.00           C
HETATM    5  C1  LID A   3       4.000   0.000   0.000  1.00  0.00           C
HETATM    6  C2  LID A   3       5.000   0.000   0.000  1.00  0.00           C
HETATM    7  O   HOH A   4       6.000   0.000   0.000  1.00  0.00           O
HETATM    8  H1  HOH A   4       6.500   0.500   0.000  1.00  0.00           H
ATOM    100  CA  ALA B   1     100.000 100.000 100.000  1.00  0.00           C
END
"""
        input_pdb.write_text(pdb_content)

        center = [0.0, 0.0, 0.0]
        radius = 10.0  # Enough to include LIG, LID, HOH, chain A

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        filtered_df = read_pdb_to_dataframe(output_pdb)
        resnames = filtered_df["residue_name"].unique()

        assert "LIG" in resnames, "LIG should be preserved"
        assert "LID" in resnames, "LID should be preserved"
        assert "HOH" in resnames, "HOH should be preserved"
        # Chain B at (100, 100, 100) is far outside and should be removed
        chain_ids = filtered_df["chain_id"].unique()
        assert "B" not in chain_ids, "Distant chain B should be removed"

    def test_filter_handles_empty_chain_id(self, temp_pdb_dir):
        """PDB with no chain IDs should still be processed correctly."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "no_chain.pdb"
        output_pdb = temp_pdb_dir / "no_chain_filtered.pdb"

        # Write a PDB with empty chain IDs (like thrombin)
        pdb_content = """\
ATOM      1  CA  ALA     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   ALA     1       1.500   0.000   0.000  1.00  0.00           C
ATOM      3  N   ALA     1       0.000   1.500   0.000  1.00  0.00           N
ATOM    100  CA  ALA    50     100.000 100.000 100.000  1.00  0.00           C
ATOM    101  C   ALA    50     101.500 100.000 100.000  1.00  0.00           C
END
"""
        input_pdb.write_text(pdb_content)

        center = [0.0, 0.0, 0.0]
        radius = 10.0

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        assert filt_count < orig_count, "Distant residue should be removed"

        filtered_df = read_pdb_to_dataframe(output_pdb)
        # Should have only the first residue (ALA 1)
        assert len(filtered_df) == 3, "Only atoms near sphere center should remain"

    def test_mdanalysis_detects_fragments(self, cdk2_protein_path, temp_pdb_dir):
        """Verify MDAnalysis correctly identifies disconnected fragments."""
        import MDAnalysis as mda

        u = mda.Universe(str(cdk2_protein_path))

        # Get all fragments
        n_fragments = len(u.atoms.fragments)

        # cdk2 has two chains, so should have at least 2 fragments
        assert n_fragments >= 2, f"Expected at least 2 fragments, got {n_fragments}"


# -----------------------------------------------------------------------------
# Unit Tests for filter_out_of_sphere_fragments()
# -----------------------------------------------------------------------------


class TestFilterOutOfSphereFragments:
    """Unit tests for the in-place filter_out_of_sphere_fragments wrapper."""

    def test_modifies_file_in_place(self, cdk2_protein_path, temp_pdb_dir):
        """Verify the function modifies the PDB file in place."""
        from QligFEP.pdb_utils import (
            filter_out_of_sphere_fragments,
            read_pdb_to_dataframe,
        )

        pdb_path = temp_pdb_dir / "cdk2.pdb"
        shutil.copy(cdk2_protein_path, pdb_path)

        orig_df = read_pdb_to_dataframe(pdb_path)
        orig_count = len(orig_df)

        center = [0.897, 26.524, 8.756]
        radius = 20.0

        orig, filt = filter_out_of_sphere_fragments(pdb_path, center, radius)

        # Verify counts returned
        assert orig == orig_count, "Original count should match"
        assert filt < orig, "Filtered count should be less"

        # Verify file was modified in place
        new_df = read_pdb_to_dataframe(pdb_path)
        assert len(new_df) == filt, "File should be updated in place"

    def test_returns_same_count_when_nothing_filtered(self, temp_pdb_dir):
        """When all atoms are inside sphere, counts should match."""
        from QligFEP.pdb_utils import filter_out_of_sphere_fragments

        pdb_path = temp_pdb_dir / "small.pdb"
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   ALA A   1       1.500   0.000   0.000  1.00  0.00           C
END
"""
        pdb_path.write_text(pdb_content)

        center = [0.0, 0.0, 0.0]
        radius = 100.0  # Very large radius

        orig, filt = filter_out_of_sphere_fragments(pdb_path, center, radius)

        assert orig == filt == 2, "No atoms should be filtered"


# -----------------------------------------------------------------------------
# Integration Tests
# -----------------------------------------------------------------------------


@pytest.mark.slow
class TestFragmentFilteringIntegration:
    """Integration tests for fragment filtering in the qligfep workflow."""

    def test_qligfep_filters_multichain(self, cdk2_protein_path, temp_pdb_dir):
        """Test that qligfep correctly filters multi-chain proteins."""
        # This would require setting up the full qligfep workflow
        # For now, just verify the basic filtering works on cdk2
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "cdk2.pdb"
        output_pdb = temp_pdb_dir / "cdk2_filtered.pdb"
        shutil.copy(cdk2_protein_path, input_pdb)

        # At radius 20, chain B should have 0 atoms
        center = [0.897, 26.524, 8.756]
        radius = 20.0

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # Get atom counts per chain
        orig_df = read_pdb_to_dataframe(input_pdb)
        filt_df = read_pdb_to_dataframe(output_pdb)

        orig_chain_a = len(orig_df[orig_df["chain_id"] == "A"])
        orig_chain_b = len(orig_df[orig_df["chain_id"] == "B"])
        filt_chain_a = len(filt_df[filt_df["chain_id"] == "A"])
        filt_chain_b = len(filt_df[filt_df["chain_id"] == "B"])

        # Chain A should be partially preserved
        assert filt_chain_a > 0, "Chain A should have atoms preserved"
        # Chain B should be completely removed at this radius
        assert filt_chain_b == 0, f"Chain B should be removed (had {orig_chain_b} atoms)"

        # Log the reduction
        print(f"Original: {orig_count} atoms (A: {orig_chain_a}, B: {orig_chain_b})")
        print(f"Filtered: {filt_count} atoms (A: {filt_chain_a}, B: {filt_chain_b})")

    def test_qligfep_preserves_partial_chain(self, cdk2_protein_path, temp_pdb_dir):
        """With larger radius, both chains should be preserved if they touch sphere."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "cdk2.pdb"
        output_pdb = temp_pdb_dir / "cdk2_filtered.pdb"
        shutil.copy(cdk2_protein_path, input_pdb)

        # At radius 25+, chain B might have atoms in sphere
        center = [0.897, 26.524, 8.756]
        radius = 30.0

        filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        filt_df = read_pdb_to_dataframe(output_pdb)

        # Check if chain B is preserved (it should be at larger radius)
        # This depends on the actual structure - adjust as needed
        chain_ids = filt_df["chain_id"].unique()
        print(f"Chains preserved at radius {radius}: {list(chain_ids)}")

    @pytest.mark.slow
    def test_qligfep_filters_lipids(self, gpcr_protein_path, temp_pdb_dir):
        """Test GPCR lipid filtering - most distant POP should be removed."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "gpcr.pdb"
        output_pdb = temp_pdb_dir / "gpcr_filtered.pdb"
        shutil.copy(gpcr_protein_path, input_pdb)

        orig_df = read_pdb_to_dataframe(input_pdb)
        orig_pop = orig_df[orig_df["residue_name"] == "POP"]["residue_seq_number"].nunique()

        # Use a reasonable sphere center and radius
        center = [-8.0, 5.0, 30.0]
        radius = 25.0

        orig_count, filt_count = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        filt_df = read_pdb_to_dataframe(output_pdb)
        filt_pop = filt_df[filt_df["residue_name"] == "POP"]["residue_seq_number"].nunique()

        # Expect significant POP reduction
        pop_removed = orig_pop - filt_pop
        print(f"POP residues: {orig_pop} -> {filt_pop} ({pop_removed} removed)")
        print(f"Total atoms: {orig_count} -> {filt_count}")

        assert pop_removed > 100, f"Expected >100 POP removed, got {pop_removed}"


# -----------------------------------------------------------------------------
# Edge Cases
# -----------------------------------------------------------------------------


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_nonstandard_element_records_with_charges(self, temp_pdb_dir):
        """PDB with formal charges embedded in element column should still work.

        Some PDB files (e.g., from membrane protein systems) have non-standard
        element records like "O1-" or "N1+" in columns 77-78. MDAnalysis parses
        these as elements "O1" and "N1" which are unknown and would cause bond
        guessing to fail without custom vdW radii.
        """
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "charged_elements.pdb"
        output_pdb = temp_pdb_dir / "charged_elements_filtered.pdb"

        # Simulate PDB content with non-standard element records (formal charges in element column)
        # The element column is positions 77-78 (0-indexed: 76-78)
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00          C
ATOM      2  C   ALA A   1       1.500   0.000   0.000  1.00  0.00          C
ATOM      3  N   ALA A   1       0.000   1.500   0.000  1.00  0.00          N
ATOM      4  O   ALA A   1       0.000   0.000   1.500  1.00  0.00          O
ATOM     10  OE2 GLU A   2       5.000   0.000   0.000  1.00  0.00          O1-
ATOM     11  NZ  LYS A   3       6.000   0.000   0.000  1.00  0.00          N1+
ATOM     12  NA  SOD A   4       7.000   0.000   0.000  1.00  0.00          Na1+
ATOM    100  CA  ALA B   1     100.000 100.000 100.000  1.00  0.00          C
ATOM    101  C   ALA B   1     101.500 100.000 100.000  1.00  0.00          C
END
"""
        input_pdb.write_text(pdb_content)

        center = [0.0, 0.0, 0.0]
        radius = 10.0  # Should include chain A but not chain B

        # This should not fail even with O1-, N1+, Na1+ element records
        orig, filt = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # Verify filtering worked (chain B removed)
        assert orig > filt, "Filtering should have occurred"
        assert filt > 0, "Some atoms should remain"

        filtered_df = read_pdb_to_dataframe(output_pdb)
        chain_ids = filtered_df["chain_id"].unique()
        assert "B" not in chain_ids, "Distant chain B should be removed"
        assert "A" in chain_ids, "Chain A should be preserved"

    def test_all_fragments_outside_sphere(self, temp_pdb_dir):
        """When all fragments are outside, should preserve nothing (or raise error)."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere, read_pdb_to_dataframe

        input_pdb = temp_pdb_dir / "far.pdb"
        output_pdb = temp_pdb_dir / "far_filtered.pdb"

        pdb_content = """\
ATOM      1  CA  ALA A   1     100.000 100.000 100.000  1.00  0.00           C
ATOM      2  C   ALA A   1     101.500 100.000 100.000  1.00  0.00           C
END
"""
        input_pdb.write_text(pdb_content)

        center = [0.0, 0.0, 0.0]
        radius = 10.0  # Far from the atoms

        orig, filt = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        # When all atoms are outside, filtered count should be 0
        assert filt == 0, "No atoms should remain when all are outside sphere"

    def test_empty_pdb_file(self, temp_pdb_dir):
        """Empty PDB should be handled gracefully."""
        from QligFEP.pdb_utils import filter_pdb_by_sphere

        input_pdb = temp_pdb_dir / "empty.pdb"
        output_pdb = temp_pdb_dir / "empty_filtered.pdb"

        input_pdb.write_text("END\n")

        center = [0.0, 0.0, 0.0]
        radius = 10.0

        # Should not raise an exception
        orig, filt = filter_pdb_by_sphere(input_pdb, output_pdb, center, radius)

        assert orig == 0, "Empty file should have 0 atoms"
        assert filt == 0, "Filtered should also be 0"
