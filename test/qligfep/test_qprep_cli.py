"""Tests for the QligFEP.CLI.qprep_cli module."""

import argparse
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import numpy as np
import pandas as pd
import pytest

from QligFEP.CLI.qprep_cli import Neutralizer


def _make_atom(record_type, serial, atom_name, residue_name, chain, resnum, x, y, z, element="C"):
    """Build a single-row dict matching the read_pdb_to_dataframe column layout."""
    return {
        "record_type": record_type,
        "atom_serial_number": serial,
        "atom_name": atom_name,
        "alt_loc": "",
        "residue_name": residue_name,
        "chain_id": chain,
        "residue_seq_number": resnum,
        "insertion_code": "",
        "x": float(x),
        "y": float(y),
        "z": float(z),
        "occupancy": 1.0,
        "temp_factor": 0.0,
        "segment_id": "",
        "element_symbol": element,
        "charge": "",
    }


def _build_dna_residue(resname, chain, resnum, x, y, z):
    """Build a minimal DNA residue with P, OP1, OP2, C1', O5' atoms at (x, y, z)."""
    atoms = []
    serial_base = resnum * 100
    for i, (aname, elem) in enumerate([("P", "P"), ("OP1", "O"), ("OP2", "O"), ("O5'", "O"), ("C1'", "C")]):
        atoms.append(_make_atom("ATOM", serial_base + i, aname, resname, chain, resnum, x, y, z, elem))
    return atoms


def _build_protein_residue(resname, chain, resnum, x, y, z, key_atom="CA"):
    """Build a minimal protein residue with key sidechain atom at (x, y, z)."""
    atoms = []
    serial_base = resnum * 100
    # Add backbone atoms
    for i, (aname, elem) in enumerate([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]):
        atoms.append(_make_atom("ATOM", serial_base + i, aname, resname, chain, resnum, x, y, z, elem))
    # Add the sidechain key atom if different from backbone
    key_atom_map = {"GLU": "CD", "ASP": "CG", "ARG": "CZ", "LYS": "NZ", "HIP": "CG"}
    if resname in key_atom_map:
        ka = key_atom_map[resname]
        atoms.append(_make_atom("ATOM", serial_base + 10, ka, resname, chain, resnum, x, y, z, ka[0]))
    return atoms


class TestNeutralizerDNA:
    """Tests for DNA nucleotide handling in the Neutralizer.

    DNA nucleotides whose C1' atom is outside the sphere radius are removed.
    DNA in the restrained shell (between rest_bound and radius) is kept with
    native charges, consistent with protein handling. The last inside residue
    adjacent to the removed region gets a 5' terminal form.
    """

    def _make_neutralizer(self, center=(0, 0, 0), radius=25.0, offset=3.0):
        return Neutralizer(center, radius, offset)

    def test_outside_dna_removed(self):
        """DA nucleotide fully outside sphere should be removed."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0))   # inside
        rows.extend(_build_dna_residue("DA", "E", 2, 30.0, 0.0, 0.0))  # outside
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert len(result_df[result_df["residue_seq_number"] == 2]) == 0

    def test_5prime_cap_on_inside_boundary(self):
        """First inside residue with removed upstream DNA gets 5' terminal form."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 30.0, 0.0, 0.0))  # outside, removed
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))   # inside, gets 5' cap
        rows.extend(_build_dna_residue("DC", "E", 3, 5.0, 0.0, 0.0))   # inside, unchanged
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        # Res 1: removed
        assert len(result_df[result_df["residue_seq_number"] == 1]) == 0
        # Res 2: 5' terminal (P removed)
        res2 = result_df[result_df["residue_seq_number"] == 2]
        assert (res2["residue_name"] == "DG5").all()
        assert "P" not in res2["atom_name"].values
        # Res 3: unchanged
        assert (result_df[result_df["residue_seq_number"] == 3]["residue_name"] == "DC").all()

    def test_dna_inside_sphere_unchanged(self):
        """DA nucleotide inside sphere should remain DA with all atoms."""
        rows = _build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df["residue_name"] == "DA").all()
        assert "P" in result_df["atom_name"].values

    def test_multiple_outside_removed_single_cap(self):
        """All outside DNA removed; only one 5' cap on the inside boundary."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 40.0, 0.0, 0.0))  # outside
        rows.extend(_build_dna_residue("DG", "E", 2, 35.0, 0.0, 0.0))  # outside
        rows.extend(_build_dna_residue("DC", "E", 3, 30.0, 0.0, 0.0))  # outside
        rows.extend(_build_dna_residue("DT", "E", 4, 5.0, 0.0, 0.0))   # inside, gets cap
        rows.extend(_build_dna_residue("DA", "E", 5, 5.0, 0.0, 0.0))   # inside
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        for rn in [1, 2, 3]:
            assert len(result_df[result_df["residue_seq_number"] == rn]) == 0
        res4 = result_df[result_df["residue_seq_number"] == 4]
        assert (res4["residue_name"] == "DT5").all()
        assert "P" not in res4["atom_name"].values
        assert (result_df[result_df["residue_seq_number"] == 5]["residue_name"] == "DA").all()

    def test_dna_with_c1_inside_sphere_unchanged(self):
        """DNA with C1' inside sphere radius should keep its charged form, even if other atoms are outside."""
        rows = []
        serial = 100
        # P and OP1/OP2 outside sphere, but C1' at 24 Å (inside radius=25)
        positions = [(26.0, 0, 0), (27.0, 0, 0), (27.0, 0, 0), (26.0, 0, 0), (24.0, 0, 0)]
        for (aname, elem), (px, py, pz) in zip(
            [("P", "P"), ("OP1", "O"), ("OP2", "O"), ("O5'", "O"), ("C1'", "C")],
            positions,
        ):
            rows.append(_make_atom("ATOM", serial, aname, "DA", "E", 1, px, py, pz, elem))
            serial += 1
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df["residue_name"] == "DA").all()

    def test_mixed_protein_and_dna(self):
        """Protein residues neutralized by rest_bound; DNA removed by same boundary."""
        rows = []
        rows.extend(_build_protein_residue("GLU", "A", 1, 23.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DA", "E", 10, 5.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DA", "E", 11, 30.0, 0.0, 0.0))
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df[result_df["residue_seq_number"] == 1]["residue_name"] == "GLH").all()
        assert (result_df[result_df["residue_seq_number"] == 10]["residue_name"] == "DA").all()
        assert len(result_df[result_df["residue_seq_number"] == 11]) == 0

    def test_downstream_removal_no_5prime_cap(self):
        """When removed DNA is only downstream (3' side), no 5' cap is needed."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0))   # inside
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))   # inside
        rows.extend(_build_dna_residue("DC", "E", 3, 30.0, 0.0, 0.0))  # outside
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        # Res 3: removed
        assert len(result_df[result_df["residue_seq_number"] == 3]) == 0
        # Res 1-2: inside, unchanged (no terminal conversion)
        assert (result_df[result_df["residue_seq_number"] == 1]["residue_name"] == "DA").all()
        assert (result_df[result_df["residue_seq_number"] == 2]["residue_name"] == "DG").all()


def _build_nterm_residue(resname, chain, resnum, x, y, z):
    """Build a minimal N-terminal residue with N, H1, H2, H3, CA atoms."""
    atoms = []
    serial_base = resnum * 100
    for i, (aname, elem) in enumerate([("N", "N"), ("H1", "H"), ("H2", "H"), ("H3", "H"), ("CA", "C")]):
        atoms.append(_make_atom("ATOM", serial_base + i, aname, resname, chain, resnum, x, y, z, elem))
    return atoms


class TestNeutralizerNTerminals:
    """Tests for N-terminal residue neutralization in the Neutralizer.

    N-terminal residues (NILE, NGLY, etc.) carry +1 charge from the protonated
    NH3+ group. When outside the neutralization boundary, they are converted to
    their internal forms (ILE, GLY, etc.) with H1→H rename and H2, H3 removal.
    """

    def _make_neutralizer(self, center=(0, 0, 0), radius=25.0, offset=3.0):
        return Neutralizer(center, radius, offset)

    def test_nterm_outside_boundary_neutralized(self):
        """NILE outside rest_bound should become ILE with H2/H3 removed."""
        rows = _build_nterm_residue("NILE", "A", 1, 30.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "ILE").all()
        assert "H" in res["atom_name"].values  # H1 renamed to H
        assert "H2" not in res["atom_name"].values
        assert "H3" not in res["atom_name"].values
        assert len(res) == 3  # N, H, CA (H2 and H3 removed)

    def test_nterm_inside_boundary_unchanged(self):
        """NILE inside rest_bound should remain NILE with all atoms."""
        rows = _build_nterm_residue("NILE", "A", 1, 5.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df["residue_name"] == "NILE").all()
        assert "H1" in result_df["atom_name"].values
        assert "H2" in result_df["atom_name"].values
        assert "H3" in result_df["atom_name"].values

    def test_ngly_outside_boundary_neutralized(self):
        """NGLY outside rest_bound should become GLY."""
        rows = _build_nterm_residue("NGLY", "B", 1, 30.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "GLY").all()
        assert "H2" not in res["atom_name"].values
        assert "H3" not in res["atom_name"].values


def _build_cterm_residue(resname, chain, resnum, x, y, z):
    """Build a minimal C-terminal residue with N, CA, C, O, OXT atoms."""
    atoms = []
    serial_base = resnum * 100
    for i, (aname, elem) in enumerate([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("OXT", "O")]):
        atoms.append(_make_atom("ATOM", serial_base + i, aname, resname, chain, resnum, x, y, z, elem))
    return atoms


class TestNeutralizerCTerminals:
    """Tests for C-terminal residue neutralization in the Neutralizer.

    C-terminal residues (CALA, CGLY, etc.) carry -1 charge from the deprotonated
    COO- group. When outside the neutralization boundary, they are converted to
    their internal forms (ALA, GLY, etc.) with OXT removal.
    """

    def _make_neutralizer(self, center=(0, 0, 0), radius=25.0, offset=3.0):
        return Neutralizer(center, radius, offset)

    def test_cterm_outside_boundary_neutralized(self):
        """CALA outside rest_bound should become ALA with OXT removed."""
        rows = _build_cterm_residue("CALA", "A", 1, 30.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "ALA").all()
        assert "OXT" not in res["atom_name"].values
        assert "O" in res["atom_name"].values
        assert len(res) == 4  # N, CA, C, O (OXT removed)

    def test_cterm_inside_boundary_unchanged(self):
        """CALA inside rest_bound should remain CALA with all atoms."""
        rows = _build_cterm_residue("CALA", "A", 1, 5.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df["residue_name"] == "CALA").all()
        assert "OXT" in result_df["atom_name"].values

    def test_cgly_outside_boundary_neutralized(self):
        """CGLY outside rest_bound should become GLY."""
        rows = _build_cterm_residue("CGLY", "B", 1, 30.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "GLY").all()
        assert "OXT" not in res["atom_name"].values


class TestFragmentFilteringInQprep:
    """Tests for fragment filtering during protein preparation."""

    def create_mock_args(self, skip_fragment_filter=False):
        """Create args with COG at origin and radius 10 (chain B at 100,100,100 is outside)."""
        return argparse.Namespace(
            input_pdb_file="protein.pdb",
            FF="AMBER14sb",
            cog=["0.000", "0.000", "0.000"],
            sphereradius=10.0,
            cysbond="auto",
            solvent_pack=3.0,
            log_level="info",
            cofactors=None,
            skip_neutralization=True,
            neutralize_boundary_offset=3.0,
            salt_bridge_cutoff=4.0,
            skip_fragment_filter=skip_fragment_filter,
        )

    TWO_CHAIN_PDB = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM    100  N   ALA B   1     100.000 100.000 100.000  1.00  0.00           N\n"
        "ATOM    101  CA  ALA B   1     101.458 100.000 100.000  1.00  0.00           C\n"
        "ATOM    102  C   ALA B   1     102.009 101.420 100.000  1.00  0.00           C\n"
    )

    MINIMAL_COMPLEX_PDB = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  O   HOH W 100       5.000   5.000   5.000  1.00  0.00           O\n"
    )

    def _run_qprep_main(self, tmp_path, monkeypatch, args):
        """Run qprep_cli.main() with mocked external dependencies."""
        monkeypatch.chdir(tmp_path)

        def mock_run_qprep(*a, **kw):
            (tmp_path / "complexnotexcluded.pdb").write_text(self.MINIMAL_COMPLEX_PDB)
            return MagicMock(returncode=0)

        with patch("QligFEP.CLI.qprep_cli.run_qprep", side_effect=mock_run_qprep):
            with patch(
                "QligFEP.CLI.qprep_cli.get_force_field_paths",
                return_value=("/path/to/ff.lib", "/path/to/ff.prm"),
            ):
                with patch("QligFEP.CLI.qprep_cli.CONFIGS", {"QPREP": "/path/to/qprep"}):
                    with patch("QligFEP.CLI.qprep_cli.handle_cysbonds", return_value=""):
                        from QligFEP.CLI.qprep_cli import main

                        main(args)

    def test_distant_chain_removed_during_protein_prep(self, tmp_path, monkeypatch):
        """Distant chain B should be removed from the protein PDB before qprep runs."""
        from QligFEP.pdb_utils import read_pdb_to_dataframe

        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text(self.TWO_CHAIN_PDB)

        args = self.create_mock_args(skip_fragment_filter=False)
        self._run_qprep_main(tmp_path, monkeypatch, args)

        # The processed PDB (protein.pdb, modified in-place) should have chain B removed
        df = read_pdb_to_dataframe(pdb_file)
        assert "B" not in df["chain_id"].values, "Chain B should be filtered out"
        assert "A" in df["chain_id"].values, "Chain A should be preserved"

    def test_skip_fragment_filter_preserves_distant_chain(self, tmp_path, monkeypatch):
        """When --skip-fragment-filter is set, distant chains should be preserved."""
        from QligFEP.pdb_utils import read_pdb_to_dataframe

        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text(self.TWO_CHAIN_PDB)

        args = self.create_mock_args(skip_fragment_filter=True)
        self._run_qprep_main(tmp_path, monkeypatch, args)

        # Chain B should still be present when filtering is skipped
        df = read_pdb_to_dataframe(pdb_file)
        assert "B" in df["chain_id"].values, "Chain B should be preserved when filtering is skipped"
