"""Tests for DNA hydrogen renaming and ion renaming in pdb2amber."""

import tempfile
from pathlib import Path

import pytest

from QligFEP.CLI.pdb_to_amber import fix_pdb
from QligFEP.pdb_utils import read_pdb_to_dataframe


def _write_pdb(tmp_path, lines):
    """Write PDB lines to a temp file and return the path."""
    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text("".join(lines))
    return pdb_file


def _get_atom_names(pdb_path, residue_name=None):
    """Read PDB and return atom names, optionally filtered by residue."""
    df = read_pdb_to_dataframe(pdb_path)
    if residue_name:
        df = df[df["residue_name"] == residue_name]
    return list(df["atom_name"])


# Minimal DNA nucleotide PDB lines with Maestro-style hydrogen naming
DA_MAESTRO = [
    "ATOM      1  P    DA E   5     200.000 200.000 200.000  1.00  0.00           P\n",
    "ATOM      2  OP1  DA E   5     200.000 200.000 200.000  1.00  0.00           O\n",
    "ATOM      3  OP2  DA E   5     200.000 200.000 200.000  1.00  0.00           O\n",
    "ATOM      4  O5'  DA E   5     200.000 200.000 200.000  1.00  0.00           O\n",
    "ATOM      5  C5'  DA E   5     200.000 200.000 200.000  1.00  0.00           C\n",
    "ATOM      6  C1'  DA E   5     200.000 200.000 200.000  1.00  0.00           C\n",
    "ATOM      7 H5'1  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      8 H5'2  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      9 HC4'  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     10 HC3'  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     11 H2'1  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     12 H2'2  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     13 HC1'  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     14  HC8  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     15 H6_1  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     16 H6_2  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     17  HC2  DA E   5     200.000 200.000 200.000  1.00  0.00           H\n",
]

DT_MAESTRO = [
    "ATOM      1  P    DT E  10     200.000 200.000 200.000  1.00  0.00           P\n",
    "ATOM      2  C1'  DT E  10     200.000 200.000 200.000  1.00  0.00           C\n",
    "ATOM      3 H5'1  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      4 H5'2  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      5 HC4'  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      6 HC3'  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      7 H2'1  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      8 H2'2  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      9 HC1'  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     10  HN3  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     11 H7_1  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     12 H7_2  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     13 H7_3  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     14  HC6  DT E  10     200.000 200.000 200.000  1.00  0.00           H\n",
]

DC_MAESTRO = [
    "ATOM      1  P    DC E  15     200.000 200.000 200.000  1.00  0.00           P\n",
    "ATOM      2  C1'  DC E  15     200.000 200.000 200.000  1.00  0.00           C\n",
    "ATOM      3 H5'1  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      4 H5'2  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      5 HC4'  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      6 HC3'  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      7 H2'1  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      8 H2'2  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      9 HC1'  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     10 H4_1  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     11 H4_2  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     12  HC5  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     13  HC6  DC E  15     200.000 200.000 200.000  1.00  0.00           H\n",
]

DG_MAESTRO = [
    "ATOM      1  P    DG E  20     200.000 200.000 200.000  1.00  0.00           P\n",
    "ATOM      2  C1'  DG E  20     200.000 200.000 200.000  1.00  0.00           C\n",
    "ATOM      3 H5'1  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      4 H5'2  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      5 HC4'  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      6 HC3'  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      7 H2'1  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      8 H2'2  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM      9 HC1'  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     10  HC8  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     11  HN1  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     12 H2_1  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
    "ATOM     13 H2_2  DG E  20     200.000 200.000 200.000  1.00  0.00           H\n",
]

MG_PDB = [
    "HETATM  100 MG    MG B 901     222.053 200.827 188.112  1.00 27.46          Mg2+\n",
]


def _line(serial, name, resname, chain, resnum, occ=1.00, alt=" ", x=0.0, y=0.0, z=0.0, elem=""):
    """Build a single PDB ATOM line with exact column placement."""
    name4 = name.ljust(4) if len(name) >= 4 else (" " + name).ljust(4)
    res4 = resname.ljust(4)
    return (
        f"ATOM  {serial:>5} {name4}{alt}{res4}{chain}{resnum:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}  0.00          {elem:>2}\n"
    )


class TestDNAHydrogenRenaming:
    """Tests for DNA hydrogen renaming in pdb2amber."""

    def test_da_hydrogen_renaming(self, tmp_path):
        pdb_file = _write_pdb(tmp_path, DA_MAESTRO)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file)
        # Sugar hydrogens
        assert "H5'" in atoms
        assert "H5''" in atoms
        assert "H4'" in atoms
        assert "H3'" in atoms
        assert "H2'" in atoms
        assert "H2''" in atoms
        assert "H1'" in atoms
        # Base hydrogens
        assert "H8" in atoms
        assert "H61" in atoms
        assert "H62" in atoms
        assert "H2" in atoms
        # Old names should be gone
        for old in ["H5'1", "H5'2", "HC4'", "HC3'", "H2'1", "H2'2", "HC1'", "HC8", "H6_1", "H6_2", "HC2"]:
            assert old not in atoms, f"{old} should have been renamed"

    def test_dt_hydrogen_renaming(self, tmp_path):
        pdb_file = _write_pdb(tmp_path, DT_MAESTRO)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file)
        assert "H3" in atoms
        assert "H71" in atoms
        assert "H72" in atoms
        assert "H73" in atoms
        assert "H6" in atoms
        for old in ["HN3", "H7_1", "H7_2", "H7_3", "HC6"]:
            assert old not in atoms, f"{old} should have been renamed"

    def test_dc_hydrogen_renaming(self, tmp_path):
        pdb_file = _write_pdb(tmp_path, DC_MAESTRO)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file)
        assert "H41" in atoms
        assert "H42" in atoms
        assert "H5" in atoms
        assert "H6" in atoms
        for old in ["H4_1", "H4_2", "HC5", "HC6"]:
            assert old not in atoms, f"{old} should have been renamed"

    def test_dg_hydrogen_renaming(self, tmp_path):
        pdb_file = _write_pdb(tmp_path, DG_MAESTRO)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file)
        assert "H8" in atoms
        assert "H1" in atoms
        assert "H21" in atoms
        assert "H22" in atoms
        for old in ["HC8", "HN1", "H2_1", "H2_2"]:
            assert old not in atoms, f"{old} should have been renamed"


class TestNTerminalDetection:
    """Tests for N-terminal detection with Maestro-style atom naming."""

    def test_n_terminal_detected_by_h1_h2(self, tmp_path):
        """N-terminal ILE with H1+H2 (no H3) should be labeled NILE."""
        lines = [
            "ATOM      1  N   ILE A   8       0.000   0.000   0.000  1.00  0.00           N\n",
            "ATOM      2  CA  ILE A   8       1.458   0.000   0.000  1.00  0.00           C\n",
            "ATOM      3  C   ILE A   8       2.009   1.420   0.000  1.00  0.00           C\n",
            "ATOM      4  O   ILE A   8       1.246   2.390   0.000  1.00  0.00           O\n",
            "ATOM      5  CB  ILE A   8       2.000   0.000   1.500  1.00  0.00           C\n",
            "ATOM      6  H1  ILE A   8      -0.500   0.800   0.000  1.00  0.00           H\n",
            "ATOM      7  H2  ILE A   8      -0.500  -0.800   0.000  1.00  0.00           H\n",
            "ATOM      8  HA  ILE A   8       1.800   0.000  -1.000  1.00  0.00           H\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "NILE").all()

    def test_n_terminal_gly_duplicate_h(self, tmp_path):
        """N-terminal GLY with duplicate H atoms should be labeled NGLY with H1,H2."""
        lines = [
            "ATOM      1  N   GLY B 402       0.000   0.000   0.000  1.00  0.00           N\n",
            "ATOM      2  CA  GLY B 402       1.458   0.000   0.000  1.00  0.00           C\n",
            "ATOM      3  C   GLY B 402       2.009   1.420   0.000  1.00  0.00           C\n",
            "ATOM      4  O   GLY B 402       1.246   2.390   0.000  1.00  0.00           O\n",
            "ATOM      5  H   GLY B 402      -0.500   0.800   0.000  1.00  0.00           H\n",
            "ATOM      6  H   GLY B 402      -0.500  -0.800   0.000  1.00  0.00           H\n",
            "ATOM      7  HA3 GLY B 402       1.800   0.000  -1.000  1.00  0.00           H\n",
            "ATOM      8  HA2 GLY B 402       1.800   0.000   1.000  1.00  0.00           H\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "NGLY").all()
        atom_names = list(df["atom_name"])
        assert "H1" in atom_names
        assert "H2" in atom_names


class TestCTerminalDetection:
    """Tests for C-terminal detection with Maestro-style atom naming."""

    def test_hxt_removed_residue_stays_standard(self, tmp_path):
        """C-terminal ASN with HXT: HXT removed, residue stays ASN (not CASN).

        Maestro places HXT but not OXT. Since qprep can't add missing heavy atoms,
        we remove HXT and keep the standard residue name to avoid OXT errors.
        """
        lines = [
            "ATOM      1  N   ASN A 524       0.000   0.000   0.000  1.00  0.00           N\n",
            "ATOM      2  CA  ASN A 524       1.458   0.000   0.000  1.00  0.00           C\n",
            "ATOM      3  C   ASN A 524       2.009   1.420   0.000  1.00  0.00           C\n",
            "ATOM      4  O   ASN A 524       1.246   2.390   0.000  1.00  0.00           O\n",
            "ATOM      5  H   ASN A 524      -0.500   0.800   0.000  1.00  0.00           H\n",
            "ATOM      6  HA  ASN A 524       1.800   0.000  -1.000  1.00  0.00           H\n",
            "ATOM      7  HXT ASN A 524       3.000   1.500   0.000  1.00  0.00           H\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "ASN").all()
        assert "HXT" not in list(df["atom_name"])


class TestMGRenaming:
    """Tests for MG → MAG renaming in pdb2amber."""

    def test_mg_renamed_to_mag(self, tmp_path):
        # Add a minimal protein residue so fix_pdb has something to process
        lines = [
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n",
        ] + MG_PDB
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "output.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        mg_rows = df[df["residue_name"].isin(["MG", "MAG"])]
        assert len(mg_rows) == 1
        assert mg_rows.iloc[0]["residue_name"] == "MAG"
        assert mg_rows.iloc[0]["atom_name"] == "MAG"


class TestNTerminalChargedResidues:
    """N-terminal residues with H1/H2/H3 must keep all three protons (not be corrupted)."""

    def _nterm_arg(self):
        # Free N-terminal ARG with NH3+ protons H1, H2, H3 (standard naming)
        return [
            _line(1, "N", "ARG", "A", 1, elem="N"),
            _line(2, "CA", "ARG", "A", 1, elem="C"),
            _line(3, "C", "ARG", "A", 1, elem="C"),
            _line(4, "O", "ARG", "A", 1, elem="O"),
            _line(5, "CB", "ARG", "A", 1, elem="C"),
            _line(6, "CG", "ARG", "A", 1, elem="C"),
            _line(7, "CD", "ARG", "A", 1, elem="C"),
            _line(8, "NE", "ARG", "A", 1, elem="N"),
            _line(9, "CZ", "ARG", "A", 1, elem="C"),
            _line(10, "NH1", "ARG", "A", 1, elem="N"),
            _line(11, "NH2", "ARG", "A", 1, elem="N"),
            _line(12, "H1", "ARG", "A", 1, elem="H"),
            _line(13, "H2", "ARG", "A", 1, elem="H"),
            _line(14, "H3", "ARG", "A", 1, elem="H"),
            _line(15, "HA", "ARG", "A", 1, elem="H"),
        ]

    def test_nterminal_arg_keeps_three_protons(self, tmp_path):
        pdb_file = _write_pdb(tmp_path, self._nterm_arg())
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "NARG").all()
        atoms = list(df["atom_name"])
        for h in ("H1", "H2", "H3"):
            assert h in atoms, f"N-terminal proton {h} was lost/renamed"
        assert atoms.count("H") == 0, "backbone amide H2 was wrongly renamed to H"

    def test_nterminal_gln_keeps_three_protons(self, tmp_path):
        lines = [
            _line(1, "N", "GLN", "A", 1, elem="N"),
            _line(2, "CA", "GLN", "A", 1, elem="C"),
            _line(3, "C", "GLN", "A", 1, elem="C"),
            _line(4, "O", "GLN", "A", 1, elem="O"),
            _line(5, "CB", "GLN", "A", 1, elem="C"),
            _line(6, "H1", "GLN", "A", 1, elem="H"),
            _line(7, "H2", "GLN", "A", 1, elem="H"),
            _line(8, "H3", "GLN", "A", 1, elem="H"),
            _line(9, "HA", "GLN", "A", 1, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "NGLN").all()
        atoms = list(df["atom_name"])
        for h in ("H1", "H2", "H3"):
            assert h in atoms
        assert atoms.count("H") == 0


class TestCappedBackboneH:
    """A residue following an ACE cap has a single backbone amide H misnamed 'H1' -> 'H'."""

    def test_lone_h1_renamed_to_h(self, tmp_path):
        # ACE cap then internal SER whose single amide proton is named H1
        lines = [
            _line(1, "CH3", "ACE", "A", 0, elem="C"),
            _line(2, "C", "ACE", "A", 0, elem="C"),
            _line(3, "O", "ACE", "A", 0, elem="O"),
            _line(4, "H1", "ACE", "A", 0, elem="H"),
            _line(5, "H2", "ACE", "A", 0, elem="H"),
            _line(6, "H3", "ACE", "A", 0, elem="H"),
            _line(7, "N", "SER", "A", 1, elem="N"),
            _line(8, "CA", "SER", "A", 1, elem="C"),
            _line(9, "C", "SER", "A", 1, elem="C"),
            _line(10, "O", "SER", "A", 1, elem="O"),
            _line(11, "CB", "SER", "A", 1, elem="C"),
            _line(12, "OG", "SER", "A", 1, elem="O"),
            _line(13, "H1", "SER", "A", 1, elem="H"),
            _line(14, "HA", "SER", "A", 1, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        ser = df[df["residue_name"] == "SER"]
        assert len(ser) == 8, "SER should stay an internal residue"
        atoms = list(ser["atom_name"])
        assert "H" in atoms, "lone backbone H1 should be renamed to H"
        assert "H1" not in atoms


class TestIonRenaming:
    """Monatomic ions must be renamed to AMBER14sb library names."""

    @pytest.mark.parametrize(
        "in_res,in_atom,out_name",
        [("ZN", "ZN", "ZIN"), ("NA", "NA", "SOD"), ("CL", "CL", "CHL")],
    )
    def test_ion_renamed(self, tmp_path, in_res, in_atom, out_name):
        lines = [
            _line(1, "N", "ALA", "A", 1, elem="N"),
            "HETATM  100 "
            + in_atom.ljust(4)
            + " "
            + in_res.ljust(4)
            + "B 901     222.053 200.827 188.112  1.00 27.46          XX\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        ion = df[df["residue_name"] == out_name]
        assert len(ion) == 1, f"{in_res} not renamed to {out_name}"
        assert ion.iloc[0]["atom_name"] == out_name

    def test_calcium_mixed_case_renamed(self, tmp_path):
        # factor_xa ships calcium as residue "Ca" / atom "CA1"; the lookup must be
        # case-insensitive and the atom name must follow the AMBER14sb convention.
        lines = [
            _line(1, "N", "ALA", "A", 1, elem="N"),
            "HETATM  100 CA1  Ca  X   1       8.846   2.038  34.536  1.00 62.09          Ca\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        ion = df[df["residue_name"] == "CAL"]
        assert len(ion) == 1, "calcium 'Ca' not renamed to CAL"
        assert ion.iloc[0]["atom_name"] == "CAL"


class TestNMECap:
    """NME C-terminal cap: CH3 -> C and methyl hydrogens -> H1/H2/H3."""

    def test_nme_methyl_renamed(self, tmp_path):
        lines = [
            _line(1, "N", "NME", "A", 224, elem="N"),
            _line(2, "CH3", "NME", "A", 224, elem="C"),
            _line(3, "H", "NME", "A", 224, elem="H"),
            "ATOM      4 1HH3 NME A 224     200.000 200.000 200.000  1.00  0.00           H\n",
            "ATOM      5 2HH3 NME A 224     200.000 200.000 200.000  1.00  0.00           H\n",
            "ATOM      6 3HH3 NME A 224     200.000 200.000 200.000  1.00  0.00           H\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file, residue_name="NME")
        assert set(atoms) == {"N", "C", "H", "H1", "H2", "H3"}, atoms

    def test_nme_ca_methyl_renamed(self, tmp_path):
        # NME cap using CA + 1HA/2HA/3HA methyl naming (common in benchmark structures)
        lines = [
            _line(1, "N", "NME", "A", 167, elem="N"),
            _line(2, "CA", "NME", "A", 167, elem="C"),
            _line(3, "H", "NME", "A", 167, elem="H"),
            "ATOM      4 1HA  NME A 167     200.000 200.000 200.000  1.00  0.00           H\n",
            "ATOM      5 2HA  NME A 167     200.000 200.000 200.000  1.00  0.00           H\n",
            "ATOM      6 3HA  NME A 167     200.000 200.000 200.000  1.00  0.00           H\n",
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = sorted(_get_atom_names(out_file, residue_name="NME"))
        assert atoms == ["C", "H", "H1", "H2", "H3", "N"], atoms


class TestAltloc:
    """Alternate conformations must be collapsed to the highest-occupancy altloc."""

    def test_altloc_collapsed_keeps_highest_occupancy(self, tmp_path):
        lines = [
            _line(1, "N", "CYS", "A", 93, elem="N"),
            _line(2, "CA", "CYS", "A", 93, elem="C"),
            _line(3, "C", "CYS", "A", 93, elem="C"),
            _line(4, "O", "CYS", "A", 93, elem="O"),
            _line(5, "CB", "CYS", "A", 93, occ=0.60, alt="A", x=1.0, elem="C"),
            _line(6, "CB", "CYS", "A", 93, occ=0.40, alt="B", x=2.0, elem="C"),
            _line(7, "SG", "CYS", "A", 93, occ=0.60, alt="A", x=1.0, elem="S"),
            _line(8, "SG", "CYS", "A", 93, occ=0.40, alt="B", x=2.0, elem="S"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        assert list(df["atom_name"]).count("CB") == 1, "duplicate altloc CB not collapsed"
        assert list(df["atom_name"]).count("SG") == 1
        cb = df[df["atom_name"] == "CB"].iloc[0]
        assert abs(cb["x"] - 1.0) < 1e-6, "should keep altloc A (occ 0.60)"
        assert (df["alt_loc"] == "").all(), "altloc indicator should be cleared"


class TestCTerminalOxygens:
    """C-terminal carboxylate oxygens named O1/O2 must map to O/OXT."""

    def test_o1_o2_renamed_and_cterminal(self, tmp_path):
        # C-terminal ASP with carboxylate oxygens named O1/O2 (no OXT)
        lines = [
            _line(1, "N", "ASP", "A", 298, elem="N"),
            _line(2, "CA", "ASP", "A", 298, elem="C"),
            _line(3, "C", "ASP", "A", 298, elem="C"),
            _line(4, "O2", "ASP", "A", 298, elem="O"),
            _line(5, "CB", "ASP", "A", 298, elem="C"),
            _line(6, "CG", "ASP", "A", 298, elem="C"),
            _line(7, "OD1", "ASP", "A", 298, elem="O"),
            _line(8, "OD2", "ASP", "A", 298, elem="O"),
            _line(9, "O1", "ASP", "A", 298, elem="O"),
            _line(10, "H", "ASP", "A", 298, elem="H"),
            _line(11, "HA", "ASP", "A", 298, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        df = read_pdb_to_dataframe(out_file)
        atoms = list(df["atom_name"])
        assert "O" in atoms and "OXT" in atoms
        assert "O1" not in atoms and "O2" not in atoms
        assert (df["residue_name"] == "CASP").all(), "should be detected as C-terminal"

    def test_existing_o_plus_protonated_o1_becomes_charged_cterminus(self, tmp_path):
        lines = [
            _line(1, "N", "VAL", "A", 633, y=4.0, elem="N"),
            _line(2, "CA", "VAL", "A", 633, y=3.0, elem="C"),
            _line(3, "C", "VAL", "A", 633, elem="C"),
            _line(4, "O", "VAL", "A", 633, x=1.2, elem="O"),
            _line(5, "O1", "VAL", "A", 633, x=-1.2, elem="O"),
            _line(6, "H", "VAL", "A", 633, y=5.0, elem="H"),
            _line(7, "H1", "VAL", "A", 633, x=-2.15, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        df = read_pdb_to_dataframe(out_file)
        atoms = set(df["atom_name"])
        assert {"O", "OXT", "H"}.issubset(atoms)
        assert "O1" not in atoms and "H1" not in atoms
        assert (df["residue_name"] == "CVAL").all()

    def test_o1_with_existing_oxt_renamed_to_o(self, tmp_path):
        # C-terminal ILE that already has OXT plus a carbonyl oxygen named O1
        lines = [
            _line(1, "N", "ILE", "A", 615, elem="N"),
            _line(2, "CA", "ILE", "A", 615, elem="C"),
            _line(3, "C", "ILE", "A", 615, elem="C"),
            _line(4, "OXT", "ILE", "A", 615, elem="O"),
            _line(5, "CB", "ILE", "A", 615, elem="C"),
            _line(6, "O1", "ILE", "A", 615, elem="O"),
            _line(7, "H", "ILE", "A", 615, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file)
        assert "O" in atoms and "OXT" in atoms and "O1" not in atoms


class TestASNSidechainAndBackbone:
    """ASN sidechain amide HD1/HD2 -> HD21/HD22 and lone backbone H11 -> H."""

    def test_asn_hd1_hd2_renamed(self, tmp_path):
        lines = [
            _line(1, "N", "ASN", "A", 109, elem="N"),
            _line(2, "CA", "ASN", "A", 109, elem="C"),
            _line(3, "C", "ASN", "A", 109, elem="C"),
            _line(4, "O", "ASN", "A", 109, elem="O"),
            _line(5, "CB", "ASN", "A", 109, elem="C"),
            _line(6, "CG", "ASN", "A", 109, elem="C"),
            _line(7, "OD1", "ASN", "A", 109, elem="O"),
            _line(8, "ND2", "ASN", "A", 109, elem="N"),
            _line(9, "H", "ASN", "A", 109, elem="H"),
            _line(10, "HD1", "ASN", "A", 109, elem="H"),
            _line(11, "HD2", "ASN", "A", 109, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file, residue_name="ASN")
        assert "HD21" in atoms and "HD22" in atoms
        assert "HD1" not in atoms and "HD2" not in atoms

    def test_asn_lone_h11_backbone_renamed(self, tmp_path):
        # ASN after an ACE cap whose lone backbone amide proton is named H11
        lines = [
            _line(1, "N", "ASN", "A", 917, elem="N"),
            _line(2, "CA", "ASN", "A", 917, elem="C"),
            _line(3, "C", "ASN", "A", 917, elem="C"),
            _line(4, "O", "ASN", "A", 917, elem="O"),
            _line(5, "CB", "ASN", "A", 917, elem="C"),
            _line(6, "CG", "ASN", "A", 917, elem="C"),
            _line(7, "OD1", "ASN", "A", 917, elem="O"),
            _line(8, "ND2", "ASN", "A", 917, elem="N"),
            _line(9, "H11", "ASN", "A", 917, elem="H"),
            _line(10, "HA", "ASN", "A", 917, elem="H"),
            _line(11, "HD21", "ASN", "A", 917, elem="H"),
            _line(12, "HD22", "ASN", "A", 917, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"
        fix_pdb(pdb_file, out_name=out_file)
        atoms = _get_atom_names(out_file, residue_name="ASN")
        assert "H" in atoms and "H11" not in atoms


def _atom_x(pdb_path, atom_name):
    df = read_pdb_to_dataframe(pdb_path)
    return float(df.loc[df["atom_name"] == atom_name, "x"].iloc[0])


class TestProtonatedCarboxylateNormalization:
    """Acidic proton geometry and AMBER14sb oxygen identities must agree."""

    @staticmethod
    def _acid_lines(resname, atom_names, hydrogen_name, hydrogen_x):
        carbon, oxygen1, oxygen2 = atom_names
        return [
            _line(1, carbon, resname, "A", 167, x=-1.2, elem="C"),
            _line(2, oxygen1, resname, "A", 167, x=0.0, elem="O"),
            _line(3, oxygen2, resname, "A", 167, x=2.4, elem="O"),
            _line(4, hydrogen_name, resname, "A", 167, x=hydrogen_x, elem="H"),
        ]

    def test_glu_proton_on_oe1_swaps_oxygens_and_becomes_glh(self, tmp_path):
        lines = self._acid_lines("GLU", ("CD", "OE1", "OE2"), "HE1", 0.96)
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "GLH").all()
        assert "HE2" in set(df["atom_name"])
        assert "HE1" not in set(df["atom_name"])
        assert _atom_x(out_file, "OE2") == pytest.approx(0.0)
        assert _atom_x(out_file, "OE1") == pytest.approx(2.4)

    def test_correct_glh_oe2_geometry_is_unchanged(self, tmp_path):
        lines = self._acid_lines("GLH", ("CD", "OE1", "OE2"), "HE2", 3.36)
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        assert _atom_x(out_file, "OE1") == pytest.approx(0.0)
        assert _atom_x(out_file, "OE2") == pytest.approx(2.4)

    def test_asp_hd2_on_od1_swaps_oxygens_and_becomes_ash(self, tmp_path):
        lines = self._acid_lines("ASP", ("CG", "OD1", "OD2"), "HD2", 0.96)
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        df = read_pdb_to_dataframe(out_file)
        assert (df["residue_name"] == "ASH").all()
        assert _atom_x(out_file, "OD2") == pytest.approx(0.0)
        assert _atom_x(out_file, "OD1") == pytest.approx(2.4)

    def test_acidic_hydrogen_far_from_both_oxygens_is_rejected(self, tmp_path):
        lines = self._acid_lines("GLH", ("CD", "OE1", "OE2"), "HE2", 10.0)
        pdb_file = _write_pdb(tmp_path, lines)

        with pytest.raises(ValueError, match="not bonded to either side-chain oxygen"):
            fix_pdb(pdb_file, out_name=tmp_path / "out.pdb")


class TestNumberedMethylHydrogenNormalization:
    @pytest.mark.parametrize(
        ("resname", "heavy_atom", "source_names", "expected_names"),
        [
            ("ALA", "CB", ("1HB", "2HB", "3HB"), ("HB1", "HB2", "HB3")),
            ("VAL", "CG1", ("1HG1", "2HG1", "3HG1"), ("HG11", "HG12", "HG13")),
            ("VAL", "CG2", ("1HG2", "2HG2", "3HG2"), ("HG21", "HG22", "HG23")),
            ("MET", "CE", ("1HE", "2HE", "3HE"), ("HE1", "HE2", "HE3")),
        ],
    )
    def test_three_methyl_hydrogens_keep_unique_names(
        self, tmp_path, resname, heavy_atom, source_names, expected_names
    ):
        lines = [_line(1, heavy_atom, resname, "A", 1, elem="C")]
        lines.extend(
            _line(index, name, resname, "A", 1, x=float(index), elem="H")
            for index, name in enumerate(source_names, 2)
        )
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        atoms = _get_atom_names(out_file, residue_name=resname)
        assert set(expected_names).issubset(atoms)
        assert len(atoms) == len(set(atoms))


class TestDuplicateAtomValidation:
    def test_unhandled_duplicate_atom_name_is_rejected(self, tmp_path):
        lines = [
            _line(1, "CA", "ALA", "A", 1, elem="C"),
            _line(2, "HA", "ALA", "A", 1, x=1.0, elem="H"),
            _line(3, "HA", "ALA", "A", 1, x=-1.0, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)

        with pytest.raises(ValueError, match="Duplicate atom name HA"):
            fix_pdb(pdb_file, out_name=tmp_path / "out.pdb")

    def test_alias_normalization_cannot_create_duplicate_atom_names(self, tmp_path):
        lines = [
            _line(1, "CB", "ALA", "A", 1, elem="C"),
            _line(2, "1HB", "ALA", "A", 1, x=1.0, elem="H"),
            _line(3, "HB1", "ALA", "A", 1, x=-1.0, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)

        with pytest.raises(ValueError, match="Duplicate atom name HB1"):
            fix_pdb(pdb_file, out_name=tmp_path / "out.pdb")


class TestBackboneAmideGeometryValidation:
    def test_stretched_backbone_hydrogen_is_rejected(self, tmp_path):
        lines = [
            _line(1, "N", "MET", "A", 111, elem="N"),
            _line(2, "H", "MET", "A", 111, x=1.57, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)

        with pytest.raises(ValueError, match="backbone H-N distance"):
            fix_pdb(pdb_file, out_name=tmp_path / "out.pdb")


class TestNeutralLysineHydrogenNormalization:
    def test_hz1_hz2_source_convention_becomes_hz2_hz3(self, tmp_path):
        lines = [
            _line(1, "NZ", "LYN", "A", 11, elem="N"),
            _line(2, "HZ1", "LYN", "A", 11, x=0.9, elem="H"),
            _line(3, "HZ2", "LYN", "A", 11, x=-0.9, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        atoms = _get_atom_names(out_file, residue_name="LYN")
        assert {"HZ2", "HZ3"}.issubset(atoms)
        assert "HZ1" not in atoms
        assert len(atoms) == len(set(atoms))


class TestNeutralArginineGeometryNormalization:
    def test_incomplete_hydrogen_set_is_left_for_qprep_to_complete(self, tmp_path):
        lines = [
            _line(1, "NH1", "ARN", "A", 42, x=0.0, elem="N"),
            _line(2, "NH2", "ARN", "A", 42, x=3.0, elem="N"),
            _line(3, "HH11", "ARN", "A", 42, x=0.9, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        atoms = _get_atom_names(out_file, residue_name="ARN")
        assert {"NH1", "NH2", "HH11"}.issubset(atoms)

    def test_opposite_nitrogen_numbering_is_exchanged(self, tmp_path):
        lines = [
            _line(1, "NH1", "ARN", "A", 42, x=0.0, elem="N"),
            _line(2, "NH2", "ARN", "A", 42, x=3.0, elem="N"),
            _line(3, "HH11", "ARN", "A", 42, x=3.9, y=0.4, elem="H"),
            _line(4, "HH12", "ARN", "A", 42, x=3.9, y=-0.4, elem="H"),
            _line(5, "HH21", "ARN", "A", 42, x=0.9, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        assert _atom_x(out_file, "NH1") == pytest.approx(3.0)
        assert _atom_x(out_file, "NH2") == pytest.approx(0.0)

    def test_duplicate_source_hydrogen_names_are_normalized_before_nesting(self, tmp_path):
        lines = [
            _line(1, "NH1", "ARN", "A", 20, x=0.0, elem="N"),
            _line(2, "NH2", "ARN", "A", 20, x=3.0, elem="N"),
            _line(3, "HH11", "ARN", "A", 20, x=0.9, elem="H"),
            _line(4, "HH21", "ARN", "A", 20, x=3.9, y=0.4, elem="H"),
            _line(5, "HH21", "ARN", "A", 20, x=3.9, y=-0.4, elem="H"),
        ]
        pdb_file = _write_pdb(tmp_path, lines)
        out_file = tmp_path / "out.pdb"

        fix_pdb(pdb_file, out_name=out_file)

        atoms = _get_atom_names(out_file, residue_name="ARN")
        assert {"HH11", "HH12", "HH21"}.issubset(atoms)
        assert len(atoms) == len(set(atoms))
        assert _atom_x(out_file, "NH1") == pytest.approx(3.0)
        assert _atom_x(out_file, "NH2") == pytest.approx(0.0)
