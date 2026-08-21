"""Tests for the QligFEP.CLI.qprep_cli module."""

import argparse
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from QligFEP.CLI.qprep_cli import Neutralizer, parse_arguments


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


class TestNeutralizationDefaults:
    def test_default_is_outer_three_angstrom_scaas_layer(self):
        result, _ = Neutralizer((0, 0, 0), radius=25.0).neutralize_outside_residues_dataframe(
            pd.DataFrame(_build_protein_residue("GLU", "A", 1, 23.0, 0.0, 0.0))
        )
        assert (result["residue_name"] == "GLH").all()

    def test_cli_default_protocol_layer_is_three_angstrom(self, monkeypatch):
        monkeypatch.setattr("sys.argv", ["qprep_prot", "-i", "protein.pdb"])
        assert parse_arguments().neutralize_boundary_offset == 3.0

    def test_q_classifies_exclusion_then_historical_anchor_classifies_included_group(self):
        neutralizer = Neutralizer(
            (0, 0, 0), radius=25.0, boundary_offset=3.0, force_field="AMBER14sb"
        )
        q_group = neutralizer._formal_charge_group("ARG", +1)

        def arg_at(*, q_x, anchor_x, residue_number):
            other_q_x = (q_x * len(q_group) - anchor_x) / (len(q_group) - 1)
            return [
                _make_atom(
                    "ATOM",
                    residue_number * 100 + index,
                    atom["name"],
                    "ARG",
                    "A",
                    residue_number,
                    anchor_x
                    if atom["name"] == "CZ"
                    else other_q_x
                    if atom["name"] in q_group
                    else 0.0,
                    0.0,
                    0.0,
                    atom["name"][0],
                )
                for index, atom in enumerate(neutralizer.residue_library["ARG"]["atoms"], 1)
            ]

        # Q includes residues 1 and 2. Their CZ anchors, not charge-group centers,
        # determine which one occupies the 22--25 A protocol layer. Q excludes
        # residue 3 even though its CZ anchor lies inside the protocol boundary.
        rows = arg_at(q_x=21.0, anchor_x=23.0, residue_number=1)
        rows += arg_at(q_x=23.0, anchor_x=21.0, residue_number=2)
        rows += arg_at(q_x=26.0, anchor_x=21.0, residue_number=3)
        result, stats = neutralizer.neutralize_outside_residues_dataframe(pd.DataFrame(rows), 0.0)

        names = result.groupby("residue_seq_number")["residue_name"].first().to_dict()
        assert names == {1: "ARN", 2: "ARG", 3: "ARN"}
        assert stats["modifications"]["A:1"]["boundary_reference"] == "historical anchor CZ"
        assert stats["modifications"]["A:3"]["boundary_reference"] == "Q charge-group center"

    def test_legacy_parameter_file_defaults_to_switch_atoms(self, monkeypatch):
        monkeypatch.setattr("QligFEP.CLI.qprep_cli.parse_lib", lambda _: {})
        monkeypatch.setattr("QligFEP.CLI.qprep_cli.parse_prm_options", lambda _: {})

        assert Neutralizer((0, 0, 0), force_field="legacy").switch_atoms is True


class TestQChargeGroupNeutralization:
    def test_uses_q_default_whole_residue_charge_group_center(self):
        glu = _build_protein_residue("GLU", "A", 1, 30.0, 0.0, 0.0)
        # Q uses the centroid of the complete default charge group, not CD alone.
        for atom in glu:
            if atom["atom_name"] == "CD":
                atom["x"] = 5.0
        df = pd.DataFrame(glu)

        result, _ = Neutralizer(
            (0, 0, 0), radius=25.0, boundary_offset=0.0
        ).neutralize_outside_residues_dataframe(df)

        assert (result["residue_name"] == "GLH").all()

    def test_salt_bridge_uses_charged_heavy_atom_distance(self):
        rows = []
        glu = _build_protein_residue("GLU", "A", 1, 30.0, 0.0, 0.0)
        glu.append(_make_atom("ATOM", 111, "OE1", "GLU", "A", 1, 28.0, 0.0, 0.0, "O"))
        lys = _build_protein_residue("LYS", "A", 2, 24.0, 0.0, 0.0)
        rows.extend(glu)
        rows.extend(lys)
        df = pd.DataFrame(rows)

        result, stats = Neutralizer(
            (0, 0, 0), radius=25.0, boundary_offset=0.0
        ).neutralize_outside_residues_dataframe(df, salt_bridge_cutoff=4.0)

        assert (result[result["residue_seq_number"] == 1]["residue_name"] == "GLH").all()
        assert (result[result["residue_seq_number"] == 2]["residue_name"] == "LYN").all()
        assert stats["salt_bridges_neutralized"] == 1


class TestOPLSChargeGroupNeutralization:
    @pytest.mark.parametrize(
        ("force_field", "switch_atom"),
        [("OPLS2005", "CZ"), ("OPLS2015", "CG")],
    )
    def test_uses_force_field_formal_charge_group_switch_atom(self, force_field, switch_atom):
        arg = _build_protein_residue("ARG", "A", 1, 0.0, 0.0, 0.0)
        if switch_atom != "CZ":
            arg.append(_make_atom("ATOM", 111, switch_atom, "ARG", "A", 1, 26.0, 0.0, 0.0))
        else:
            for atom in arg:
                if atom["atom_name"] == switch_atom:
                    atom["x"] = 26.0
        result, stats = Neutralizer(
            (0, 0, 0), radius=25.0, force_field=force_field
        ).neutralize_outside_residues_dataframe(pd.DataFrame(arg))

        assert (result["residue_name"] == "ARN").all()
        assert stats["modifications"]["A:1"]["boundary_reference"] == f"Q switch atom {switch_atom}"

    def test_does_not_use_whole_residue_centroid_with_switch_atoms(self):
        arg = _build_protein_residue("ARG", "A", 1, 30.0, 0.0, 0.0)
        for atom in arg:
            if atom["atom_name"] == "CZ":
                atom["x"] = 21.0

        result, _ = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="OPLS2005"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(arg))

        assert (result["residue_name"] == "ARG").all()

    def test_atom_removal_comes_from_selected_force_field(self):
        opls = Neutralizer((0, 0, 0), force_field="OPLS2005")
        amber = Neutralizer((0, 0, 0), force_field="AMBER14sb")
        assert opls._get_atoms_to_remove("LYS", "LYN") == ["HZ3"]
        assert opls._get_atoms_to_remove("ARG", "ARN") == ["HH12"]
        assert amber._get_atoms_to_remove("LYS", "LYN") == ["HZ1"]
        assert amber._get_atoms_to_remove("ARG", "ARN") == ["HH22"]

    def test_opls2005_nonstandard_terminal_names_preserve_sidechain_state(self):
        neutralizer = Neutralizer((0, 0, 0), force_field="OPLS2005")
        assert neutralizer._terminal_neutral_form("NAR+", terminal_charge=+1.0) == "ARG"
        assert neutralizer._terminal_neutral_form("NARG", terminal_charge=+1.0) == "ARN"
        assert neutralizer._terminal_neutral_form("CAR+", terminal_charge=-1.0) == "ARG"
        assert neutralizer._terminal_neutral_form("CARG", terminal_charge=-1.0) == "ARN"

    @pytest.mark.parametrize(
        ("charged_name", "neutral_name", "base_name", "charge"),
        [
            ("NAR+", "NARG", "ARG", +1),
            ("CAR+", "CARG", "ARG", +1),
            ("NLY+", "NLYS", "LYS", +1),
            ("CLY+", "CLYS", "LYS", +1),
        ],
    )
    def test_opls2005_charged_terminal_sidechain_uses_its_library_alias(
        self, charged_name, neutral_name, base_name, charge
    ):
        neutralizer = Neutralizer((0, 0, 0), radius=25.0, force_field="OPLS2005")
        markers = neutralizer._charged_heavy_atoms[base_name]
        switch_atom = neutralizer._formal_charge_group(charged_name, charge, markers)[0]
        rows = [
            _make_atom(
                "ATOM",
                index,
                atom["name"],
                charged_name,
                "A",
                1,
                26.0 if atom["name"] == switch_atom else 10.0,
                0.0,
                0.0,
                atom["name"][0],
            )
            for index, atom in enumerate(neutralizer.residue_library[charged_name]["atoms"], 1)
        ]

        result, stats = neutralizer.neutralize_outside_residues_dataframe(pd.DataFrame(rows), 0.0)

        assert (result["residue_name"] == neutral_name).all()
        assert stats["final_total_charge"] == round(neutralizer._template_charge(neutral_name))
        assert stats["modifications"]["A:1"]["boundary_reference"] == (
            f"Q switch atom {switch_atom}"
        )

    @pytest.mark.parametrize(
        ("terminal_name", "neutral_name"),
        [
            ("NGLU", "GLH"),
            ("CGLU", "GLH"),
            ("NASP", "ASH"),
            ("CASP", "ASH"),
            ("NARG", "ARN"),
            ("CARG", "ARN"),
            ("NLYS", "LYN"),
            ("CLYS", "LYN"),
        ],
    )
    def test_opls2005_terminal_templates_without_sidechain_charge_are_not_classified(
        self, terminal_name, neutral_name
    ):
        neutralizer = Neutralizer((0, 0, 0), force_field="OPLS2005")
        rows = [
            _make_atom("ATOM", index, atom["name"], terminal_name, "A", 1, 30.0, 0.0, 0.0)
            for index, atom in enumerate(neutralizer.residue_library[terminal_name]["atoms"], 1)
        ]

        result, _ = neutralizer.neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == neutral_name).all()


class TestNeutralizerDNA:
    """Tests for DNA nucleotide handling in the Neutralizer.

    DNA nucleotides whose C1' atom is outside the sphere radius are removed.
    DNA in the outer SCAAS layer is retained. Each retained fragment is capped
    with complementary 5' and 3' terminal templates at cut boundaries.
    """

    def _make_neutralizer(self, center=(0, 0, 0), radius=25.0, offset=3.0):
        return Neutralizer(center, radius, offset)

    def test_unsupported_force_field_reports_dna_limitation(self):
        neutralizer = Neutralizer((0, 0, 0), force_field="CHARMM36")
        df = pd.DataFrame(_build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0))

        with pytest.raises(ValueError, match="CHARMM36 does not provide DNA residue templates"):
            neutralizer.neutralize_outside_residues_dataframe(df)

    def test_outside_dna_removed(self):
        """DA nucleotide fully outside sphere should be removed."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0))  # inside
        rows.extend(_build_dna_residue("DA", "E", 2, 30.0, 0.0, 0.0))  # outside
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert len(result_df[result_df["residue_seq_number"] == 2]) == 0

    def test_5prime_cap_on_inside_boundary(self):
        """First inside residue with removed upstream DNA gets 5' terminal form."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 30.0, 0.0, 0.0))  # outside, removed
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))  # inside, gets 5' cap
        rows.extend(_build_dna_residue("DC", "E", 3, 5.0, 0.0, 0.0))  # inside, unchanged
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
        rows.extend(_build_dna_residue("DT", "E", 4, 5.0, 0.0, 0.0))  # inside, gets cap
        rows.extend(_build_dna_residue("DA", "E", 5, 5.0, 0.0, 0.0))  # inside
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
        """Protein uses rest_bound while DNA is cut and capped at the sphere."""
        rows = []
        rows.extend(_build_protein_residue("GLU", "A", 1, 23.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DA", "E", 10, 5.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DA", "E", 11, 30.0, 0.0, 0.0))
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        assert (result_df[result_df["residue_seq_number"] == 1]["residue_name"] == "GLH").all()
        assert (result_df[result_df["residue_seq_number"] == 10]["residue_name"] == "DA3").all()
        assert len(result_df[result_df["residue_seq_number"] == 11]) == 0

    def test_downstream_removal_adds_3prime_cap(self):
        """DNA next to a removed downstream segment gets a 3' cap."""
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 5.0, 0.0, 0.0))  # inside
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))  # inside
        rows.extend(_build_dna_residue("DC", "E", 3, 30.0, 0.0, 0.0))  # outside
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        # Res 3: removed
        assert len(result_df[result_df["residue_seq_number"] == 3]) == 0
        assert (result_df[result_df["residue_seq_number"] == 1]["residue_name"] == "DA").all()
        assert (result_df[result_df["residue_seq_number"] == 2]["residue_name"] == "DG3").all()

    def test_cut_fragment_has_integral_template_charge(self):
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 30.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DC", "E", 3, 5.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DT", "E", 4, 30.0, 0.0, 0.0))

        neutralizer = self._make_neutralizer()
        result_df, _ = neutralizer.neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        names = result_df.drop_duplicates("residue_seq_number")["residue_name"]
        charge = sum(neutralizer._template_charge(name) for name in names)
        assert set(names) == {"DG5", "DC3"}
        assert charge == pytest.approx(-1.0)

    def test_single_retained_nucleotide_uses_neutral_template(self):
        rows = []
        rows.extend(_build_dna_residue("DA", "E", 1, 30.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DG", "E", 2, 5.0, 0.0, 0.0))
        rows.extend(_build_dna_residue("DC", "E", 3, 30.0, 0.0, 0.0))

        neutralizer = self._make_neutralizer()
        result_df, _ = neutralizer.neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result_df["residue_name"] == "DGN").all()
        assert "P" not in result_df["atom_name"].values
        assert neutralizer._template_charge("DGN") == pytest.approx(0.0)


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

    def test_npro_outside_boundary_drops_backbone_h(self):
        """NPRO outside rest_bound should become PRO with no backbone amide H.

        Proline's nitrogen is tertiary (bonded to CA, CD, and previous-C only).
        The standard PRO library entry in AMBER14sb has no N-H atom at all.
        Naively renaming H1→H (as for NILE etc.) leaves an orphan H atom and
        qprep aborts with "Too many atoms in residue PRO ...".
        """
        atoms = []
        serial_base = 100
        for i, (aname, elem) in enumerate(
            [
                ("N", "N"),
                ("CA", "C"),
                ("C", "C"),
                ("O", "O"),
                ("CB", "C"),
                ("CG", "C"),
                ("CD", "C"),
                ("H1", "H"),
                ("H2", "H"),
            ]
        ):
            atoms.append(_make_atom("ATOM", serial_base + i, aname, "NPRO", "F", 1, 30.0, 0.0, 0.0, elem))
        df = pd.DataFrame(atoms)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "PRO").all()
        for forbidden in ("H", "H1", "H2", "H3"):
            assert (
                forbidden not in res["atom_name"].values
            ), f"PRO library has no backbone amide H but found {forbidden!r}"


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

    def test_cterm_protocol_layer_preserves_historical_n_anchor(self):
        rows = _build_cterm_residue("CALA", "A", 1, 23.0, 0.0, 0.0)
        for atom in rows:
            if atom["atom_name"] == "N":
                atom["x"] = 21.0

        result, _ = self._make_neutralizer().neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "CALA").all()

    def test_q_excluded_cterm_is_neutralized_even_when_historical_anchor_is_inside(self):
        rows = _build_cterm_residue("CALA", "A", 1, 21.0, 0.0, 0.0)
        for atom in rows:
            if atom["atom_name"] != "N":
                atom["x"] = 27.0

        result, stats = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="AMBER14sb"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "ALA").all()
        assert stats["modifications"]["A:1"]["boundary_reference"] == "Q charge-group center"

    def test_cgly_outside_boundary_neutralized(self):
        """CGLY outside rest_bound should become GLY."""
        rows = _build_cterm_residue("CGLY", "B", 1, 30.0, 0.0, 0.0)
        df = pd.DataFrame(rows)

        neutralizer = self._make_neutralizer()
        result_df, stats = neutralizer.neutralize_outside_residues_dataframe(df)

        res = result_df[result_df["residue_seq_number"] == 1]
        assert (res["residue_name"] == "GLY").all()
        assert "OXT" not in res["atom_name"].values


class TestCHARMMNeutralization:
    @staticmethod
    def _template_residue(
        neutralizer: Neutralizer,
        residue_name: str,
        group_positions: dict[int, float],
    ) -> pd.DataFrame:
        entry = neutralizer.residue_library[residue_name]
        group_by_atom = {
            atom_name: group_number
            for group_number, group in enumerate(entry["charge_groups"], 1)
            for atom_name in group
        }
        rows = []
        for serial, atom in enumerate(entry["atoms"], 1):
            x = group_positions[group_by_atom[atom["name"]]]
            rows.append(
                _make_atom(
                    "ATOM",
                    serial,
                    atom["name"],
                    residue_name,
                    "A",
                    1,
                    x,
                    0.0,
                    0.0,
                    atom["name"][0],
                )
            )
        return pd.DataFrame(rows)

    def test_standard_terminal_names_preserve_residue_identity(self):
        neutralizer = Neutralizer((0, 0, 0), force_field="CHARMM36")

        assert neutralizer._terminal_neutral_form("NSER", terminal_charge=+1.0) == "SER"
        assert neutralizer._terminal_neutral_form("CSER", terminal_charge=-1.0) == "SER"

    def test_charged_sidechain_uses_charmm_charge_groups_and_atom_names(self):
        rows = _build_protein_residue("LYS", "A", 1, 23.0, 0.0, 0.0)
        rows.extend(
            _make_atom("ATOM", 110 + index, name, "LYS", "A", 1, 23.0, 0.0, 0.0, "H")
            for index, name in enumerate(("HZ1", "HZ2", "HZ3"), 1)
        )

        result, stats = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="CHARMM36"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "LYN").all()
        assert "HZ3" not in result["atom_name"].values
        assert stats["modifications"]["A:1"]["boundary_reference"] == "historical anchor NZ"

    def test_hip_uses_merged_integer_charge_group(self):
        first_group = ("CB", "HB1", "HB2", "CD2", "HD2", "CG")
        second_group = ("NE2", "HE2", "ND1", "HD1", "CE1", "HE1")
        rows = [
            _make_atom("ATOM", index, name, "HIP", "A", 1, 21.0, 0.0, 0.0, name[0])
            for index, name in enumerate(first_group, 1)
        ]
        rows.extend(
            _make_atom("ATOM", index, name, "HIP", "A", 1, 31.0, 0.0, 0.0, name[0])
            for index, name in enumerate(second_group, len(rows) + 1)
        )

        result, stats = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="CHARMM36"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "HID").all()
        assert stats["modifications"]["A:1"]["boundary_reference"] == "Q charge-group center"

    def test_terminal_sidechain_is_classified_before_terminal_charge(self):
        neutralizer = Neutralizer((0, 0, 0), radius=25.0, force_field="CHARMM36")
        rows = self._template_residue(
            neutralizer,
            "NGLU",
            {1: 10.0, 2: 10.0, 3: 23.0, 4: 10.0},
        )

        result, stats = neutralizer.neutralize_outside_residues_dataframe(rows)

        assert (result["residue_name"] == "NGLH").all()
        assert stats["terminals_neutralized"] == 0
        assert stats["modifications"]["A:1"]["original"] == "NGLU"
        assert stats["modifications"]["A:1"]["modified"] == "NGLH"

    def test_sidechain_and_terminal_conversions_count_residue_once(self):
        neutralizer = Neutralizer((0, 0, 0), radius=25.0, force_field="CHARMM36")
        rows = self._template_residue(
            neutralizer,
            "NGLU",
            {1: 23.0, 2: 23.0, 3: 23.0, 4: 23.0},
        )

        result, stats = neutralizer.neutralize_outside_residues_dataframe(rows)

        assert (result["residue_name"] == "GLH").all()
        assert stats["terminals_neutralized"] == 1
        assert stats["residues_neutralized"] == 1
        assert stats["modifications"]["A:1"]["original"] == "NGLU"
        assert stats["modifications"]["A:1"]["modified"] == "GLH"

    def test_terminal_sidechain_group_is_disambiguated_from_same_sign_terminal_group(self):
        neutralizer = Neutralizer((0, 0, 0), radius=25.0, force_field="CHARMM36")
        rows = self._template_residue(
            neutralizer,
            "NARG",
            {1: 10.0, 2: 10.0, 3: 10.0, 4: 23.0, 5: 10.0},
        )

        result, stats = neutralizer.neutralize_outside_residues_dataframe(rows)

        # CHARMM has no NARN template, so conservatively neutralize both sites.
        assert (result["residue_name"] == "ARN").all()
        assert stats["residues_neutralized"] == 1

    def test_fractional_prefixed_neutral_template_uses_conservative_fallback(self):
        neutralizer = Neutralizer((0, 0, 0), radius=25.0, force_field="OPLS2015")
        entry = neutralizer.residue_library["CARG"]
        charges = {atom["name"]: atom["charge"] for atom in entry["atoms"]}
        positions = {
            group_number: (23.0 if "CZ" in group else 10.0)
            for group_number, group in enumerate(entry["charge_groups"], 1)
        }
        assert any(
            abs(sum(charges[atom] for atom in group) - 1.0) <= 1e-6
            for group in entry["charge_groups"]
            if "CZ" in group
        )
        rows = self._template_residue(neutralizer, "CARG", positions)

        result, _ = neutralizer.neutralize_outside_residues_dataframe(rows)

        assert (result["residue_name"] == "ARN").all()

    def test_nterminal_ht1_coordinate_is_preserved_as_internal_h(self):
        rows = [
            _make_atom("ATOM", 1, name, "NALA", "A", 1, 30.0, 0.0, 0.0, element)
            for name, element in (("N", "N"), ("HT1", "H"), ("HT2", "H"), ("HT3", "H"), ("CA", "C"))
        ]

        result, _ = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="CHARMM36"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "ALA").all()
        assert "H" in result["atom_name"].values
        assert not {"HT1", "HT2", "HT3"} & set(result["atom_name"])

    def test_cterminal_ot1_coordinate_is_preserved_as_internal_o(self):
        rows = [
            _make_atom("ATOM", 1, name, "CALA", "A", 1, 30.0, 0.0, 0.0, element)
            for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("OT1", "O"), ("OT2", "O"))
        ]

        result, _ = Neutralizer(
            (0, 0, 0), radius=25.0, force_field="CHARMM36"
        ).neutralize_outside_residues_dataframe(pd.DataFrame(rows))

        assert (result["residue_name"] == "ALA").all()
        assert "O" in result["atom_name"].values
        assert not {"OT1", "OT2"} & set(result["atom_name"])


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


class TestNeutralizerInsertionCode:
    """Residues sharing a (chain, seq_number) via insertion codes must not collide.

    A C-terminal NME cap is commonly numbered the same as the preceding residue
    with an insertion code (e.g. GLU 167 and NME 167A). Neutralizing the GLU must
    not rewrite the NME cap that shares chain+seq_number.
    """

    def _make_neutralizer(self, center=(0, 0, 0), radius=25.0, offset=3.0):
        return Neutralizer(center, radius, offset)

    def test_neutralizing_glu_does_not_touch_capped_nme(self):
        # GLU 167 outside boundary (x=30 > rest_bound 22) -> should become GLH.
        glu = _build_protein_residue("GLU", "A", 167, 30.0, 0.0, 0.0)
        # NME cap sharing chain A + seq 167 but insertion code "A" -> must be untouched.
        nme = []
        for i, aname in enumerate(["N", "C", "H", "H1", "H2", "H3"]):
            row = _make_atom("ATOM", 99000 + i, aname, "NME", "A", 167, 30.0, 0.0, 0.0, "H")
            row["insertion_code"] = "A"
            nme.append(row)
        df = pd.DataFrame(glu + nme)

        neutralizer = self._make_neutralizer()
        out, _ = neutralizer.neutralize_outside_residues_dataframe(df)

        nme_out = out[out["insertion_code"] == "A"]
        assert (nme_out["residue_name"] == "NME").all(), "NME cap was corrupted by GLU neutralization"
        assert sorted(nme_out["atom_name"]) == ["C", "H", "H1", "H2", "H3", "N"]
        glu_out = out[(out["insertion_code"] == "") & (out["residue_name"].isin(["GLU", "GLH"]))]
        assert (glu_out["residue_name"] == "GLH").all(), "GLU outside boundary should be neutralized to GLH"
