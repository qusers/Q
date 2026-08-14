"""Tests for the amino-acid mutation FEP setup (QresFEP).

The end-to-end tests run the real ``qprep`` binary on T4 lysozyme and check the
files it produces; nothing here is mocked. The reference values they assert
against come from the published QresFEP implementation and from the worked
example shipped with its tutorial.
"""

import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from QligFEP import amino_acids
from QligFEP.analyze_resfep import (
    LegResult,
    MutationResult,
    collect,
    read_qfep_stage,
    to_frame,
)
from QligFEP.functions import resfep_lambda_ladder
from QligFEP.peptide_caps import PeptideBuildError, build_reference_peptide
from QligFEP.qresfep import MutationError, QresFEP, parse_mutation
from QligFEP.settings.settings import BIN
from QligFEP.sphere_prep import (
    ResidueMapping,
    SpherePrep,
    build_residue_map,
    parse_disulfide_pairs,
)

PROJECT_ROOT = Path(__file__).parent.parent.parent
T4L_RESOURCES = Path(__file__).parent / "resources" / "t4l"
T4L_TUTORIAL = PROJECT_ROOT / "tutorials" / "T4L"

#: CB of LEU39 chain A -- the sphere centre the QresFEP tutorial uses.
LEU39_CB = ("41.088", "31.308", "21.828")

qprep_required = pytest.mark.skipif(
    not (BIN / "qprep").exists(),
    reason="qprep binary not compiled (cd src/q6 && make all)",
)


# ----------------------------------------------------------------------
# Mutation strings
# ----------------------------------------------------------------------


class TestParseMutation:
    def test_three_letter_codes(self):
        assert parse_mutation("LEU39ALA") == ("LEU", 39, "ALA")

    def test_one_letter_codes(self):
        assert parse_mutation("L39A") == ("LEU", 39, "ALA")

    def test_mixed_codes(self):
        assert parse_mutation("LEU39A") == ("LEU", 39, "ALA")

    def test_lowercase_is_accepted(self):
        assert parse_mutation("leu39ala") == ("LEU", 39, "ALA")

    def test_protonation_variants_are_distinct_residues(self):
        """A neutralised residue is a different wild type, not the same one."""
        assert parse_mutation("ASH92ASN") == ("ASH", 92, "ASN")

    @pytest.mark.parametrize("mutation", ["LEU39", "39ALA", "LEUALA", "", "LEU39XYZ"])
    def test_unreadable_mutations_raise(self, mutation):
        with pytest.raises(MutationError):
            parse_mutation(mutation)


# ----------------------------------------------------------------------
# Residue chemistry
# ----------------------------------------------------------------------


class TestAtomNaming:
    """PyMOL numbers branched hydrogens the opposite way from the Q libraries."""

    def test_leading_index_moves_to_the_end(self):
        line = "ATOM      8 1HB  ALA A  39      41.547  31.477  20.801  1.00  0.00           H"
        assert amino_acids.normalize_pdb_atom_name(line)[12:16] == " HB1"

    def test_four_character_name_fills_the_field(self):
        line = "ATOM     12 1HG1 ILE A   3      41.547  31.477  20.801  1.00  0.00           H"
        assert amino_acids.normalize_pdb_atom_name(line)[12:16] == "HG11"

    def test_glycine_alpha_hydrogen_is_numbered(self):
        line = "ATOM      7  HA  GLY A  25      39.755  29.732  21.169  1.00  0.00           H"
        assert amino_acids.normalize_pdb_atom_name(line)[12:16] == " HA2"

    def test_non_glycine_alpha_hydrogen_is_left_alone(self):
        line = "ATOM      7  HA  ALA A  39      39.755  29.732  21.169  1.00  0.00           H"
        assert amino_acids.normalize_pdb_atom_name(line)[12:16] == " HA "

    def test_conforming_names_are_unchanged(self):
        line = "ATOM      5  CB  ALA A  39      41.123  31.295  21.807  1.00  0.00           C"
        assert amino_acids.normalize_pdb_atom_name(line) == line

    def test_non_atom_lines_are_unchanged(self):
        for line in ("TER   ", "END", ""):
            assert amino_acids.normalize_pdb_atom_name(line) == line

    def test_the_rest_of_the_line_survives(self):
        line = "ATOM      8 1HB  ALA A  39      41.547  31.477  20.801  1.00  0.00           H"
        result = amino_acids.normalize_pdb_atom_name(line)
        assert result[:12] == line[:12]
        assert result[16:] == line[16:]


class TestBackboneAtoms:
    def test_hа_is_shared_for_an_ordinary_mutation(self):
        assert "HA" in amino_acids.backbone_atoms(from_gly=False, to_gly=False)

    @pytest.mark.parametrize("from_gly,to_gly", [(True, False), (False, True)])
    def test_ha_is_perturbed_for_glycine_mutations(self, from_gly, to_gly):
        """Glycine changes the alpha carbon's hydrogen count, so HA cannot be shared."""
        assert "HA" not in amino_acids.backbone_atoms(from_gly, to_gly)

    def test_the_rest_of_the_backbone_is_always_shared(self):
        for from_gly, to_gly in [(False, False), (True, False), (False, True)]:
            backbone = amino_acids.backbone_atoms(from_gly, to_gly)
            assert {"N", "H", "C", "O", "CA"} <= set(backbone)


class TestSideChains:
    def test_mutant_atoms_are_lower_cased(self):
        wild_type, mutant = amino_acids.side_chain_pair("LEU", "ALA")
        assert all(atom.isupper() or atom.isdigit() for atom in "".join(wild_type))
        assert mutant == [a.lower() for a in amino_acids.SIDE_CHAINS["ALA"]]

    def test_glycine_has_no_side_chain(self):
        _, mutant = amino_acids.side_chain_pair("LEU", "GLY")
        assert mutant == []

    def test_every_side_chain_residue_is_in_the_atom_name_table(self):
        """The two tables have to agree, or atoms silently drop out of the topology."""
        known = set(amino_acids.SIDE_CHAIN_ATOM_NAMES)
        for residue, atoms in amino_acids.SIDE_CHAINS.items():
            missing = set(atoms) - known
            assert not missing, f"{residue} has atoms absent from SIDE_CHAIN_ATOM_NAMES: {missing}"

    def test_chi1_atoms_belong_to_the_side_chain(self):
        for residue in amino_acids.SIDE_CHAINS:
            for atom in amino_acids.chi1_atoms(residue):
                assert atom in amino_acids.SIDE_CHAINS[residue], f"{residue}: {atom}"

    @pytest.mark.parametrize(
        "variant,charged_form",
        [("ASH", "ASP"), ("GLH", "GLU"), ("HIP", "HID"), ("ARN", "ARG"), ("LYN", "LYS")],
    )
    def test_protonation_variants_keep_the_parent_chi1(self, variant, charged_form):
        assert amino_acids.chi1_atoms(variant) == amino_acids.chi1_atoms(charged_form)

    def test_glycine_has_no_chi1(self):
        assert amino_acids.chi1_atoms("GLY") == []


class TestZeroTorsions:
    @pytest.mark.parametrize("variant", ["ASH", "GLH", "HIP", "ARN", "LYN", "PHE"])
    @pytest.mark.parametrize("as_mutant", [False, True])
    def test_variants_and_phe_zero_cross_topology_chi1_torsions(self, variant, as_mutant):
        mutation = f"ALA1{variant}" if as_mutant else f"{variant}1ALA"
        run = QresFEP(
            mutation=mutation,
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            replicates=1,
        )
        atom_names = (
            "HB1",
            "HB2",
            "HB3",
            "CG",
            "CB",
            "CA",
            "cb",
            "cg",
            "hb1",
            "hb2",
            "hb3",
        )
        run.atom_ids = {name: str(index) for index, name in enumerate(atom_names, 1)}

        torsions = run._zero_torsions()

        assert len(torsions) == 6  # three for each side-chain topology
        assert all(line.split()[-2:] == ["0", "0"] for line in torsions)


class TestHeavyAtomMatching:
    @staticmethod
    def _atom(name, xyz):
        # Positional record as produced by pdb_utils.pdb_parse_in.
        record = [None] * 14
        record[2] = name
        record[8], record[9], record[10] = xyz
        return record

    def test_superimposed_equivalents_are_paired(self):
        hybrid = [self._atom("CB", (0.0, 0.0, 0.0)), self._atom("cb", (0.1, 0.0, 0.0))]
        assert amino_acids.match_heavy_atoms(["CB"], ["cb"], hybrid) == [("CB", "cb")]

    def test_distant_atoms_are_not_paired(self):
        hybrid = [self._atom("CB", (0.0, 0.0, 0.0)), self._atom("cb", (5.0, 0.0, 0.0))]
        assert amino_acids.match_heavy_atoms(["CB"], ["cb"], hybrid) == []

    def test_hydrogens_are_excluded(self):
        """Restraining hydrogens would fight the force field's own geometry."""
        hybrid = [self._atom("HB2", (0.0, 0.0, 0.0)), self._atom("hb2", (0.05, 0.0, 0.0))]
        assert amino_acids.match_heavy_atoms(["HB2"], ["hb2"], hybrid) == []

    def test_chemically_unrelated_overlaps_are_not_paired(self):
        """Only atoms whose names start with the same letter can correspond."""
        hybrid = [self._atom("CB", (0.0, 0.0, 0.0)), self._atom("og", (0.1, 0.0, 0.0))]
        assert amino_acids.match_heavy_atoms(["CB"], ["og"], hybrid) == []


class TestProteinComplexCoordinates:
    def test_crystallographic_waters_are_left_for_the_solvent_sphere(self, tmp_path):
        """The solute complex must not protect a crystal water that clashes with the hybrid."""
        protein = tmp_path / "protein_processed.pdb"
        protein.write_text(
            "ATOM      1  N   LEU A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  LEU A   1       1.450   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   LEU A   1       2.000   1.400   0.000  1.00  0.00           C\n"
            "ATOM      4  O   LEU A   1       1.400   2.400   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  LEU A   1       1.900  -0.800   1.200  1.00  0.00           C\n"
            "HETATM    6  O   HOH W   2       2.500  -1.000   1.500  1.00 20.00           O\n"
            "HETATM    7  H1  HOH W   2       3.300  -1.000   1.500  1.00 20.00           H\n"
            "HETATM    8  H2  HOH W   2       2.200  -0.250   1.500  1.00 20.00           H\n"
        )
        (tmp_path / "ALA1.pdb").write_text(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       1.400   2.400   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA A   1       1.900  -0.800   1.200  1.00  0.00           C\n"
        )
        run = QresFEP(
            mutation="LEU1ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=tmp_path,
        )
        run.protein_pdb = protein
        run.q_position = 1

        run.read_pdb()

        assert all(atom[4] != "HOH" for residue in run.pdb.values() for atom in residue)
        assert set(run.pdb) == {1}


class TestProductionSteps:
    def test_custom_step_count_is_stored(self, tmp_path):
        run = QresFEP(
            mutation="LEU1ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            production_steps=25_000,
            workdir=tmp_path,
        )

        assert run.production_steps == 25_000

    def test_custom_step_count_is_written_to_each_window(self, tmp_path):
        run = QresFEP(
            mutation="LEU1ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            windows=1,
            production_steps=25_000,
            workdir=tmp_path,
        )
        run.inputfiles.mkdir(parents=True)
        run.shell_radius = 25
        run.system_size = 0
        run.counter_ions = 0
        run.anchor_atom = None
        run._distance_restraints = lambda: ""
        run.lambda_schedule = lambda: [["1.000"], ["1.000"]]

        stage_files = run.write_production()

        for names in stage_files:
            for name in names:
                assert "steps                     25000" in (run.inputfiles / f"{name}.inp").read_text()

    @pytest.mark.parametrize("steps", [0, -1])
    def test_step_count_must_be_positive(self, tmp_path, steps):
        with pytest.raises(MutationError, match="Production steps must be a positive integer"):
            QresFEP(
                mutation="LEU1ALA",
                chain="A",
                system="protein",
                force_field="OPLSAAM",
                cluster="SNELLIUS",
                production_steps=steps,
                workdir=tmp_path,
            )


class TestReplicateSeeds:
    def test_seed_count_must_match_replicates(self, tmp_path):
        with pytest.raises(MutationError, match="Expected 3 random seeds, received 2"):
            QresFEP(
                mutation="LEU1ALA",
                chain="A",
                system="protein",
                force_field="OPLSAAM",
                cluster="SNELLIUS",
                replicates=3,
                seeds=[971, 2856],
                workdir=tmp_path,
            )

    @pytest.mark.parametrize("seed", [0, 10000, 1.5, "971"])
    def test_seed_values_follow_qs_valid_range(self, tmp_path, seed):
        with pytest.raises(MutationError, match="integers from 1 through 9999"):
            QresFEP(
                mutation="LEU1ALA",
                chain="A",
                system="protein",
                force_field="OPLSAAM",
                cluster="SNELLIUS",
                replicates=1,
                seeds=[seed],
                workdir=tmp_path,
            )


# ----------------------------------------------------------------------
# Lambda ladders
# ----------------------------------------------------------------------


class TestLambdaLadder:
    @pytest.mark.parametrize("sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"])
    def test_runs_from_one_to_zero(self, sampling):
        ladder = resfep_lambda_ladder(25, sampling)
        assert ladder[0] == "1.000"
        assert ladder[-1] == "0.000"

    @pytest.mark.parametrize("sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"])
    def test_length_is_exactly_the_window_count(self, sampling):
        """One value per window -- unlike the ligand ladder, which returns windows + 1."""
        assert len(resfep_lambda_ladder(25, sampling)) == 25

    @pytest.mark.parametrize("sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"])
    def test_decreases_monotonically(self, sampling):
        values = [float(v) for v in resfep_lambda_ladder(25, sampling)]
        assert all(a >= b for a, b in zip(values, values[1:]))

    def test_no_negative_zero(self):
        """`-0.000` would end up in file names and in Q's [lambdas] section."""
        for sampling in ("linear", "sigmoidal", "exponential", "reverse_exponential"):
            assert "-0.000" not in resfep_lambda_ladder(25, sampling)

    def test_exponential_pair_is_mirrored(self):
        """The two stages' ladders mirror each other, which is why they pair up."""
        forward = [float(v) for v in resfep_lambda_ladder(25, "exponential")]
        reverse = [float(v) for v in resfep_lambda_ladder(25, "reverse_exponential")]
        for f, r in zip(forward, reversed(reverse)):
            assert f == pytest.approx(1.0 - r, abs=1e-3)

    def test_matches_the_published_exponential_ladder(self):
        """Values from the reference QresFEP implementation, 25 windows."""
        expected = [
            "1.000", "0.993", "0.985", "0.976", "0.966", "0.955", "0.941", "0.927",
            "0.910", "0.891", "0.870", "0.845", "0.818", "0.786", "0.751", "0.711",
            "0.665", "0.614", "0.555", "0.489", "0.414", "0.329", "0.233", "0.124",
            "0.000",
        ] # fmt: skip
        assert resfep_lambda_ladder(25, "exponential") == expected

    def test_unknown_scheme_raises(self):
        with pytest.raises(ValueError, match="Unknown lambda sampling"):
            resfep_lambda_ladder(25, "quadratic")

    @pytest.mark.parametrize("windows", [0, 1, -5])
    def test_too_few_windows_raises(self, windows):
        with pytest.raises(ValueError, match="at least 2 windows"):
            resfep_lambda_ladder(windows, "linear")


# ----------------------------------------------------------------------
# Sphere preparation record
# ----------------------------------------------------------------------


class TestSpherePrep:
    @staticmethod
    def _record():
        return SpherePrep(
            input_pdb="protein.pdb",
            prepared_pdb="protein_processed.pdb",
            force_field="OPLSAAM",
            center=[1.0, 2.0, 3.0],
            radius=25.0,
            total_charge=2,
            disulfides=[(30, 97)],
            residues=[
                ResidueMapping(1, 10, "A", "MET"),
                ResidueMapping(2, 11, "A", "LEU"),
                ResidueMapping(3, 10, "B", "GLU"),
            ],
        )

    def test_round_trip_through_json(self, tmp_path):
        original = self._record()
        original.write(tmp_path)
        loaded = SpherePrep.read(tmp_path)
        assert loaded == original

    def test_read_accepts_a_directory_or_a_file(self, tmp_path):
        path = self._record().write(tmp_path)
        assert SpherePrep.read(tmp_path) == SpherePrep.read(path)

    def test_missing_file_names_the_command_that_writes_it(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="qprep_prot"):
            SpherePrep.read(tmp_path)

    def test_translates_pdb_numbering_to_q_numbering(self):
        assert self._record().q_number(11, "A") == 2

    def test_ambiguous_residue_number_requires_a_chain(self):
        """Residue 10 exists in both chains, so the answer is not well defined."""
        with pytest.raises(KeyError, match="chains"):
            self._record().q_number(10)

    def test_chain_disambiguates(self):
        record = self._record()
        assert record.q_number(10, "A") == 1
        assert record.q_number(10, "B") == 3

    def test_absent_residue_raises(self):
        with pytest.raises(KeyError, match="not in the prepared sphere"):
            self._record().q_number(999, "A")

    def test_residue_name_follows_the_same_lookup(self):
        assert self._record().residue_name(11, "A") == "LEU"


class TestResidueMap:
    PDB = (
        "ATOM      1  N   MET A 100       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  MET A 100       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   LEU A 101       2.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      4  O   HOH W 500       9.000   9.000   9.000  1.00  0.00           O\n"
    )

    def test_q_numbering_is_positional(self, tmp_path):
        pdb = tmp_path / "prepared.pdb"
        pdb.write_text(self.PDB)
        mapping = build_residue_map(pdb)
        assert [entry.q_number for entry in mapping] == [1, 2]

    def test_water_is_excluded(self, tmp_path):
        """qprep adds the sphere's solvent itself, after the solute."""
        pdb = tmp_path / "prepared.pdb"
        pdb.write_text(self.PDB)
        assert all(entry.name != "HOH" for entry in build_residue_map(pdb))

    def test_original_numbering_is_reported_when_known(self, tmp_path):
        """Reindexing destroys the user's numbering; only the returned map has it."""
        pdb = tmp_path / "prepared.pdb"
        pdb.write_text(self.PDB)
        original = {100: (7, "A", ""), 101: (8, "A", "")}
        mapping = build_residue_map(pdb, original_numbering=original)
        assert [(e.q_number, e.pdb_number) for e in mapping] == [(1, 7), (2, 8)]

    def test_falls_back_to_the_prepared_numbering(self, tmp_path):
        pdb = tmp_path / "prepared.pdb"
        pdb.write_text(self.PDB)
        mapping = build_residue_map(pdb)
        assert [e.pdb_number for e in mapping] == [100, 101]


class TestDisulfideTranslation:
    def test_prepared_numbers_become_q_numbers(self):
        lines = "addbond 30:SG 97:SG y\n"
        assert parse_disulfide_pairs(lines, {30: 3, 97: 9}) == [(3, 9)]

    def test_commented_bonds_are_still_read(self):
        """qprep_prot writes them commented out; they are still the detected bonds."""
        lines = "!addbond 30:SG 97:SG y\n"
        assert parse_disulfide_pairs(lines, {30: 3, 97: 9}) == [(3, 9)]

    def test_empty_input_gives_no_pairs(self):
        assert parse_disulfide_pairs("", {}) == []

    def test_unknown_partner_is_kept_unchanged(self):
        assert parse_disulfide_pairs("addbond 30:SG 97:SG y\n", {30: 3}) == [(3, 97)]


# ----------------------------------------------------------------------
# Reading results
# ----------------------------------------------------------------------

REFERENCE_EXAMPLE = Path("/Users/davidararipe/projects/qligfep/tutorials/2.QresFEP_T4L/FEP_example")


class TestQfepReading:
    """Read the worked example shipped with the reference tutorial.

    Its expected values are printed in that tutorial's README, so these assert
    against numbers published independently of this implementation.
    """

    @pytest.fixture
    def stage_one(self):
        qfep_out = REFERENCE_EXAMPLE / "FEP1" / "298" / "1" / "qfep.out"
        if not qfep_out.exists():
            pytest.skip("reference QresFEP tutorial output not available")
        return qfep_out

    def test_bar_matches_the_published_value(self, stage_one):
        assert read_qfep_stage(stage_one, reverse_direction=False)["dG_bar"] == pytest.approx(65.81, abs=0.01)

    def test_all_estimators_are_read(self, stage_one):
        values = read_qfep_stage(stage_one, reverse_direction=False)
        assert not any(math.isnan(v) for v in values.values())

    def test_reversing_the_direction_reads_the_other_end_and_negates(self, tmp_path):
        """A mutation out of glycine is run from lambda 0 to 1, so qfep accumulates
        towards the opposite row and the sign of the result is flipped."""
        qfep_out = tmp_path / "qfep.out"
        qfep_out.write_text(
            "# Part 1: Free energy perturbation summary:\n"
            "# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>\n"
            "   1.000000    0.000    7.000   -0.375   -7.000    7.000\n"
            "   0.000000    0.388    3.000    0.000    0.000    3.000\n"
            "# Part 2: Reaction free energy summary:\n"
            "# Part 6: BAR Bennet:\n"
            "   1.000000    0.000    7.000\n"
            "   0.000000    0.388    3.000\n"
        )
        forward = read_qfep_stage(qfep_out, reverse_direction=False)
        reverse = read_qfep_stage(qfep_out, reverse_direction=True)

        assert forward["dG_bar"] == pytest.approx(3.0)  # read at lambda 0
        assert reverse["dG_bar"] == pytest.approx(-7.0)  # read at lambda 1, negated

    def test_qfep_errors_are_raised(self, tmp_path):
        broken = tmp_path / "qfep.out"
        broken.write_text("ERROR: Too few data points in file.\n")
        with pytest.raises(OSError, match="qfep reported an error"):
            read_qfep_stage(broken, reverse_direction=False)


class TestLegAggregation:
    def test_mean_and_sem_over_replicates(self):
        leg = LegResult(leg="protein", replicates={"dG_bar": [10.0, 12.0, 14.0]})
        assert leg.mean() == pytest.approx(12.0)
        assert leg.sem() == pytest.approx(2.0 / math.sqrt(3))

    def test_sem_needs_two_replicates(self):
        leg = LegResult(leg="protein", replicates={"dG_bar": [10.0]})
        assert math.isnan(leg.sem())

    def test_failed_replicates_are_ignored_not_counted_as_zero(self):
        leg = LegResult(leg="protein", replicates={"dG_bar": [10.0, math.nan, 14.0]})
        assert leg.mean() == pytest.approx(12.0)
        assert leg.n_replicates == 2


class TestFoldingFreeEnergy:
    @staticmethod
    def _result(protein, tripeptide):
        record = MutationResult("LEU39ALA", "LEU", 39, "ALA")
        if protein is not None:
            record.legs["protein"] = LegResult("protein", {"dG_bar": protein})
        if tripeptide is not None:
            record.legs["tripeptide"] = LegResult("tripeptide", {"dG_bar": tripeptide})
        return record

    def test_ddG_is_the_difference_between_the_legs(self):
        record = self._result([10.0, 12.0], [4.0, 6.0])
        assert record.ddG() == pytest.approx(6.0)

    def test_error_is_propagated_from_both_legs(self):
        record = self._result([10.0, 12.0], [4.0, 6.0])
        expected = math.sqrt(record.legs["protein"].sem() ** 2 + record.legs["tripeptide"].sem() ** 2)
        assert record.ddG_sem() == pytest.approx(expected)

    @pytest.mark.parametrize("protein,tripeptide", [(None, [4.0, 6.0]), ([10.0, 12.0], None), (None, None)])
    def test_a_missing_leg_gives_no_answer(self, protein, tripeptide):
        """One leg is not a folding free energy, so it must not be reported as one."""
        record = self._result(protein, tripeptide)
        assert math.isnan(record.ddG())
        assert math.isnan(record.ddG_sem())


class TestResultDiscovery:
    """Mutations are found from directory names -- there is no mapping file."""

    @staticmethod
    def _make_leg(root: Path, mutation: str, dG_per_stage: float, replicates: int = 2):
        qfep_text = _qfep_output(dG_per_stage)
        for stage in (1, 2):
            for replicate in range(1, replicates + 1):
                run = root / f"FEP_{mutation}" / f"FEP{stage}" / "298" / str(replicate)
                run.mkdir(parents=True)
                (run / "qfep.out").write_text(qfep_text)

    def test_both_legs_are_matched_by_name(self, tmp_path):
        protein, tripeptide = tmp_path / "protein", tmp_path / "tripeptide"
        self._make_leg(protein, "LEU39ALA", 5.0)
        self._make_leg(tripeptide, "LEU39ALA", 2.0)

        results = collect(protein, tripeptide, "298")
        assert len(results) == 1
        record = results[0]
        assert (record.wild_type, record.position, record.mutant) == ("LEU", 39, "ALA")
        # Two stages summed per leg: 10 - 4.
        assert record.ddG() == pytest.approx(6.0)

    def test_results_are_ordered_by_position(self, tmp_path):
        protein, tripeptide = tmp_path / "protein", tmp_path / "tripeptide"
        for mutation in ("TYR25GLY", "LEU39ALA", "ILE3THR"):
            self._make_leg(protein, mutation, 5.0)
            self._make_leg(tripeptide, mutation, 2.0)
        assert [r.position for r in collect(protein, tripeptide, "298")] == [3, 25, 39]

    def test_a_mutation_with_only_one_leg_is_still_listed(self, tmp_path):
        """An incomplete set has to be visible, not silently dropped."""
        protein, tripeptide = tmp_path / "protein", tmp_path / "tripeptide"
        self._make_leg(protein, "LEU39ALA", 5.0)
        tripeptide.mkdir()

        frame = to_frame(collect(protein, tripeptide, "298"))
        assert list(frame["mutation"]) == ["LEU39ALA"]
        assert math.isnan(frame.loc[0, "ddG_fold"])
        assert frame.loc[0, "n_protein"] == 2

    def test_directories_that_are_not_mutations_are_skipped(self, tmp_path):
        protein, tripeptide = tmp_path / "protein", tmp_path / "tripeptide"
        self._make_leg(protein, "LEU39ALA", 5.0)
        (protein / "FEP_scratch").mkdir()
        tripeptide.mkdir()
        assert [r.mutation for r in collect(protein, tripeptide, "298")] == ["LEU39ALA"]

    def test_a_replicate_missing_a_stage_is_recorded_as_a_crash(self, tmp_path):
        """Half a transformation is not a free energy."""
        protein, tripeptide = tmp_path / "protein", tmp_path / "tripeptide"
        self._make_leg(protein, "LEU39ALA", 5.0, replicates=2)
        (protein / "FEP_LEU39ALA" / "FEP2" / "298" / "2" / "qfep.out").unlink()
        tripeptide.mkdir()

        record = collect(protein, tripeptide, "298")[0]
        assert record.legs["protein"].n_replicates == 1
        assert len(record.legs["protein"].crashed) == 1


def _qfep_output(dG: float) -> str:
    """A qfep.out carrying `dG`, with the section headers qfep really writes.

    The reader keys on the fourth whitespace-separated token of each header line,
    so the headers have to be spelled as qfep spells them (including its
    "Termodynamic").
    """
    return (
        "# Part 1: Free energy perturbation summary:\n"
        "# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>\n"
        f"   1.000000    0.000    0.000   -0.375 {-dG:9.3f}    0.000\n"
        f"   0.000000    0.388 {dG:9.3f}    0.000    0.000 {dG:9.3f}\n"
        "# Part 2: Reaction free energy summary:\n"
        "# Part 6: BAR Bennet:\n"
        "   1.000000    0.000    0.000\n"
        f"   0.000000    0.388 {dG:9.3f}\n"
    )


# ----------------------------------------------------------------------
# Run-script safety (does not require compiled Q binaries)
# ----------------------------------------------------------------------


class TestRunScriptSafety:
    @pytest.fixture(scope="class")
    def script(self, tmp_path_factory):
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="HABROK",
            windows=2,
            temperature="298",
            replicates=1,
            seeds=[971],
            to_clean=["en", "re", "inp", "top", "dcd", "log"],
            write_trajectories=False,
            separate_scaling=False,
            workdir=tmp_path_factory.mktemp("run-script"),
        )
        run.inputfiles.mkdir(parents=True)
        run.write_runfile(
            [
                ["md_1_1000_0000", "md_1_0000_1000"],
                ["md_2_1000_0000", "md_2_0000_1000"],
            ]
        )
        return run.inputfiles / "runHABROK.sh"

    def test_generated_script_is_valid_bash(self, script):
        subprocess.run(["bash", "-n", script], check=True)

    def test_eq_logs_are_recorded_and_removed_before_production(self, script):
        content = script.read_text()
        status = content.index("> equilibration.status")
        cleanup = content.index('rm -f -- "${eq_logs[@]}"', status)
        production = content.index("echo md_1_", cleanup)
        assert status < cleanup < production

    def test_qfep_is_validated_before_stage_cleanup(self, script):
        content = script.read_text()
        production_validation = content.index('"FEP$stage replicate $run_num')
        qfep = content.index("timeout 3m", production_validation)
        part_three = content.index("qfep Part 3 is incomplete", qfep)
        status = content.index("> stage_validation.status", part_three)
        cleanup = content.index("rm -f -- *en *re *inp *top *dcd *log", status)
        loop_end = content.index("\ndone", cleanup)
        assert production_validation < qfep < part_three < status < cleanup < loop_end


# ----------------------------------------------------------------------
# End-to-end setup
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def prepared_t4l(tmp_path_factory):
    """Prepare a T4L sphere centred on LEU39, using the real qprep binary."""
    if not (BIN / "qprep").exists():
        pytest.skip("qprep binary not compiled")
    structure = T4L_TUTORIAL / "2LZM_prep.pdb"
    if not structure.exists():
        pytest.skip(f"tutorial structure not found: {structure}")

    work = tmp_path_factory.mktemp("t4l")
    shutil.copy(structure, work / "protein.pdb")
    shutil.copy(T4L_RESOURCES / "ALA39.pdb", work / "ALA39.pdb")

    result = subprocess.run(
        ["qprep_prot", "-i", "protein.pdb", "-FF", "OPLSAAM", "-r", "25", "-cog", *LEU39_CB],
        cwd=work,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, f"qprep_prot failed:\n{result.stderr}"
    return work


@qprep_required
class TestPreparationRecord:
    def test_prep_json_is_written(self, prepared_t4l):
        assert (prepared_t4l / "prep.json").exists()

    def test_it_records_the_sphere(self, prepared_t4l):
        prep = SpherePrep.read(prepared_t4l)
        assert prep.radius == 25.0
        assert prep.center == pytest.approx([float(c) for c in LEU39_CB])
        assert prep.force_field == "OPLSAAM"

    def test_t4l_has_all_its_residues_and_no_disulfides(self, prepared_t4l):
        """2LZM has 164 residues and two free cysteines that are not bridged."""
        prep = SpherePrep.read(prepared_t4l)
        assert len(prep.residues) == 164
        assert prep.disulfides == []

    def test_residue_numbering_agrees_with_the_topology_qprep_wrote(self, prepared_t4l):
        """top_p.pdb carries Q's own numbering, so it is the authority."""
        from QligFEP.sphere_prep import verify_residue_map

        prep = SpherePrep.read(prepared_t4l)
        assert verify_residue_map(prep.residues, prepared_t4l / "top_p.pdb")

    def test_leu39_is_found_under_its_pdb_number(self, prepared_t4l):
        prep = SpherePrep.read(prepared_t4l)
        assert prep.residue_name(39, "A") == "LEU"


@qprep_required
class TestProteinLegSetup:
    @pytest.fixture(scope="class")
    def fep_dir(self, prepared_t4l):
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            windows=25,
            sampling="exponential",
            to_clean=["inp", "re", "top", "dcd"],
            write_trajectories=False,
            workdir=prepared_t4l,
        )
        return run.run() / "inputfiles"

    def test_the_expected_files_are_written(self, fep_dir):
        for name in (
            "L2A.lib",
            "FEP1.fep",
            "FEP2.fep",
            "complex.pdb",
            "dualtop.top",
            "top_p.pdb",
            "qfep.inp",
            "qprep.inp",
            "runSNELLIUS.sh",
            "OPLSAAM_merged.prm",
            "eq1.inp",
            "eq5.inp",
        ):
            assert (fep_dir / name).exists(), f"{name} was not written"

    def test_two_stages_of_production_files(self, fep_dir):
        for stage in (1, 2):
            assert len(list(fep_dir.glob(f"md_{stage}_*.inp"))) == 25

    def test_qprep_accepted_the_topology(self, fep_dir):
        """A hybrid residue that does not match its library entry is rejected here."""
        assert "ERROR" not in (fep_dir / "qprep.out").read_text()

    def test_the_solute_restraint_spans_the_topology_not_the_input_pdb(self, fep_dir):
        """qprep builds atoms the input PDB does not carry -- the proton of a
        residue prepared in a neutral form, a missing N-terminal hydrogen. A
        restraint sized from the PDB leaves the last of the solute free while the
        rest of the protein is held: 2626 atoms restrained out of 2638 for T4L."""
        solute = sum(
            1
            for line in (fep_dir / "top_p.pdb").read_text().splitlines()
            if line.startswith("ATOM") and line[17:20].strip() not in ("HOH", "SOL", "WAT")
        )
        assert solute > _atom_records(fep_dir / "complex.pdb"), "qprep added no atoms here"

        restraint = _sequence_restraint(fep_dir / "eq1.inp")
        assert restraint == (1, solute)

    def test_hybrid_library_holds_both_side_chains(self, fep_dir):
        library = (fep_dir / "L2A.lib").read_text()
        assert "{L2A}" in library
        assert "CB" in library  # wild-type leucine
        assert "cb" in library  # mutant alanine, lower-cased

    def test_hybrid_library_bonds_the_mutant_to_the_shared_backbone(self, fep_dir):
        assert "CA    cb" in (fep_dir / "L2A.lib").read_text()

    def test_stage_one_discharges_the_wild_type(self, fep_dir):
        """CG belongs to leucine only, so stage 1 must take it to zero charge."""
        charges = _fep_section(fep_dir / "FEP1.fep", "change_charges")
        atoms = _fep_atom_names(fep_dir / "FEP1.fep")
        index = atoms.index("CG")
        assert float(charges[index].split()[2]) == 0.0

    def test_stage_one_keeps_the_mutant_uncharged(self, fep_dir):
        charges = _fep_section(fep_dir / "FEP1.fep", "change_charges")
        for line, name in zip(charges, _library_atom_names(fep_dir / "L2A.lib")):
            if name.islower():
                assert [float(v) for v in line.split()[1:3]] == [0.0, 0.0]

    def test_stage_two_grows_the_mutant_in(self, fep_dir):
        """The mutant side chain starts as a dummy type and ends as its real type."""
        types = _fep_section(fep_dir / "FEP2.fep", "change_atoms")
        names = _library_atom_names(fep_dir / "L2A.lib")
        mutant = [line for line, name in zip(types, names) if name.islower()]
        assert mutant and all(line.split()[1] != "DUM" for line in mutant)

    def test_stage_two_turns_the_wild_type_into_a_dummy(self, fep_dir):
        types = _fep_section(fep_dir / "FEP2.fep", "change_atoms")
        names = _library_atom_names(fep_dir / "L2A.lib")
        backbone = amino_acids.backbone_atoms(False, False)
        wild_type = [line for line, name in zip(types, names) if not name.islower() and name not in backbone]
        assert wild_type and all(line.split()[2] == "DUM" for line in wild_type)

    def test_the_backbone_is_never_softcored(self, fep_dir):
        softcore = _fep_section(fep_dir / "FEP1.fep", "softcore")
        names = _library_atom_names(fep_dir / "L2A.lib")
        backbone = amino_acids.backbone_atoms(False, False)
        for line, name in zip(softcore, names):
            if name in backbone:
                assert [v for v in line.split()[1:3]] == ["0", "0"], name

    def test_the_two_side_chains_cannot_see_each_other(self, fep_dir):
        """Both occupy the pocket at once, so every cross pair must be excluded."""
        excluded = _fep_section(fep_dir / "FEP1.fep", "excluded_pairs")
        names = _fep_atom_names(fep_dir / "FEP1.fep")
        backbone = amino_acids.backbone_atoms(False, False)
        wild_type = [n for n in names if not n.islower() and n not in backbone]
        mutant = [n for n in names if n.islower()]
        assert (len(wild_type), len(mutant)) == (13, 4)  # leucine vs alanine
        assert len(excluded) == len(wild_type) * len(mutant)

    def test_chi1_torsions_across_the_topologies_are_zeroed(self, fep_dir):
        torsions = _fep_section(fep_dir / "FEP1.fep", "change_torsions")
        assert torsions
        for line in torsions:
            assert line.split()[-2:] == ["0", "0"]

    def test_the_run_script_chains_the_stages(self, fep_dir):
        script = (fep_dir / "runSNELLIUS.sh").read_text()
        assert "fepfiles=(FEP1.fep FEP2.fep)" in script
        # Stage 2 restarts from stage 1's last window.
        assert "restartfile=md_1_0000_1000.re" in script

    def test_each_stage_cleans_its_own_files(self, fep_dir):
        script = (fep_dir / "runSNELLIUS.sh").read_text()
        cleanup_line = "rm -f -- *inp *re *top *dcd"
        cleanup = script.index(cleanup_line)
        stage_loop_end = script.index("\ndone", cleanup)

        assert cleanup < stage_loop_end
        assert script.count(cleanup_line) == 1

    def test_stage_one_restart_survives_until_stage_two(self, fep_dir):
        script = (fep_dir / "runSNELLIUS.sh").read_text()
        bridge_copy = script.index('cp -- "$restartfile" "$bridge_restart"')
        cleanup = script.index("rm -f -- *inp *re *top *dcd", bridge_copy)
        bridge_move = script.index('mv -- "$bridge_restart" "$rundir/eq5.re"')

        assert bridge_copy < cleanup
        assert bridge_move < bridge_copy  # text lives in the earlier stage-selection branch

    def test_generated_inputs_do_not_request_trajectories(self, fep_dir):
        for input_file in fep_dir.glob("*.inp"):
            content = input_file.read_text()
            if "[intervals]" not in content:
                continue
            intervals = content.split("[intervals]", 1)[1].split("[files]", 1)[0]
            files = content.split("[files]", 1)[1].split("[", 1)[0]
            assert "trajectory                0" in intervals
            assert "trajectory" not in files

    def test_the_protein_leg_has_no_anchor_restraint(self, fep_dir):
        """The protein holds the residue in place; an anchor would be redundant."""
        content = (fep_dir / "eq5.inp").read_text()
        section = content.split("[sequence_restraints]")[1].split("[")[0]
        assert not section.strip()

    def test_distance_restraints_hold_the_side_chains_together(self, fep_dir):
        content = (fep_dir / "eq5.inp").read_text()
        section = content.split("[distance_restraints]")[1].split("[")[0]
        assert section.strip(), "no distance restraints were written"


@qprep_required
class TestReferencePeptideSetup:
    @pytest.fixture(scope="class")
    def fep_dir(self, prepared_t4l, tmp_path_factory):
        work = tmp_path_factory.mktemp("t4l_tripeptide")
        for name in ("protein_processed.pdb", "water.pdb", "prep.json", "ALA39.pdb"):
            shutil.copy(prepared_t4l / name, work / name)
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="tripeptide",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            windows=25,
            tripeptide_flanks="A",
            workdir=work,
        )
        return run.run() / "inputfiles"

    def test_the_peptide_is_capped_and_flanked(self, fep_dir):
        """ACE-ALA-<hybrid>-ALA-NMA: the mutated residue with Ala flanks and caps."""
        residues = []
        for line in (fep_dir / "complex.pdb").read_text().splitlines():
            if line.startswith("ATOM"):
                name = line[17:20].strip()
                if not residues or residues[-1] != name:
                    residues.append(name)
        assert residues == ["ACE", "ALA", "L2A", "ALA", "NMA"]

    def test_the_cap_grows_off_the_backbone_not_the_side_chain(self, fep_dir):
        """The C-terminal cap attaches to the carbonyl carbon; anchoring it on CA or
        CB puts the rest of the peptide on the side chain."""
        atoms = {}
        for line in (fep_dir / "complex.pdb").read_text().splitlines():
            if line.startswith("ATOM"):
                atoms[(int(line[22:26]), line[12:16].strip())] = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
        carbonyl = atoms[(3, "C")]
        cap_nitrogen = atoms[(4, "N")]
        beta_carbon = atoms[(3, "CB")]
        assert math.dist(carbonyl, cap_nitrogen) < 1.6  # a real peptide bond
        assert math.dist(beta_carbon, cap_nitrogen) > 2.0

    def test_the_peptide_is_anchored_in_its_sphere(self, fep_dir):
        """Without protein around it, the peptide would drift to the boundary."""
        content = (fep_dir / "eq5.inp").read_text()
        section = content.split("[sequence_restraints]")[1].split("[")[0].strip()
        assert section, "the reference peptide has no anchor"
        first, last = section.split()[:2]
        assert first == last, "the anchor should pin a single atom"

    def test_qprep_accepted_the_topology(self, fep_dir):
        assert "ERROR" not in (fep_dir / "qprep.out").read_text()

    def test_the_same_hybrid_residue_is_perturbed_as_in_the_protein_leg(self, fep_dir, prepared_t4l):
        """Both legs perturb the same residue, so the atoms and their charge and type
        changes must agree. The topology *indices* differ, because a capped peptide is
        a far smaller system than the protein."""
        protein_fep = prepared_t4l / "FEP_LEU39ALA" / "inputfiles" / "FEP1.fep"
        if not protein_fep.exists():
            pytest.skip("protein leg was not set up in this session")

        assert _fep_atom_names(fep_dir / "FEP1.fep") == _fep_atom_names(protein_fep)
        for section in ("change_charges", "change_atoms", "softcore"):
            assert _fep_section(fep_dir / "FEP1.fep", section) == _fep_section(protein_fep, section), section


@qprep_required
class TestResidueMustBeInTheSphere:
    """A residue beyond the sphere is not simulated, whatever its name.

    This is why a set of mutations needs one prepared sphere per mutation: the
    sphere has to enclose the residue being mutated, and each centre also decides
    which other charges get neutralised.
    """

    def test_a_distant_neutral_residue_is_refused(self, prepared_t4l):
        """The case a name-based check cannot see: nothing renames a far-away PHE,
        so only geometry catches it."""
        from QligFEP.sphere_prep import residue_distance_from_center

        prep = SpherePrep.read(prepared_t4l)
        distant = [
            r
            for r in prep.residues
            if residue_distance_from_center(prepared_t4l / prep.prepared_pdb, r.q_number, prep.center)
            > prep.radius
            and r.name not in ("GLH", "ASH", "LYN", "ARN")
        ]
        if not distant:
            pytest.skip("no neutral residue lies outside this sphere")
        residue = distant[0]

        run = QresFEP(
            mutation=f"{residue.name}{residue.pdb_number}ALA",
            chain=residue.chain,
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=prepared_t4l,
        )
        with pytest.raises(MutationError, match="not part of the simulated system"):
            run.read_prep()

    def test_the_residue_at_the_centre_is_accepted(self, prepared_t4l):
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=prepared_t4l,
        )
        run.read_prep()  # must not raise
        assert run.q_position == 39

    def test_a_residue_in_the_boundary_shell_warns_but_proceeds(self, prepared_t4l):
        """It is simulated, so it is allowed -- but its environment is the boundary."""
        from QligFEP.logger import logger
        from QligFEP.sphere_prep import residue_distance_from_center

        prep = SpherePrep.read(prepared_t4l)
        in_shell = [
            r
            for r in prep.residues
            if prep.boundary_radius
            < residue_distance_from_center(prepared_t4l / prep.prepared_pdb, r.q_number, prep.center)
            <= prep.radius
            and r.name in amino_acids.SIDE_CHAINS
        ]
        if not in_shell:
            pytest.skip("no residue lies in the boundary shell of this sphere")
        residue = in_shell[0]

        run = QresFEP(
            mutation=f"{residue.name}{residue.pdb_number}ALA",
            chain=residue.chain,
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=prepared_t4l,
        )
        # loguru does not feed pytest's caplog, so capture it with a sink.
        messages: list[str] = []
        sink = logger.add(messages.append, level="WARNING")
        try:
            run.read_prep()
        except FileNotFoundError:
            pass  # the mutant PDB is not needed to reach the warning
        finally:
            logger.remove(sink)

        assert any("restrained boundary shell" in m for m in messages), messages

    def test_the_boundary_radius_comes_from_the_preparation(self, prepared_t4l):
        prep = SpherePrep.read(prepared_t4l)
        assert prep.boundary_radius == prep.radius - prep.neutralization_offset
        assert prep.neutralization_offset == 3.0  # qprep_prot's default


@qprep_required
class TestSetupGuards:
    def test_every_neutralised_residue_is_rejected_with_the_right_reason(self, prepared_t4l):
        """Naming the charged form of a neutralised residue must never proceed.

        Which message it gets depends on why: beyond the sphere it is not
        simulated at all, inside it the charge was merely stripped. Both refuse,
        but only the second has a same-sphere alternative.
        """
        from QligFEP.sphere_prep import residue_distance_from_center

        prep = SpherePrep.read(prepared_t4l)
        charged_form = {"GLH": "GLU", "ASH": "ASP", "LYN": "LYS", "ARN": "ARG"}
        neutralised = [r for r in prep.residues if r.name in charged_form]
        if not neutralised:
            pytest.skip("no residue was neutralised in this sphere")

        seen_outside = seen_inside = False
        for residue in neutralised:
            run = QresFEP(
                mutation=f"{charged_form[residue.name]}{residue.pdb_number}ALA",
                chain=residue.chain,
                system="protein",
                force_field="OPLSAAM",
                cluster="SNELLIUS",
                workdir=prepared_t4l,
            )
            with pytest.raises(MutationError) as raised:
                run.read_prep()

            distance = residue_distance_from_center(
                prepared_t4l / prep.prepared_pdb, residue.q_number, prep.center
            )
            message = str(raised.value)
            if distance > prep.radius:
                assert "not part of the simulated system" in message
                seen_outside = True
            else:
                assert f"is {residue.name}" in message
                assert "different perturbation" in message
                seen_inside = True

        # T4L with a 25 A sphere on residue 39 exercises both branches.
        assert seen_outside, "expected at least one neutralised residue beyond the sphere"
        assert seen_inside, "expected at least one neutralised residue inside the sphere"

    def test_a_missing_mutant_pdb_is_reported(self, prepared_t4l):
        run = QresFEP(
            mutation="LEU39TRP",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=prepared_t4l,
        )
        with pytest.raises(FileNotFoundError, match="TRP39.pdb"):
            run.read_prep()

    def test_local_execution_is_refused_rather_than_written_wrong(self, prepared_t4l):
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="LOCAL",
            workdir=prepared_t4l,
        )
        with pytest.raises(MutationError, match="no local run script"):
            run.write_runfile([["md_1_1000_0000"], ["md_2_1000_0000"]])


@qprep_required
class TestConfigRecord:
    def test_the_setup_is_recorded(self, prepared_t4l, tmp_path_factory):
        """resfep_config.json is what the analysis reads the run direction from."""
        from QligFEP.CLI import qresfep_cli

        work = tmp_path_factory.mktemp("t4l_config")
        for name in ("protein_processed.pdb", "water.pdb", "prep.json", "ALA39.pdb"):
            shutil.copy(prepared_t4l / name, work / name)

        directory = qresfep_cli.main(
            mutation="LEU39ALA",
            chain="A",
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            windows=25,
            workdir=work,
        )
        config = json.loads((directory / "inputfiles" / "resfep_config.json").read_text())
        assert config["system"] == "protein"
        assert config["q_position"] == 39
        assert config["hybrid_residue"] == "L2A"
        assert config["start"] == "1"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _fep_section(fep_file: Path, section: str) -> list[str]:
    """Return the non-empty lines of one section of a FEP file."""
    lines = []
    inside = False
    for line in fep_file.read_text().splitlines():
        stripped = line.strip()
        if stripped == f"[{section}]":
            inside = True
            continue
        if inside and stripped.startswith("["):
            break
        if inside and stripped:
            lines.append(stripped)
    return lines


def _fep_atom_names(fep_file: Path) -> list[str]:
    """Return the atom names commented into a FEP file's [atoms] section."""
    return [line.split("!")[1].strip() for line in _fep_section(fep_file, "atoms")]


class TestReferencePeptideGeometry:
    """The capped peptide is built from a stored fragment, not by a molecular editor.

    What matters is that the junctions are real peptide bonds, the built residues
    keep L chirality, and the anchoring residue's carbonyl and amide hydrogen are
    moved to match the idealised conformation the caps are built in.
    """

    #: An alanine lifted out of T4 lysozyme, as (name, x, y, z).
    ALANINE = [
        ("N", 38.865, 31.620, 20.999),
        ("CA", 39.713, 30.684, 21.730),
        ("C", 39.089, 30.259, 23.055),
        ("O", 39.611, 29.376, 23.738),
        ("CB", 41.089, 31.295, 21.987),
        ("H", 38.021, 31.170, 20.669),
        ("HA", 39.999, 29.799, 21.160),
    ]

    def residue(self, name="ALA", atoms=None):
        records = []
        for atom, x, y, z in atoms or self.ALANINE:
            record = ["ATOM  ", 0, atom, " ", name, "A", 39, "", x, y, z, 1.0, 0.0, "  ", "  "]
            records.append(record)
        return records

    @staticmethod
    def _by_name(peptide):
        return {(atom[6], atom[2]): np.array(atom[8:11]) for atom in peptide}

    @pytest.mark.parametrize(
        "flanks,expected",
        [
            ("A", ["ACE", "ALA", "ALA", "ALA", "NMA"]),
            ("G", ["ACE", "GLY", "ALA", "GLY", "NMA"]),
            ("X", ["ACE", "ALA", "NMA"]),
        ],
    )
    def test_the_peptide_is_flanked_and_capped_as_asked(self, flanks, expected):
        peptide = build_reference_peptide([self.residue()], flanks)
        residues = []
        for atom in peptide:
            if not residues or residues[-1][0] != atom[6]:
                residues.append((atom[6], atom[4]))
        assert [name for _, name in residues] == expected
        assert [number for number, _ in residues] == sorted(number for number, _ in residues)

    def test_the_junctions_are_peptide_bonds(self):
        """1.29 A into the nitrogen, 1.47 A out of the carbonyl carbon -- the
        lengths the capping builds with."""
        atoms = self._by_name(build_reference_peptide([self.residue()], "A"))
        assert np.linalg.norm(atoms[(2, "C")] - atoms[(3, "N")]) == pytest.approx(1.29, abs=0.02)
        assert np.linalg.norm(atoms[(3, "C")] - atoms[(4, "N")]) == pytest.approx(1.47, abs=0.02)

    def test_the_built_residues_keep_l_chirality(self):
        """Kabsch superposition is free to return a reflection; a mirrored flank
        would be a D-amino acid the force field has no parameters for."""
        atoms = self._by_name(build_reference_peptide([self.residue()], "A"))
        for residue in (2, 4):
            n, ca, c, cb = (atoms[(residue, name)] for name in ("N", "CA", "C", "CB"))
            assert np.dot(np.cross(n - ca, c - ca), cb - ca) > 0

    def test_the_anchoring_carbonyl_and_amide_hydrogen_are_moved(self):
        """Capping sets the residue's own psi and phi, which swings O and H."""
        original = {name: np.array([x, y, z]) for name, x, y, z in self.ALANINE}
        atoms = self._by_name(build_reference_peptide([self.residue()], "A"))
        assert np.linalg.norm(atoms[(3, "O")] - original["O"]) > 0.5
        assert np.linalg.norm(atoms[(3, "H")] - original["H"]) > 0.1
        # The frame itself must not move: it is what everything was placed against.
        for name in ("N", "CA", "C", "CB"):
            assert atoms[(3, name)] == pytest.approx(original[name])

    def test_native_flanks_cap_the_outer_residues(self):
        """With -t Z the neighbours come from the protein, so the caps go on them."""
        peptide = build_reference_peptide([self.residue()] * 3, "Z")
        names = {atom[4] for atom in peptide}
        assert names == {"ACE", "ALA", "NMA"}
        assert sum(1 for atom in peptide if atom[4] == "ALA") == 3 * len(self.ALANINE)

    def test_an_unknown_flank_choice_is_refused(self):
        with pytest.raises(PeptideBuildError, match="Unknown flanking"):
            build_reference_peptide([self.residue()], "Q")

    def test_a_residue_without_a_backbone_cannot_be_capped(self):
        side_chain_only = self.residue(atoms=[("CB", 41.089, 31.295, 21.987)])
        with pytest.raises(PeptideBuildError, match="no N, CA, C atom"):
            build_reference_peptide([side_chain_only], "A")

    def test_proline_keeps_its_nitrogen_free(self):
        """Proline has no amide hydrogen; putting one there is not a residue any
        library has an entry for."""
        proline = [a for a in self.ALANINE if a[0] != "H"]
        peptide = build_reference_peptide([self.residue("PRO", proline)], "X")
        assert not any(atom[2] == "H" and atom[4] == "PRO" for atom in peptide)


class TestSoluteSizeComesFromTheTopology:
    """What eq1--eq4 restrain has to be what qprep built, not what it was given."""

    def _run(self, tmp_path, topology: str, counter_ions: int) -> QresFEP:
        run = QresFEP(
            mutation="LEU39ALA",
            chain="A",
            system="tripeptide",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=tmp_path,
        )
        run.inputfiles = tmp_path
        run.radius = 25.0
        (tmp_path / "top_p.pdb").write_text(topology)
        run.counter_ions = counter_ions
        run.system_size = 0
        run._read_topology_size()
        return run

    @staticmethod
    def _topology(*residues: tuple[str, int]) -> str:
        lines, serial = [], 1
        for name, count in residues:
            for _ in range(count):
                lines.append(f"ATOM  {serial:5d}  X   {name}     1       0.000   0.000   0.000")
                serial += 1
        return "\n".join(lines) + "\n"

    def test_water_is_not_part_of_the_solute(self, tmp_path):
        topology = self._topology(("ALA", 10), ("HOH", 30))
        assert self._run(tmp_path, topology, counter_ions=0).system_size == 10

    def test_counter_ions_are_left_out(self, tmp_path):
        """They sit past the restrained solute and are held by a wall restraint,
        which addresses them as the atoms just past system_size."""
        topology = self._topology(("ALA", 10), ("SOD", 2), ("HOH", 30))
        run = self._run(tmp_path, topology, counter_ions=2)
        assert run.system_size == 10
        assert "11 12" in run._wall_restraints().replace("  ", " ")


def _atom_records(pdb_file: Path) -> int:
    """Return how many atoms a PDB holds."""
    return sum(1 for line in pdb_file.read_text().splitlines() if line.startswith(("ATOM", "HETATM")))


def _sequence_restraint(md_input: Path) -> tuple[int, int]:
    """Return the atom range of an md input's first sequence restraint."""
    lines = md_input.read_text().splitlines()
    start = lines.index("[sequence_restraints]")
    first, last = lines[start + 1].split()[:2]
    return int(first), int(last)


def _library_atom_names(lib_file: Path) -> list[str]:
    """Return the atom names of a hybrid residue library, in order."""
    names = []
    inside = False
    for line in lib_file.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[atoms]"):
            inside = True
            continue
        if inside and stripped.startswith("["):
            break
        if inside and stripped:
            names.append(stripped.split()[1])
    return names
