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
from QligFEP.qresfep import MutationError, QresFEP, parse_mutation
from QligFEP.sphere_prep import (
    ResidueMapping,
    SpherePrep,
    build_residue_map,
    parse_disulfide_pairs,
)

PROJECT_ROOT = Path(__file__).parent.parent.parent
Q6_BIN = PROJECT_ROOT / "src" / "q6" / "bin" / "q6"
T4L_RESOURCES = Path(__file__).parent / "resources" / "t4l"
T4L_TUTORIAL = PROJECT_ROOT / "tutorials" / "T4L"

#: CB of LEU39 chain A -- the sphere centre the QresFEP tutorial uses.
LEU39_CB = ("41.088", "31.308", "21.828")

pymol_required = pytest.mark.skipif(
    shutil.which("pymol") is None,
    reason="PyMOL builds and caps the reference peptide",
)
qprep_required = pytest.mark.skipif(
    not (Q6_BIN / "qprep").exists(),
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

    def test_glycine_has_no_chi1(self):
        assert amino_acids.chi1_atoms("GLY") == []


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


# ----------------------------------------------------------------------
# Lambda ladders
# ----------------------------------------------------------------------


class TestLambdaLadder:
    @pytest.mark.parametrize(
        "sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"]
    )
    def test_runs_from_one_to_zero(self, sampling):
        ladder = resfep_lambda_ladder(25, sampling)
        assert ladder[0] == "1.000"
        assert ladder[-1] == "0.000"

    @pytest.mark.parametrize(
        "sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"]
    )
    def test_length_is_exactly_the_window_count(self, sampling):
        """One value per window -- unlike the ligand ladder, which returns windows + 1."""
        assert len(resfep_lambda_ladder(25, sampling)) == 25

    @pytest.mark.parametrize(
        "sampling", ["linear", "sigmoidal", "exponential", "reverse_exponential"]
    )
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
        ]
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

REFERENCE_EXAMPLE = (
    Path("/Users/davidararipe/projects/qligfep/tutorials/2.QresFEP_T4L/FEP_example")
)


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
        assert read_qfep_stage(stage_one, reverse_direction=False)["dG_bar"] == pytest.approx(
            65.81, abs=0.01
        )

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

    @pytest.mark.parametrize(
        "protein,tripeptide", [(None, [4.0, 6.0]), ([10.0, 12.0], None), (None, None)]
    )
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
# End-to-end setup
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def prepared_t4l(tmp_path_factory):
    """Prepare a T4L sphere centred on LEU39, using the real qprep binary."""
    if not (Q6_BIN / "qprep").exists():
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
            workdir=prepared_t4l,
        )
        return run.run() / "inputfiles"

    def test_the_expected_files_are_written(self, fep_dir):
        for name in (
            "L2A.lib", "FEP1.fep", "FEP2.fep", "complex.pdb", "dualtop.top",
            "top_p.pdb", "qfep.inp", "qprep.inp", "runSNELLIUS.sh",
            "OPLSAAM_merged.prm", "eq1.inp", "eq5.inp",
        ):
            assert (fep_dir / name).exists(), f"{name} was not written"

    def test_two_stages_of_production_files(self, fep_dir):
        for stage in (1, 2):
            assert len(list(fep_dir.glob(f"md_{stage}_*.inp"))) == 25

    def test_qprep_accepted_the_topology(self, fep_dir):
        """A hybrid residue that does not match its library entry is rejected here."""
        assert "ERROR" not in (fep_dir / "qprep.out").read_text()

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
        wild_type = [
            line for line, name in zip(types, names)
            if not name.islower() and name not in backbone
        ]
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
@pymol_required
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
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
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
            assert _fep_section(fep_dir / "FEP1.fep", section) == _fep_section(
                protein_fep, section
            ), section


@qprep_required
class TestSetupGuards:
    def test_a_neutralised_wild_type_is_rejected(self, prepared_t4l):
        """qprep_prot prepares out-of-sphere GLU as GLH; mutating "GLU" would then
        build the hybrid from the wrong library entry."""
        prep = SpherePrep.read(prepared_t4l)
        neutralised = next(
            (r for r in prep.residues if r.name in ("GLH", "ASH", "LYN", "ARN")), None
        )
        if neutralised is None:
            pytest.skip("no residue was neutralised in this sphere")
        charged = {"GLH": "GLU", "ASH": "ASP", "LYN": "LYS", "ARN": "ARG"}[neutralised.name]

        run = QresFEP(
            mutation=f"{charged}{neutralised.pdb_number}ALA",
            chain=neutralised.chain,
            system="protein",
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=prepared_t4l,
        )
        with pytest.raises(MutationError, match="neutralised"):
            run.read_prep()

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
