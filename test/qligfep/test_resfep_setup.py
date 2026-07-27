"""Tests for the mutation-series setup (``setup_resFEP``).

The end-to-end tests run the real ``qprep`` binary and PyMOL on T4 lysozyme.
Nothing here is mocked.

They exercise the installed ``qprep_prot``/``qresfep`` console scripts, because
that is how the driver runs them; with an editable install that is this checkout.
"""

import shutil
from pathlib import Path

import pytest

from QligFEP.pdb_utils import read_pdb_to_dataframe
from QligFEP.resfep_setup import (
    Mutation,
    MutationSeries,
    SetupError,
    library_residues,
    read_mutations,
    residue_center,
    residue_reach,
    unknown_residues,
    write_mutagenesis_script,
    validate_series,
)
from QligFEP.sphere_prep import SpherePrep

PROJECT_ROOT = Path(__file__).parent.parent.parent
Q6_BIN = PROJECT_ROOT / "src" / "q6" / "bin" / "q6"
T4L_TUTORIAL = PROJECT_ROOT / "tutorials" / "T4L"
T4L_STRUCTURE = T4L_TUTORIAL / "2LZM_prep.pdb"

#: CB of LEU39 chain A -- the sphere centre the QresFEP tutorial uses.
LEU39_CB = [41.088, 31.308, 21.828]

qprep_required = pytest.mark.skipif(
    not (Q6_BIN / "qprep").exists(),
    reason="qprep binary not compiled (cd src/q6 && make all)",
)
pymol_required = pytest.mark.skipif(
    shutil.which("pymol") is None,
    reason="PyMOL builds the mutant side chains",
)
structure_required = pytest.mark.skipif(
    not T4L_STRUCTURE.exists(), reason="T4L tutorial structure not found"
)


@pytest.fixture(scope="module")
def t4l():
    if not T4L_STRUCTURE.exists():
        pytest.skip("T4L tutorial structure not found")
    return read_pdb_to_dataframe(str(T4L_STRUCTURE))


def mutations(*entries: str, chain: str = "A") -> list[Mutation]:
    return [Mutation.from_string(entry, chain) for entry in entries]


# ----------------------------------------------------------------------
# Reading the request
# ----------------------------------------------------------------------


class TestMutation:
    def test_three_letter_codes(self):
        mutation = Mutation.from_string("LEU39ALA", "A")
        assert (mutation.wild_type, mutation.position, mutation.mutant) == ("LEU", 39, "ALA")
        assert mutation.name == "LEU39ALA"
        assert mutation.mutant_pdb_name == "ALA39.pdb"

    def test_one_letter_codes_mean_the_same_mutation(self):
        assert Mutation.from_string("L39A", "A") == Mutation.from_string("LEU39ALA", "A")

    def test_nonsense_is_refused(self):
        with pytest.raises(SetupError):
            Mutation.from_string("not-a-mutation", "A")


class TestMutationList:
    def test_comments_and_blank_lines_are_ignored(self, tmp_path):
        listing = tmp_path / "mutations.txt"
        listing.write_text("# neutral set\nLEU39ALA\n\n  VAL149ALA  # keep\n")
        assert [m.name for m in read_mutations(listing, "A")] == ["LEU39ALA", "VAL149ALA"]

    def test_an_empty_list_is_refused(self, tmp_path):
        listing = tmp_path / "mutations.txt"
        listing.write_text("# nothing here\n")
        with pytest.raises(SetupError, match="lists no mutations"):
            read_mutations(listing, "A")


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------


@structure_required
class TestSphereCentre:
    def test_the_centre_is_the_cb_of_the_mutated_residue(self, t4l):
        """The centre the published protocol uses, and the one the tutorial records."""
        assert residue_center(t4l, 39, "A") == pytest.approx(LEU39_CB)

    def test_glycine_falls_back_to_ha3(self, t4l):
        """Glycine has no CB, so there is nothing else to centre on."""
        centre = residue_center(t4l, 12, "A")
        ha3 = t4l[
            (t4l["residue_seq_number"] == 12) & (t4l["atom_name"].str.strip() == "HA3")
        ]
        assert centre == pytest.approx([float(ha3.iloc[0][axis]) for axis in ("x", "y", "z")])

    def test_a_residue_that_is_not_there_is_refused(self, t4l):
        with pytest.raises(SetupError, match="not in the structure"):
            residue_center(t4l, 9999, "A")

    def test_reach_is_measured_to_the_farthest_atom(self, t4l):
        """The side chain is what is perturbed, so its tip decides whether the
        residue is inside the sphere."""
        assert residue_reach(t4l, 39, "A", LEU39_CB) == pytest.approx(3.5, abs=0.1)
        assert residue_reach(t4l, 149, "A", LEU39_CB) == pytest.approx(34.0, abs=0.1)


# ----------------------------------------------------------------------
# Force field coverage
# ----------------------------------------------------------------------


class TestForceFieldCoverage:
    def test_library_residues_are_read_from_the_lib_file(self, tmp_path):
        library = tmp_path / "cofactor.lib"
        library.write_text("{LIG}   ! a ligand\n    [atoms]\n        1 C1 CA -0.1\n{HEM}\n")
        assert library_residues(library) == {"LIG", "HEM"}

    @structure_required
    def test_a_prepared_structure_is_fully_covered(self, t4l):
        """2LZM_prep.pdb is prepared for OPLSAAM, so nothing should be missing."""
        from QligFEP.IO import parse_lib

        assert unknown_residues(t4l, set(parse_lib("OPLSAAM"))) == {}

    def test_water_is_never_reported(self, tmp_path):
        """qprep adds the sphere's own solvent, so the library need not define it."""
        path = tmp_path / "tiny.pdb"
        path.write_text(
            "ATOM      1  O   HOH A   1      0.000   0.000   0.000  1.00  0.00           O\n"
            "ATOM      2  CA  LIG A   2      1.000   0.000   0.000  1.00  0.00           C\n"
        )
        missing = unknown_residues(read_pdb_to_dataframe(str(path)), {"ALA"})
        assert set(missing) == {"LIG"}
        assert missing["LIG"] == ["A:2"]


# ----------------------------------------------------------------------
# Validating the series
# ----------------------------------------------------------------------


@structure_required
class TestValidation:
    def test_a_sound_series_passes(self, t4l):
        validate_series(t4l, mutations("LEU39ALA", "VAL149ALA"), "OPLSAAM", 25.0)

    def test_the_wild_type_has_to_be_the_residue_that_is_there(self, t4l):
        with pytest.raises(SetupError, match="is LEU, not TYR"):
            validate_series(t4l, mutations("TYR39ALA"), "OPLSAAM", 25.0)

    def test_a_residue_that_is_not_in_the_structure_is_refused(self, t4l):
        with pytest.raises(SetupError, match="not in the structure"):
            validate_series(t4l, mutations("LEU999ALA"), "OPLSAAM", 25.0)

    def test_a_repeated_mutation_is_refused(self, t4l):
        """Both would be written to the same directory, so one would be lost."""
        with pytest.raises(SetupError, match="listed more than once"):
            validate_series(t4l, mutations("LEU39ALA", "L39A"), "OPLSAAM", 25.0)

    def test_every_problem_is_reported_at_once(self, t4l):
        """A series runs unattended; finding out about the last mutation after the
        others have been set up is not worth the time."""
        with pytest.raises(SetupError) as error:
            validate_series(t4l, mutations("TYR39ALA", "LEU999ALA"), "OPLSAAM", 25.0)
        assert "2 problem(s)" in str(error.value)
        assert "TYR39ALA" in str(error.value) and "LEU999ALA" in str(error.value)

    def test_a_mutant_pymol_cannot_build_needs_a_ready_made_pdb(self, t4l):
        """Protonation variants are Q's names, not PyMOL's."""
        with pytest.raises(SetupError, match="mutagenesis wizard cannot build"):
            validate_series(t4l, mutations("LEU39ASH"), "OPLSAAM", 25.0)

    def test_supplied_mutant_pdbs_are_checked_for(self, t4l, tmp_path):
        with pytest.raises(SetupError, match="ALA39.pdb not found"):
            validate_series(
                t4l, mutations("LEU39ALA"), "OPLSAAM", 25.0, mutant_pdb_dir=tmp_path
            )
        (tmp_path / "ALA39.pdb").write_text("")
        validate_series(t4l, mutations("LEU39ALA"), "OPLSAAM", 25.0, mutant_pdb_dir=tmp_path)


@structure_required
class TestFixedCentreIsCheckedAgainstTheSphere:
    """With one centre for the whole series, nothing guarantees a residue is
    simulated -- so every mutation is checked against the sphere up front."""

    def test_a_residue_at_the_centre_passes(self, t4l):
        validate_series(t4l, mutations("LEU39ALA"), "OPLSAAM", 25.0, center=LEU39_CB)

    def test_a_residue_beyond_the_radius_is_refused(self, t4l):
        with pytest.raises(SetupError, match="not in the simulated system"):
            validate_series(t4l, mutations("VAL149ALA"), "OPLSAAM", 25.0, center=LEU39_CB)

    def test_a_residue_reaching_into_the_restrained_shell_is_refused(self, t4l):
        """ILE29 reaches 22.6 A, inside a 25 A sphere but into the 3 A shell where
        `qprep_prot` neutralises charges."""
        with pytest.raises(SetupError, match="restrained shell"):
            validate_series(t4l, mutations("ILE29ALA"), "OPLSAAM", 25.0, center=LEU39_CB)

    def test_the_shell_follows_the_neutralisation_offset(self, t4l):
        """The same residue is fine once nothing is neutralised around it."""
        validate_series(
            t4l,
            mutations("ILE29ALA"),
            "OPLSAAM",
            25.0,
            center=LEU39_CB,
            neutralization_offset=0.0,
        )

    def test_the_refusal_says_what_to_do_instead(self, t4l):
        with pytest.raises(SetupError, match="centre each sphere on its own mutated residue"):
            validate_series(t4l, mutations("VAL149ALA"), "OPLSAAM", 25.0, center=LEU39_CB)

    def test_centring_on_each_residue_never_hits_either_case(self, t4l):
        """The point of the default: the mutation is at the centre by construction."""
        validate_series(t4l, mutations("LEU39ALA", "VAL149ALA"), "OPLSAAM", 25.0)


# ----------------------------------------------------------------------
# Mutant residues
# ----------------------------------------------------------------------


class TestMutagenesisScript:
    def test_one_block_per_mutation_selecting_the_right_chain(self, tmp_path):
        script = write_mutagenesis_script(
            Path("2LZM_prep.pdb"), mutations("LEU39ALA", "VAL149GLY", chain="B"), tmp_path / "m.pml"
        )
        text = script.read_text()
        assert text.count("cmd.wizard('mutagenesis')") == 2
        assert "cmd.get_wizard().set_mode('ALA')" in text
        assert "chain B and resi 39" in text
        assert "save ALA39.pdb, chain B and resi 39" in text

    def test_mutations_that_build_the_same_residue_are_built_once(self, tmp_path):
        """LEU39ALA and ILE39ALA would both write ALA39.pdb."""
        script = write_mutagenesis_script(
            Path("p.pdb"), mutations("LEU39ALA", "LEU39ALA"), tmp_path / "m.pml"
        )
        assert script.read_text().count("cmd.wizard('mutagenesis')") == 1


# ----------------------------------------------------------------------
# The commands the series runs
# ----------------------------------------------------------------------


@structure_required
class TestSeriesCommands:
    def series(self, tmp_path, **kwargs):
        kwargs.setdefault("mutations", mutations("LEU39ALA"))
        return MutationSeries(
            structure=T4L_STRUCTURE,
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            workdir=tmp_path,
            **kwargs,
        )

    def test_each_mutation_is_centred_on_its_own_residue(self, tmp_path):
        series = self.series(tmp_path, mutations=mutations("LEU39ALA", "VAL149ALA"))
        centres = [series.center_for(m) for m in series.mutations]
        assert centres[0] == pytest.approx(LEU39_CB)
        assert centres[0] != centres[1]

    def test_a_fixed_centre_is_used_for_every_mutation(self, tmp_path):
        series = self.series(
            tmp_path, mutations=mutations("LEU39ALA", "ASN40ALA"), center=LEU39_CB
        )
        assert [series.center_for(m) for m in series.mutations] == [LEU39_CB, LEU39_CB]

    def test_the_preparation_carries_the_centre_radius_and_offset(self, tmp_path):
        series = self.series(tmp_path, radius=18.0, neutralization_offset=2.0)
        command = series._qprep_command(LEU39_CB)
        assert command[:3] == ["qprep_prot", "-i", "protein.pdb"]
        assert command[command.index("-r") + 1] == "18.0"
        assert command[command.index("-nbo") + 1] == "2.0"
        assert command[command.index("-cog") + 1 :] == ["41.088", "31.308", "21.828"]

    def test_options_are_passed_through_to_qresfep(self, tmp_path):
        series = self.series(tmp_path, qresfep_options=["-w", "5", "-r", "2"])
        command = series._qresfep_command(series.mutations[0], "tripeptide")
        assert command[:9] == [
            "qresfep", "-m", "LEU39ALA", "-mc", "A", "-S", "tripeptide", "-FF", "OPLSAAM",
        ]
        assert command[-4:] == ["-w", "5", "-r", "2"]

    def test_cofactors_reach_both_steps_in_the_form_each_expects(self, tmp_path):
        series = self.series(tmp_path, cofactors=["NAD"])
        assert "-cof" in series._qprep_command(LEU39_CB)
        assert series._qprep_command(LEU39_CB)[-1] == "NAD.pdb"
        assert series._qresfep_command(series.mutations[0], "protein")[-1] == "NAD"

    def test_a_missing_cofactor_file_stops_the_series(self, tmp_path):
        series = self.series(tmp_path, cofactors=["NAD"])
        with pytest.raises(SetupError, match="Missing cofactor file"):
            series.validate()


# ----------------------------------------------------------------------
# End to end
# ----------------------------------------------------------------------


@qprep_required
@pymol_required
@structure_required
class TestSeriesEndToEnd:
    @pytest.fixture(scope="class")
    def series(self, tmp_path_factory):
        work = tmp_path_factory.mktemp("series")
        series = MutationSeries(
            structure=T4L_STRUCTURE,
            mutations=mutations("LEU39ALA", "VAL149ALA"),
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            radius=25.0,
            workdir=work,
            qresfep_options=["-w", "5", "-r", "2"],
        )
        outcomes = series.run()
        return series, outcomes

    def test_every_mutation_is_set_up(self, series):
        _, outcomes = series
        assert [o.mutation.name for o in outcomes] == ["LEU39ALA", "VAL149ALA"]
        assert all(o.ok for o in outcomes), [o.error for o in outcomes]
        assert all(o.legs == ["protein", "tripeptide"] for o in outcomes)

    def test_both_legs_land_where_the_analysis_looks_for_them(self, series):
        setup, _ = series
        for leg in ("protein", "tripeptide"):
            for name in ("LEU39ALA", "VAL149ALA"):
                assert (setup.workdir / leg / f"FEP_{name}" / "inputfiles" / "FEP1.fep").exists()

    def test_each_mutation_got_a_sphere_of_its_own(self, series):
        """The reason the series is not one prepared sphere: a different centre
        also means a different set of neutralised charges."""
        setup, _ = series
        prepared = [SpherePrep.read(setup.workdir / "work" / name)
                    for name in ("LEU39ALA", "VAL149ALA")]
        assert prepared[0].center == pytest.approx(LEU39_CB)
        assert prepared[0].center != prepared[1].center

        neutralised = [
            {r.pdb_number for r in prep.residues if r.name in ("GLH", "ASH", "LYN", "ARN")}
            for prep in prepared
        ]
        assert neutralised[0] != neutralised[1]

    def test_the_mutated_residue_keeps_its_charge_in_its_own_sphere(self, series):
        """Which is the whole point of re-centring: a neutralised residue is no
        longer the residue the mutation names."""
        setup, _ = series
        prep = SpherePrep.read(setup.workdir / "work" / "LEU39ALA")
        assert prep.residue_name(39, "A") == "LEU"

    def test_the_mutant_residues_were_built(self, series):
        setup, _ = series
        assert (setup.workdir / "mutants" / "ALA39.pdb").exists()
        assert (setup.workdir / "mutants" / "ALA149.pdb").exists()


@qprep_required
@structure_required
class TestFixedCentreSeriesRefusesBeforeWriting:
    def test_nothing_is_written_when_a_mutation_is_out_of_the_sphere(self, tmp_path):
        series = MutationSeries(
            structure=T4L_STRUCTURE,
            mutations=mutations("LEU39ALA", "VAL149ALA"),
            force_field="OPLSAAM",
            cluster="SNELLIUS",
            center=LEU39_CB,
            workdir=tmp_path,
        )
        with pytest.raises(SetupError, match="VAL149ALA"):
            series.run()
        assert not (tmp_path / "work").exists()
        assert not (tmp_path / "protein").exists()
