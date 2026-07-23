"""Golden file tests for FEP input file generation.

Validates that generated .inp files (equilibration and production MD)
match golden reference files. This enables safe refactoring of the
template-based input file generation system.

The golden files were generated using:
    qligfep -l1 'ejm_31' -l2 'ejm_45' -FF AMBER14sb -s protein -c SNELLIUS \
      -R 10 -S sigmoidal -r 25 -l 0.5 -w 100 -T 298 -ts 2fs \
      -rest heavyatom_p -drf 0.5 -log info -b auto -clean dcd -rs 42
"""

import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import (
    GOLDEN_FEP_SETUP_DIR,
    TUTORIALS_DIR,
    InputFileGoldenManager,
    parse_inp_file,
)

GOLDEN_DIR = GOLDEN_FEP_SETUP_DIR / "tyk2_ejm_31_ejm_45"

# Keys to ignore during comparison (values that vary between runs)
IGNORE_KEYS = [
    "MD.random_seed",  # SEED_VAR placeholder
    "MD.temperature",  # T_VAR placeholder
    "files.fep",  # FEP_VAR placeholder
]


@pytest.fixture(scope="module")
def generated_fep_dir(tmp_path_factory) -> Path:
    """Run qligfep with reference parameters and return output dir.

    Uses the pre-prepared ligand files from tutorials/Tyk2/setupFEP/
    to avoid running qparams (which is slow).
    """
    work_dir = tmp_path_factory.mktemp("fep_test")

    # Copy required files from tutorials/Tyk2/setupFEP/
    src_dir = TUTORIALS_DIR / "Tyk2" / "setupFEP"

    # Copy ligand files for ejm_31 and ejm_45
    for lig in ["ejm_31", "ejm_45"]:
        for ext in [".lib", ".prm", ".sdf", ".pdb"]:
            shutil.copy(src_dir / f"{lig}{ext}", work_dir / f"{lig}{ext}")

    # Copy protein.pdb and water.pdb for protein system
    shutil.copy(src_dir / "protein.pdb", work_dir / "protein.pdb")
    shutil.copy(src_dir / "water.pdb", work_dir / "water.pdb")

    # Run qligfep with the same parameters as the golden files
    cmd = [
        "qligfep",
        "-l1",
        "ejm_31",
        "-l2",
        "ejm_45",
        "-FF",
        "AMBER14sb",
        "-s",
        "protein",
        "-c",
        "SNELLIUS",
        "-R",
        "10",
        "-S",
        "sigmoidal",
        "-r",
        "25",
        "-l",
        "0.5",
        "-w",
        "100",
        "-T",
        "298",
        "-ts",
        "2fs",
        "-rest",
        "heavyatom_p",
        "-drf",
        "0.5",
        "-log",
        "info",
        "-b",
        "auto",
        "-clean",
        "dcd",
        "-rs",
        "42",
    ]

    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)  # noqa: S603, S607

    if result.returncode != 0:
        pytest.fail(f"qligfep failed:\nstdout: {result.stdout}\nstderr: {result.stderr}")

    inputfiles_dir = work_dir / "FEP_ejm_31_ejm_45" / "inputfiles"

    if not inputfiles_dir.exists():
        pytest.fail(f"FEP output directory not created: {inputfiles_dir}")

    return inputfiles_dir


class TestFEPInputFileGolden:
    """Golden file tests for FEP input file generation."""

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "inp_file", ["eq1.inp", "eq2.inp", "eq3.inp", "eq4.inp", "eq5.inp", "md_0500_0500.inp"]
    )
    def test_inp_file_matches_golden(self, inp_file: str, generated_fep_dir: Path):
        """Verify generated .inp file matches golden reference."""
        # Load golden file
        golden_path = GOLDEN_DIR / inp_file
        if not golden_path.exists():
            pytest.fail(f"Golden file not found: {golden_path}")

        golden_content = golden_path.read_text()
        golden_parsed = parse_inp_file(golden_content)

        # Load generated file
        actual_path = generated_fep_dir / inp_file
        if not actual_path.exists():
            pytest.fail(f"Generated file not found: {actual_path}")

        actual_content = actual_path.read_text()
        actual_parsed = parse_inp_file(actual_content)

        # Compare using InputFileGoldenManager
        manager = InputFileGoldenManager(GOLDEN_DIR)
        match, message = manager.compare_inp_files(actual_parsed, golden_parsed, ignore_keys=IGNORE_KEYS)

        assert match, f"{inp_file} mismatch: {message}"

    @pytest.mark.slow
    def test_distance_restraints_match_golden(self, generated_fep_dir: Path):
        """Verify distance restraints section matches exactly."""
        golden_content = (GOLDEN_DIR / "md_0500_0500.inp").read_text()
        golden_parsed = parse_inp_file(golden_content)

        actual_content = (generated_fep_dir / "md_0500_0500.inp").read_text()
        actual_parsed = parse_inp_file(actual_content)

        golden_restraints = golden_parsed.get("distance_restraints", [])
        actual_restraints = actual_parsed.get("distance_restraints", [])

        assert len(actual_restraints) == len(
            golden_restraints
        ), f"Distance restraints count mismatch: expected {len(golden_restraints)}, got {len(actual_restraints)}"

        # Compare each restraint line
        for i, (golden_line, actual_line) in enumerate(zip(golden_restraints, actual_restraints)):
            assert (
                actual_line == golden_line
            ), f"Distance restraint {i} mismatch:\n  expected: {golden_line}\n  actual: {actual_line}"

    @pytest.mark.slow
    def test_sequence_restraints_empty_for_protein(self, generated_fep_dir: Path):
        """Verify sequence restraints section is empty for protein system."""
        actual_content = (generated_fep_dir / "md_0500_0500.inp").read_text()
        actual_parsed = parse_inp_file(actual_content)

        seq_restraints = actual_parsed.get("sequence_restraints", [])

        # For protein system, sequence restraints should be empty
        assert (
            len(seq_restraints) == 0
        ), f"Expected empty sequence_restraints for protein system, got {len(seq_restraints)} entries"

    @pytest.mark.slow
    def test_eq1_has_sequence_restraints(self, generated_fep_dir: Path):
        """Verify eq1.inp has sequence restraints (ligand restraints during equilibration)."""
        actual_content = (generated_fep_dir / "eq1.inp").read_text()
        actual_parsed = parse_inp_file(actual_content)

        golden_content = (GOLDEN_DIR / "eq1.inp").read_text()
        golden_parsed = parse_inp_file(golden_content)

        actual_seq = actual_parsed.get("sequence_restraints", [])
        golden_seq = golden_parsed.get("sequence_restraints", [])

        # eq1 should have sequence restraints for ligand position
        assert len(actual_seq) == len(
            golden_seq
        ), f"eq1 sequence restraints count mismatch: expected {len(golden_seq)}, got {len(actual_seq)}"

        if golden_seq:
            assert actual_seq == golden_seq, "eq1 sequence restraints mismatch"

    @pytest.mark.slow
    def test_md_section_settings(self, generated_fep_dir: Path):
        """Verify MD section has correct settings for 2fs timestep."""
        actual_content = (generated_fep_dir / "md_0500_0500.inp").read_text()
        actual_parsed = parse_inp_file(actual_content)

        md_section = actual_parsed.get("MD", {})

        # 2fs timestep should have stepsize=2.0
        assert (
            md_section.get("stepsize") == "2.0"
        ), f"Expected stepsize=2.0 for 2fs, got {md_section.get('stepsize')}"

        # Shake should be on for hydrogens
        assert md_section.get("shake_hydrogens") == "on", "Expected shake_hydrogens=on for 2fs timestep"

    @pytest.mark.slow
    def test_sphere_settings(self, generated_fep_dir: Path):
        """Verify sphere settings match (radius 25)."""
        actual_content = (generated_fep_dir / "md_0500_0500.inp").read_text()
        actual_parsed = parse_inp_file(actual_content)

        sphere_section = actual_parsed.get("sphere", {})

        assert (
            sphere_section.get("shell_radius") == "25"
        ), f"Expected shell_radius=25, got {sphere_section.get('shell_radius')}"
