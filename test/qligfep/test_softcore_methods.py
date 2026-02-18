"""Tests for configurable soft-core methods in FEP file generation.

Verifies that the --softcore-method CLI argument is correctly parsed and
that write_FEP_file() generates FEP files with the appropriate keywords
for the Fortran parser in qatom.f90.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

# Path setup
TEST_DIR = Path(__file__).parent.parent
PROJECT_ROOT = TEST_DIR.parent
TUTORIALS_DIR = PROJECT_ROOT / "tutorials"
SETUP_FEP_DIR = TUTORIALS_DIR / "Tyk2" / "setupFEP"


def parse_fep_section(fep_content: str) -> dict[str, str]:
    """Parse the [FEP] section from a .fep file into key-value pairs."""
    result = {}
    in_section = False

    for line in fep_content.splitlines():
        stripped = line.strip()

        if stripped == "[FEP]":
            in_section = True
            continue
        if in_section and stripped.startswith("["):
            break
        if in_section and stripped and not stripped.startswith("!"):
            parts = stripped.split(None, 1)
            if len(parts) == 2:
                result[parts[0]] = parts[1]

    return result


@pytest.fixture(scope="module")
def fep_work_dir(tmp_path_factory) -> Path:
    """Create a temporary working directory with ligand and protein files."""
    work_dir = tmp_path_factory.mktemp("softcore_test")

    # Copy ligand files for ejm_31 and ejm_45
    for lig in ["ejm_31", "ejm_45"]:
        for ext in [".lib", ".prm", ".sdf", ".pdb"]:
            shutil.copy(SETUP_FEP_DIR / f"{lig}{ext}", work_dir / f"{lig}{ext}")

    # Copy protein/water PDB files
    amber_dir = SETUP_FEP_DIR / "amber"
    shutil.copy(amber_dir / "protein.pdb", work_dir / "protein.pdb")
    shutil.copy(amber_dir / "water.pdb", work_dir / "water.pdb")

    return work_dir


def run_qligfep(work_dir: Path, softcore_method: str | None = None) -> Path:
    """Run qligfep and return path to the FEP1.fep file."""
    cmd = [
        "qligfep",
        "-l1", "ejm_31",
        "-l2", "ejm_45",
        "-FF", "AMBER14sb",
        "-s", "water",
        "-c", "LOCAL",
        "-r", "25",
        "-w", "3",
        "-rs", "42",
        "-log", "error",
    ]
    if softcore_method is not None:
        cmd.extend(["-sc", softcore_method])

    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)  # noqa: S603, S607
    if result.returncode != 0:
        pytest.fail(f"qligfep failed:\nstdout: {result.stdout}\nstderr: {result.stderr}")

    fep_file = work_dir / "FEP_ejm_31_ejm_45" / "inputfiles" / "FEP1.fep"
    if not fep_file.exists():
        pytest.fail(f"FEP file not created: {fep_file}")
    return fep_file


def cleanup_fep_dir(work_dir: Path):
    """Remove generated FEP directory for re-use."""
    fep_dir = work_dir / "FEP_ejm_31_ejm_45"
    if fep_dir.exists():
        shutil.rmtree(fep_dir)


class TestSoftcoreMethodCLI:
    """Test CLI argument validation for --softcore-method."""

    def test_invalid_softcore_method_rejected(self, fep_work_dir):
        """Invalid softcore methods should cause argparse to exit with error."""
        cmd = [
            "qligfep",
            "-l1", "ejm_31",
            "-l2", "ejm_45",
            "-FF", "AMBER14sb",
            "-s", "water",
            "-c", "LOCAL",
            "-sc", "invalid_method",
        ]
        result = subprocess.run(cmd, cwd=fep_work_dir, capture_output=True, text=True)  # noqa: S603, S607
        assert result.returncode != 0, "Expected non-zero exit code for invalid softcore method"
        assert "invalid choice" in result.stderr.lower()


class TestSoftcoreFEPFileGeneration:
    """Test that FEP files contain correct softcore_method keywords."""

    @pytest.mark.slow
    def test_standard_omits_softcore_method(self, fep_work_dir):
        """Default (standard) should NOT write softcore_method to FEP file."""
        try:
            fep_file = run_qligfep(fep_work_dir)
            fep_section = parse_fep_section(fep_file.read_text())

            assert "softcore_method" not in fep_section, (
                "FEP file should NOT contain softcore_method when using standard (default)"
            )
            assert fep_section.get("states") == "2"
            assert fep_section.get("softcore_use_max_potential") == "on"
        finally:
            cleanup_fep_dir(fep_work_dir)

    @pytest.mark.slow
    def test_beutler_coul_in_fep_file(self, fep_work_dir):
        """beutler_coul should write softcore_method to FEP file."""
        try:
            fep_file = run_qligfep(fep_work_dir, softcore_method="beutler_coul")
            fep_section = parse_fep_section(fep_file.read_text())

            assert fep_section.get("softcore_method") == "beutler_coul", (
                f"Expected softcore_method=beutler_coul, got: {fep_section.get('softcore_method')}"
            )
        finally:
            cleanup_fep_dir(fep_work_dir)

    @pytest.mark.slow
    def test_gapsys_in_fep_file(self, fep_work_dir):
        """gapsys should write softcore_method to FEP file."""
        try:
            fep_file = run_qligfep(fep_work_dir, softcore_method="gapsys")
            fep_section = parse_fep_section(fep_file.read_text())

            assert fep_section.get("softcore_method") == "gapsys", (
                f"Expected softcore_method=gapsys, got: {fep_section.get('softcore_method')}"
            )
        finally:
            cleanup_fep_dir(fep_work_dir)

    @pytest.mark.slow
    def test_explicit_standard_omits_softcore_method(self, fep_work_dir):
        """Explicitly passing standard should also omit softcore_method from FEP file."""
        try:
            fep_file = run_qligfep(fep_work_dir, softcore_method="standard")
            fep_section = parse_fep_section(fep_file.read_text())

            assert "softcore_method" not in fep_section, (
                "FEP file should NOT contain softcore_method when explicitly using standard"
            )
        finally:
            cleanup_fep_dir(fep_work_dir)

    @pytest.mark.slow
    def test_softcore_method_in_fep_config_json(self, fep_work_dir):
        """softcore_method should be recorded in fep_config.json."""
        try:
            import json

            run_qligfep(fep_work_dir, softcore_method="beutler_coul")
            config_path = fep_work_dir / "FEP_ejm_31_ejm_45" / "inputfiles" / "fep_config.json"
            config = json.loads(config_path.read_text())

            assert config.get("softcore_method") == "beutler_coul", (
                f"Expected softcore_method=beutler_coul in config, got: {config.get('softcore_method')}"
            )
        finally:
            cleanup_fep_dir(fep_work_dir)


class TestSoftcoreFEPFileFormat:
    """Test that FEP file format is compatible with Fortran parser expectations."""

    @pytest.mark.slow
    def test_fep_file_sections_present(self, fep_work_dir):
        """Verify all required FEP sections are present regardless of softcore method."""
        required_sections = ["FEP", "atoms", "change_charges", "atom_types", "softcore", "change_atoms"]

        for method in [None, "beutler_coul", "gapsys"]:
            try:
                fep_file = run_qligfep(fep_work_dir, softcore_method=method)
                content = fep_file.read_text()

                for section in required_sections:
                    assert f"[{section}]" in content, (
                        f"Missing [{section}] in FEP file with softcore_method={method}"
                    )
            finally:
                cleanup_fep_dir(fep_work_dir)

    @pytest.mark.slow
    def test_softcore_section_has_alpha_values(self, fep_work_dir):
        """Verify [softcore] section has alpha values for all atoms."""
        try:
            fep_file = run_qligfep(fep_work_dir)
            content = fep_file.read_text()

            in_softcore = False
            softcore_lines = []
            for line in content.splitlines():
                stripped = line.strip()
                if stripped == "[softcore]":
                    in_softcore = True
                    continue
                if in_softcore and stripped.startswith("["):
                    break
                if in_softcore and stripped:
                    softcore_lines.append(stripped)

            assert len(softcore_lines) > 0, "No entries in [softcore] section"

            # Each line should have: atom_index alpha_state1 alpha_state2
            for line in softcore_lines:
                parts = line.split()
                assert len(parts) == 3, f"Expected 3 columns in softcore line, got: {line}"
                # Alpha values should be 0 or 20 (disappearing/appearing atoms)
                alpha1, alpha2 = int(parts[1]), int(parts[2])
                assert {alpha1, alpha2} == {0, 20}, f"Unexpected alpha values: {alpha1}, {alpha2}"
        finally:
            cleanup_fep_dir(fep_work_dir)
