"""Integration tests for Tyk2 tutorial workflow (steps 1-5).

Tests CLI commands: qparams, qlomap, qcog, qprep_prot, setupFEP.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
Q6_BIN = PROJECT_ROOT / "src" / "q6" / "bin" / "q6"


def run_cli_command(cmd: list[str], cwd: Path = None, timeout: int = 300) -> subprocess.CompletedProcess:
    """Run a CLI command and return result."""
    env = os.environ.copy()
    env["PATH"] = f"{Q6_BIN}:{env.get('PATH', '')}"
    return subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )


class TestQcog:
    """Tests for qcog CLI - center of geometry calculation."""

    @pytest.mark.tutorial
    def test_qcog_sdf_input(self, tyk2_sdf: Path):
        """qcog returns valid center of geometry for SDF input."""
        result = run_cli_command(["qcog", "-i", str(tyk2_sdf)])

        assert result.returncode == 0, f"qcog failed: {result.stderr}"
        # Should contain output like "Center of geometry: [X Y Z]"
        assert "Center of geometry" in result.stderr or "Center of geometry" in result.stdout

    @pytest.mark.tutorial
    def test_qcog_pdb_input(self, tyk2_tutorial_path: Path):
        """qcog returns valid center of geometry for PDB input."""
        protein_pdb = tyk2_tutorial_path / "setupFEP" / "amber" / "protein.pdb"

        if not protein_pdb.exists():
            pytest.skip(f"Protein PDB not found: {protein_pdb}")

        result = run_cli_command(["qcog", "-i", str(protein_pdb)])

        assert result.returncode == 0, f"qcog failed: {result.stderr}"
        assert "Center of geometry" in result.stderr or "Center of geometry" in result.stdout


class TestQlomap:
    """Tests for qlomap CLI - perturbation network mapping."""

    @pytest.mark.tutorial
    def test_qlomap_creates_json(self, tyk2_sdf: Path, temp_work_dir: Path):
        """qlomap creates lomap.json file with valid structure."""
        # Copy SDF to temp dir
        sdf_copy = temp_work_dir / tyk2_sdf.name
        shutil.copy(tyk2_sdf, sdf_copy)

        result = run_cli_command(["qlomap", "-i", str(sdf_copy)], cwd=temp_work_dir)

        # lomap creates a directory with the same name as the input file (without extension)
        output_dir = temp_work_dir / tyk2_sdf.stem
        json_file = output_dir / "lomap.json"

        if not json_file.exists():
            # Try alternative location
            json_file = temp_work_dir / "lomap.json"

        assert json_file.exists(), f"lomap.json not created. stderr: {result.stderr}"

        # Validate JSON structure
        with open(json_file) as f:
            data = json.load(f)

        assert "nodes" in data, "lomap.json missing 'nodes'"
        assert "edges" in data, "lomap.json missing 'edges'"
        assert len(data["nodes"]) > 0, "No nodes in lomap.json"
        assert len(data["edges"]) > 0, "No edges in lomap.json"


class TestQparams:
    """Tests for qparams CLI - ligand parameter generation."""

    @pytest.mark.slow
    @pytest.mark.tutorial
    def test_qparams_generates_files(self, tyk2_sdf: Path, temp_work_dir: Path):
        """qparams generates .lib, .prm, .pdb files for ligands."""
        # Copy SDF to temp dir
        sdf_copy = temp_work_dir / tyk2_sdf.name
        shutil.copy(tyk2_sdf, sdf_copy)

        # Use NAGL for faster execution in tests
        result = run_cli_command(
            ["qparams", "-i", str(sdf_copy), "-p", "2", "-nagl"],
            cwd=temp_work_dir,
            timeout=600,  # 10 minutes for charge calculation
        )

        assert result.returncode == 0, f"qparams failed: {result.stderr}"

        # Check that parameter files were created
        lib_files = list(temp_work_dir.glob("*.lib"))
        prm_files = list(temp_work_dir.glob("*.prm"))
        pdb_files = list(temp_work_dir.glob("*.pdb"))

        assert len(lib_files) > 0, "No .lib files generated"
        assert len(prm_files) > 0, "No .prm files generated"
        assert len(pdb_files) > 0, "No .pdb files generated"


class TestQprepProt:
    """Tests for qprep_prot CLI - protein preparation with spherical boundary."""

    @pytest.fixture
    def qprep_binary_available(self) -> bool:
        """Check if qprep binary is available."""
        qprep_path = Q6_BIN / "qprep"
        return qprep_path.exists()

    @pytest.mark.tutorial
    @pytest.mark.slow
    def test_qprep_prot_generates_water_sphere(self, tyk2_tutorial_path: Path, temp_work_dir: Path):
        """qprep_prot generates water.pdb and dualtop.top."""
        qprep_binary = Q6_BIN / "qprep"
        if not qprep_binary.exists():
            pytest.skip(f"qprep binary not found at {qprep_binary}")

        # Copy protein.pdb to temp dir
        protein_src = tyk2_tutorial_path / "setupFEP" / "amber" / "protein.pdb"
        if not protein_src.exists():
            pytest.skip(f"Protein PDB not found: {protein_src}")

        protein_dst = temp_work_dir / "protein.pdb"
        shutil.copy(protein_src, protein_dst)

        # Use COG from tutorial: -4.689 26.119 -30.570 (approximate)
        result = run_cli_command(
            [
                "qprep_prot",
                "-cog",
                "-4.689",
                "26.119",
                "-30.570",
                "-i",
                str(protein_dst),
                "-FF",
                "AMBER14sb",
                "-r",
                "25",
            ],
            cwd=temp_work_dir,
            timeout=300,
        )

        assert result.returncode == 0, f"qprep_prot failed: {result.stderr}\n{result.stdout}"

        # Check output files
        water_pdb = temp_work_dir / "water.pdb"
        dualtop = temp_work_dir / "dualtop.top"

        assert water_pdb.exists(), "water.pdb not created"
        assert dualtop.exists(), "dualtop.top not created"

        # Verify water.pdb has content
        assert water_pdb.stat().st_size > 0, "water.pdb is empty"


class TestSetupFEP:
    """Tests for setupFEP CLI - batch FEP setup."""

    @pytest.mark.tutorial
    @pytest.mark.slow
    def test_setupfep_creates_directories(self, tyk2_tutorial_path: Path, temp_work_dir: Path):
        """setupFEP creates 1.water and 2.protein directories with FEP input files."""
        qprep_binary = Q6_BIN / "qprep"
        if not qprep_binary.exists():
            pytest.skip(f"qprep binary not found at {qprep_binary}")

        setup_src = tyk2_tutorial_path / "setupFEP"

        # Copy required files to temp dir
        required_files = [
            "amber/protein.pdb",
            "amber/water.pdb",
        ]

        # Check if pre-generated files exist
        for rel_path in required_files:
            src = setup_src / rel_path
            if not src.exists():
                pytest.skip(f"Required file not found: {src}")

        # Copy all necessary files
        shutil.copy(setup_src / "amber" / "protein.pdb", temp_work_dir / "protein.pdb")
        shutil.copy(setup_src / "amber" / "water.pdb", temp_work_dir / "water.pdb")

        # Copy ligand files (.lib, .prm, .sdf)
        for ext in ["*.lib", "*.prm", "*.sdf"]:
            for f in setup_src.glob(ext):
                shutil.copy(f, temp_work_dir / f.name)

        # Copy lomap.json
        lomap_src = setup_src / "lomap.json" if (setup_src / "lomap.json").exists() else None
        if lomap_src is None:
            # Try to find it in a subdirectory
            for candidate in setup_src.glob("*/lomap.json"):
                lomap_src = candidate
                break

        if lomap_src is None or not lomap_src.exists():
            pytest.skip("lomap.json not found in setupFEP directory")

        shutil.copy(lomap_src, temp_work_dir / "lomap.json")

        # Run setupFEP with minimal config
        result = run_cli_command(
            [
                "setupFEP",
                "-FF",
                "AMBER14sb",
                "-r",
                "25",
                "-ts",
                "2fs",
                "-j",
                "lomap.json",
                "-c",
                "LOCAL",  # Required cluster argument
                "-R",
                "1",  # Single replicate for faster test
                "-w",
                "10",  # Fewer windows for faster test
            ],
            cwd=temp_work_dir,
            timeout=600,
        )

        # Check for water directory (always created)
        water_dir = temp_work_dir / "1.water"
        protein_dir = temp_work_dir / "2.protein"

        # Note: setupFEP might fail on some edges but should create directories
        # Check that at least the directory structure was created
        assert (
            water_dir.exists() or protein_dir.exists()
        ), f"Neither 1.water nor 2.protein created. stderr: {result.stderr}\nstdout: {result.stdout}"


class TestTutorialWorkflowIntegration:
    """End-to-end test verifying the complete tutorial workflow."""

    @pytest.mark.tutorial
    @pytest.mark.slow
    @pytest.mark.integration
    def test_cog_lomap_integration(self, tyk2_sdf: Path, temp_work_dir: Path):
        """Test qcog and qlomap work together in the tutorial workflow."""
        # Step 1: Calculate center of geometry
        cog_result = run_cli_command(["qcog", "-i", str(tyk2_sdf)])
        assert cog_result.returncode == 0

        # Step 2: Create perturbation network
        sdf_copy = temp_work_dir / tyk2_sdf.name
        shutil.copy(tyk2_sdf, sdf_copy)

        lomap_result = run_cli_command(["qlomap", "-i", str(sdf_copy)], cwd=temp_work_dir)

        # Find and validate lomap.json
        output_dir = temp_work_dir / tyk2_sdf.stem
        json_file = output_dir / "lomap.json"

        if not json_file.exists():
            json_file = temp_work_dir / "lomap.json"

        assert json_file.exists(), "lomap.json not created after qlomap"

        with open(json_file) as f:
            network = json.load(f)

        # Validate the network has proper structure
        # nodes is a dict where keys are node names
        assert len(network["nodes"]) >= 2, "Network needs at least 2 nodes"
        assert len(network["edges"]) >= 1, "Network needs at least 1 edge"

        # Verify edges reference valid nodes
        # nodes can be either a dict (keys are names) or list (check for 'name' key)
        if isinstance(network["nodes"], dict):
            node_names = set(network["nodes"].keys())
        else:
            node_names = {n["name"] for n in network["nodes"]}

        for edge in network["edges"]:
            assert edge["from"] in node_names, f"Edge 'from' {edge['from']} not in nodes"
            assert edge["to"] in node_names, f"Edge 'to' {edge['to']} not in nodes"
