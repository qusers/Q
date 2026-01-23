"""Golden file regression tests for FEP setup consistency.

Validates that FEP setup output maintains consistent structure across code changes.
Tests validate the "shape" of generated files rather than exact values that may
change with OpenFF versions.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
from rdkit import Chem

from .fep_parsers import (
    get_lig_start_indices,
    normalize_restraints,
    parse_distance_restraints,
    parse_fep_atom_count,
    parse_sequence_restraints,
    parse_topology_header,
)


# Test pairs from lomap.json (sorted by weight)
TEST_PAIRS = [
    ("ejm_31", "ejm_42"),  # weight 0.90
    ("ejm_31", "ejm_50"),  # weight 0.90
    ("ejm_31", "ejm_45"),  # weight 0.74
]

SYSTEMS = ["water", "protein"]

# Path to Tyk2 tutorial
TEST_DIR = Path(__file__).parent.parent
PROJECT_ROOT = TEST_DIR.parent
TYK2_TUTORIAL_PATH = PROJECT_ROOT / "tutorials" / "Tyk2"


def write_mol_to_sdf(mol: Chem.Mol, path: Path) -> Path:
    """Write a single molecule to an SDF file."""
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


def extract_ligand_by_name(sdf_path: Path, ligand_name: str, dest_dir: Path) -> Path:
    """Extract a specific ligand from a multi-molecule SDF file."""
    supplier = Chem.SDMolSupplier(str(sdf_path))
    for mol in supplier:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if name == ligand_name:
            out_path = dest_dir / f"{ligand_name}.sdf"
            return write_mol_to_sdf(mol, out_path)
    raise ValueError(f"Ligand {ligand_name} not found in {sdf_path}")


class TestFEPSetupConsistencyGolden:
    """Golden file regression tests for FEP setup output structure."""

    @pytest.fixture(scope="class")
    def fep_work_dir(self) -> Path:
        """Create a temporary working directory with required source files."""
        temp_dir = Path(tempfile.mkdtemp())

        # Copy source files
        src_sdf = TYK2_TUTORIAL_PATH / "ligprep" / "tyk2_ligands.sdf"
        src_protein = TYK2_TUTORIAL_PATH / "setupFEP" / "amber" / "protein_neutralized.pdb"
        src_water = TYK2_TUTORIAL_PATH / "setupFEP" / "amber" / "water.pdb"

        shutil.copy(src_water, temp_dir / "water.pdb")
        shutil.copy(src_protein, temp_dir / "protein.pdb")
        shutil.copy(src_sdf, temp_dir / "tyk2_ligands.sdf")

        yield temp_dir

        # Cleanup
        shutil.rmtree(temp_dir, ignore_errors=True)

    @pytest.fixture(scope="class")
    def prepared_ligands(self, fep_work_dir: Path) -> dict[str, Path]:
        """Extract and prepare ligands needed for test pairs."""
        src_sdf = fep_work_dir / "tyk2_ligands.sdf"
        required_ligands = {"ejm_31", "ejm_42", "ejm_45", "ejm_50"}

        # Extract each ligand
        extracted = {}
        for lig_name in required_ligands:
            sdf_path = extract_ligand_by_name(src_sdf, lig_name, fep_work_dir)
            extracted[lig_name] = sdf_path

        # Run qparams for each ligand
        for lig_name, sdf_path in extracted.items():
            cmd = [
                "micromamba",
                "run",
                "-n",
                "qligfep_new",
                "qparams",
                "-i",
                str(sdf_path),
                "-nagl",
                "-log",
                "error",
            ]
            result = subprocess.run(cmd, cwd=fep_work_dir, capture_output=True, text=True)
            if result.returncode != 0:
                pytest.fail(f"qparams failed for {lig_name}: {result.stderr}")

        return extracted

    def run_qligfep(
        self, work_dir: Path, lig1: str, lig2: str, system: str, random_state: int = 42
    ) -> Path:
        """Run qligfep and return path to inputfiles directory.

        Note: Using windows=3 because windows=1 causes division by zero with sigmoidal sampling.
        """
        cmd = [
            "micromamba",
            "run",
            "-n",
            "qligfep_new",
            "qligfep",
            "-l1",
            lig1,
            "-l2",
            lig2,
            "-FF",
            "AMBER14sb",
            "-s",
            system,
            "-c",
            "LOCAL",
            "-r",
            "25",
            "-w",
            "3",
            "-rs",
            str(random_state),
            "-log",
            "error",
        ]
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            pytest.fail(f"qligfep failed for {lig1}_{lig2} ({system}): {result.stderr}")

        return work_dir / f"FEP_{lig1}_{lig2}" / "inputfiles"

    @pytest.mark.golden
    @pytest.mark.slow
    @pytest.mark.fep_setup
    @pytest.mark.parametrize("lig1,lig2", TEST_PAIRS)
    @pytest.mark.parametrize("system", SYSTEMS)
    def test_fep_setup_structure(
        self,
        fep_work_dir: Path,
        prepared_ligands: dict[str, Path],
        fep_golden_manager,
        lig1: str,
        lig2: str,
        system: str,
    ):
        """Test that FEP setup output matches golden file expectations."""
        # Run qligfep
        inputfiles_dir = self.run_qligfep(fep_work_dir, lig1, lig2, system)

        try:
            # Read output files
            md_file = inputfiles_dir / "md_0500_0500.inp"
            fep_file = inputfiles_dir / "FEP1.fep"
            top_file = inputfiles_dir / "dualtop.top"
            pdb_file = inputfiles_dir / f"{lig1}_{lig2}.pdb"

            md_content = md_file.read_text()
            fep_content = fep_file.read_text()
            top_content = top_file.read_text()
            pdb_content = pdb_file.read_text()

            # Parse actual output
            restraints = parse_distance_restraints(md_content)
            lig1_start, lig2_start = get_lig_start_indices(pdb_content)
            lig1_indices, lig2_indices = normalize_restraints(restraints, lig1_start, lig2_start)
            has_seq_restraints = parse_sequence_restraints(md_content)
            _, solute_atoms = parse_topology_header(top_content)
            fep_atom_count = parse_fep_atom_count(fep_content)

            # Load golden expectations
            pair_key = f"{lig1}_{lig2}"
            golden = fep_golden_manager.load_golden("tyk2", pair_key, system)
            if golden is None:
                pytest.skip(f"No golden data for {pair_key}/{system}")

            # Compare distance restraints
            match, msg = fep_golden_manager.compare_distance_restraints(
                lig1_indices, lig2_indices, golden["distance_restraints"]
            )
            assert match, f"Distance restraints mismatch: {msg}"

            # Compare sequence restraints presence
            expected_seq = golden["sequence_restraints_present"]
            assert has_seq_restraints == expected_seq, (
                f"Sequence restraints presence mismatch: {has_seq_restraints} != {expected_seq}"
            )

            # Compare topology
            match, msg = fep_golden_manager.compare_topology(solute_atoms, golden["topology"])
            assert match, f"Topology mismatch: {msg}"

            # Compare FEP atom count
            match, msg = fep_golden_manager.compare_fep_atoms(fep_atom_count, golden["fep_file"])
            assert match, f"FEP file mismatch: {msg}"

        finally:
            # Clean up FEP directory for next test
            fep_dir = fep_work_dir / f"FEP_{lig1}_{lig2}"
            if fep_dir.exists():
                shutil.rmtree(fep_dir)


class TestFEPOutputParsers:
    """Unit tests for FEP output parser functions."""

    def test_parse_distance_restraints(self):
        """Test parsing distance restraints from MD input file."""
        content = """[MD]
steps                     5000

[sequence_restraints]
1      67      1.0 0 1

[distance_restraints]
1 33 0.0 0.1 0.5 0
2 34 0.0 0.1 0.5 0
3 35 0.0 0.1 0.5 0

[lambdas]
0.500 0.500
"""
        restraints = parse_distance_restraints(content)
        assert len(restraints) == 3
        assert restraints[0] == (1, 33, 0.5)
        assert restraints[1] == (2, 34, 0.5)
        assert restraints[2] == (3, 35, 0.5)

    def test_parse_sequence_restraints_present(self):
        """Test detecting presence of sequence restraints."""
        content_with = """[sequence_restraints]
1      67      1.0 0 1

[distance_restraints]
"""
        content_without = """[sequence_restraints]

[distance_restraints]
1 33 0.0 0.1 0.5 0
"""
        assert parse_sequence_restraints(content_with) is True
        assert parse_sequence_restraints(content_without) is False

    def test_parse_topology_header(self):
        """Test parsing topology header for atom counts."""
        content = """Q topology file
TITLE      test
    7603    4753 = Total no. of atoms, no. of solute atoms.
  -12.791    13.140   -47.770
"""
        total, solute = parse_topology_header(content)
        assert total == 7603
        assert solute == 4753

    def test_parse_fep_atom_count(self):
        """Test counting atoms in FEP file."""
        content = """[FEP]
states 2

[atoms]
1    4687
2    4688
3    4689

[change_charges]
1        -0.136     0.000
"""
        count = parse_fep_atom_count(content)
        assert count == 3

    def test_get_lig_start_indices(self):
        """Test finding ligand start indices in PDB."""
        content = """ATOM      1 C    PRO A   1       0.0   0.0   0.0  1.00  0.00           C
ATOM      2 N    PRO A   1       1.0   0.0   0.0  1.00  0.00           N
HETATM    3 C1   LIG A   2       2.0   0.0   0.0  1.00  0.00           C
HETATM    4 C2   LIG A   2       3.0   0.0   0.0  1.00  0.00           C
HETATM    5 C1   LID A   3       4.0   0.0   0.0  1.00  0.00           C
HETATM    6 C2   LID A   3       5.0   0.0   0.0  1.00  0.00           C
END
"""
        lig_start, lid_start = get_lig_start_indices(content)
        assert lig_start == 3
        assert lid_start == 5

    def test_normalize_restraints(self):
        """Test converting absolute indices to relative."""
        restraints = [(103, 203, 0.5), (104, 204, 0.5), (105, 205, 0.5)]
        lig1_indices, lig2_indices = normalize_restraints(restraints, 100, 200)
        assert lig1_indices == [3, 4, 5]
        assert lig2_indices == [3, 4, 5]
