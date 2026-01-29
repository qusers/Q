"""Shared pytest fixtures for QligFEP test suite.

Provides paths to test resources, temporary directories, golden file helpers,
FEP output parsers, and golden file generation functions.

To generate golden files, run pytest with the following options:
    pytest --generate-restraint-golden
    pytest --generate-fep-golden
"""

import json
import re
import shutil
import subprocess
import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest

# -----------------------------------------------------------------------------
# Path Configuration
# -----------------------------------------------------------------------------

TEST_DIR = Path(__file__).parent
PROJECT_ROOT = TEST_DIR.parent
QLIGFEP_TEST_DIR = TEST_DIR / "qligfep"
RESOURCES_DIR = QLIGFEP_TEST_DIR / "resources"
GOLDEN_RESTRAINTS_DIR = QLIGFEP_TEST_DIR / "golden_restraints"
GOLDEN_FEP_SETUP_DIR = QLIGFEP_TEST_DIR / "golden_fep_setup"
TUTORIALS_DIR = PROJECT_ROOT / "tutorials"


# -----------------------------------------------------------------------------
# FEP Output Parsers (inlined from fep_parsers.py)
# -----------------------------------------------------------------------------


def parse_distance_restraints(md_content: str) -> list[tuple[int, int, float]]:
    """Extract (lig1_atom, lig2_atom, force) from [distance_restraints] section."""
    restraints = []
    in_section = False

    for line in md_content.splitlines():
        line = line.strip()

        if line.startswith("[distance_restraints]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            parts = line.split()
            if len(parts) >= 5:
                lig1_atom = int(parts[0])
                lig2_atom = int(parts[1])
                force = float(parts[4])
                restraints.append((lig1_atom, lig2_atom, force))

    return restraints


def parse_sequence_restraints(md_content: str) -> bool:
    """Check if [sequence_restraints] section has content."""
    in_section = False

    for line in md_content.splitlines():
        line = line.strip()

        if line.startswith("[sequence_restraints]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            return True

    return False


def parse_topology_header(top_content: str) -> tuple[int, int]:
    """Extract (total_atoms, solute_atoms) from topology header.

    The header line format is:
    '    7603    4753 = Total no. of atoms, no. of solute atoms.'
    """
    for line in top_content.splitlines():
        match = re.match(r"\s*(\d+)\s+(\d+)\s*=\s*Total no\. of atoms", line)
        if match:
            return int(match.group(1)), int(match.group(2))

    raise ValueError("Could not parse topology header")


def parse_fep_atom_count(fep_content: str) -> int:
    """Count entries in [atoms] section of FEP file."""
    count = 0
    in_section = False

    for line in fep_content.splitlines():
        line = line.strip()

        if line.startswith("[atoms]"):
            in_section = True
            continue
        if in_section and line.startswith("["):
            break
        if in_section and line:
            count += 1

    return count


def parse_inp_file(content: str) -> dict:
    """Parse a Q .inp file into a structured dict with all sections.

    Returns a dict with section names as keys. For key-value sections (like [MD]),
    the value is a dict of settings. For list sections (like [distance_restraints]),
    the value is a list of parsed entries.
    """
    result = {}
    current_section = None
    current_data: dict | list = {}

    for line in content.splitlines():
        line_stripped = line.strip()

        # Skip empty lines
        if not line_stripped:
            continue

        # Check for section header
        if line_stripped.startswith("[") and line_stripped.endswith("]"):
            # Save previous section
            if current_section is not None:
                result[current_section] = current_data

            current_section = line_stripped[1:-1]

            # Determine section type (list vs dict)
            if current_section in ("distance_restraints", "sequence_restraints", "trajectory_atoms"):
                current_data = []
            else:
                current_data = {}
            continue

        # Parse content based on section type
        if current_section is None:
            continue

        if isinstance(current_data, list):
            # List section - just store the line
            current_data.append(line_stripped)
        else:
            # Key-value section - parse "key value" format
            parts = line_stripped.split(None, 1)
            if len(parts) == 2:
                current_data[parts[0]] = parts[1]
            elif len(parts) == 1:
                # Handle lines with just a value (like lambdas section content)
                current_data["_value"] = parts[0]

    # Save last section
    if current_section is not None:
        result[current_section] = current_data

    return result


def parse_submission_script(content: str) -> dict:
    """Extract key settings from a SLURM submission script.

    Returns dict with:
    - seeds: list of seed values
    - temperatures: list of temperature values
    - fepfiles: list of FEP files
    - array_size: SBATCH array size
    - partition: SBATCH partition
    """
    result = {
        "seeds": [],
        "temperatures": [],
        "fepfiles": [],
        "array_size": None,
        "partition": None,
    }

    for line in content.splitlines():
        line = line.strip()

        # Parse seeds array
        if line.startswith("seeds=("):
            # Extract values from seeds=(1 2 3) format
            match = re.search(r"seeds=\(([^)]+)\)", line)
            if match:
                result["seeds"] = [int(x) for x in match.group(1).split()]

        # Parse temperatures array
        elif line.startswith("temperatures=("):
            match = re.search(r"temperatures=\(([^)]+)\)", line)
            if match:
                result["temperatures"] = [int(x) for x in match.group(1).split()]

        # Parse fepfiles array
        elif line.startswith("fepfiles=("):
            match = re.search(r"fepfiles=\(([^)]+)\)", line)
            if match:
                result["fepfiles"] = match.group(1).split()

        # Parse SBATCH array
        elif line.startswith("#SBATCH --array="):
            match = re.search(r"--array=(\d+)-(\d+)", line)
            if match:
                result["array_size"] = int(match.group(2)) - int(match.group(1)) + 1

        # Parse SBATCH partition
        elif line.startswith("#SBATCH -p "):
            result["partition"] = line.split()[-1]

    return result


def get_lig_start_indices(pdb_content: str) -> tuple[int, int]:
    """Find first atom index for LIG and LID residues in combined PDB."""
    lig_start = None
    lid_start = None

    for line in pdb_content.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        atom_num = int(line[6:11].strip())
        resname = line[17:20].strip()

        if resname == "LIG" and lig_start is None:
            lig_start = atom_num
        elif resname == "LID" and lid_start is None:
            lid_start = atom_num

        if lig_start is not None and lid_start is not None:
            break

    if lig_start is None or lid_start is None:
        raise ValueError("Could not find LIG and LID residues in PDB")

    return lig_start, lid_start


def normalize_restraints(
    restraints: list[tuple[int, int, float]], lig1_start: int, lig2_start: int
) -> tuple[list[int], list[int]]:
    """Convert absolute indices to relative (0-based from each ligand start).

    Returns (lig1_relative_indices, lig2_relative_indices).
    """
    lig1_indices = [r[0] - lig1_start for r in restraints]
    lig2_indices = [r[1] - lig2_start for r in restraints]
    return lig1_indices, lig2_indices


# -----------------------------------------------------------------------------
# Pytest Fixtures
# -----------------------------------------------------------------------------


@pytest.fixture
def test_resources_path() -> Path:
    """Path to test/qligfep/resources/ containing SDF datasets."""
    return RESOURCES_DIR


@pytest.fixture
def golden_restraints_path() -> Path:
    """Path to golden restraint files directory."""
    return GOLDEN_RESTRAINTS_DIR


@pytest.fixture
def all_sdf_files(test_resources_path: Path) -> list[Path]:
    """All SDF files available in test resources."""
    return sorted(test_resources_path.glob("*.sdf"))


@pytest.fixture
def tyk2_sdf(test_resources_path: Path) -> Path:
    """Path to Tyk2 ligands SDF file."""
    return test_resources_path / "tyk2_ligands.sdf"


@pytest.fixture
def tutorials_path() -> Path:
    """Path to tutorials directory."""
    return TUTORIALS_DIR


@pytest.fixture
def tyk2_tutorial_path(tutorials_path: Path) -> Path:
    """Path to Tyk2 tutorial directory."""
    return tutorials_path / "Tyk2"


@pytest.fixture
def temp_work_dir() -> Generator[Path, None, None]:
    """Create a temporary working directory that is cleaned up after the test."""
    temp_dir = Path(tempfile.mkdtemp())
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# -----------------------------------------------------------------------------
# Golden File Manager for Restraint Tests
# -----------------------------------------------------------------------------


class GoldenFileManager:
    """Helper class for golden file (snapshot) testing.

    Handles loading and comparing golden files for regression tests.
    Golden files are stored as one JSON file per dataset with structure:
    {
        "lig1_lig2": {
            "element_p": {"0": 0, ...},
            "element_strict": {...},
            ...
        },
        ...
    }
    """

    def __init__(self, golden_dir: Path):
        self.golden_dir = golden_dir
        self._cache: dict[str, dict] = {}

    def _load_dataset(self, dataset: str) -> dict | None:
        """Load and cache a dataset's golden file."""
        if dataset in self._cache:
            return self._cache[dataset]

        golden_path = self.golden_dir / f"{dataset}.json"
        if golden_path.exists():
            with open(golden_path) as f:
                self._cache[dataset] = json.load(f)
            return self._cache[dataset]
        return None

    def load_golden(self, dataset: str, lig1: str, lig2: str, method: str) -> dict | None:
        """Load golden data for a specific ligand pair and method."""
        dataset_data = self._load_dataset(dataset)
        if dataset_data is None:
            return None

        pair_key = f"{lig1}_{lig2}"
        pair_data = dataset_data.get(pair_key)
        if pair_data is None:
            return None

        return pair_data.get(method)

    def compare(self, actual: dict, dataset: str, lig1: str, lig2: str, method: str) -> tuple[bool, str]:
        """Compare actual output to golden file.

        Returns:
            tuple of (match: bool, message: str)
        """
        expected = self.load_golden(dataset, lig1, lig2, method)
        if expected is None:
            return False, f"Golden data not found for {dataset}/{lig1}_{lig2}/{method}"

        # Convert keys to strings for comparison (JSON loads int keys as strings)
        actual_str_keys = {str(k): v for k, v in actual.items()}
        expected_str_keys = {str(k): v for k, v in expected.items()}

        if actual_str_keys == expected_str_keys:
            return True, "Match"

        # Build detailed diff message
        missing_in_actual = set(expected_str_keys.keys()) - set(actual_str_keys.keys())
        extra_in_actual = set(actual_str_keys.keys()) - set(expected_str_keys.keys())
        value_mismatches = {
            k: (expected_str_keys[k], actual_str_keys[k])
            for k in set(expected_str_keys.keys()) & set(actual_str_keys.keys())
            if expected_str_keys[k] != actual_str_keys[k]
        }

        msg_parts = []
        if missing_in_actual:
            msg_parts.append(f"Missing keys: {missing_in_actual}")
        if extra_in_actual:
            msg_parts.append(f"Extra keys: {extra_in_actual}")
        if value_mismatches:
            msg_parts.append(f"Value mismatches: {value_mismatches}")

        return False, "; ".join(msg_parts)


@pytest.fixture
def golden_manager(golden_restraints_path: Path) -> GoldenFileManager:
    """GoldenFileManager instance for snapshot testing."""
    return GoldenFileManager(golden_restraints_path)


def get_ligand_names_from_sdf(sdf_path: Path) -> list[str]:
    """Extract ligand names from a multi-molecule SDF file."""
    from rdkit import Chem

    names = []
    supplier = Chem.SDMolSupplier(str(sdf_path))
    for mol in supplier:
        if mol is not None:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{len(names)}"
            names.append(name)
    return names


def get_ligand_pairs_from_sdf(sdf_path: Path, max_pairs: int | None = None) -> list[tuple[str, str]]:
    """Get all unique ligand pairs from an SDF file.

    Returns list of (ligand1_name, ligand2_name) tuples.
    """
    names = get_ligand_names_from_sdf(sdf_path)
    pairs = []
    for i, name1 in enumerate(names):
        for name2 in names[i + 1 :]:
            pairs.append((name1, name2))
            if max_pairs and len(pairs) >= max_pairs:
                return pairs
    return pairs


# -----------------------------------------------------------------------------
# FEP Golden File Manager
# -----------------------------------------------------------------------------


class FEPGoldenFileManager:
    """Helper class for FEP setup golden file (snapshot) testing.

    Validates FEP setup output structure rather than exact values.
    """

    def __init__(self, golden_dir: Path):
        self.golden_dir = golden_dir
        self._cache: dict[str, dict] = {}

    def _load_dataset(self, dataset: str) -> dict | None:
        """Load and cache a dataset's golden file."""
        if dataset in self._cache:
            return self._cache[dataset]

        golden_path = self.golden_dir / f"{dataset}.json"
        if golden_path.exists():
            with open(golden_path) as f:
                self._cache[dataset] = json.load(f)
            return self._cache[dataset]
        return None

    def load_golden(self, dataset: str, pair_key: str, system: str) -> dict | None:
        """Load golden data for a specific ligand pair and system (water/protein)."""
        dataset_data = self._load_dataset(dataset)
        if dataset_data is None:
            return None

        pair_data = dataset_data.get(pair_key)
        if pair_data is None:
            return None

        return pair_data.get(system)

    def compare_distance_restraints(
        self, actual_lig1: list[int], actual_lig2: list[int], expected: dict
    ) -> tuple[bool, str]:
        """Compare distance restraints against golden expectations."""
        expected_count = expected.get("count", 0)
        expected_lig1 = expected.get("lig1_indices", [])
        expected_lig2 = expected.get("lig2_indices", [])

        if len(actual_lig1) != expected_count:
            return False, f"Restraint count mismatch: {len(actual_lig1)} != {expected_count}"

        if actual_lig1 != expected_lig1:
            return False, f"LIG1 indices mismatch: {actual_lig1} != {expected_lig1}"

        if actual_lig2 != expected_lig2:
            return False, f"LIG2 indices mismatch: {actual_lig2} != {expected_lig2}"

        return True, "Match"

    def compare_topology(self, actual_solute_atoms: int, expected: dict) -> tuple[bool, str]:
        """Compare topology solute atom count against golden expectations."""
        if "solute_atoms" in expected:
            if actual_solute_atoms != expected["solute_atoms"]:
                return False, f"Solute atoms mismatch: {actual_solute_atoms} != {expected['solute_atoms']}"
        elif "solute_atoms_min" in expected and actual_solute_atoms < expected["solute_atoms_min"]:
            return (
                False,
                f"Solute atoms below minimum: {actual_solute_atoms} < {expected['solute_atoms_min']}",
            )
        return True, "Match"

    def compare_fep_atoms(self, actual_count: int, expected: dict) -> tuple[bool, str]:
        """Compare FEP atom count against golden expectations."""
        expected_count = expected.get("total_fep_atoms", 0)
        if actual_count != expected_count:
            return False, f"FEP atom count mismatch: {actual_count} != {expected_count}"
        return True, "Match"


@pytest.fixture
def fep_golden_manager() -> FEPGoldenFileManager:
    """FEPGoldenFileManager instance for FEP setup snapshot testing."""
    return FEPGoldenFileManager(GOLDEN_FEP_SETUP_DIR)


# -----------------------------------------------------------------------------
# Input File Golden Manager
# -----------------------------------------------------------------------------


class InputFileGoldenManager:
    """Compare generated .inp files against golden references.

    Validates that FEP input files match expected structure and values.
    Golden files are stored as raw .inp files that can be parsed and compared.
    """

    def __init__(self, golden_dir: Path):
        self.golden_dir = golden_dir
        self._cache: dict[str, dict] = {}

    def load_golden_inp(self, filename: str) -> dict:
        """Load and parse a golden .inp file.

        Args:
            filename: Name of the .inp file (e.g., "eq1.inp")

        Returns:
            Parsed inp file as a dict of sections
        """
        cache_key = str(self.golden_dir / filename)
        if cache_key in self._cache:
            return self._cache[cache_key]

        golden_path = self.golden_dir / filename
        if not golden_path.exists():
            raise FileNotFoundError(f"Golden file not found: {golden_path}")

        content = golden_path.read_text()
        parsed = parse_inp_file(content)
        self._cache[cache_key] = parsed
        return parsed

    def load_golden_submission_script(self) -> dict:
        """Load and parse the golden submission script."""
        # Look for any run*.sh file
        scripts = list(self.golden_dir.glob("run*.sh"))
        if not scripts:
            raise FileNotFoundError(f"No submission script found in {self.golden_dir}")

        content = scripts[0].read_text()
        return parse_submission_script(content)

    def compare_inp_files(self, actual: dict, expected: dict, ignore_keys: list[str] | None = None) -> tuple[bool, str]:
        """Compare parsed .inp file structures.

        Args:
            actual: Parsed actual .inp file
            expected: Parsed expected (golden) .inp file
            ignore_keys: List of section/key combinations to ignore (e.g., ["MD.random_seed"])

        Returns:
            Tuple of (match: bool, message: str)
        """
        ignore_keys = ignore_keys or []
        differences = []

        # Check all expected sections exist
        for section in expected:
            if section not in actual:
                differences.append(f"Missing section: [{section}]")
                continue

            expected_section = expected[section]
            actual_section = actual[section]

            # Compare section contents
            if isinstance(expected_section, dict):
                for key, expected_value in expected_section.items():
                    ignore_path = f"{section}.{key}"
                    if ignore_path in ignore_keys:
                        continue

                    if key not in actual_section:
                        differences.append(f"[{section}] missing key: {key}")
                    elif actual_section[key] != expected_value:
                        differences.append(
                            f"[{section}].{key}: expected '{expected_value}', got '{actual_section[key]}'"
                        )

                # Check for extra keys in actual
                for key in actual_section:
                    if key not in expected_section:
                        ignore_path = f"{section}.{key}"
                        if ignore_path not in ignore_keys:
                            differences.append(f"[{section}] extra key: {key}")

            elif isinstance(expected_section, list):
                if len(actual_section) != len(expected_section):
                    differences.append(
                        f"[{section}] length mismatch: expected {len(expected_section)}, got {len(actual_section)}"
                    )
                else:
                    for i, (exp_line, act_line) in enumerate(zip(expected_section, actual_section)):
                        if exp_line != act_line:
                            differences.append(f"[{section}][{i}]: expected '{exp_line}', got '{act_line}'")

        # Check for extra sections in actual
        for section in actual:
            if section not in expected:
                differences.append(f"Extra section: [{section}]")

        if differences:
            return False, "; ".join(differences[:5])  # Limit to first 5 differences
        return True, "Match"

    def compare_submission_scripts(self, actual: dict, expected: dict) -> tuple[bool, str]:
        """Compare submission script settings.

        Args:
            actual: Parsed actual script
            expected: Parsed expected (golden) script

        Returns:
            Tuple of (match: bool, message: str)
        """
        differences = []

        # Compare seeds count (not values, as they're random)
        if len(actual.get("seeds", [])) != len(expected.get("seeds", [])):
            differences.append(
                f"Seeds count mismatch: expected {len(expected.get('seeds', []))}, got {len(actual.get('seeds', []))}"
            )

        # Compare temperatures
        if actual.get("temperatures") != expected.get("temperatures"):
            differences.append(
                f"Temperatures mismatch: expected {expected.get('temperatures')}, got {actual.get('temperatures')}"
            )

        # Compare fepfiles
        if actual.get("fepfiles") != expected.get("fepfiles"):
            differences.append(
                f"FEP files mismatch: expected {expected.get('fepfiles')}, got {actual.get('fepfiles')}"
            )

        # Compare array size
        if actual.get("array_size") != expected.get("array_size"):
            differences.append(
                f"Array size mismatch: expected {expected.get('array_size')}, got {actual.get('array_size')}"
            )

        if differences:
            return False, "; ".join(differences)
        return True, "Match"


@pytest.fixture
def input_file_golden_manager() -> InputFileGoldenManager:
    """InputFileGoldenManager instance for FEP input file snapshot testing."""
    return InputFileGoldenManager(GOLDEN_FEP_SETUP_DIR / "tyk2_ejm_31_ejm_45")


# -----------------------------------------------------------------------------
# Golden File Generation - Restraints
# -----------------------------------------------------------------------------

# Methods to generate golden files for
RESTRAINT_METHODS = [
    {
        "name": "element_p",
        "atom_compare_method": "element",
        "strict_surround": False,
        "ignore_surround_atom_type": False,
    },
    {
        "name": "element_strict",
        "atom_compare_method": "element",
        "strict_surround": True,
        "ignore_surround_atom_type": False,
    },
    {
        "name": "heavyatom_p",
        "atom_compare_method": "heavyatom",
        "strict_surround": False,
        "ignore_surround_atom_type": False,
    },
    {
        "name": "heavyatom_strict",
        "atom_compare_method": "heavyatom",
        "strict_surround": True,
        "ignore_surround_atom_type": False,
    },
]

# Maximum pairs per dataset to avoid combinatorial explosion
MAX_PAIRS_PER_DATASET = 5


def _write_mol_to_sdf(mol, path: Path) -> Path:
    """Write a single molecule to an SDF file."""
    from rdkit import Chem

    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


def _get_restraints_for_pair(mol1_path: Path, mol2_path: Path, method_config: dict) -> dict:
    """Get restraints for a molecule pair with given method configuration."""
    from QligFEP.restraints.restraint_setter import RestraintSetter

    rsetter = RestraintSetter(mol1_path, mol2_path)
    restraints = rsetter.set_restraints(
        atom_compare_method=method_config["atom_compare_method"],
        strict_surround=method_config["strict_surround"],
        ignore_surround_atom_type=method_config["ignore_surround_atom_type"],
    )
    # Convert int keys to strings for JSON serialization
    return {str(k): v for k, v in restraints.items()}


def _generate_golden_for_dataset(sdf_path: Path, max_pairs: int = MAX_PAIRS_PER_DATASET) -> dict:
    """Generate golden data for a single dataset.

    Returns a dict with structure:
    {
        "pair1_pair2": {
            "element_p": {...},
            "element_strict": {...},
            ...
        },
        ...
    }
    """
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(sdf_path))
    mols = [
        (m, m.GetProp("_Name") if m.HasProp("_Name") else f"mol_{i}")
        for i, m in enumerate(supplier)
        if m is not None
    ]

    if len(mols) < 2:
        print(f"  Skipping {sdf_path.name}: fewer than 2 molecules")
        return {}

    dataset_data = {}
    pairs_processed = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for i, (mol1, name1) in enumerate(mols):
            for mol2, name2 in mols[i + 1 :]:
                if pairs_processed >= max_pairs:
                    break

                mol1_path = tmpdir / "mol1.sdf"
                mol2_path = tmpdir / "mol2.sdf"
                _write_mol_to_sdf(mol1, mol1_path)
                _write_mol_to_sdf(mol2, mol2_path)

                pair_key = f"{name1}_{name2}"
                pair_data = {}

                for method in RESTRAINT_METHODS:
                    method_name = method["name"]
                    try:
                        restraints = _get_restraints_for_pair(mol1_path, mol2_path, method)
                        pair_data[method_name] = restraints
                    except Exception as e:
                        print(f"  Error for {pair_key}/{method_name}: {e}")

                if pair_data:
                    dataset_data[pair_key] = pair_data

                pairs_processed += 1

            if pairs_processed >= max_pairs:
                break

    return dataset_data


def generate_all_restraint_golden_files(max_pairs: int = MAX_PAIRS_PER_DATASET, dataset: str | None = None):
    """Entry point for pytest --generate-restraint-golden command."""
    from tqdm import tqdm

    GOLDEN_RESTRAINTS_DIR.mkdir(parents=True, exist_ok=True)

    sdf_files = sorted(RESOURCES_DIR.glob("*.sdf"))

    if dataset:
        sdf_files = [f for f in sdf_files if dataset in f.stem]

    print(f"Generating golden files for {len(sdf_files)} datasets...")
    print(f"Methods: {[m['name'] for m in RESTRAINT_METHODS]}")
    print(f"Max pairs per dataset: {max_pairs}")
    print()

    total_pairs = 0
    for sdf_path in tqdm(sdf_files, desc="Datasets"):
        dataset_name = sdf_path.stem.replace("_ligands", "")
        dataset_data = _generate_golden_for_dataset(sdf_path, max_pairs)

        if dataset_data:
            golden_file = GOLDEN_RESTRAINTS_DIR / f"{dataset_name}.json"
            with open(golden_file, "w") as f:
                json.dump(dataset_data, f, indent=2, sort_keys=True)
            total_pairs += len(dataset_data)

    print(f"\nGenerated {len(sdf_files)} golden files ({total_pairs} pairs) in {GOLDEN_RESTRAINTS_DIR}")


# -----------------------------------------------------------------------------
# Golden File Generation - FEP Setup
# -----------------------------------------------------------------------------

# Test pairs from lomap.json (sorted by weight)
FEP_TEST_PAIRS = [
    ("ejm_31", "ejm_42"),  # weight 0.90
    ("ejm_31", "ejm_50"),  # weight 0.90
    ("ejm_31", "ejm_45"),  # weight 0.74
]

FEP_SYSTEMS = ["water", "protein"]

# Ligands needed for test pairs
FEP_REQUIRED_LIGANDS = ["ejm_31", "ejm_42", "ejm_45", "ejm_50"]


def _extract_ligands_from_sdf(src_sdf: Path, dest_dir: Path, ligand_names: list[str]) -> dict[str, Path]:
    """Extract specific ligands from multi-molecule SDF to individual files."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(src_sdf))
    extracted = {}

    for mol in supplier:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if name in ligand_names:
            out_path = dest_dir / f"{name}.sdf"
            writer = Chem.SDWriter(str(out_path))
            writer.write(mol)
            writer.close()
            extracted[name] = out_path

    missing = set(ligand_names) - set(extracted.keys())
    if missing:
        raise ValueError(f"Could not find ligands: {missing}")

    return extracted


def _run_qparams(work_dir: Path, sdf_path: Path) -> None:
    """Run qparams to generate .lib/.prm files for a ligand."""
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
        "warning",
    ]
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"qparams failed for {sdf_path.name}:")
        print(result.stderr)
        raise RuntimeError(f"qparams failed for {sdf_path.name}")


def _run_qligfep(work_dir: Path, lig1: str, lig2: str, system: str, random_state: int = 42) -> Path:
    """Run qligfep to set up FEP for a ligand pair.

    Returns the path to the FEP output directory.
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
        "warning",
    ]
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"qligfep failed for {lig1}_{lig2} ({system}):")
        print(result.stderr)
        raise RuntimeError(f"qligfep failed for {lig1}_{lig2} ({system})")

    return work_dir / f"FEP_{lig1}_{lig2}" / "inputfiles"


def _parse_fep_output(inputfiles_dir: Path, lig1: str, lig2: str) -> dict:
    """Parse FEP output files and extract structural data."""
    md_file = inputfiles_dir / "md_0500_0500.inp"
    fep_file = inputfiles_dir / "FEP1.fep"
    top_file = inputfiles_dir / "dualtop.top"
    pdb_file = inputfiles_dir / f"{lig1}_{lig2}.pdb"

    md_content = md_file.read_text()
    fep_content = fep_file.read_text()
    top_content = top_file.read_text()
    pdb_content = pdb_file.read_text()

    # Parse distance restraints
    restraints = parse_distance_restraints(md_content)
    lig1_start, lig2_start = get_lig_start_indices(pdb_content)
    lig1_indices, lig2_indices = normalize_restraints(restraints, lig1_start, lig2_start)

    # Get force from first restraint (should be same for all)
    force = restraints[0][2] if restraints else 0.5

    # Parse other sections
    has_seq_restraints = parse_sequence_restraints(md_content)
    total_atoms, solute_atoms = parse_topology_header(top_content)
    fep_atom_count = parse_fep_atom_count(fep_content)

    return {
        "distance_restraints": {
            "count": len(restraints),
            "lig1_indices": lig1_indices,
            "lig2_indices": lig2_indices,
            "force": force,
        },
        "sequence_restraints_present": has_seq_restraints,
        "topology": {
            "solute_atoms": solute_atoms,
        },
        "fep_file": {
            "total_fep_atoms": fep_atom_count,
        },
    }


def _generate_fep_golden_data() -> dict:
    """Run FEP setup for all test pairs and systems, collect golden data."""
    golden_data = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Copy source files
        src_sdf = TUTORIALS_DIR / "Tyk2" / "ligprep" / "tyk2_ligands.sdf"
        src_protein = TUTORIALS_DIR / "Tyk2" / "setupFEP" / "amber" / "protein_neutralized.pdb"
        src_water = TUTORIALS_DIR / "Tyk2" / "setupFEP" / "amber" / "water.pdb"

        # Copy water.pdb and protein.pdb to work directory
        shutil.copy(src_water, work_dir / "water.pdb")
        shutil.copy(src_protein, work_dir / "protein.pdb")

        print(f"Working in: {work_dir}")
        print("Extracting ligands from SDF...")

        # Extract ligands
        ligand_paths = _extract_ligands_from_sdf(src_sdf, work_dir, FEP_REQUIRED_LIGANDS)

        # Run qparams for each ligand
        print("Running qparams for each ligand...")
        for lig_name, sdf_path in ligand_paths.items():
            print(f"  {lig_name}...")
            _run_qparams(work_dir, sdf_path)

        # Run qligfep for each pair and system
        for lig1, lig2 in FEP_TEST_PAIRS:
            pair_key = f"{lig1}_{lig2}"
            golden_data[pair_key] = {}

            for system in FEP_SYSTEMS:
                print(f"Running qligfep for {pair_key} ({system})...")
                inputfiles_dir = _run_qligfep(work_dir, lig1, lig2, system)

                # Parse output
                data = _parse_fep_output(inputfiles_dir, lig1, lig2)
                golden_data[pair_key][system] = data

                # Clean up FEP directory for next run
                fep_dir = work_dir / f"FEP_{lig1}_{lig2}"
                shutil.rmtree(fep_dir)

    return golden_data


def generate_all_fep_golden_files():
    """Entry point for pytest --generate-fep-golden command."""
    print("Generating FEP setup golden files...")

    golden_data = _generate_fep_golden_data()

    # Ensure golden directory exists
    GOLDEN_FEP_SETUP_DIR.mkdir(parents=True, exist_ok=True)

    # Save golden data
    golden_path = GOLDEN_FEP_SETUP_DIR / "tyk2.json"
    with open(golden_path, "w") as f:
        json.dump(golden_data, f, indent=2)

    print(f"\nGolden data saved to: {golden_path}")

    # Print summary
    print("\nSummary:")
    for pair_key, pair_data in golden_data.items():
        print(f"\n{pair_key}:")
        for system, data in pair_data.items():
            dr = data["distance_restraints"]
            print(f"  {system}:")
            print(f"    distance_restraints: {dr['count']} pairs")
            print(f"    sequence_restraints_present: {data['sequence_restraints_present']}")
            print(f"    solute_atoms: {data['topology']['solute_atoms']}")
            print(f"    fep_atoms: {data['fep_file']['total_fep_atoms']}")


# -----------------------------------------------------------------------------
# Pytest hooks for golden file generation
# -----------------------------------------------------------------------------


def pytest_addoption(parser):
    """Register command-line options for golden file generation."""
    parser.addoption(
        "--generate-restraint-golden",
        action="store_true",
        default=False,
        help="Generate golden files for RestraintSetter tests",
    )
    parser.addoption(
        "--generate-fep-golden",
        action="store_true",
        default=False,
        help="Generate golden files for FEP setup tests",
    )


def pytest_configure(config):
    """Run golden file generation if requested."""
    if config.getoption("--generate-restraint-golden"):
        generate_all_restraint_golden_files()
        pytest.exit("Golden restraint files generated successfully", returncode=0)

    if config.getoption("--generate-fep-golden"):
        generate_all_fep_golden_files()
        pytest.exit("Golden FEP setup files generated successfully", returncode=0)
