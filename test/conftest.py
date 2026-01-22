"""Shared pytest fixtures for QligFEP test suite.

Provides paths to test resources, temporary directories, and golden file helpers.
"""

import json
import shutil
import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest

# Root paths
TEST_DIR = Path(__file__).parent
PROJECT_ROOT = TEST_DIR.parent
QLIGFEP_TEST_DIR = TEST_DIR / "QligFEP"
RESOURCES_DIR = QLIGFEP_TEST_DIR / "resources"
GOLDEN_DIR = QLIGFEP_TEST_DIR / "golden_restraints"
TUTORIALS_DIR = PROJECT_ROOT / "tutorials"


@pytest.fixture
def test_resources_path() -> Path:
    """Path to test/QligFEP/resources/ containing SDF datasets."""
    return RESOURCES_DIR


@pytest.fixture
def golden_restraints_path() -> Path:
    """Path to golden restraint files directory."""
    return GOLDEN_DIR


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
