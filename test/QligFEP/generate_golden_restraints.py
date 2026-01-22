"""Generate golden files for RestraintSetter tests.

Usage:
    python test/QligFEP/generate_golden_restraints.py

This script generates expected restraint mappings for all ligand pairs in the test resources directory.
Run this script when:
- Setting up the test suite for the first time
- The RestraintSetter algorithm intentionally changes
"""

import argparse
import json
import tempfile
from pathlib import Path

from rdkit import Chem
from tqdm import tqdm

from QligFEP.restraints.restraint_setter import RestraintSetter

# Configuration
SCRIPT_DIR = Path(__file__).parent
RESOURCES_DIR = SCRIPT_DIR / "resources"
GOLDEN_DIR = SCRIPT_DIR / "golden_restraints"

# Methods to generate golden files for
METHODS = [
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


def write_mol_to_sdf(mol: Chem.Mol, path: Path) -> Path:
    """Write a single molecule to an SDF file."""
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


def get_restraints_for_pair(mol1_path: Path, mol2_path: Path, method_config: dict) -> dict:
    """Get restraints for a molecule pair with given method configuration."""
    rsetter = RestraintSetter(mol1_path, mol2_path)
    restraints = rsetter.set_restraints(
        atom_compare_method=method_config["atom_compare_method"],
        strict_surround=method_config["strict_surround"],
        ignore_surround_atom_type=method_config["ignore_surround_atom_type"],
    )
    # Convert int keys to strings for JSON serialization
    return {str(k): v for k, v in restraints.items()}


def generate_golden_for_dataset(
    sdf_path: Path, output_dir: Path, max_pairs: int = MAX_PAIRS_PER_DATASET
) -> int:
    """Generate golden files for a single dataset."""
    dataset_name = sdf_path.stem.replace("_ligands", "")
    dataset_output = output_dir / dataset_name
    dataset_output.mkdir(parents=True, exist_ok=True)

    supplier = Chem.SDMolSupplier(str(sdf_path))
    mols = [
        (m, m.GetProp("_Name") if m.HasProp("_Name") else f"mol_{i}")
        for i, m in enumerate(supplier)
        if m is not None
    ]

    if len(mols) < 2:
        print(f"  Skipping {sdf_path.name}: fewer than 2 molecules")
        return 0

    count = 0
    pairs_processed = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for i, (mol1, name1) in enumerate(mols):
            for j, (mol2, name2) in enumerate(mols[i + 1 :], start=i + 1):
                if pairs_processed >= max_pairs:
                    break

                mol1_path = tmpdir / "mol1.sdf"
                mol2_path = tmpdir / "mol2.sdf"
                write_mol_to_sdf(mol1, mol1_path)
                write_mol_to_sdf(mol2, mol2_path)

                for method in METHODS:
                    method_name = method["name"]
                    golden_file = dataset_output / f"{name1}_{name2}_{method_name}.json"

                    try:
                        restraints = get_restraints_for_pair(mol1_path, mol2_path, method)
                        with open(golden_file, "w") as f:
                            json.dump(restraints, f, indent=2, sort_keys=True)
                        count += 1
                    except Exception as e:
                        print(f"  Error for {name1}_{name2}_{method_name}: {e}")

                pairs_processed += 1

            if pairs_processed >= max_pairs:
                break

    return count


def main():
    parser = argparse.ArgumentParser(description="Generate golden files for RestraintSetter tests")
    parser.add_argument(
        "--max-pairs", type=int, default=MAX_PAIRS_PER_DATASET, help="Maximum pairs per dataset"
    )
    parser.add_argument("--dataset", type=str, default=None, help="Only generate for specific dataset")
    args = parser.parse_args()

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    sdf_files = sorted(RESOURCES_DIR.glob("*.sdf"))

    if args.dataset:
        sdf_files = [f for f in sdf_files if args.dataset in f.stem]

    print(f"Generating golden files for {len(sdf_files)} datasets...")
    print(f"Methods: {[m['name'] for m in METHODS]}")
    print(f"Max pairs per dataset: {args.max_pairs}")
    print()

    total_count = 0
    for sdf_path in tqdm(sdf_files, desc="Datasets"):
        count = generate_golden_for_dataset(sdf_path, GOLDEN_DIR, args.max_pairs)
        total_count += count

    print(f"\nGenerated {total_count} golden files in {GOLDEN_DIR}")


if __name__ == "__main__":
    main()
