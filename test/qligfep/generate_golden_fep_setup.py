"""Generate golden files for FEP setup consistency tests.

This module runs the full FEP setup workflow for selected ligand pairs
and captures structural data to use as golden file expectations.

Usage: pytest --generate-fep-golden
"""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from fep_parsers import (
    get_lig_start_indices,
    normalize_restraints,
    parse_distance_restraints,
    parse_fep_atom_count,
    parse_sequence_restraints,
    parse_topology_header,
)

# Paths
TEST_DIR = Path(__file__).parent.parent
PROJECT_ROOT = TEST_DIR.parent
TUTORIALS_DIR = PROJECT_ROOT / "tutorials"
GOLDEN_FEP_SETUP_DIR = Path(__file__).parent / "golden_fep_setup"

# Test pairs from lomap.json (sorted by weight)
TEST_PAIRS = [
    ("ejm_31", "ejm_42"),  # weight 0.90
    ("ejm_31", "ejm_50"),  # weight 0.90
    ("ejm_31", "ejm_45"),  # weight 0.74
]

SYSTEMS = ["water", "protein"]

# Ligands needed for test pairs
REQUIRED_LIGANDS = ["ejm_31", "ejm_42", "ejm_45", "ejm_50"]


def extract_ligands_from_sdf(src_sdf: Path, dest_dir: Path, ligand_names: list[str]) -> dict[str, Path]:
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


def run_qparams(work_dir: Path, sdf_path: Path) -> None:
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


def run_qligfep(work_dir: Path, lig1: str, lig2: str, system: str, random_state: int = 42) -> Path:
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


def parse_fep_output(inputfiles_dir: Path, lig1: str, lig2: str) -> dict:
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


def generate_golden_data() -> dict:
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
        ligand_paths = extract_ligands_from_sdf(src_sdf, work_dir, REQUIRED_LIGANDS)

        # Run qparams for each ligand
        print("Running qparams for each ligand...")
        for lig_name, sdf_path in ligand_paths.items():
            print(f"  {lig_name}...")
            run_qparams(work_dir, sdf_path)

        # Run qligfep for each pair and system
        for lig1, lig2 in TEST_PAIRS:
            pair_key = f"{lig1}_{lig2}"
            golden_data[pair_key] = {}

            for system in SYSTEMS:
                print(f"Running qligfep for {pair_key} ({system})...")
                inputfiles_dir = run_qligfep(work_dir, lig1, lig2, system)

                # Parse output
                data = parse_fep_output(inputfiles_dir, lig1, lig2)
                golden_data[pair_key][system] = data

                # Clean up FEP directory for next run
                fep_dir = work_dir / f"FEP_{lig1}_{lig2}"
                shutil.rmtree(fep_dir)

    return golden_data


def generate_all_fep_golden_files():
    """Entry point for pytest --generate-fep-golden command."""
    print("Generating FEP setup golden files...")

    golden_data = generate_golden_data()

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
