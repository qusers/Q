"""Golden file regression tests for RestraintSetter consistency.

Ensures restraint mappings remain stable across code changes.
"""

from pathlib import Path

import pytest
from rdkit import Chem

from QligFEP.restraints.restraint_setter import RestraintSetter


def write_mol_to_sdf(mol: Chem.Mol, path: Path) -> Path:
    """Write a single molecule to an SDF file."""
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


# Method configurations matching the golden file generation script
METHOD_CONFIGS = {
    "element_p": {
        "atom_compare_method": "element",
        "strict_surround": False,
        "ignore_surround_atom_type": False,
    },
    "element_strict": {
        "atom_compare_method": "element",
        "strict_surround": True,
        "ignore_surround_atom_type": False,
    },
    "heavyatom_p": {
        "atom_compare_method": "heavyatom",
        "strict_surround": False,
        "ignore_surround_atom_type": False,
    },
}


class TestRestraintConsistencyGolden:
    """Golden file regression tests for RestraintSetter."""

    @pytest.mark.golden
    @pytest.mark.parametrize("method_name", ["element_p", "element_strict", "heavyatom_p"])
    def test_tyk2_golden_consistency(
        self,
        tyk2_sdf: Path,
        golden_manager,
        temp_work_dir: Path,
        method_name: str,
    ):
        """Test Tyk2 dataset restraints match golden files."""
        supplier = Chem.SDMolSupplier(str(tyk2_sdf))
        mols = [
            (m, m.GetProp("_Name") if m.HasProp("_Name") else f"mol_{i}")
            for i, m in enumerate(supplier)
            if m is not None
        ]

        method_config = METHOD_CONFIGS[method_name]

        # Test first 3 pairs (matching golden file generation)
        pairs_tested = 0
        max_pairs = 3

        for i, (mol1, name1) in enumerate(mols):
            for mol2, name2 in mols[i + 1 :]:
                if pairs_tested >= max_pairs:
                    break

                mol1_path = temp_work_dir / "mol1.sdf"
                mol2_path = temp_work_dir / "mol2.sdf"
                write_mol_to_sdf(mol1, mol1_path)
                write_mol_to_sdf(mol2, mol2_path)

                rsetter = RestraintSetter(mol1_path, mol2_path)
                restraints = rsetter.set_restraints(
                    atom_compare_method=method_config["atom_compare_method"],
                    strict_surround=method_config["strict_surround"],
                    ignore_surround_atom_type=method_config["ignore_surround_atom_type"],
                )

                match, msg = golden_manager.compare(restraints, "tyk2", name1, name2, method_name)
                assert match, f"Mismatch for {name1}_{name2}_{method_name}: {msg}"

                pairs_tested += 1

            if pairs_tested >= max_pairs:
                break

    @pytest.mark.golden
    @pytest.mark.slow
    def test_all_datasets_golden_consistency(
        self,
        all_sdf_files: list[Path],
        golden_manager,
        temp_work_dir: Path,
    ):
        """Test all datasets' restraints match golden files."""
        for sdf_path in all_sdf_files:
            dataset_name = sdf_path.stem.replace("_ligands", "")

            supplier = Chem.SDMolSupplier(str(sdf_path))
            mols = [
                (m, m.GetProp("_Name") if m.HasProp("_Name") else f"mol_{i}")
                for i, m in enumerate(supplier)
                if m is not None
            ]

            if len(mols) < 2:
                continue

            pairs_tested = 0
            max_pairs = 3

            for i, (mol1, name1) in enumerate(mols):
                for mol2, name2 in mols[i + 1 :]:
                    if pairs_tested >= max_pairs:
                        break

                    mol1_path = temp_work_dir / f"{dataset_name}_mol1.sdf"
                    mol2_path = temp_work_dir / f"{dataset_name}_mol2.sdf"
                    write_mol_to_sdf(mol1, mol1_path)
                    write_mol_to_sdf(mol2, mol2_path)

                    for method_name, method_config in METHOD_CONFIGS.items():
                        rsetter = RestraintSetter(mol1_path, mol2_path)
                        restraints = rsetter.set_restraints(
                            atom_compare_method=method_config["atom_compare_method"],
                            strict_surround=method_config["strict_surround"],
                            ignore_surround_atom_type=method_config["ignore_surround_atom_type"],
                        )

                        match, msg = golden_manager.compare(
                            restraints, dataset_name, name1, name2, method_name
                        )
                        assert match, f"Mismatch for {dataset_name}/{name1}_{name2}_{method_name}: {msg}"

                    pairs_tested += 1

                if pairs_tested >= max_pairs:
                    break


class TestRestraintMethodDifferences:
    """Tests verifying that different methods produce different results when expected."""

    @pytest.mark.golden
    def test_strict_vs_permissive_differ_appropriately(
        self,
        tyk2_sdf: Path,
        temp_work_dir: Path,
    ):
        """Strict mode should produce same or fewer restraints than permissive."""
        supplier = Chem.SDMolSupplier(str(tyk2_sdf))
        mols = [m for m in supplier if m is not None][:2]

        mol1_path = temp_work_dir / "mol1.sdf"
        mol2_path = temp_work_dir / "mol2.sdf"
        write_mol_to_sdf(mols[0], mol1_path)
        write_mol_to_sdf(mols[1], mol2_path)

        rsetter = RestraintSetter(mol1_path, mol2_path)

        permissive = rsetter.set_restraints(
            atom_compare_method="element",
            strict_surround=False,
        )

        strict = rsetter.set_restraints(
            atom_compare_method="element",
            strict_surround=True,
            ignore_surround_atom_type=False,
        )

        # Strict should have same or fewer restraints
        assert len(strict) <= len(permissive)

    @pytest.mark.golden
    def test_heavyatom_more_permissive_than_element(
        self,
        tyk2_sdf: Path,
        temp_work_dir: Path,
    ):
        """Heavyatom method should produce same or more restraints than element."""
        supplier = Chem.SDMolSupplier(str(tyk2_sdf))
        mols = [m for m in supplier if m is not None][:2]

        mol1_path = temp_work_dir / "mol1.sdf"
        mol2_path = temp_work_dir / "mol2.sdf"
        write_mol_to_sdf(mols[0], mol1_path)
        write_mol_to_sdf(mols[1], mol2_path)

        rsetter = RestraintSetter(mol1_path, mol2_path)

        element = rsetter.set_restraints(
            atom_compare_method="element",
            strict_surround=False,
        )

        heavyatom = rsetter.set_restraints(
            atom_compare_method="heavyatom",
            strict_surround=False,
        )

        # Heavyatom treats all non-H as equivalent, should have same or more restraints
        assert len(heavyatom) >= len(element)
