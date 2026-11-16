"""Unit tests for RestraintSetter class - atom mapping for FEP restraints.

Test initialization, comparison methods, surround modes, and smoke tests.
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


def get_mol_pair_from_sdf(sdf_path: Path, idx1: int = 0, idx2: int = 1) -> tuple[Chem.Mol, Chem.Mol]:
    """Get two molecules from an SDF file by index."""
    supplier = Chem.SDMolSupplier(str(sdf_path))
    mols = [m for m in supplier if m is not None]
    return mols[idx1], mols[idx2]


@pytest.fixture
def tyk2_mol_pair(tyk2_sdf: Path, temp_work_dir: Path) -> tuple[Path, Path]:
    """Get first two molecules from Tyk2 SDF as separate files."""
    mol1, mol2 = get_mol_pair_from_sdf(tyk2_sdf)

    mol1_path = temp_work_dir / "mol1.sdf"
    mol2_path = temp_work_dir / "mol2.sdf"

    write_mol_to_sdf(mol1, mol1_path)
    write_mol_to_sdf(mol2, mol2_path)

    return mol1_path, mol2_path


class TestRestraintSetterInit:
    """Tests for RestraintSetter initialization."""

    def test_init_with_valid_sdf_paths(self, tyk2_mol_pair: tuple[Path, Path]):
        """RestraintSetter initializes correctly with valid SDF paths."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(str(mol1_path), str(mol2_path))

        assert rsetter.molA is not None
        assert rsetter.molB is not None
        assert rsetter.kartograf_mapping is not None
        assert isinstance(rsetter.atom_mapping, dict)

    def test_init_with_path_objects(self, tyk2_mol_pair: tuple[Path, Path]):
        """RestraintSetter accepts Path objects."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)

        assert rsetter.molA is not None

    def test_init_with_custom_max_distance(self, tyk2_mol_pair: tuple[Path, Path]):
        """RestraintSetter respects custom kartograf_max_atom_distance."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path, kartograf_max_atom_distance=1.5)

        assert rsetter.atom_mapping is not None

    def test_init_rejects_pdb_files(self, temp_work_dir: Path):
        """RestraintSetter raises ValueError for PDB files."""
        pdb_path = temp_work_dir / "test.pdb"
        pdb_path.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N")

        with pytest.raises(ValueError, match="RDKit can not safely read PDBs"):
            RestraintSetter(str(pdb_path), str(pdb_path))


class TestSetRestraints:
    """Tests for set_restraints method with different parameters."""

    @pytest.mark.parametrize(
        "atom_compare_method",
        ["element", "hybridization", "aromaticity", "heavyatom"],
    )
    def test_all_atom_compare_methods(self, tyk2_mol_pair: tuple[Path, Path], atom_compare_method: str):
        """set_restraints works with all atom comparison methods."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        restraints = rsetter.set_restraints(
            atom_compare_method=atom_compare_method,
            strict_surround=False,
            ignore_surround_atom_type=False,
        )

        assert isinstance(restraints, dict)
        assert len(restraints) > 0
        # All keys and values should be non-negative integers
        for k, v in restraints.items():
            assert isinstance(k, int) and k >= 0
            assert isinstance(v, int) and v >= 0

    @pytest.mark.parametrize(
        "strict_surround,ignore_surround_atom_type",
        [
            (False, False),  # permissive
            (True, True),  # less_strict
            (True, False),  # strict
        ],
    )
    def test_all_surround_modes(
        self, tyk2_mol_pair: tuple[Path, Path], strict_surround: bool, ignore_surround_atom_type: bool
    ):
        """set_restraints works with all surround comparison modes."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        restraints = rsetter.set_restraints(
            atom_compare_method="element",
            strict_surround=strict_surround,
            ignore_surround_atom_type=ignore_surround_atom_type,
        )

        assert isinstance(restraints, dict)
        assert len(restraints) > 0

    def test_kartograf_native_mode(self, tyk2_mol_pair: tuple[Path, Path]):
        """set_restraints returns kartograf mapping when kartograf_native=True."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        native_restraints = rsetter.set_restraints(kartograf_native=True)

        # Should return the raw kartograf mapping
        assert native_restraints == rsetter.atom_mapping

    def test_restraints_stored_as_attribute(self, tyk2_mol_pair: tuple[Path, Path]):
        """Restraints are stored as an attribute after calling set_restraints."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        restraints = rsetter.set_restraints()

        assert hasattr(rsetter, "restraints")
        assert rsetter.restraints == restraints


class TestRestraintSetterAllDatasets:
    """Smoke tests running RestraintSetter on all available datasets."""

    @pytest.mark.slow
    def test_all_datasets_no_crash(self, all_sdf_files: list[Path], temp_work_dir: Path):
        """RestraintSetter runs without crashing on all dataset first pairs."""
        for sdf_path in all_sdf_files:
            supplier = Chem.SDMolSupplier(str(sdf_path))
            mols = [m for m in supplier if m is not None]

            if len(mols) < 2:
                continue

            mol1_path = temp_work_dir / f"{sdf_path.stem}_mol1.sdf"
            mol2_path = temp_work_dir / f"{sdf_path.stem}_mol2.sdf"

            write_mol_to_sdf(mols[0], mol1_path)
            write_mol_to_sdf(mols[1], mol2_path)

            rsetter = RestraintSetter(mol1_path, mol2_path)
            restraints = rsetter.set_restraints()

            assert isinstance(restraints, dict), f"Failed for {sdf_path.name}"


class TestAreAtomsEquivalent:
    """Tests for the static method are_atoms_equivalent."""

    def test_element_same_atomic_number(self, tyk2_mol_pair: tuple[Path, Path]):
        """are_atoms_equivalent returns True for same atomic number."""
        mol1_path, mol2_path = tyk2_mol_pair
        rsetter = RestraintSetter(mol1_path, mol2_path)

        mol = rsetter.molA.to_rdkit()
        # Get two carbon atoms (should be equivalent by element)
        carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
        if len(carbons) >= 2:
            result = rsetter.are_atoms_equivalent(carbons[0], carbons[1], compare_method="element")
            assert result is True

    def test_element_different_atomic_number(self, tyk2_mol_pair: tuple[Path, Path]):
        """are_atoms_equivalent returns False for different atomic numbers."""
        mol1_path, mol2_path = tyk2_mol_pair
        rsetter = RestraintSetter(mol1_path, mol2_path)

        mol = rsetter.molA.to_rdkit()
        carbon = next((a for a in mol.GetAtoms() if a.GetAtomicNum() == 6), None)
        nitrogen = next((a for a in mol.GetAtoms() if a.GetAtomicNum() == 7), None)

        if carbon and nitrogen:
            result = rsetter.are_atoms_equivalent(carbon, nitrogen, compare_method="element")
            assert result is False

    def test_heavyatom_all_heavy_atoms_equivalent(self, tyk2_mol_pair: tuple[Path, Path]):
        """heavyatom method treats all non-H atoms as equivalent."""
        mol1_path, mol2_path = tyk2_mol_pair
        rsetter = RestraintSetter(mol1_path, mol2_path)

        mol = rsetter.molA.to_rdkit()
        carbon = next((a for a in mol.GetAtoms() if a.GetAtomicNum() == 6), None)
        nitrogen = next((a for a in mol.GetAtoms() if a.GetAtomicNum() == 7), None)

        if carbon and nitrogen:
            result = rsetter.are_atoms_equivalent(carbon, nitrogen, compare_method="heavyatom")
            assert result is True

    def test_invalid_compare_method_raises(self, tyk2_mol_pair: tuple[Path, Path]):
        """are_atoms_equivalent raises ValueError for invalid compare method."""
        mol1_path, mol2_path = tyk2_mol_pair
        rsetter = RestraintSetter(mol1_path, mol2_path)

        mol = rsetter.molA.to_rdkit()
        atom = mol.GetAtomWithIdx(0)

        with pytest.raises(ValueError, match="Invalid compare method"):
            rsetter.are_atoms_equivalent(atom, atom, compare_method="invalid")


class TestRestraintConsistency:
    """Tests ensuring restraints are consistent and valid."""

    def test_restraints_map_valid_atoms(self, tyk2_mol_pair: tuple[Path, Path]):
        """All restraint indices should be valid atom indices."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        restraints = rsetter.set_restraints()

        molA = rsetter.molA.to_rdkit()
        molB = rsetter.molB.to_rdkit()
        num_atoms_A = molA.GetNumAtoms()
        num_atoms_B = molB.GetNumAtoms()

        for k, v in restraints.items():
            assert 0 <= k < num_atoms_A, f"Key {k} out of range for molA ({num_atoms_A} atoms)"
            assert 0 <= v < num_atoms_B, f"Value {v} out of range for molB ({num_atoms_B} atoms)"

    def test_restraints_subset_of_kartograf_mapping(self, tyk2_mol_pair: tuple[Path, Path]):
        """Restraints should be a subset of the kartograf mapping."""
        mol1_path, mol2_path = tyk2_mol_pair

        rsetter = RestraintSetter(mol1_path, mol2_path)
        restraints = rsetter.set_restraints()

        # Our restraints should be a subset of what kartograf mapped
        for k, v in restraints.items():
            assert k in rsetter.atom_mapping
            assert rsetter.atom_mapping[k] == v
