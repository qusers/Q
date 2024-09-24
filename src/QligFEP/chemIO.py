from io import StringIO
from itertools import product
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from openff.toolkit import Molecule
from openff.toolkit.utils import UndefinedStereochemistryError
from rdkit import Chem

from .logger import logger
from .pdb_utils import read_pdb_to_dataframe, write_dataframe_to_pdb
from .restraints.hydrogen_utils import are_hydrogens_at_end, reindex_hydrogens_to_end


class MoleculeIO:
    """A class to handle the input/output of molecules. Ligands are usually initialized
    from (.sdf files | directories containing .sdf files) and then processed into individual
    molecules.

    Attributes:
        self.lig: the input ligand `.sdf` file or directory
        self.lig_files: a list of ligand files if the input is a directory
        self.molecules: a list of Molecule objects
        self.lig_names: a list with the names of the ligands
        self.sdf_contents: a dictionary of the sdf contents for each ligand
    """

    def __init__(self, lig, pattern: str = "*.sdf", reindex_hydrogens: bool = True):
        """Initialize a Molecule Input/Output object. This helper class has a base functionality
        used for handling `.sdf`, like reading it, outputting separate `.sdf` files (required by lomap),
        and storing the molecules & their names into a single object.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
            reindex_hydrogens: If True, loading molecules will assert that hydrogen atoms are at the end
                of the atom list and reindex them if they are not (needed by restraint setting algorithm).
                If False, the molecules will be loaded as is. Defaults to True.
        """
        self._reindex_hydrogens = reindex_hydrogens
        self.lig = lig
        self.setup_mols_and_names(self.lig, pattern)
        self.parse_sdf_contents()  # add the sdf content to the dictionary

    def _force_H_reindexing(self, mol: Molecule) -> Molecule:
        rdkit_mol = mol.to_rdkit()
        if not are_hydrogens_at_end(rdkit_mol):
            rdkit_mol = reindex_hydrogens_to_end(rdkit_mol)
            logger.warning(f"Hydrogens not at the end of the atom list for molecule {mol.name}. Reindexed.")
        return Molecule.from_rdkit(rdkit_mol, hydrogens_are_explicit=True)

    def _parse_mol(self, ligpath: Union[Path, str]) -> tuple[list[Molecule], list[str]]:
        """Function to parse a .sdf file into a list of Molecule objects and their names.

        Args:
            ligpath: path to the .sdf file to be processed.

        Returns:
            Tuple[List[Molecule], List[str]]: a tuple containing the molecules and their names.
        """
        if isinstance(ligpath, Path):
            ligpath = str(ligpath)
        try:
            mols = Molecule.from_file(ligpath)
        except UndefinedStereochemistryError:
            logger.warning(
                "Undefined stereochemistry in the input file!! Will try to process the ligands anyway."
            )
            mols = Molecule.from_file(ligpath, allow_undefined_stereo=True)
        except Exception as e:
            logger.error(f"Error processing file {ligpath}: {e}")
            return [], []

        mols = [mols] if not isinstance(mols, list) else mols
        lig_names = [mol.name if mol.name else f"lig_{idx}" for idx, mol in enumerate(mols)]
        if self._reindex_hydrogens:
            mols = [self._force_H_reindexing(mol) for mol in mols]
        for mol, name in zip(mols, lig_names):  # update name property
            mol.name = name
        return mols, lig_names

    def setup_mols_and_names(self, lig: str, pattern: str = "*.sdf"):
        """Function to setup the molecules and their names from the input ligand file (lig).
        Running this function after the object is already initialized will reset the `molecules`
        and `lig_names` attributes.

        Args:
            lig: the input ligand `.sdf` file or a directory containing `.sdf` files.
            pattern: pattern for finding the `.sdf` files in case `lig` is a directory. This
                argument is ignored if `lig` is a file path. Defaults to '*.sdf'.

        Raises:
            ValueError: if the input ligand file does not exist.
            ValueError: if the input ligand file is not a `.sdf` file or a directory containing

        Returns:
            Tuple[List[Molecule], List[str]]: a tuple containing the molecules and their names.
                Both of which are saved as the attributes `self.molecules` and `self.lig_names`.
        """
        self.lig = lig
        self.lig_files = []
        self.molecules = []
        self.lig_names = []
        if not Path(lig).exists():
            raise ValueError(f"File | directory: {lig} does not exist.")
        elif Path(lig).is_dir():
            lig = Path(lig)
            if pattern is not None:
                logger.debug(f"Searched for {pattern} files in {Path(lig).absolute()}")
                self.lig_files = sorted(list(lig.glob(pattern)))  # glob doesn't return sorted list
            mols_and_names = [self._parse_mol(lig) for lig in self.lig_files]
        elif Path(lig).exists():
            mols_and_names = [self._parse_mol(lig)]
        else:
            raise ValueError(
                f"No ligands found in {lig}. Make sure it's a `.sdf` file or a directory containing those."
            )
        for mols, names in mols_and_names:
            self.molecules.extend(mols)
            self.lig_names.extend(names)
        return self.molecules, self.lig_names

    def parse_sdf_contents(self):
        """Parse the SDF content into individual entries and have them saved in a dictionary."""
        self.sdf_contents = {}
        for mol, name in zip(self.molecules, self.lig_names):
            string_buffer = StringIO()
            mol.to_file(string_buffer, file_format="sdf")
            self.sdf_contents.update({name: string_buffer.getvalue().splitlines()})

    def write_sdf_separate(self, output_dir, molecules: Optional[list[Molecule]] = None) -> None:
        """Function to write the separate multiple molecules within a sdf file into their own
        .sdf, placed under `output_dir`.

        Args:
            output_dir: the directory to save the single sdf molecules.
            molecules: a list of Molecule objects to write to the output directory. If None, the
                `self.molecules` will be written. Defaults to None.

        Raises:
            ValueError: if `output_dir` is neither a `str` or a `pathlib.Path` object.
        """
        if isinstance(output_dir, str):
            output_dir = Path(output_dir)
        elif not isinstance(output_dir, Path):
            raise ValueError("output_dir must be either a str or pathlib.Path object")

        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        if molecules is None:
            for idx, lname in enumerate(self.lig_names):
                fpath = str(output_dir / f"{lname}.sdf")
                self.molecules[idx].to_file(file_path=fpath, file_format="sdf")
        else:
            for mol in molecules:
                mol.to_file(file_path=f"{mol.name}.sdf", file_format="sdf")

    def write_to_single_sdf(self, output_name: str, molecules: Optional[list[Molecule]] = None) -> None:
        """Writes all `self.molecules` to a single `.sdf` file.

        Args:
            output_name: name of the output file to write the aligned ligands to.
            molecules: a list of Molecule objects to write to the output directory. If None, the
                `self.molecules` will be written. Defaults to None.
        """
        writer = Chem.SDWriter(output_name)
        for mol in self.molecules if molecules is None else molecules:
            writer.write(mol.to_rdkit())
        writer.close()
        logger.info(f"{'`self.molecules`' if molecules is None else 'molecules'} written to {output_name}")

    def write_to_single_pdb(self, output_name: str, init_offset: int = 0) -> pd.DataFrame:
        """Writes all `self.molecules` to a single `.pdb` file.

        Args:
            output_name: name of the output file to write the aligned ligands to.
        """
        init_offset = 0
        ligands_dfs = []
        lig_resn = ["LI", "LG", "LH"]
        last_lig_resn = [d for d in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
        lig_resnames = ["".join(i) for i in product(lig_resn, last_lig_resn)]
        for mol, resn in zip(self.molecules, lig_resnames):
            # write the molecule pdb lines in memory and convert them to pd.DataFrame
            output = StringIO()
            mol.to_file(output, file_format="pdb")
            pdb_data = output.getvalue()
            pdb_lines = pdb_data.splitlines()
            output.close()  # Close the StringIO instance when done
            lig_pdb = read_pdb_to_dataframe(pdb_lines)
            # update the init_offset for the atom_serial_number & assign unique residue names
            lig_pdb = lig_pdb.assign(
                atom_serial_number=lambda x: x["atom_serial_number"].astype(int) + init_offset,  # noqa: B023
                residue_name=resn,
            )
            init_offset = lig_pdb["atom_serial_number"].max() + 1
            logger.trace(f"offset is now {init_offset}")
            ligands_dfs.append(lig_pdb)
        ligands_dfs = pd.concat(ligands_dfs, ignore_index=True)
        write_dataframe_to_pdb(ligands_dfs, output_name)
        return ligands_dfs
