from pathlib import Path
from typing import Union

from openff.toolkit import Molecule
from openff.toolkit.utils import UndefinedStereochemistryError
from rdkit import Chem

from .logger import logger


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

    def __init__(self, lig, pattern: str = "*.sdf", force_rdkit: bool = False):
        """Initialize a Molecule Input/Output object. This helper class has a base functionality
        used for handling `.sdf`, like reading it, outputting separate `.sdf` files (required by lomap),
        and storing the molecules & their names into a single object.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
        """
        self.lig = lig
        self.force_rdkit = force_rdkit
        self.setup_mols_and_names(self.lig, pattern)
        self.parse_sdf_contents()  # add the sdf content to the dictionary

    def _parse_mol(self, ligpath: Union[Path, str]) -> tuple[list[Molecule], list[str]]:
        """Function to parse a .sdf file into a list of Molecule objects and their names.

        Args:
            ligpath: path to the .sdf file to be processed.

        Returns:
            Tuple[List[Molecule], List[str]]: a tuple containing the molecules and their names.
        """
        if isinstance(ligpath, Path):
            ligpath = str(ligpath)
        if self.force_rdkit:  # temporary fix for processing aligned molecules
            # For details on why is this here: https://github.com/openforcefield/openff-toolkit/issues/1872
            rdmols = [
                mol
                for mol in Chem.SDMolSupplier(
                    str(ligpath),
                    sanitize=False,
                )
                if mol is not None
            ]
            mols = []
            for mol in rdmols:
                try:
                    mols.append(Molecule.from_rdkit(mol))
                except UndefinedStereochemistryError:
                    logger.warning(
                        "Undefined stereochemistry in the input file!! Will try to process the ligands anyway."
                    )
                    mols.append(Molecule.from_rdkit(mol, allow_undefined_stereo=True))
                except Exception as e:
                    logger.error(f"Error processing file {ligpath}: {e}")
                    return [], []
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
                logger.debug(f"Searched for {pattern} files in {Path(lig).absolute}")
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
        if self.lig_files != []:
            ligands = []
            for ligfile in self.lig_files:
                ligands.extend(self._parse_sdf(ligfile))
        else:
            ligands = self._parse_sdf(self.lig)

        for lname, sdf_content in zip(self.lig_names, ligands):
            self.sdf_contents.update({lname: sdf_content})

    def _parse_sdf(self, lig: str) -> list[list[str]]:
        """Reads a `.sdf` file and returns a list containing the lines for each ligand in the file.

        Args:
            lig: the path to the `.sdf` file.

        Returns:
            list[list[str]]: a list of ligands, where each ligand is a list of lines.
        """
        with open(lig) as infile:
            content = infile.readlines()
            ligands = []
            current_ligand = []

            for line in content:
                if line.strip() == "$$$$":  # End of a ligand entry
                    ligands.append(current_ligand)
                    current_ligand = []  # Start a new ligand entry
                else:
                    current_ligand.append(line.rstrip())
        return ligands

    def write_sdf_separate(self, output_dir):
        """Function to write the separate multiple molecules within a sdf file into their own
        .sdf, placed under `output_dir`.

        Args:
            output_dir: the directory to save the single sdf molecules.

        Raises:
            ValueError: if `output_dir` is neither a `str` or a `pathlib.Path` object.
        """
        if isinstance(output_dir, str):
            output_dir = Path(output_dir)
        elif not isinstance(output_dir, Path):
            raise ValueError("output_dir must be either a str or pathlib.Path object")

        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        for idx, lname in enumerate(self.lig_names):
            fpath = str(output_dir / f"{lname}.sdf")
            self.molecules[idx].to_file(file_path=fpath, file_format="sdf")

    def write_to_single_sdf(self, output_name: str) -> None:
        """Writes all `self.molecules` to a single `.sdf` file.

        Args:
            output_name: name of the output file to write the aligned ligands to.
        """
        writer = Chem.SDWriter(output_name)
        for mol in self.molecules:
            writer.write(mol.to_rdkit())
        writer.close()
        logger.info(f"`self.molecules` written to {output_name}")
