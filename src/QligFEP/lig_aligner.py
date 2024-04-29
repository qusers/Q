"""Module with the LigandAligner class to align ligands to a reference ligand using kcombu."""

import shutil
import subprocess
from functools import partial
from pathlib import Path

from joblib import Parallel, delayed, parallel_config
from openff.toolkit import Molecule
from rdkit import Chem
from tqdm import tqdm

from . import SRC
from .chemIO import MoleculeIO
from .logger import logger


class LigandAligner(MoleculeIO):

    def __init__(
        self,
        lig,
        pattern: str = "*.sdf",
        n_threads: int = 1,
        temp_ligalign_dir: str = "to_align_ligands",
    ):
        """initalize the ligand aligner. This class inherits from MoleculeIO and adds the
        functionality to align the ligands to a reference using kcombu.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
            n_threads: Number of threads to create for the ligand alignment part. Defaults to 1.
            temp_ligalign_dir: name for the temporary directory to store the separate sdf files.
                Defaults to "to_align_ligands".

        Raises:
            FileNotFoundError: If the kcombu executable is not found.
        """
        super().__init__(lig, pattern)
        self.lig_is_dir = Path(lig).is_dir()
        self._setup_tempdir(temp_ligalign_dir)
        self.kcombu_exe = str(SRC / "kcombu/fkcombu")
        self.n_threads = n_threads
        self.reference_mol: Molecule = None
        if not Path(self.kcombu_exe).exists():
            raise FileNotFoundError(
                f"Could not find kcombu executable at {self.kcombu_exe} make sure it is installed."
            )

    def _setup_tempdir(self, tempdir: Path):
        """If the input is not a directory, create a temporary directory to store the
        separate sdf files and write them there. If the input is a directory, the temporary
        directory is the same as the input directory.

        Args:
            tempdir: Path to the temporary directory.
        """
        tempdir = Path(tempdir)
        if self.lig_files != []:  # created when the input is a directory
            self.tempdir = tempdir
            self.write_sdf_separate(self.tempdir)
        else:
            self.tempdir = Path(self.lig)

    def _transfer_sdf_metadata(self, original_file, aligned_file):
        """
        Copies metadata from the original SDF file to the aligned SDF file and adds
        hydrogens to the aligned molecules.

        Args:
        original_file (str or Path): Path to the original SDF file.
        aligned_file (str or Path): Path to the aligned SDF file after transformation.
        """
        original_supplier = Chem.SDMolSupplier(str(original_file))
        aligned_supplier = Chem.SDMolSupplier(str(aligned_file))

        aligned_mols = []
        # Iterate over molecules from both suppliers simultaneously
        for original_mol, aligned_mol in zip(original_supplier, aligned_supplier):
            if original_mol is not None and aligned_mol is not None:
                # Copy properties from the original molecule to the aligned molecule
                aligned_mols.append(aligned_mol)
                for prop_name in original_mol.GetPropNames():
                    prop_value = original_mol.GetProp(prop_name)
                    aligned_mol.SetProp(prop_name, prop_value)
                    aligned_mol = Chem.AddHs(aligned_mol, addCoords=True)
        aligned_writer = Chem.SDWriter(str(aligned_file))
        for mol in aligned_mols:
            aligned_writer.write(aligned_mol)
        aligned_writer.close()

    def _update_aligned_sdf_files(self, ligand_names, reference, ligpath):
        """Update the aligned SDF files with metadata from the original SDF files.

        Args:
            ligand_names: List of ligand names to process.
            reference: Name of the reference ligand.
            ligpath: Path to the directory containing the SDF files.
        """
        modification = f"_{reference}_aligned.sdf"

        for name in ligand_names:
            if name != reference:
                original_sdf = ligpath / f"{name}.sdf"
                aligned_sdf = ligpath / f"{name}{modification}"
                if aligned_sdf.exists():
                    self._transfer_sdf_metadata(original_sdf, aligned_sdf)
                    logger.debug(f"Updated metadata for {aligned_sdf}")
                else:
                    logger.error(f"Aligned file {aligned_sdf} not found.")

    def kcombu_align(self, reference: str) -> None:
        """Aligns the ligands to a reference ligand using kcombu and saves the path
        to the aligned ligand in `self.reference`.

        Args:
            reference: The name of the ligand to align the other ligands to.
        """
        cwd = Path.cwd()
        ligpath = cwd / self.tempdir
        self.refname = reference
        self.reference_path = ligpath / f"{reference}.sdf"
        self.suffix = f"_{reference}_aligned"

        to_align_sdfs = [
            ligpath / f"{name}.sdf" for name in self.lig_names if name != reference
        ]

        commands = []
        logger.info(f"Aligning ligands to ref `{reference}` with kcombu...")
        for lig in to_align_sdfs:
            kcombu_options = (
                f"-T {lig} -R {self.reference_path} "
                f"-osdfT {lig.with_stem(lig.stem + self.suffix)} -E 'V'"
            )
            commands.append(f"{self.kcombu_exe} {kcombu_options}")

        commands = tqdm(commands)
        partial_func = partial(subprocess.run, capture_output=True, text=True)
        with parallel_config(backend="threading", n_jobs=self.n_threads):
            Parallel()(delayed(partial_func)(cmd.split()) for cmd in commands)
        self._update_aligned_sdf_files(self.lig_names, reference, ligpath)

    def output_aligned_molecules(self, output_name: str) -> None:
        """Writes the aligned molecules to a single `.sdf` file.

        Args:
            output_name: name of the output file to write the aligned ligands to.
        """
        # make a temporary copy of the reference ligand with f"_{reference}_aligned"
        temp_ref_stem = f"{self.refname}_tempRef_{self.suffix}"
        shutil.copy(self.reference_path, self.reference_path.with_stem(temp_ref_stem))
        molio = MoleculeIO(self.tempdir, pattern=f"*{self.suffix}.sdf")
        self.reference_path.with_stem(temp_ref_stem).unlink()
        molio.write_to_single_sdf(output_name)
