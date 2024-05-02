"""Module with the LigandAligner class to align ligands to a reference ligand using kcombu."""

import shutil
import subprocess
from functools import partial
from pathlib import Path

from joblib import Parallel, delayed, parallel_config
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import rdFMCS
from tqdm import tqdm

from . import SRC
from .chemIO import MoleculeIO
from .logger import logger


class LigandAligner(MoleculeIO):

    # TODO: in the future it would be good to make this more flexible so that it
    # could also perform 1:1 alignment between ligands. E.g.: FEP_lig1_lig2 would
    # have lig 1 aligned to lig 2 or vice-versa.

    def __init__(
        self,
        lig,
        pattern: str = "*.sdf",
        n_threads: int = 1,
        tempdir: str = "to_align_ligands",
        delete_tempdir: bool = True,
    ):
        """initalize the ligand aligner. This class inherits from MoleculeIO and adds the
        functionality to align the ligands to a reference using kcombu.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
            n_threads: Number of threads to create for the ligand alignment part. Defaults to 1.
            tempdir: name for the temporary directory to store the separate sdf files.
                Defaults to "to_align_ligands".
            delete_tempdir: If True, the temporary directory will be deleted upon calling
                `output_aligned_molecules`.

        Raises:
            FileNotFoundError: If the kcombu executable is not found.
        """
        super().__init__(lig, pattern)
        self._setup_tempdir(tempdir, delete_tempdir=delete_tempdir)
        self.kcombu_exe = str(SRC / "kcombu/fkcombu")
        self.n_threads = n_threads
        self.reference_mol: Molecule = None
        if not Path(self.kcombu_exe).exists():
            raise FileNotFoundError(
                f"Could not find kcombu executable at {self.kcombu_exe} make sure it is installed."
            )

    def _setup_tempdir(self, tempdir: Path, delete_tempdir):
        """If the input is not a directory, create a temporary directory to store the
        separate sdf files and write them there. If the input is a directory, the temporary
        directory is the same as the input directory.

        Args:
            tempdir: Path to the temporary directory.
            delete_tempdir: If True, the temporary directory will be deleted upon calling
                `output_aligned_molecules`.
        """
        tempdir = Path(tempdir)
        if self.lig_files != []:  # created when the input is a directory
            self.tempdir = Path(self.lig)
            logger.warning(
                "Input directory is used as tempdir for the aligned ligands " "and won't be deleted"
            )
            self.delete_tempdir = False  # never delete the input directory
        elif Path(self.lig).suffix == ".sdf":
            self.tempdir = tempdir
            self.delete_tempdir = delete_tempdir
        self.write_sdf_separate(self.tempdir)

    @staticmethod
    def _transfer_charges_metadata(molA, molB):
        """Method to transfer the charges from one molecule to another based on the MCS.
        This is needed because the aligned molecules created by kcombu do not have the
        charges & hydrogens from the reference molecule. MCS is used to make sure the
        atoms are mapped correctly.

        Args:
            molA: molecule A to transfer the charges from.
            molB: molecule B to transfer the charges to.

        Returns:
            molA, molB: Molecules with the same charges.
        """
        mcs_result = rdFMCS.FindMCS(
            [molA, molB],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            completeRingsOnly=True,
        )
        mcs_smarts = mcs_result.smartsString  # SMARTS pattern of the MCS
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)  # Convert SMARTS to an RDKit Mol object

        # Map atoms based on MCS
        matchA = molA.GetSubstructMatch(mcs_mol)
        matchB = molB.GetSubstructMatch(mcs_mol)

        # Optional: Print mappings to debug
        logger.trace("Mapping of atoms:")
        for a, b in zip(matchA, matchB):
            logger.trace(f"MolA atom {a} maps to MolB atom {b}")

        # Example application: transferring formal charges based on the mapping
        for a, b in zip(matchA, matchB):
            atomA = molA.GetAtomWithIdx(a)
            atomB = molB.GetAtomWithIdx(b)
            formal_charge = atomA.GetFormalCharge()
            atomB.SetFormalCharge(formal_charge)
            if atomA.HasProp("_GasteigerCharge"):
                gaister_charge = atomA.GetProp("_GasteigerCharge")
                atomB.SetProp("_GasteigerCharge", gaister_charge)
            # if atom is oxygen, check for valence 3 and remove hydrogens bound to it
            if atomA.GetAtomicNum() == 8 and atomA.GetTotalDegree() == 3:
                # forcefully remove the hydrogen bound to it
                atomB.SetNumExplicitHs(0)
                atomB.UpdatePropertyCache()
        return molA, molB

    def _transfer_sdf_metadata(self, original_file, aligned_file):
        """
        Copies metadata from the original SDF file to the aligned SDF file and adds
        hydrogens to the aligned ligands.

        Args:
        original_file (str or Path): Path to the original SDF file.
        aligned_file (str or Path): Path to the aligned SDF file after transformation.
        """
        original_supplier = Chem.SDMolSupplier(str(original_file), removeHs=True)
        aligned_supplier = Chem.SDMolSupplier(str(aligned_file), removeHs=True, sanitize=False)

        aligned_mols = []
        # Iterate over molecules from both suppliers simultaneously
        for original_mol, aligned_mol in zip(original_supplier, aligned_supplier):
            if original_mol is not None and aligned_mol is not None:
                Chem.SanitizeMol(aligned_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
                # Copy all molecular properties
                for prop_name in original_mol.GetPropNames():
                    prop_value = original_mol.GetProp(prop_name)
                    aligned_mol.SetProp(prop_name, prop_value)

                # remove the H's to make it easier on the MCS & transfer the charges
                original_mol = Chem.RemoveHs(original_mol)
                logger.debug(f"Transferring metadata to {aligned_file}")
                original_mol, aligned_mol = self._transfer_charges_metadata(original_mol, aligned_mol)
                aligned_mols.append(aligned_mol)

        aligned_writer = Chem.SDWriter(str(aligned_file))
        for mol in aligned_mols:
            aligned_writer.write(mol)
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

        to_align_sdfs = [ligpath / f"{name}.sdf" for name in self.lig_names if name != reference]

        commands = []
        logger.info(f"Aligning ligands to ref `{reference}` with kcombu...")
        for lig in to_align_sdfs:
            kcombu_options = (
                f"-T {lig} -R {self.reference_path} " f"-osdfT {lig.with_stem(lig.stem + self.suffix)} -E 'V'"
            )
            commands.append(f"{self.kcombu_exe} {kcombu_options}")

        commands = tqdm(commands)
        partial_func = partial(subprocess.run, capture_output=True, text=True)
        with parallel_config(backend="threading", n_jobs=self.n_threads):
            Parallel()(delayed(partial_func)(cmd.split()) for cmd in commands)
        self._update_aligned_sdf_files(self.lig_names, reference, ligpath)

    def output_aligned_ligands(self, output_name: str) -> None:
        """Writes the aligned molecules to a single `.sdf` file.

        Args:
            output_name: name of the output file to write the aligned ligands to.
        """
        # temporary copy of the reference ligand with f"_{reference}_aligned" suffix so we get it with glob
        temp_ref_stem = f"{self.refname}_tempRef_{self.suffix}"
        shutil.copy(self.reference_path, self.reference_path.with_stem(temp_ref_stem))
        molio = MoleculeIO(self.tempdir, pattern=f"*{self.suffix}.sdf")
        self.reference_path.with_stem(temp_ref_stem).unlink()
        molio.write_to_single_sdf(output_name)
        if self.delete_tempdir:
            logger.info(f"Deleting temporary directory with aligned ligands: {self.tempdir}")
            shutil.rmtree(self.tempdir)
