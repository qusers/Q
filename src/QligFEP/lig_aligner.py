from copy import deepcopy
from rdkit import Chem
from tqdm import tqdm
from pathlib import Path
import subprocess
from joblib import delayed, Parallel, parallel_config
from functools import partial

from . import SRC
from .chemIO import MoleculeIO
from .logger import logger

class LigandAligner:
    def __init__(self, molio: MoleculeIO, n_threads: int = 1):
        """initalize the ligand aligner class from a MoleculeIO object.

        Args:
            molio: MoleculeIO object (input/output). See `chemIO.py` for more details.
            n_threads: Number of threads to create for the ligand alignment part. Defaults to 1.

        Raises:
            FileNotFoundError: If the kcombu executable is not found.
        """        
        self.molecules = [deepcopy(m.to_rdkit()) for m in molio.molecules]
        self.lig_names = molio.lig_names
        self.kcombu_exe = str(SRC  / "kcombu/fkcombu")
        if not Path(self.kcombu_exe).exists():
            raise FileNotFoundError(f"Could not find kcombu executable at {self.kcombu_exe} make sure it is installed.")
        molio.parse_sdf_contents()
        molio.write_sdf_separate("to_align_ligands")
        
    def _transfer_sdf_metadata(self, original_file, aligned_file):
        """
        Copies metadata from the original SDF file to the aligned SDF file.

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

        aligned_writer = Chem.SDWriter(str(aligned_file))
        for mol in aligned_mols:
            aligned_writer.write(aligned_mol)
        aligned_writer.close()

    def _update_aligned_sdf_files(self, ligand_names, reference, ligpath):
        """
        Update the aligned SDF files with metadata from the original SDF files.

        Args:
        ligand_names (list): List of ligand names to process.
        reference (str): Name of the reference ligand.
        ligpath (Path): Path to the directory containing the SDF files.
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

    def kcombu_align(self, reference:str) -> None:
        """Aligns the ligands to a reference ligand using kcombu.

        Args:
            reference: The name of the ligand to align the other ligands to.
        """
        cwd = Path.cwd()
        ligpath = cwd / "to_align_ligands"
        ref_sdf = ligpath / f"{reference}.sdf"
        suffix = f"_{reference}_aligned"
        
        self.reference = ref_sdf
        to_align_sdfs = [ligpath / f"{name}.sdf" for name in self.lig_names if name != reference]
        
        commands = []
        logger.info(f'Aligning ligands to ref `{reference}` with kcombu...')
        for lig in to_align_sdfs:
            kcombu_options = (
                f"-T {lig} -R {ref_sdf} -osdfT {lig.with_stem(lig.stem + suffix)} -E 'V'"
            )
            commands.append(f"{self.kcombu_exe} {kcombu_options}")
        
        commands = tqdm(commands)
        partial_func = partial(subprocess.run, capture_output=True, text=True)
        with parallel_config(backend='threading', n_jobs=2):
            Parallel()(delayed(partial_func)(cmd.split()) for cmd in commands)
        self._update_aligned_sdf_files(self.lig_names, reference, ligpath)
            
    def aligned_ligands_to_sigle_file(self,output_name: str, align_dir= None) -> None:
        """Writes the aligned ligands to a single file.

        Args:
            align_dir: _description_
            output_name: _description_
        """        
        if align_dir is None:
            align_dir = Path.cwd() / "to_align_ligands"
        else: 
            align_dir = Path(align_dir)
        sdf_files = [lig for lig in align_dir.glob("*_aligned.sdf")] + [self.reference]

        molecules = []
        for sdf_file in sdf_files:
            # Use the supplier to load molecules from each file
            suppl = Chem.SDMolSupplier(str(sdf_file))
            for mol in suppl:
                if mol is not None:
                    molecules.append(mol)

        writer = Chem.SDWriter(output_name)
        for mol in molecules:
            writer.write(mol)
        writer.close()
        
        logger.info(f"Aligned ligands written to {output_name}")
        