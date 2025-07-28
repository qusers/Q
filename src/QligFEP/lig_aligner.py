import atexit
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from typing import Any, Optional, Union

from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import rdFMCS
from tqdm import tqdm

from .chemIO import MoleculeIO
from .logger import logger


class GlobalLigandAligner(MoleculeIO):
    """Align ligands based on three-dimensional coordinates using the fkcombu program.

        For more information on fkcombu, see their docs:
            https://pdbj.org/kcombu/doc/README_fkcombu.html
        For information on the `connectivity`, `top_constraint_tol` parameters, see:
            https://pdbj.org/kcombu/doc/README_pkcombu.html

    Attributes:
        SDF_EXTENSION (str): Default file extension for SDF files.
        ALIGNED_SUFFIX (str): Suffix added to aligned molecule files.
        kcombu_exe (Path): Path to the fkcombu executable.
        n_threads (int): Number of threads to use for parallel alignment.
        reference_mol (str): Name of the molecule used as reference for the last alignment performed.
        aligned_molecules (dict): Aligned molecules, stored with [key] as name of the molecule.
        temp_dir (Path): Temporary directory for alignment operations.
        fkparams (dict): Dictionary of fkcombu parameters.
        fkcombu_command (str): Base fkcombu command string.
    """

    SDF_EXTENSION = ".sdf"
    ALIGNED_SUFFIX = "_aligned"

    def __init__(
        self,
        lig,
        pattern: str = f"*{SDF_EXTENSION}",
        reindex_hydrogens: bool = True,
        n_threads: int = 1,
        protein: Optional[str] = None,
        energy: str = "a",
        search: str = "f",
        steep_descend: bool = True,
        connectivity: str = "t",
        top_constraint_tol: Optional[int] = None,
        **fkcombu_params,
    ):
        """
        Initialize a GlobalLigandAligner object for aligning ligands based on three-dimensional coordinates
        using the fkcombu program. Additional parameters can be passed to fkcombu through keyword arguments.

        For more information on fkcombu, see their docs:
            https://pdbj.org/kcombu/doc/README_fkcombu.html
        For information on the `connectivity`, `top_constraint_tol` parameters, see:
            https://pdbj.org/kcombu/doc/README_pkcombu.html

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to "*.sdf".
            reindex_hydrogens: If True, loading molecules will assert that hydrogen atoms are at the end
                of the atom list and reindex them if they are not (needed by restraint setting algorithm).
                If False, the molecules will be loaded as is. Defaults to True.
            n_threads: Number of threads to use for parallel alignment. Defaults to 1.
            energy: fkcombu energy calculation method. `a` for atom-match, `v` for volume-overlap. Defaults to `a`.
            search: fkcombu search method `f` for flexible, `r` for rigid, `n` for nothing. Defaults to `f`.
            steep_descend: fkcombu perform Gradient-based Steepest Descent fitting. Defaults to True.
            connectivity: fkcombu connectivity method for finding the MCS. `c` for connected,
                `s` for substructure, `i` for isomorphic, `t` for topo_constrained_disconnected. If more
                flexible correspondences are needed, use `t` together with the `top_constraint_tol` parameter.
            top_constraint_tol: the maximum number of bonds (shortest path) allowed as a tolerance for
                not breaking the connectivity of the MCS. Only used if `connectivity` is set to `t`.

        Raises:
            FileNotFoundError: If the kcombu executable is not found.
        """
        super().__init__(lig, pattern=pattern, reindex_hydrogens=reindex_hydrogens)
        self.kcombu_exe = self._set_fkcombu_exe()
        self.n_threads = n_threads
        self.reference_mol: Optional[Molecule] = None
        self.aligned_molecules: dict[str, Molecule] = {}
        self.temp_dir: Optional[tempfile.TemporaryDirectory] = None
        self.fkparams = self._process_fkparams(
            {"P": protein, "E": energy.upper(), "S": search.upper(), "SD": steep_descend, **fkcombu_params},
            connectivity.upper(),
            top_constraint_tol,
        )
        self.fkcombu_command = self._build_fkcombu_command()
        atexit.register(self.cleanup)  # Clean up temp directory on exit of the program

        if not Path(self.kcombu_exe).exists():
            raise FileNotFoundError(
                f"Could not find kcombu executable at {self.kcombu_exe}. Make sure it is installed."
            )

    def _set_fkcombu_exe(self) -> str:
        if shutil.which("fkcombu") is not None:
            fkcombu_exe = Path(shutil.which("fkcombu"))
            logger.debug(f"Using fkcombu executable from {fkcombu_exe}")
        else:
            raise ImportError(
                "fkcombu not installed. To install, run the command `micromamba install michellab::fkcombu`."
            )
        return str(fkcombu_exe)

    def _process_fkparams(self, fkparams: dict[str, Any], connectivity, top_constraint_tol) -> dict[str, str]:
        """Process and validate fkcombu parameters, called during initialization."""
        valid_params = {
            "P": None,  # Protein file
            "E": ("A", "V"),  # Energy calculation method; [A]tom-match, [V]olume-overlap
            "S": ("F", "R", "N"),  # Search method; [F]lexible, [R]igid, [N]othing
            "SD": ("T", "F"),  # Perform Gradient-based Steepest Descent fitting
        }
        processed_params = {}

        connectivity_supported = ("C", "S", "I", "T")
        if connectivity not in connectivity_supported:
            raise ValueError(
                f"Invalid connectivity method '{connectivity}'. Supported methods: {connectivity_supported}"
            )
        else:
            processed_params["con"] = connectivity
            if connectivity == "T" and top_constraint_tol is not None:
                processed_params["mtd"] = str(top_constraint_tol)

        for key, value in fkparams.items():
            if key in valid_params:

                if key == "P":
                    if value is None:
                        continue
                    elif not Path(value).exists():
                        raise FileNotFoundError(f"Protein file '{value}' not found.")
                    processed_params[key] = str(value) if isinstance(value, Path) else value
                elif value in valid_params[key]:
                    processed_params[key] = value
                elif isinstance(value, bool):
                    processed_params[key] = "T" if value else "F"
                else:
                    raise ValueError(
                        f"Invalid value '{value}' for parameter '{key}'. Valid values: {valid_params[key]}"
                    )
            else:
                logger.info(f"Parsing keyword parameter: {key}; {value}")
                if isinstance(value, bool):
                    processed_params[key] = "T" if value else "F"
                else:
                    processed_params[key] = str(value)

        return processed_params

    def _build_fkcombu_command(self) -> str:
        """Build the base fkcombu command string."""
        command = [self.kcombu_exe]
        for key, value in self.fkparams.items():
            if key in ("E", "S"):
                command.extend([f"-{key}", value])
            else:
                command.extend([f"-{key}", value])
        return command

    def _setup_temp_dir(self) -> Path:
        """Set up a temporary directory for alignment operations."""
        self.temp_dir = tempfile.TemporaryDirectory(prefix="ligand_alignment_")
        temp_path = Path(self.temp_dir.name)
        logger.debug(f"Temporary directory created for the alignment at {temp_path}")
        self.write_sdf_separate(temp_path)
        return temp_path

    def _run_kcombu(self, mol_path, reference_path, output_file) -> None:
        """Run fkcombu to align a molecule to a reference and write the aligned molecule
        to an output file.

        Args:
            mol_path: path to the molecule to be aligned.
            reference_path: path to the reference molecule.
            output_file: path to save the aligned molecule to.
        """

        command = self.fkcombu_command + [
            "-T",
            str(mol_path),
            "-R",
            str(reference_path),
            "-osdfT",
            str(output_file),
        ]
        logger.debug(f"Running kcombu: {' '.join(command)}")
        try:
            subprocess.run(command, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"kcombu process error: {e}")
        except Exception as e:
            logger.error(f"Unexpected error running kcombu: {e}")

    def _load_aligned_molecules(self, temp_path: Path) -> list[Molecule]:
        """Load aligned molecules from temporary files into memory.

        Args:
            temp_path: Path to the temporary directory containing the aligned molecules.

        Returns:
            List of aligned Molecule objects
        """
        for name in self.lig_names:
            if name != self.refname:
                aligned_file = temp_path / f"{name}{self.ALIGNED_SUFFIX}{self.SDF_EXTENSION}"
                original_file = temp_path / f"{name}{self.SDF_EXTENSION}"
                if aligned_file.exists():
                    self._transfer_sdf_metadata(original_file, aligned_file)
                    aligned_mol = Molecule.from_file(str(aligned_file))

                    self.aligned_molecules[name] = aligned_mol
                else:
                    logger.warning(f"Aligned file not found for {name}")

        # Add the reference molecule to aligned_molecules
        self.aligned_molecules[self.refname] = self.reference_mol
        return self.aligned_molecules

    @staticmethod
    def _transfer_charges_metadata(molA, molB):
        """Transfer the charges from one molecule to another based on the MCS. Necessary
        step because the aligned molecules created by kcombu don't have the charges' metadata
        from the reference molecule, causing the implicit hydrogens to be incorrectly assinged.
        MCS is used to make sure the atoms are mapped correctly.

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
        mcs_smarts = mcs_result.smartsString
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)

        matchA = molA.GetSubstructMatch(mcs_mol)
        matchB = molB.GetSubstructMatch(mcs_mol)

        logger.trace("Mapping of atoms:")
        for a, b in zip(matchA, matchB):  # trace the mapping; used for debugging
            logger.trace(f"MolA atom {a} maps to MolB atom {b}")

        for a, b in zip(matchA, matchB):  # iterate atoms and transfer charges
            atomA = molA.GetAtomWithIdx(a)
            atomB = molB.GetAtomWithIdx(b)
            formal_charge = atomA.GetFormalCharge()
            atomB.SetFormalCharge(formal_charge)
            if atomA.HasProp("_GasteigerCharge"):
                gaister_charge = atomA.GetProp("_GasteigerCharge")
                atomB.SetProp("_GasteigerCharge", gaister_charge)
            # if atom is oxygen, check for valence 3 and remove bound hydrogens
            if atomA.GetAtomicNum() == 8 and atomA.GetTotalDegree() == 3:
                atomB.SetNumExplicitHs(0)
                atomB.UpdatePropertyCache()
        return molA, molB

    def _transfer_sdf_metadata(self, original_file: Path, aligned_file: Path):
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
        for original_mol, aligned_mol in zip(original_supplier, aligned_supplier):
            if original_mol is not None and aligned_mol is not None:
                Chem.SanitizeMol(aligned_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)

                aligned_mol.SetProp("_Name", original_mol.GetProp("_Name"))
                for prop_name in original_mol.GetPropNames():  # copy all properties
                    prop_value = original_mol.GetProp(prop_name)
                    aligned_mol.SetProp(prop_name, prop_value)

                original_mol = Chem.RemoveHs(original_mol)  # rm H's to transfer charges transfer the charges
                logger.debug(f"Transferring metadata from {original_file} to {aligned_file}")
                original_mol, aligned_mol = self._transfer_charges_metadata(original_mol, aligned_mol)
                aligned_mols.append(aligned_mol)

        aligned_writer = Chem.SDWriter(str(aligned_file))
        for mol in aligned_mols:
            aligned_writer.write(mol)
        aligned_writer.close()

    def align_single_molecule(
        self, molecule: Union[str, Molecule], reference: Union[str, Molecule]
    ) -> Molecule:
        """
        Align a single molecule to a reference molecule.

        Args:
            molecule: The molecule to align (either a Molecule object or a name of a molecule in self.molecules)
            reference: The reference molecule (either a Molecule object or a name of a molecule in self.molecules)

        Returns:
            The aligned Molecule object
        """
        temp_path = self._setup_temp_dir()

        if isinstance(molecule, str):
            molecule = next((mol for mol in self.molecules if mol.name == molecule), None)
            if molecule is None:
                raise ValueError(f"Molecule {molecule} not found in self.molecules")
        if isinstance(reference, str):
            reference = next((mol for mol in self.molecules if mol.name == reference), None)
            if reference is None:
                raise ValueError(f"Reference molecule {reference} not found in self.molecules")

        molecule_path = temp_path / f"to_align{self.SDF_EXTENSION}"
        reference_path = temp_path / f"reference{self.SDF_EXTENSION}"
        molecule.to_file(molecule_path, file_format="sdf")
        reference.to_file(reference_path, file_format="sdf")

        output_path = temp_path / f"aligned{self.SDF_EXTENSION}"

        self._run_kcombu(molecule_path, reference_path, output_path)

        self._transfer_sdf_metadata(molecule_path, output_path)
        aligned_molecule = Molecule.from_file(str(output_path))

        if self.temp_dir:
            self.temp_dir.cleanup()

        return aligned_molecule

    def kcombu_align(
        self, reference: Union[str, Molecule], molecules_to_align: Optional[list[Union[str, Molecule]]] = None
    ) -> list[Molecule]:
        """
        Aligns the specified molecules to a reference molecule using kcombu. The aligned molecules returned
        by this function are stored in the `aligned_molecules` attribute for convenience. Upon using a reference
        molecule, the original reference molecule is stored in the `aligned_molecules` attribute as well.

        Args:
            reference: The reference molecule (either a Molecule object or a name of a molecule
                in self.molecules)
            molecules_to_align: List of molecules to align (either Molecule objects or names of
                molecules in self.molecules). If None, aligns all molecules except the reference.

        Returns:
            List of aligned Molecule objects
        """
        temp_path = self._setup_temp_dir()

        if isinstance(reference, str):  # prepare the reference
            self.refname = reference
            self.reference_mol = next((mol for mol in self.molecules if mol.name == reference), None)
            if self.reference_mol is None:
                raise ValueError(f"Reference molecule {reference} not found in self.molecules")
        else:
            self.refname = reference.name
            self.reference_mol = reference

        self.reference_path = temp_path / f"{self.refname}{self.SDF_EXTENSION}"
        self.reference_mol.to_file(str(self.reference_path), file_format="sdf")

        if molecules_to_align is None:  # prepare molecules to align
            molecules_to_align = [mol for mol in self.molecules if mol.name != self.refname]
        else:
            molecules_to_align = [
                next((mol for mol in self.molecules if mol.name == m), m) if isinstance(m, str) else m
                for m in molecules_to_align
            ]

        logger.info(f"Aligning {len(molecules_to_align)} molecules to ref `{self.refname}` with kcombu...")

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = []
            for mol in molecules_to_align:
                mol_path = temp_path / f"{mol.name}{self.SDF_EXTENSION}"
                mol.to_file(str(mol_path), file_format="sdf")
                output_file = mol_path.with_stem(f"{mol.name}{self.ALIGNED_SUFFIX}")
                partial_func = partial(
                    self._run_kcombu,
                    mol_path=mol_path,
                    reference_path=self.reference_path,
                    output_file=output_file,
                )
                futures.append(executor.submit(partial_func))

            for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning ligands"):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error in kcombu alignment: {e}")

        aligned_ligands = self._load_aligned_molecules(temp_path)
        self.cleanup()
        return aligned_ligands

    def output_aligned_ligands(
        self, output_name: str, ref_names: Optional[Union[str, list[str]]] = None
    ) -> None:
        """
        Write the aligned molecules to a single .sdf file, optionally including the original reference ligand(s).

        Args:
            output_name (str): Name of the output .sdf file.
            ref_names: Name(s) of the reference ligand(s) to include in their original conformation.
                Note that reference ligands are aligned when using the `kcombu_align` method.
                If None, only molecules within self.aligned_molecules will be written. Defaults to None.

        Raises:
            ValueError: If a specified reference name is not found in the original molecules.
        """
        writer = Chem.SDWriter(output_name)

        # Write aligned molecules
        for _, mol in self.aligned_molecules.items():
            writer.write(mol.to_rdkit())

        # Process and write reference ligands if specified
        if ref_names:
            if isinstance(ref_names, str):
                ref_names = [ref_names]

            for ref_name in ref_names:
                original_ref = self.get_molecule(ref_name, aligned=False)
                if original_ref is None:
                    raise ValueError(f"Reference molecule '{ref_name}' not found in original molecules.")
                writer.write(original_ref.to_rdkit())
                logger.info(f"Original reference '{ref_name}' added to output file.")

        writer.close()
        logger.info(f"Aligned molecules and specified references written to {output_name}")

        if self.temp_dir:
            self.temp_dir.cleanup()
            logger.info("Temporary directory cleaned up")

    def get_molecule(self, name: str, aligned: bool = True) -> Optional[Molecule]:
        """
        Retrieve a molecule by name, either aligned or original.

        Args:
            name: The name of the molecule to retrieve.
            aligned: If True, return the aligned molecule; if False, return the original.

        Returns:
            The requested Molecule object, or None if not found.
        """
        if aligned:
            return self.aligned_molecules.get(name)
        return next((mol for mol in self.molecules if mol.name == name), None)

    def cleanup(self) -> None:
        """Clean up the temporary directory used for the ligand alignment operations."""
        logger.debug(f"Temporary directory {self.temp_dir.name} cleaned up")
        if self.temp_dir:
            self.temp_dir.cleanup()
            self.temp_dir = None
