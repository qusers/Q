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
from rdkit.Chem import rdFMCS, rdShapeAlign
from tqdm import tqdm

from .chemIO import MoleculeIO
from .logger import logger

# def find_best_reference(molecules: list[Molecule], lig_names: list[str]) -> str:
#     """Find the molecule with the highest aggregate MCS similarity to all others.

#     Computes a pairwise MCS similarity matrix and returns the name of the molecule
#     whose total similarity to all others is highest, making it a good alignment reference.

#     Args:
#         molecules: List of Molecule objects to compare.
#         lig_names: List of molecule names corresponding to `molecules`.

#     Returns:
#         Name of the best reference molecule.
#     """
#     from MolClusterkit.mcs import MCSClustering

#     smiles_list = [mol.to_smiles() for mol in molecules]
#     mcs_kwargs = {
#         "bondCompare": "CompareOrderExact",
#         "ringCompare": "StrictRingFusion",
#         "ringMatchesRingOnly": True,
#         "completeRingsOnly": True,
#         "timeout": 20,
#     }
#     mcs_cluster = MCSClustering(smiles_list, **mcs_kwargs)
#     logger.info(
#         "Scoring ligands based on Maximum Common Substructure to find the best reference."
#     )
#     _, simi_matrix = mcs_cluster.compute_similarity_matrix()

#     highest_score_idx = simi_matrix.sum(axis=1).argmax()
#     best_ref = lig_names[highest_score_idx]
#     logger.info(f"Best reference molecule: {best_ref}")
#     return best_ref


class fkcombuLigandAligner(MoleculeIO):
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
        atom_type: str = "X",
        bond_type: str = "X",
        **fkcombu_params,
    ):
        """
        Initialize a fkcombuLigandAligner object for aligning ligands based on three-dimensional coordinates
        using the fkcombu program. Additional parameters can be passed to fkcombu through keyword arguments.

        For more information on fkcombu, see their docs:
            https://pdbj.org/kcombu/doc/README_fkcombu.html
        For information on the `connectivity`, `top_constraint_tol`, `atom_type` parameters, see:
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
            atom_type: fkcombu atom type classification for MCS matching. Controls how strictly atoms
                must match during MCS calculation. Options:
                - `X`: ignore atom types entirely (topology-only matching, best for FEP)
                - `E`: element-only (C, N, O must match)
                - `T`: element + bond count + ring status
                - `K`: combu-recommend 12-type classification (fkcombu default)
                Defaults to `X` because FEP transmutes atom types — the MCS should capture the
                shared scaffold regardless of element identity.
            bond_type: fkcombu bond type classification for MCS matching. Options:
                - `X`: ignore bond types (default, good for FEP)
                - `C`: care about SDF/MOL2 bond types
                - `R`: distinguish rotatable vs non-rotatable
                Defaults to `X`.

        Raises:
            FileNotFoundError: If the kcombu executable is not found.
        """
        super().__init__(lig, pattern=pattern, reindex_hydrogens=reindex_hydrogens)
        self.kcombu_exe = self._set_fkcombu_exe()
        self.n_threads = n_threads
        self.reference_mol: Optional[Molecule] = None
        self.aligned_molecules: dict[str, Molecule] = {}
        self.alignment_scores: dict[str, dict[str, float]] = {}
        self.temp_dir: Optional[tempfile.TemporaryDirectory] = None
        self.fkparams = self._process_fkparams(
            {"P": protein, "E": energy.upper(), "S": search.upper(), "SD": steep_descend,
             "at": atom_type.upper(), "bo": bond_type.upper(), **fkcombu_params},
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

    @staticmethod
    def options() -> dict[str, str]:
        """Parse the fkcombu -h output and return a dict mapping option flags to descriptions.

        Returns:
            Dict mapping option names (e.g. "-T", "-E") to their descriptions.

        Raises:
            ImportError: If fkcombu is not installed.
        """
        fkcombu_exe = shutil.which("fkcombu")
        if fkcombu_exe is None:
            raise ImportError(
                "fkcombu not installed. To install, run the command `micromamba install michellab::fkcombu`."
            )

        result = subprocess.run([fkcombu_exe, "-h"], capture_output=True, text=True)
        help_text = result.stdout + result.stderr

        options = {}
        for line in help_text.splitlines():
            line = line.strip()
            if not line or line.startswith(">") or line.startswith("$") or line.startswith("<"):
                continue
            # Option lines start with a dash, e.g. " -T   :target molecule ..."
            if line.startswith("-"):
                parts = line.split(":", 1)
                if len(parts) == 2:
                    flag = parts[0].strip().split()[0]  # first token is the flag
                    description = parts[1].strip()
                    options[flag] = description
        return options

    def _process_fkparams(self, fkparams: dict[str, Any], connectivity, top_constraint_tol) -> dict[str, str]:
        """Process and validate fkcombu parameters, called during initialization."""
        valid_params = {
            "P": None,  # Protein file
            "E": ("A", "V"),  # Energy calculation method; [A]tom-match, [V]olume-overlap
            "S": ("F", "R", "N"),  # Search method; [F]lexible, [R]igid, [N]othing
            "SD": ("T", "F"),  # Perform Gradient-based Steepest Descent fitting
            "at": ("X", "E", "T", "K", "D"),  # Atom type classification for MCS
            "bo": ("X", "C", "R"),  # Bond type classification for MCS
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

    # Scores that are only meaningful when MCS atom matching was performed (NATOM_MATCH > 0)
    _MCS_DEPENDENT_SCORES = ("TANIMOTO", "RMSD_MATCH", "NATOM_MATCH")
    # Scores to extract from fkcombu output
    _SCORE_KEYS = ("TANIMOTO", "TANIMOTO_VOL", "RMSD_MATCH", "NATOM_MATCH", "Etotal")

    def _parse_kcombu_output(self, stdout: str) -> dict[str, float]:
        """Parse fkcombu key-value stdout lines into a scores dict.

        Only non-comment lines (not starting with #) are parsed. When NATOM_MATCH is 0
        (no MCS was computed, e.g. with -E V), MCS-dependent scores are excluded since
        their zero values are not meaningful.

        Args:
            stdout: Raw stdout string from fkcombu.

        Returns:
            Dict mapping score names to their values.
        """
        kv_pairs = {}
        for line in stdout.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 2:
                kv_pairs[parts[0]] = parts[1]

        scores = {}
        for key in self._SCORE_KEYS:
            if key in kv_pairs:
                if key == "NATOM_MATCH":
                    scores[key] = int(kv_pairs[key])
                else:
                    scores[key] = float(kv_pairs[key])

        # Exclude MCS-dependent scores when no MCS was computed
        if scores.get("NATOM_MATCH", 0) == 0:
            for key in self._MCS_DEPENDENT_SCORES:
                scores.pop(key, None)

        return scores

    def _run_kcombu(self, mol_path, reference_path, output_file) -> dict[str, float]:
        """Run fkcombu to align a molecule to a reference and write the aligned molecule
        to an output file.

        Args:
            mol_path: path to the molecule to be aligned.
            reference_path: path to the reference molecule.
            output_file: path to save the aligned molecule to.

        Returns:
            Dict of alignment scores parsed from fkcombu stdout.
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
            result = subprocess.run(command, capture_output=True, text=True)
            return self._parse_kcombu_output(result.stdout + result.stderr)
        except subprocess.CalledProcessError as e:
            logger.error(f"kcombu process error: {e}")
        except Exception as e:
            logger.error(f"Unexpected error running kcombu: {e}")
        return {}

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
    ) -> tuple[Molecule, dict[str, float]]:
        """
        Align a single molecule to a reference molecule.

        Args:
            molecule: The molecule to align (either a Molecule object or a name of a molecule in self.molecules)
            reference: The reference molecule (either a Molecule object or a name of a molecule in self.molecules)

        Returns:
            Tuple of (aligned Molecule object, dict of alignment scores).
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

        scores = self._run_kcombu(molecule_path, reference_path, output_path)

        self._transfer_sdf_metadata(molecule_path, output_path)
        aligned_molecule = Molecule.from_file(str(output_path))

        if self.temp_dir:
            self.temp_dir.cleanup()

        return aligned_molecule, scores

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
            futures = {}
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
                futures[executor.submit(partial_func)] = mol.name

            for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning ligands"):
                mol_name = futures[future]
                try:
                    scores = future.result()
                    if scores:
                        self.alignment_scores[mol_name] = scores
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

        # Write aligned molecules with alignment scores as SDF properties
        for name, mol in self.aligned_molecules.items():
            rdkit_mol = mol.to_rdkit()
            if name in self.alignment_scores:
                for score_key, score_val in self.alignment_scores[name].items():
                    rdkit_mol.SetProp(f"Align_{score_key}", f"{score_val}")
            writer.write(rdkit_mol)

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
        if self.temp_dir:
            logger.debug(f"Temporary directory {self.temp_dir.name} cleaned up")
            self.temp_dir.cleanup()
            self.temp_dir = None


class rdkitLigandAligner(MoleculeIO):
    """Align ligands based on three-dimensional shape and pharmacophore features
    using RDKit's rdShapeAlign.

    This aligner maximizes volumetric and chemical feature overlap between molecules,
    producing shape and color Tanimoto scores comparable to OpenEye ROCS.

    Attributes:
        opt_param (float): Balance between shape (1.0) and color (0.0) for optimization.
            0.5 gives equal weight to both. Default is 0.5.
        max_preiters (int): Maximum iterations on all poses during the first optimization phase.
        max_postiters (int): Maximum iterations on the best poses during the second phase.
        reference_mol (Optional[Molecule]): Current reference molecule.
        aligned_molecules (dict[str, Molecule]): Aligned molecules keyed by name.
        alignment_scores (dict[str, tuple[float, float]]): Per-molecule (shape_tani, color_tani).
    """

    def __init__(
        self,
        lig,
        pattern: str = "*.sdf",
        reindex_hydrogens: bool = True,
        opt_param: float = 0.5,
        max_preiters: int = 10,
        max_postiters: int = 30,
    ):
        """Initialize a rdkitLigandAligner for aligning ligands using RDKit shape alignment.

        Args:
            lig: SDF file or directory containing SDF files with the molecules to align.
            pattern: Glob pattern for finding SDF files if `lig` is a directory. Defaults to "*.sdf".
            reindex_hydrogens: If True, ensure hydrogen atoms are at the end of the atom list.
                Defaults to True.
            opt_param: Balance of shape and color for optimization.
                0.0 = color only, 0.5 = equal weight, 1.0 = shape only. Defaults to 0.5.
            max_preiters: Maximum iterations during the first optimization phase on all poses.
                Defaults to 10.
            max_postiters: Maximum iterations during the second optimization phase on the best poses.
                Defaults to 30.
        """
        super().__init__(lig, pattern=pattern, reindex_hydrogens=reindex_hydrogens)
        self.opt_param = opt_param
        self.max_preiters = max_preiters
        self.max_postiters = max_postiters
        self.reference_mol: Optional[Molecule] = None
        self.aligned_molecules: dict[str, Molecule] = {}
        self.alignment_scores: dict[str, tuple[float, float]] = {}

    def _resolve_molecule(self, mol_or_name: Union[str, Molecule]) -> Molecule:
        """Resolve a molecule from a name string or return the Molecule directly.

        Args:
            mol_or_name: Either a molecule name (str) or a Molecule object.

        Returns:
            The resolved Molecule object.

        Raises:
            ValueError: If the name is not found in self.molecules.
        """
        if isinstance(mol_or_name, str):
            mol = next((m for m in self.molecules if m.name == mol_or_name), None)
            if mol is None:
                raise ValueError(f"Molecule '{mol_or_name}' not found in self.molecules")
            return mol
        return mol_or_name

    def _align_rdkit_mol(self, ref_rdkit: Chem.Mol, probe_rdkit: Chem.Mol) -> tuple[float, float]:
        """Align a probe RDKit molecule to a reference using shape alignment.

        The probe molecule is modified in place.

        Args:
            ref_rdkit: Reference RDKit molecule.
            probe_rdkit: Probe RDKit molecule (will be modified in place).

        Returns:
            Tuple of (shape_tanimoto, color_tanimoto) scores.
        """
        shape_tani, color_tani = rdShapeAlign.AlignMol(
            ref_rdkit,
            probe_rdkit,
            opt_param=self.opt_param,
            max_preiters=self.max_preiters,
            max_postiters=self.max_postiters,
        )
        return shape_tani, color_tani

    def align_single_molecule(
        self, molecule: Union[str, Molecule], reference: Union[str, Molecule]
    ) -> tuple[Molecule, float, float]:
        """Align a single molecule to a reference molecule.

        Args:
            molecule: The molecule to align (name or Molecule object).
            reference: The reference molecule (name or Molecule object).

        Returns:
            Tuple of (aligned_molecule, shape_tanimoto, color_tanimoto).
        """
        molecule = self._resolve_molecule(molecule)
        reference = self._resolve_molecule(reference)

        ref_rdkit = reference.to_rdkit()
        probe_rdkit = Chem.RWMol(molecule.to_rdkit())

        shape_tani, color_tani = self._align_rdkit_mol(ref_rdkit, probe_rdkit)

        aligned_mol = Molecule.from_rdkit(probe_rdkit, allow_undefined_stereo=True)
        aligned_mol.name = molecule.name

        return aligned_mol, shape_tani, color_tani

    def align(
        self, reference: Union[str, Molecule], molecules_to_align: Optional[list[Union[str, Molecule]]] = None
    ) -> dict[str, Molecule]:
        """Align molecules to a reference using shape+color overlap.

        Aligned molecules are stored in `self.aligned_molecules` and scores in
        `self.alignment_scores`.

        Args:
            reference: The reference molecule (name or Molecule object).
            molecules_to_align: List of molecules to align. If None, aligns all
                molecules except the reference.

        Returns:
            Dictionary of aligned Molecule objects keyed by name.
        """
        reference = self._resolve_molecule(reference)
        self.reference_mol = reference
        ref_rdkit = reference.to_rdkit()

        if molecules_to_align is None:
            molecules_to_align = [mol for mol in self.molecules if mol.name != reference.name]
        else:
            molecules_to_align = [self._resolve_molecule(m) for m in molecules_to_align]

        logger.info(
            f"Aligning {len(molecules_to_align)} molecules to ref `{reference.name}` "
            f"with rdShapeAlign (opt_param={self.opt_param})..."
        )

        for mol in tqdm(molecules_to_align, desc="Aligning ligands"):
            probe_rdkit = Chem.RWMol(mol.to_rdkit())

            try:
                shape_tani, color_tani = self._align_rdkit_mol(ref_rdkit, probe_rdkit)
            except Exception as e:
                logger.error(f"Failed to align {mol.name}: {e}")
                continue

            aligned_mol = Molecule.from_rdkit(probe_rdkit, allow_undefined_stereo=True)
            aligned_mol.name = mol.name

            self.aligned_molecules[mol.name] = aligned_mol
            self.alignment_scores[mol.name] = (shape_tani, color_tani)

            logger.debug(
                f"  {mol.name}: shape={shape_tani:.4f} color={color_tani:.4f} "
                f"combo={shape_tani + color_tani:.4f}"
            )

        # Include reference in aligned_molecules
        self.aligned_molecules[reference.name] = reference

        return self.aligned_molecules

    def output_aligned_ligands(
        self, output_name: str, ref_names: Optional[Union[str, list[str]]] = None
    ) -> None:
        """Write aligned molecules to a single SDF file, including alignment scores as SD properties.

        Args:
            output_name: Path to the output SDF file.
            ref_names: Name(s) of reference ligand(s) to include in their original conformation.
                If None, only aligned molecules are written. Defaults to None.

        Raises:
            ValueError: If a specified reference name is not found in the original molecules.
        """
        writer = Chem.SDWriter(output_name)

        for name, mol in self.aligned_molecules.items():
            rdkit_mol = mol.to_rdkit()
            if name in self.alignment_scores:
                shape_tani, color_tani = self.alignment_scores[name]
                rdkit_mol.SetProp("Align_ShapeTanimoto", f"{shape_tani:.4f}")
                rdkit_mol.SetProp("Align_ColorTanimoto", f"{color_tani:.4f}")
                rdkit_mol.SetProp("Align_TanimotoCombo", f"{shape_tani + color_tani:.4f}")
            writer.write(rdkit_mol)

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
        logger.info(f"Aligned molecules written to {output_name}")

    def get_molecule(self, name: str, aligned: bool = True) -> Optional[Molecule]:
        """Retrieve a molecule by name, either aligned or original.

        Args:
            name: The name of the molecule to retrieve.
            aligned: If True, return the aligned molecule; if False, return the original.

        Returns:
            The requested Molecule object, or None if not found.
        """
        if aligned:
            return self.aligned_molecules.get(name)
        return next((mol for mol in self.molecules if mol.name == name), None)
