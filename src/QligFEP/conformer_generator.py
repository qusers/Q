"""Generate 3D conformers biased toward reference structures using constrained embedding.

For FEP ligand preparation, conformers should match the binding pocket geometry.
Standard conformer generators (ETKDGv3) produce arbitrary orientations that may
have ring flips or rotamer differences from the reference pose. This module pins
shared scaffold atoms to reference coordinates during embedding, guaranteeing
conformational consistency.
"""

from pathlib import Path
from typing import Union

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

from .chemIO import MoleculeIO
from .logger import logger


class BiasedConformerGenerator(MoleculeIO):
    """Generate 3D conformers with shared scaffolds pinned to reference coordinates.

    For each input SMILES, finds the reference molecule with the largest MCS
    (maximum common substructure) and generates a conformer where the MCS atoms
    are constrained to the reference's 3D coordinates. This prevents ring flips
    and rotamer mismatches that plague unconstrained conformer generation.

    Usage:
        >>> gen = BiasedConformerGenerator("ref_ligand.sdf")
        >>> gen.add("ligA", "c1ccccc1N", dg_value=-8.5)
        >>> gen.add("ligB", "c1ccccc1Cl", dg_value=-7.2)
        >>> gen.generate()
        >>> gen.write_to_single_sdf("ligands.sdf")

    The output is a MoleculeIO-compatible object: iterate, index by name,
    write to SDF/PDB, etc.
    """

    def __init__(
        self,
        references: Union[str, Path, MoleculeIO],
        reindex_hydrogens: bool = True,
        mcs_atom_compare: str = "any",
        mcs_bond_compare: str = "any",
        random_seed: int = 42,
        mcs_timeout: int = 30,
    ):
        """Initialize with one or more reference molecules that have 3D coordinates.

        Args:
            references: SDF file path or MoleculeIO containing reference structures
                with known 3D coordinates (e.g., crystal poses, docked poses).
            reindex_hydrogens: Reindex H atoms to end of atom list (needed for
                Q's restraint setting). Defaults to True.
            mcs_atom_compare: How strictly to match atoms when selecting the
                best reference molecule.
                'any' = ignore elements (topology only). Finds the largest
                    possible MCS, useful when atoms get transmuted in FEP.
                'elements' = match by element (C to C, N to N).
                Defaults to 'any'.
                Note: coordinate pinning during embedding always uses
                element-strict matching to avoid cross-element misplacement.
            mcs_bond_compare: How strictly to match bonds in MCS calculation.
                'any' = ignore bond types.
                'order' = match single/double/triple.
                Defaults to 'any'.
            random_seed: Seed for conformer generation reproducibility.
        """
        # Were we don't call MoleculeIO.__init__. instead, we populate molecules ourselves
        self._reindex_hydrogens = reindex_hydrogens
        self.lig = None
        self.lig_files = []
        self.molecules = []
        self.lig_names = []
        self.sdf_contents = {}
        self.random_seed = random_seed

        # Load references
        if isinstance(references, MoleculeIO):
            self._references = references
        else:
            self._references = MoleculeIO(str(references), reindex_hydrogens=False)

        # Pre-compute reference RDKit mols (without H for MCS)
        self._ref_mols_noH = {}
        self._ref_mols_withH = {}
        for name, mol in self._references:
            rdmol = mol.to_rdkit()
            self._ref_mols_withH[name] = rdmol
            self._ref_mols_noH[name] = Chem.RemoveHs(rdmol)

        # MCS comparison settings
        atom_compare_map = {
            "elements": rdFMCS.AtomCompare.CompareElements,
            "any": rdFMCS.AtomCompare.CompareAny,
        }
        bond_compare_map = {
            "any": rdFMCS.BondCompare.CompareAny,
            "order": rdFMCS.BondCompare.CompareOrder,
        }
        self._atom_compare = atom_compare_map.get(mcs_atom_compare, rdFMCS.AtomCompare.CompareElements)
        self._bond_compare = bond_compare_map.get(mcs_bond_compare, rdFMCS.BondCompare.CompareAny)
        self._mcs_timeout = mcs_timeout

        # Queue of molecules to generate
        self._queue: list[dict] = []

        logger.info(
            f"BiasedConformerGenerator initialized with {len(self._ref_mols_noH)} reference(s): "
            f"{', '.join(self._ref_mols_noH.keys())}"
        )

    def add(self, name: str, smiles: str, **properties):
        """Add a molecule to the generation queue.

        Args:
            name: Molecule name (used in output SDF).
            smiles: SMILES string for the molecule.
            **properties: Additional SDF properties to set (e.g., dg_value=-8.5).
        """
        self._queue.append({"name": name, "smiles": smiles, "properties": properties})

    def _find_best_reference(self, mol_noH: Chem.Mol) -> tuple[str, Chem.Mol, int]:
        """Find the reference with the largest MCS to the query molecule.

        Returns:
            (ref_name, mcs_result, mcs_num_atoms)
        """
        best_ref = None
        best_mcs = None
        best_size = 0

        for ref_name, ref_noH in self._ref_mols_noH.items():
            mcs = rdFMCS.FindMCS(
                [ref_noH, mol_noH],
                atomCompare=self._atom_compare,
                bondCompare=self._bond_compare,
                ringMatchesRingOnly=True,
                timeout=self._mcs_timeout,
            )
            if mcs.numAtoms > best_size:
                best_size = mcs.numAtoms
                best_mcs = mcs
                best_ref = ref_name

        return best_ref, best_mcs, best_size

    def _build_core_with_coords(self, mcs, ref_noH, ref_conf):
        """Build an RDKit molecule from an MCS with reference 3D coordinates.

        Args:
            mcs: rdFMCS result object.
            ref_noH: Reference molecule without hydrogens.
            ref_conf: Conformer of the reference molecule (with H).

        Returns:
            (core_mol, ref_match) or (None, None) if substructure match fails.
        """
        core = Chem.MolFromSmarts(mcs.smartsString)
        ref_match = ref_noH.GetSubstructMatch(core)
        if not ref_match:
            return None, None

        core_mol = Chem.RWMol(core)
        core_conf = Chem.Conformer(core_mol.GetNumAtoms())
        for i, ri in enumerate(ref_match):
            core_conf.SetAtomPosition(i, ref_conf.GetAtomPosition(ri))
        core_mol.AddConformer(core_conf)
        return core_mol.GetMol(), ref_match

    def _constrained_embed(self, smiles: str, ref_name: str, mcs) -> Chem.Mol:
        """Generate a 3D conformer with MCS atoms pinned to reference coordinates.

        Uses element-strict MCS for coordinate pinning to avoid mapping atoms
        to positions of different elements (e.g., pinning a C to where an F
        sits in the reference). The topology-only MCS passed in via `mcs` is
        only used as a fallback for the coordMap strategy.

        Strategies tried in order:
        1. ConstrainedEmbed with element-strict MCS
        2. ConstrainedEmbed with topology-only MCS (only if it differs from #1)
        3. coordMap embed with topology-only MCS positions
        4. Unconstrained embed + centroid translation
        """
        ref_withH = self._ref_mols_withH[ref_name]
        ref_noH = self._ref_mols_noH[ref_name]
        ref_conf = ref_withH.GetConformer()

        mol_noH = Chem.RemoveHs(Chem.MolFromSmiles(smiles))

        # Element-strict MCS: only pins atoms that match by element, preventing
        # cross-element coordinate pinning that silently produces bad geometry
        mcs_elements = rdFMCS.FindMCS(
            [ref_noH, mol_noH],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True,
            timeout=self._mcs_timeout,
        )

        # Strategy 1: ConstrainedEmbed with element-strict MCS
        if mcs_elements.numAtoms >= 3:
            core_mol, _ = self._build_core_with_coords(mcs_elements, ref_noH, ref_conf)
            if core_mol is not None:
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                try:
                    AllChem.ConstrainedEmbed(mol, core_mol, randomseed=self.random_seed)
                    if mcs_elements.numAtoms < mcs.numAtoms:
                        logger.info(
                            f"ConstrainedEmbed with element-strict MCS "
                            f"({mcs_elements.numAtoms}/{mcs.numAtoms} atoms pinned)"
                        )
                    return mol
                except (ValueError, RuntimeError):
                    logger.debug(
                        f"ConstrainedEmbed with element-strict MCS failed for {ref_name} "
                        f"({mcs_elements.numAtoms} atoms)"
                    )

        # Strategy 2: ConstrainedEmbed with topology-only MCS (if it covers more atoms)
        if mcs.numAtoms > mcs_elements.numAtoms:
            core_mol, ref_match = self._build_core_with_coords(mcs, ref_noH, ref_conf)
            if core_mol is not None:
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                try:
                    AllChem.ConstrainedEmbed(mol, core_mol, randomseed=self.random_seed)
                    logger.debug(
                        f"ConstrainedEmbed with topology MCS "
                        f"({mcs.numAtoms} atoms pinned)"
                    )
                    return mol
                except (ValueError, RuntimeError):
                    logger.debug(
                        f"ConstrainedEmbed with topology MCS also failed for {ref_name}"
                    )

        # Strategy 3: EmbedMolecule with coordMap (looser constraints)
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        mol_noH_fresh = Chem.RemoveHs(mol)

        core = Chem.MolFromSmarts(mcs.smartsString)
        ref_match = ref_noH.GetSubstructMatch(core)
        mol_match = mol_noH_fresh.GetSubstructMatch(core)

        if ref_match and mol_match:
            heavy_to_full = {}
            heavy_idx = 0
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 1:
                    heavy_to_full[heavy_idx] = atom.GetIdx()
                    heavy_idx += 1

            coord_map = {}
            for core_idx, mol_noH_idx in enumerate(mol_match):
                if mol_noH_idx in heavy_to_full:
                    coord_map[heavy_to_full[mol_noH_idx]] = ref_conf.GetAtomPosition(ref_match[core_idx])

            if coord_map:
                result = AllChem.EmbedMolecule(
                    mol, coordMap=coord_map, randomSeed=self.random_seed,
                    useRandomCoords=True, enforceChirality=True,
                )
                if result != -1:
                    AllChem.MMFFOptimizeMolecule(mol)
                    logger.info(
                        f"Used coordMap fallback for {ref_name} "
                        f"({len(coord_map)} pinned atoms)"
                    )
                    return mol
                logger.debug(f"coordMap embed also failed for {ref_name}")

        # Strategy 4: unconstrained embed + centroid translation
        logger.warning(
            f"All constrained strategies failed for {ref_name}, "
            f"using unconstrained embed + centroid translation"
        )
        return self._unconstrained_embed(mol, ref_conf=ref_conf)

    def _unconstrained_embed(self, mol: Chem.Mol, ref_conf=None) -> Chem.Mol:
        """Fallback: standard ETKDGv3 embedding without constraints.

        If ref_conf is provided, translates the generated conformer so its
        centroid matches the reference centroid. This ensures downstream
        alignment tools (fkcombu) start from the correct spatial region
        instead of from the origin.
        """
        params = AllChem.ETKDGv3()
        params.randomSeed = self.random_seed
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            params.useRandomCoords = True
            AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)

        # Translate to reference centroid if available
        if ref_conf is not None and mol.GetNumConformers() > 0:
            import numpy as np

            ref_cent = np.mean(
                [[ref_conf.GetAtomPosition(i).x,
                  ref_conf.GetAtomPosition(i).y,
                  ref_conf.GetAtomPosition(i).z]
                 for i in range(ref_conf.GetNumAtoms())],
                axis=0,
            )
            mol_conf = mol.GetConformer()
            mol_cent = np.mean(
                [[mol_conf.GetAtomPosition(i).x,
                  mol_conf.GetAtomPosition(i).y,
                  mol_conf.GetAtomPosition(i).z]
                 for i in range(mol.GetNumAtoms())],
                axis=0,
            )
            shift = ref_cent - mol_cent
            for i in range(mol.GetNumAtoms()):
                p = mol_conf.GetAtomPosition(i)
                mol_conf.SetAtomPosition(
                    i, Chem.rdGeometry.Point3D(p.x + shift[0], p.y + shift[1], p.z + shift[2]),
                )

        return mol

    def generate(self):
        """Generate conformers for all queued molecules.

        For each molecule, finds the reference with the largest MCS and pins
        shared atoms to reference coordinates during embedding. Results are
        stored in self.molecules and self.lig_names (MoleculeIO interface).
        """
        # Include reference molecules in the output (they already have 3D coords)
        for ref_name, ref_mol in self._references:
            rdmol = ref_mol.to_rdkit()
            off_mol = self._rdkit_to_openff(rdmol, hydrogens_are_explicit=True)
            off_mol.name = ref_name
            if self._reindex_hydrogens:
                off_mol = self._force_H_reindexing(off_mol)
                off_mol.name = ref_name
            self.molecules.append(off_mol)
            self.lig_names.append(ref_name)

        for entry in self._queue:
            name = entry["name"]
            smiles = entry["smiles"]
            properties = entry["properties"]

            mol_noH = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
            ref_name, mcs, mcs_size = self._find_best_reference(mol_noH)

            if mcs_size < 3:
                logger.warning(
                    f"{name}: MCS with best reference ({ref_name}) has only {mcs_size} atoms, using unconstrained embed"
                )
                rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                rdmol = self._unconstrained_embed(rdmol)
            else:
                logger.debug(f"{name}: best reference={ref_name}, MCS={mcs_size} atoms")
                rdmol = self._constrained_embed(smiles, ref_name, mcs)

            # Set molecule name and properties
            rdmol.SetProp("_Name", name)
            for prop_key, prop_val in properties.items():
                rdmol.SetProp(prop_key, str(prop_val))

            # Convert to OpenFF
            off_mol = self._rdkit_to_openff(rdmol, hydrogens_are_explicit=True)
            off_mol.name = name
            if self._reindex_hydrogens:
                off_mol = self._force_H_reindexing(off_mol)
                off_mol.name = name

            self.molecules.append(off_mol)
            self.lig_names.append(name)

        self.parse_sdf_contents()
        logger.info(
            f"Generated {len(self._queue)} conformers ({len(self.molecules)} total including references)"
        )
