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
            mcs_atom_compare: How strictly to match atoms in MCS calculation.
                'any' = ignore elements (topology only). Best for FEP where
                    atoms are transmuted — avoids ConstrainedEmbed failures
                    when ring sizes differ (e.g., phenyl vs naphthyl).
                'elements' = match by element (C to C, N to N).
                Defaults to 'any'.
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

    def _constrained_embed(self, smiles: str, ref_name: str, mcs) -> Chem.Mol:
        """Generate a 3D conformer with MCS atoms pinned to reference coordinates."""
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        mol_noH = Chem.RemoveHs(mol)

        ref_withH = self._ref_mols_withH[ref_name]
        ref_noH = self._ref_mols_noH[ref_name]
        ref_conf = ref_withH.GetConformer()

        core = Chem.MolFromSmarts(mcs.smartsString)
        ref_match = ref_noH.GetSubstructMatch(core)
        mol_match = mol_noH.GetSubstructMatch(core)

        if not ref_match or not mol_match:
            logger.warning(
                f"MCS substructure match failed for {ref_name}, falling back to unconstrained embed"
            )
            return self._unconstrained_embed(mol)

        # Build core molecule with reference 3D coordinates
        core_mol = Chem.RWMol(core)
        core_conf = Chem.Conformer(core_mol.GetNumAtoms())
        for i, ri in enumerate(ref_match):
            pos = ref_conf.GetAtomPosition(ri)
            core_conf.SetAtomPosition(i, pos)
        core_mol.AddConformer(core_conf)

        try:
            AllChem.ConstrainedEmbed(mol, core_mol.GetMol(), randomseed=self.random_seed)
        except ValueError:
            logger.warning(f"ConstrainedEmbed failed for {ref_name}, falling back to unconstrained embed")
            return self._unconstrained_embed(mol)

        return mol

    def _unconstrained_embed(self, mol: Chem.Mol) -> Chem.Mol:
        """Fallback: standard ETKDGv3 embedding without constraints."""
        params = AllChem.ETKDGv3()
        params.randomSeed = self.random_seed
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            params.useRandomCoords = True
            AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)
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
