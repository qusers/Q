"""Generate perturbation networks using konnektor with pluggable scorers and network topologies."""

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np
from kartograf import KartografAtomMapper, SmallMoleculeComponent
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

from ..chemIO import MoleculeIO
from ..logger import logger, setup_logger
from ..scoring import RestraintLomapScorer, RestraintScorer
from .lomap_wrap_cli import LomapWrap

NETWORK_GENERATORS = {
    "mst": "MinimalSpanningTreeNetworkGenerator",
    "rmst": "RedundantMinimalSpanningTreeNetworkGenerator",
    "star": "StarNetworkGenerator",
    "nedges": "NNodeEdgesNetworkGenerator",
    "cyclic": "CyclicNetworkGenerator",
}


class KonnektorWrap:
    """Generate perturbation networks using konnektor.

    Supports MST, redundant MST, star, n-edges, and cyclic network topologies
    with either RestraintSetter-based or kartograf volume ratio scoring.
    """

    def __init__(
        self,
        inp: str,
        out: Optional[str] = None,
        network: str = "mst",
        scorer: str = "combined",
        restraint_method: str = "heavyatom_p",
        processes: int = 1,
        log_level: str = "info",
        central_ligand: Optional[str] = None,
        n_redundancy: int = 2,
        connectivity: int = 3,
        separate_charges: bool = False,
        charge_changes_score: float = 0.0,
        exp_key: Optional[str] = None,
        self_solve: bool = False,
    ):
        self.inp = inp
        self.network_type = network
        self.scorer_type = scorer
        self.restraint_method = restraint_method
        self.processes = processes
        self.log_level = log_level
        self.central_ligand = central_ligand
        self.n_redundancy = n_redundancy
        self.connectivity = connectivity
        self.separate_charges = separate_charges
        self.charge_changes_score = charge_changes_score
        self.exp_key = exp_key
        self.self_solve = self_solve
        self.nodes = {}

        self.out = self._parse_output(out)
        self._sdf_dir = self._prepare_input()

    def _parse_output(self, output: Optional[str]) -> str:
        inpath = Path(self.inp)
        if output is None:
            if inpath.is_dir():
                return str(inpath)
            else:
                return str(inpath.parent / inpath.stem)
        return str(Path(output).absolute())

    def _prepare_input(self) -> str:
        """Ensure input is a directory of individual SDF files."""
        inpath = Path(self.inp)
        if inpath.is_dir():
            return str(inpath)
        elif inpath.suffix == ".sdf":
            outdir = Path(self.out)
            outdir.mkdir(parents=True, exist_ok=True)
            handler = MoleculeIO(self.inp)
            logger.info(f"Writing {self.inp} to separate .sdf files in {self.out}.")
            handler.write_sdf_separate(self.out)
            return self.out
        else:
            raise ValueError("Input must be an .sdf file or a directory containing .sdf files.")

    def _build_scorer(self):
        """Create the scoring function for konnektor."""
        if self.scorer_type == "kartograf":
            from kartograf.atom_mapping_scorer import MappingVolumeRatioScorer

            return MappingVolumeRatioScorer()
        elif self.scorer_type == "combined":
            return RestraintLomapScorer(
                restraint_method=self.restraint_method,
                charge_changes_score=self.charge_changes_score,
            )
        else:
            return RestraintScorer(restraint_method=self.restraint_method)

    def _build_generator(self, mapper, scorer):
        """Create the network generator based on network type."""
        from konnektor.network_planners import (
            CyclicNetworkGenerator,
            MinimalSpanningTreeNetworkGenerator,
            NNodeEdgesNetworkGenerator,
            RedundantMinimalSpanningTreeNetworkGenerator,
            StarNetworkGenerator,
        )

        common_kwargs = dict(mappers=[mapper], scorer=scorer, n_processes=self.processes)

        generators = {
            "mst": lambda: MinimalSpanningTreeNetworkGenerator(**common_kwargs),
            "rmst": lambda: RedundantMinimalSpanningTreeNetworkGenerator(
                **common_kwargs,
                n_redundancy=self.n_redundancy,
            ),
            "star": lambda: StarNetworkGenerator(**common_kwargs),
            "nedges": lambda: NNodeEdgesNetworkGenerator(
                **common_kwargs,
                target_component_connectivity=self.connectivity,
            ),
            "cyclic": lambda: CyclicNetworkGenerator(**common_kwargs),
        }
        if self.network_type not in generators:
            raise ValueError(
                f"Unknown network type '{self.network_type}'. Choose from: {list(generators.keys())}"
            )
        return generators[self.network_type]()

    def _load_components(self) -> list[SmallMoleculeComponent]:
        """Load SDF files from directory as SmallMoleculeComponents."""
        sdf_dir = Path(self._sdf_dir)
        sdf_files = sorted(sdf_dir.glob("*.sdf"))
        if len(sdf_files) < 2:
            raise ValueError(f"Need at least 2 SDF files, found {len(sdf_files)} in {sdf_dir}")

        components = []
        for sdf_file in sdf_files:
            comp = SmallMoleculeComponent.from_sdf_file(str(sdf_file))
            components.append(comp)

            # Build node data from the SDF file
            name = sdf_file.stem
            rdmol = comp.to_rdkit()
            smiles = Chem.MolToSmiles(rdmol)
            formal_charge = Chem.GetFormalCharge(rdmol)
            extra_data = LomapWrap.extract_user_defined_properties(sdf_file)
            self.nodes[name] = {"smiles": smiles, "formal_charge": formal_charge, **extra_data}

        # Rename experimental key to standardized dg_value on nodes
        if self.exp_key:
            for name, node in self.nodes.items():
                if self.exp_key not in node:
                    raise KeyError(
                        f"exp_key '{self.exp_key}' not found in node '{name}'. "
                        f"Available keys: {', '.join(k for k in node if k not in ('smiles', 'formal_charge'))}"
                    )
                node["dg_value"] = node.pop(self.exp_key)

        return components

    def _network_to_dict(self, network) -> dict:
        """Convert a konnektor LigandNetwork to the output dict format."""
        # Build name lookup from components
        comp_to_name = {}
        for comp in network.nodes:
            name = comp.name
            if not name or name not in self.nodes:
                # Fallback: match by SMILES
                rdmol = comp.to_rdkit()
                smiles = Chem.MolToSmiles(rdmol)
                for node_name, node_data in self.nodes.items():
                    if node_data["smiles"] == smiles:
                        name = node_name
                        break
            comp_to_name[comp] = name

        # Build name-to-index mapping for source/target
        node_names = sorted(self.nodes.keys())
        name_to_idx = {name: idx for idx, name in enumerate(node_names)}

        edges = []
        same_charges = []
        for mapping in network.edges:
            from_name = comp_to_name[mapping.componentA]
            to_name = comp_to_name[mapping.componentB]
            weight = mapping.annotations.get("score", 0.0)

            same_charge = self.nodes[from_name]["formal_charge"] == self.nodes[to_name]["formal_charge"]
            same_charges.append(same_charge)

            edge = {
                "weight": weight,
                "source": name_to_idx[from_name],
                "target": name_to_idx[to_name],
                "from": from_name,
                "to": to_name,
                "same_charge": same_charge,
            }

            if self.exp_key:
                edge["ddg_value"] = self.nodes[to_name]["dg_value"] - self.nodes[from_name]["dg_value"]

            edges.append(edge)

        if all(same_charges):
            logger.info("All ligands have the same formal charge.")
        else:
            logger.warning("Not all ligands have the same formal charge!!")

        return {"nodes": self.nodes, "edges": edges}

    def _generate_network(self, generator, components) -> dict:
        """Generate a single network from components and return edge list."""
        gen_kwargs = {}
        if self.central_ligand is not None:
            if self.network_type != "star":
                logger.warning(
                    f"--central_ligand is only used with star networks, ignoring for '{self.network_type}'."
                )
            else:
                central_comp = None
                for comp in components:
                    if comp.name == self.central_ligand:
                        central_comp = comp
                        break
                if central_comp is None:
                    available = [c.name for c in components]
                    raise ValueError(
                        f"Central ligand '{self.central_ligand}' not found. Available: {available}"
                    )
                gen_kwargs["central_component"] = central_comp
                logger.info(f"Using '{self.central_ligand}' as central node for star network.")

        return generator.generate_ligand_network(components, **gen_kwargs)

    def run(self) -> dict:
        """Generate the perturbation network and write mapping.json.

        Returns:
            Dictionary with 'nodes' and 'edges' keys.
        """
        mapper = KartografAtomMapper(atom_map_hydrogens=False, map_exact_ring_matches_only=True)
        scorer = self._build_scorer()
        generator = self._build_generator(mapper, scorer)

        components = self._load_components()
        logger.info(f"Generating {self.network_type} network for {len(components)} ligands...")

        if self.separate_charges:
            from konnektor.network_tools.clustering.charge_clustering import (
                ChargeClusterer,
            )

            clusters = ChargeClusterer().cluster_compounds(components)
            logger.info(
                f"Separated ligands into {len(clusters)} charge groups: "
                f"{', '.join(f'charge {c}: {len(v)} ligands' for c, v in sorted(clusters.items()))}"
            )

            all_edges = []
            for charge, cluster_components in sorted(clusters.items()):
                if len(cluster_components) < 2:
                    logger.warning(
                        f"Charge group {charge} has only {len(cluster_components)} ligand(s), skipping network generation."
                    )
                    continue
                network = self._generate_network(generator, cluster_components)
                result = self._network_to_dict(network)
                all_edges.extend(result["edges"])

            result = {"nodes": self.nodes, "edges": all_edges}
        else:
            network = self._generate_network(generator, components)
            result = self._network_to_dict(network)

        # Write JSON output
        outpath = Path(self.out) / "mapping.json"
        with outpath.open("w") as f:
            json.dump(result, f, indent=4)
        logger.info(f"Network written to {outpath} ({len(result['edges'])} edges)")

        # Validate spatial alignment of ligand pairs in the network
        sdf_path = self._resolve_sdf_path()
        if sdf_path is not None:
            fixed = validate_edge_alignment(sdf_path, result, self_solve=self.self_solve)
            if fixed:
                logger.info(
                    "Re-generating network from corrected ligand coordinates..."
                )
                self.nodes = {}
                components = self._load_components()
                network = self._generate_network(generator, components)
                result = self._network_to_dict(network)
                with outpath.open("w") as f:
                    json.dump(result, f, indent=4)
                logger.info(
                    f"Network re-written to {outpath} ({len(result['edges'])} edges)"
                )
                # Re-validate — should be clean now
                validate_edge_alignment(sdf_path, result)

        return result

    def _resolve_sdf_path(self) -> Optional[Path]:
        """Find the SDF file used as input."""
        inp = Path(self.inp)
        if inp.is_file() and inp.suffix == ".sdf":
            return inp
        if inp.is_dir():
            sdfs = list(inp.glob("*.sdf"))
            if len(sdfs) == 1:
                return sdfs[0]
        return None


def _mcs_rmsd(mol_a, mol_b, timeout=5):
    """Compute RMSD of element-strict MCS-matched heavy atoms between two molecules.

    Uses element-strict atom comparison to prevent matching atoms of
    different elements (e.g., C to F), which would give misleadingly
    low RMSD values even when the molecules are poorly aligned.

    Returns (mcs_fraction, rmsd, n_mcs_atoms) or (0, -1, 0) on failure.
    """
    ha = Chem.RemoveHs(mol_a)
    hb = Chem.RemoveHs(mol_b)

    mcs = rdFMCS.FindMCS(
        [ha, hb], timeout=timeout,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
    )
    if mcs.numAtoms == 0:
        return 0.0, -1.0, 0

    mcs_frac = mcs.numAtoms / min(ha.GetNumAtoms(), hb.GetNumAtoms())
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    match_a = ha.GetSubstructMatch(mcs_mol)
    match_b = hb.GetSubstructMatch(mcs_mol)

    if not match_a or not match_b:
        return mcs_frac, -1.0, mcs.numAtoms

    conf_a = ha.GetConformer()
    conf_b = hb.GetConformer()
    sq_dists = []
    for ia, ib in zip(match_a, match_b):
        pa = conf_a.GetAtomPosition(ia)
        pb = conf_b.GetAtomPosition(ib)
        sq_dists.append((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2)

    rmsd = float(np.sqrt(np.mean(sq_dists)))
    return mcs_frac, rmsd, mcs.numAtoms


def _realign_to_neighbors(outlier_mol, neighbor_mols):
    """Regenerate outlier from SMILES, constrained to the neighbor centroid.

    Rebuilds from SMILES using ConstrainedEmbed with element-strict MCS
    to avoid cross-element coordinate pinning. The MCS core is constructed
    with centroid coordinates averaged across all neighbors, producing a
    compromise pose that minimizes drift across the neighborhood.

    Uses the same element-strict MCS strategy as BiasedConformerGenerator:
    first tries element-strict MCS (prevents cross-element misplacement),
    falls back to topology-only MCS if element-strict fails.

    Returns a new Mol with updated coordinates, or None on failure.
    """
    ho = Chem.RemoveHs(outlier_mol)

    # Find the best anchor using element-strict MCS
    best_anchor = None
    best_mcs = None
    best_size = 0
    for nb_mol in neighbor_mols:
        hn = Chem.RemoveHs(nb_mol)
        mcs = rdFMCS.FindMCS(
            [ho, hn], timeout=5,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True,
        )
        if mcs.numAtoms > best_size:
            best_size = mcs.numAtoms
            best_mcs = mcs
            best_anchor = nb_mol

    if best_mcs is None or best_size < 3:
        return None

    core = Chem.MolFromSmarts(best_mcs.smartsString)
    match_o = ho.GetSubstructMatch(core)
    if not match_o:
        return None

    # Collect centroid positions from all neighbors using element-strict MCS
    core_positions = {}
    for nb_mol in neighbor_mols:
        hn = Chem.RemoveHs(nb_mol)
        mcs_nb = rdFMCS.FindMCS(
            [ho, hn], timeout=5,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True,
        )
        if mcs_nb.numAtoms < 3:
            continue
        core_nb = Chem.MolFromSmarts(mcs_nb.smartsString)
        match_o_nb = ho.GetSubstructMatch(core_nb)
        match_n_nb = hn.GetSubstructMatch(core_nb)
        if not match_o_nb or not match_n_nb:
            continue
        conf_n = hn.GetConformer()
        o_to_core = {oa: i for i, oa in enumerate(match_o)}
        for oa, na in zip(match_o_nb, match_n_nb):
            if oa in o_to_core:
                p = conf_n.GetAtomPosition(na)
                core_positions.setdefault(o_to_core[oa], []).append(
                    np.array([p.x, p.y, p.z])
                )

    if not core_positions:
        return None

    # Build core with centroid coordinates
    core_mol = Chem.RWMol(core)
    core_conf = Chem.Conformer(core_mol.GetNumAtoms())
    for ci, positions in core_positions.items():
        centroid = np.mean(positions, axis=0)
        core_conf.SetAtomPosition(ci, Chem.rdGeometry.Point3D(*centroid.tolist()))
    # Fill unmatched core atoms from best anchor
    ha = Chem.RemoveHs(best_anchor)
    match_a = ha.GetSubstructMatch(core)
    if match_a:
        conf_a = ha.GetConformer()
        for ci in range(core_mol.GetNumAtoms()):
            if ci not in core_positions:
                p = conf_a.GetAtomPosition(match_a[ci])
                core_conf.SetAtomPosition(ci, p)
    core_mol.AddConformer(core_conf)

    # Strategy 1: ConstrainedEmbed with element-strict core
    smiles = Chem.MolToSmiles(Chem.RemoveHs(outlier_mol))
    new_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    try:
        AllChem.ConstrainedEmbed(new_mol, core_mol.GetMol(), randomseed=42)
    except (ValueError, RuntimeError):
        # Strategy 2: coordMap fallback (looser constraints)
        mol_noH = Chem.RemoveHs(new_mol)
        mol_match = mol_noH.GetSubstructMatch(core)
        if not mol_match:
            return None

        heavy_to_full = {}
        heavy_idx = 0
        for atom in new_mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                heavy_to_full[heavy_idx] = atom.GetIdx()
                heavy_idx += 1

        coord_map = {}
        for ci in core_positions:
            mol_noH_idx = mol_match[ci]
            if mol_noH_idx in heavy_to_full:
                centroid = np.mean(core_positions[ci], axis=0)
                coord_map[heavy_to_full[mol_noH_idx]] = Chem.rdGeometry.Point3D(*centroid.tolist())

        if not coord_map:
            return None

        result = AllChem.EmbedMolecule(
            new_mol, coordMap=coord_map, randomSeed=42,
            useRandomCoords=True, enforceChirality=True,
        )
        if result == -1:
            return None
        AllChem.MMFFOptimizeMolecule(new_mol)

    # Verify improvement
    old_rmsds = [_mcs_rmsd(outlier_mol, nb)[1] for nb in neighbor_mols]
    new_rmsds = [_mcs_rmsd(new_mol, nb)[1] for nb in neighbor_mols]
    old_mean = np.mean([r for r in old_rmsds if r >= 0])
    new_mean = np.mean([r for r in new_rmsds if r >= 0])
    if new_mean >= old_mean:
        return None

    # Transfer name and properties
    new_mol.SetProp("_Name", outlier_mol.GetProp("_Name"))
    for prop, val in outlier_mol.GetPropsAsDict().items():
        if prop == "_Name":
            continue
        if isinstance(val, float):
            new_mol.SetDoubleProp(prop, val)
        elif isinstance(val, int):
            new_mol.SetIntProp(prop, val)
        elif isinstance(val, str):
            new_mol.SetProp(prop, val)

    return new_mol


def validate_edge_alignment(
    sdf_path, network_dict, warn_rmsd=0.7, flag_rmsd=1.0, flag_mcs=0.7,
    self_solve=False,
):
    """Validate spatial alignment quality using neighborhood-based analysis.

    Instead of checking edges in isolation, computes each ligand's median
    MCS RMSD to ALL its network neighbors. Ligands with consistently high
    RMSD across their neighborhood are flagged — these are the ones whose
    3D pose is inconsistent with the local network topology.

    With self_solve=True, re-embeds flagged ligands constrained to the
    centroid of their neighbors' MCS positions, producing a compromise
    pose that minimizes drift across the entire neighborhood. Returns True
    if any ligands were fixed (caller should re-generate the network).

    Args:
        sdf_path: Path to the ligands SDF file.
        network_dict: Dictionary with 'edges' list (from mapping.json).
        warn_rmsd: Median neighbor RMSD for warnings (default 0.7 A).
        flag_rmsd: Median neighbor RMSD for flags (default 1.0 A).
        flag_mcs: Minimum MCS fraction for flagging (default 0.7).
        self_solve: If True, attempt to fix flagged ligands automatically.

    Returns:
        True if ligand coordinates were modified (SDF rewritten), False otherwise.
    """
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = {}
    for mol in supplier:
        if mol is None:
            continue
        mols[mol.GetProp("_Name")] = mol

    # Build adjacency list
    neighbors = {}
    for edge in network_dict["edges"]:
        a, b = edge["from"], edge["to"]
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)

    # Compute per-ligand neighborhood alignment profile
    logger.info("Validating ligand alignment (neighborhood analysis)...")
    ligand_profiles = {}
    for name in neighbors:
        if name not in mols:
            continue
        neighbor_rmsds = []
        for nb in neighbors[name]:
            if nb not in mols:
                continue
            mcs_frac, rmsd, _ = _mcs_rmsd(mols[name], mols[nb])
            if rmsd >= 0 and mcs_frac >= flag_mcs:
                neighbor_rmsds.append((nb, rmsd))
        if neighbor_rmsds:
            rmsds_only = [r for _, r in neighbor_rmsds]
            ligand_profiles[name] = {
                "median_rmsd": float(np.median(rmsds_only)),
                "max_rmsd": float(max(rmsds_only)),
                "n_neighbors": len(neighbor_rmsds),
                "details": neighbor_rmsds,
            }

    # Flag ligands with high median neighbor RMSD
    flagged = []
    warned = []
    for name, profile in sorted(ligand_profiles.items(), key=lambda x: -x[1]["median_rmsd"]):
        med = profile["median_rmsd"]
        mx = profile["max_rmsd"]
        n = profile["n_neighbors"]
        if med > flag_rmsd:
            flagged.append(name)
            bad_nbs = ", ".join(
                f"{nb}({r:.1f}A)" for nb, r in profile["details"] if r > warn_rmsd
            )
            logger.warning(
                f"  MISALIGNED: {name} "
                f"(median={med:.1f}A, max={mx:.1f}A, {n} neighbors: {bad_nbs})"
            )
        elif med > warn_rmsd:
            warned.append(name)
            logger.warning(
                f"  HIGH_DRIFT: {name} "
                f"(median={med:.1f}A, max={mx:.1f}A, {n} neighbors)"
            )

    if flagged:
        logger.warning(
            f"*** {len(flagged)} ligand(s) with poor neighborhood alignment! ***\n"
            f"    These ligands are spatially inconsistent with their network neighbors.\n"
            f"    FEP convergence will likely be poor for edges involving these ligands."
        )
        if self_solve:
            fixed_names = set()
            for name in flagged:
                nb_mols = [mols[nb] for nb in neighbors[name] if nb not in flagged and nb in mols]
                if not nb_mols:
                    # All neighbors are also flagged; use all of them
                    nb_mols = [mols[nb] for nb in neighbors[name] if nb in mols]

                old_med = ligand_profiles[name]["median_rmsd"]
                new_mol = _realign_to_neighbors(mols[name], nb_mols)
                if new_mol is not None:
                    # Recompute median RMSD with new coordinates
                    new_rmsds = [
                        _mcs_rmsd(new_mol, mols[nb])[1]
                        for nb in neighbors[name] if nb in mols
                    ]
                    new_med = float(np.median([r for r in new_rmsds if r >= 0]))
                    mols[name] = new_mol
                    fixed_names.add(name)
                    logger.warning(
                        f"    -> FIXED: {name} re-aligned to neighborhood centroid "
                        f"(median RMSD {old_med:.1f}A -> {new_med:.1f}A)"
                    )
                else:
                    logger.warning(
                        f"    -> {name} re-alignment FAILED. "
                        f"Re-generate with BiasedConformerGenerator using a "
                        f"different reference ligand."
                    )

            # Rewrite SDF
            if fixed_names:
                import tempfile

                sdf_p = Path(sdf_path)
                all_mols = list(mols.values())
                with tempfile.NamedTemporaryFile(
                    suffix=".sdf", dir=sdf_p.parent, delete=False, mode="w",
                ) as tmp:
                    tmp_path = Path(tmp.name)
                writer = Chem.SDWriter(str(tmp_path))
                for mol in all_mols:
                    writer.write(mol)
                writer.close()
                tmp_path.replace(sdf_p)
                logger.info(
                    f"Rewrote {sdf_p.name} with {len(fixed_names)} re-aligned "
                    f"ligand(s): {', '.join(sorted(fixed_names))}"
                )
                return True
        else:
            for name in flagged:
                nb_names = [nb for nb in neighbors[name] if nb not in flagged]
                if nb_names:
                    logger.warning(
                        f"    -> {name}: re-run with --self-solve to re-align "
                        f"to neighborhood centroid of {', '.join(nb_names)}."
                    )
                else:
                    logger.warning(
                        f"    -> {name}: all neighbors also flagged. "
                        f"Re-generate with BiasedConformerGenerator using a "
                        f"different reference ligand."
                    )
    elif warned:
        logger.info(
            f"{len(warned)} ligand(s) with moderate spatial drift (median RMSD > {warn_rmsd}A). "
            f"May be acceptable if due to rotatable substituents."
        )
    else:
        logger.info("All ligands have good neighborhood alignment (median RMSD < 1.0A).")
    return False


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate perturbation networks using konnektor.")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="SDF file or directory of SDF files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output directory. Defaults to input stem or directory.",
    )
    parser.add_argument(
        "-n",
        "--network",
        type=str,
        default="mst",
        choices=list(NETWORK_GENERATORS.keys()),
        help="Network topology: mst, rmst, star, nedges, or cyclic. Default: mst.",
    )
    parser.add_argument(
        "-s",
        "--scorer",
        type=str,
        default="combined",
        choices=["restraint", "kartograf", "combined"],
        help="Scorer type. 'combined' uses lomap heuristics * restraint overlap. Default: combined.",
    )
    parser.add_argument(
        "-rest",
        "--restraint_method",
        type=str,
        default="heavyatom_p",
        help="Restraint method (when scorer=restraint). Default: heavyatom_p.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Parallel processes for scoring. Default: 1.",
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log",
        default="info",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Set the log level for the logger. Defaults to `info`.",
    )
    parser.add_argument(
        "-cl",
        "--central_ligand",
        type=str,
        default=None,
        help="Central ligand name for star networks. Must match a molecule name in the input SDF.",
    )
    parser.add_argument(
        "--n_redundancy",
        type=int,
        default=2,
        help="Number of redundant MST iterations (only for rmst network). Default: 2.",
    )
    parser.add_argument(
        "--connectivity",
        type=int,
        default=3,
        help="Minimum edges per node (only for nedges network). Default: 3.",
    )
    parser.add_argument(
        "--separate_charges",
        action="store_true",
        default=False,
        help="Generate separate networks per charge group, avoiding cross-charge perturbations.",
    )
    parser.add_argument(
        "--charge_changes_score",
        type=float,
        default=0.0,
        help="LOMAP penalty for charge-changing edges (only with scorer=combined). "
        "0.0 blocks them, 0.5 halves their score, 1.0 no penalty. Default: 0.0.",
    )
    parser.add_argument(
        "--exp-key",
        "-exp",
        type=str,
        default=None,
        help="SDF property containing experimental dG. Stored as dg_value on nodes, ddg_value on edges.",
    )
    parser.add_argument(
        "--self-solve",
        action="store_true",
        default=False,
        help="Automatically fix misaligned ligand pairs by re-embedding the outlier "
        "constrained to its partner's MCS coordinates. Rewrites ligands.sdf in place.",
    )
    return parser.parse_args()


def main(args):
    setup_logger(level=args.log)
    kw = KonnektorWrap(
        inp=args.input,
        out=args.output,
        network=args.network,
        scorer=args.scorer,
        restraint_method=args.restraint_method,
        processes=args.processes,
        log_level=args.log,
        central_ligand=args.central_ligand,
        n_redundancy=args.n_redundancy,
        connectivity=args.connectivity,
        separate_charges=args.separate_charges,
        charge_changes_score=args.charge_changes_score,
        exp_key=args.exp_key,
        self_solve=args.self_solve,
    )
    kw.run()


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
