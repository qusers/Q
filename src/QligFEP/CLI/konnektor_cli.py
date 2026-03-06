"""Generate perturbation networks using konnektor with pluggable scorers and network topologies."""

import argparse
import json
from pathlib import Path
from typing import Optional

from kartograf import KartografAtomMapper, SmallMoleculeComponent
from rdkit import Chem

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

        # Find which node properties are numerical (for delta computation)
        extra_numerical_keys = set()
        for node_data in self.nodes.values():
            for k, v in node_data.items():
                if k not in ("smiles", "formal_charge") and isinstance(v, (int, float)):
                    extra_numerical_keys.add(k)

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

            # Compute delta for numerical properties
            for key in extra_numerical_keys:
                try:
                    delta = self.nodes[from_name][key] - self.nodes[to_name][key]
                    edge[f"delta_{key}"] = delta
                except (TypeError, KeyError):
                    logger.info(f"Could not compute delta_{key} for {from_name} | {to_name}")

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

        return result


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
    )
    kw.run()


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
