"""Generate perturbation networks using konnektor with pluggable scorers and network topologies."""

import argparse
import json
from pathlib import Path
from typing import Optional

from kartograf import KartografAtomMapper, SmallMoleculeComponent
from rdkit import Chem

from ..chemIO import MoleculeIO
from ..logger import logger, setup_logger
from ..scoring import RestraintScorer
from .lomap_wrap_cli import LomapWrap

NETWORK_GENERATORS = {
    "mst": "MinimalSpanningTreeNetworkGenerator",
    "star": "StarNetworkGenerator",
    "cyclic": "CyclicNetworkGenerator",
}


class KonnektorWrap:
    """Generate perturbation networks using konnektor.

    Supports MST, star, and cyclic network topologies with either
    RestraintSetter-based or kartograf volume ratio scoring.
    """

    def __init__(
        self,
        inp: str,
        out: Optional[str] = None,
        network: str = "mst",
        scorer: str = "restraint",
        restraint_method: str = "heavyatom_p",
        processes: int = 1,
        verbose: str = "info",
        central_ligand: Optional[str] = None,
    ):
        self.inp = inp
        self.network_type = network
        self.scorer_type = scorer
        self.restraint_method = restraint_method
        self.processes = processes
        self.verbose = verbose
        self.central_ligand = central_ligand
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
        else:
            return RestraintScorer(restraint_method=self.restraint_method)

    def _build_generator(self, mapper, scorer):
        """Create the network generator based on network type."""
        from konnektor.network_planners import (
            CyclicNetworkGenerator,
            MinimalSpanningTreeNetworkGenerator,
            StarNetworkGenerator,
        )

        generators = {
            "mst": MinimalSpanningTreeNetworkGenerator,
            "star": StarNetworkGenerator,
            "cyclic": CyclicNetworkGenerator,
        }
        if self.network_type not in generators:
            raise ValueError(f"Unknown network type '{self.network_type}'. Choose from: {list(generators.keys())}")
        return generators[self.network_type](
            mappers=[mapper],
            scorer=scorer,
            n_processes=self.processes,
        )

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

    def run(self) -> dict:
        """Generate the perturbation network and write mapping.json.

        Returns:
            Dictionary with 'nodes' and 'edges' keys.
        """
        mapper = KartografAtomMapper(atom_map_hydrogens=False)
        scorer = self._build_scorer()
        generator = self._build_generator(mapper, scorer)

        components = self._load_components()
        logger.info(f"Generating {self.network_type} network for {len(components)} ligands...")

        # Resolve central ligand for star networks
        gen_kwargs = {}
        if self.central_ligand is not None:
            if self.network_type != "star":
                logger.warning(f"--central_ligand is only used with star networks, ignoring for '{self.network_type}'.")
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

        network = generator.generate_ligand_network(components, **gen_kwargs)

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
        "-i", "--input", type=str, required=True,
        help="SDF file or directory of SDF files.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output directory. Defaults to input stem or directory.",
    )
    parser.add_argument(
        "-n", "--network", type=str, default="mst", choices=list(NETWORK_GENERATORS.keys()),
        help="Network topology: mst, star, or cyclic. Default: mst.",
    )
    parser.add_argument(
        "-s", "--scorer", type=str, default="restraint", choices=["restraint", "kartograf"],
        help="Scorer type. Default: restraint.",
    )
    parser.add_argument(
        "-rest", "--restraint_method", type=str, default="heavyatom_p",
        help="Restraint method (when scorer=restraint). Default: heavyatom_p.",
    )
    parser.add_argument(
        "-p", "--processes", type=int, default=1,
        help="Parallel processes for scoring. Default: 1.",
    )
    parser.add_argument(
        "-v", "--verbose", type=str, default="info",
        help="Verbosity level. Default: info.",
    )
    parser.add_argument(
        "-cl", "--central_ligand", type=str, default=None,
        help="Central ligand name for star networks. Must match a molecule name in the input SDF.",
    )
    return parser.parse_args()


def main(args):
    setup_logger(level=args.verbose)
    kw = KonnektorWrap(
        inp=args.input,
        out=args.output,
        network=args.network,
        scorer=args.scorer,
        restraint_method=args.restraint_method,
        processes=args.processes,
        verbose=args.verbose,
        central_ligand=args.central_ligand,
    )
    kw.run()


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
