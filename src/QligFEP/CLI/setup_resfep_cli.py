"""Command line interface for setting up a series of amino-acid mutations."""

import argparse
import sys
from pathlib import Path
from typing import Optional

from ..logger import logger, setup_logger
from ..resfep_protocols import DEFAULT_PRODUCTION_STEPS, apply_manuscript_settings
from ..resfep_setup import Mutation, MutationSeries, SetupError, read_mutations
from ..settings.settings import CLUSTER_DICT


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="setup_resFEP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Set up an amino-acid mutation FEP series: one prepared sphere per mutation, both "
            "legs of the thermodynamic cycle, from a single protein structure.\n\n"
            "Unlike a ligand series, mutations do not share a binding site, so each one gets "
            "its own `qprep_prot` run centred on the residue being mutated. That also "
            "re-derives which other charges are neutralised as out-of-sphere.\n\n"
            "Pass -cog to keep one centre for the whole series instead -- for scanning "
            "residues around a bound ligand. Mutations that would then reach outside the "
            "sphere, or into the restrained shell where charges are neutralised, are "
            "refused.\n\n"
            "Produces:\n"
            "  work/<MUTATION>/       the prepared sphere for each mutation\n"
            "  protein/FEP_<MUTATION>/\n"
            "  tripeptide/FEP_<MUTATION>/\n\n"
            "Then analyse with: qresfep_analyze -p protein -t tripeptide -T 298"
        ),
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        "--input_pdb_file",
        dest="structure",
        required=True,
        help="Prepared protein PDB. Mutations are numbered against this file.",
    )
    mutations = required.add_mutually_exclusive_group(required=True)
    mutations.add_argument(
        "-m",
        "--mutations",
        dest="mutations",
        nargs="+",
        help="Mutations as <wild-type><position><mutant> (e.g. LEU39ALA L39G).",
    )
    mutations.add_argument(
        "-M",
        "--mutations-file",
        dest="mutations_file",
        help="File listing one mutation per line; blank lines and `#` comments are ignored.",
    )
    required.add_argument(
        "-mc",
        "--mutchain",
        dest="chain",
        required=True,
        help="Chain the mutated residues belong to, as in the input PDB.",
    )
    required.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        required=True,
        choices=list(CLUSTER_DICT.keys()),
        help="Cluster the run scripts are written for; profiles live in settings.py.",
    )

    sphere = parser.add_argument_group("sphere")
    sphere.add_argument(
        "-r",
        "--sphereradius",
        dest="radius",
        type=float,
        default=25.0,
        help="Sphere radius in Angstrom. Defaults to 25.",
    )
    sphere.add_argument(
        "-cog",
        "--center_of_geometry",
        dest="center",
        nargs=3,
        type=float,
        default=None,
        help=(
            "Keep one sphere centre for the whole series, as 'x y z' -- for mutating "
            "residues around a bound ligand. Defaults to centring each sphere on its own "
            "mutated residue (its CB, or HA3 for glycine), which is what the published "
            "protocol does."
        ),
    )
    sphere.add_argument(
        "-nbo",
        "--neutralize_boundary_offset",
        dest="neutralization_offset",
        type=float,
        default=3.0,
        help=(
            "Passed to `qprep_prot`: charged residues beyond (radius - offset) are prepared "
            "in a neutral form. Defaults to 3.0."
        ),
    )

    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "-FF",
        "--forcefield",
        dest="force_field",
        default="OPLSAAM",
        help=(
            "Force field for both preparation and setup. Either a name (OPLSAAM, OPLS2015, "
            "OPLS2005, AMBER14sb, CHARMM36) or a path without the extension. Defaults to "
            "OPLSAAM."
        ),
    )
    optional.add_argument(
        "-t",
        "--tripeptide",
        dest="tripeptide_flanks",
        default="A",
        choices=["A", "G", "X", "Z"],
        help="Residues flanking the mutated one in the reference peptide. Defaults to A.",
    )
    optional.add_argument(
        "-w",
        "--windows",
        dest="windows",
        type=int,
        default=25,
        help="Lambda windows per stage; a mutation runs two stages. Defaults to 25.",
    )
    optional.add_argument(
        "-s",
        "--sampling",
        dest="sampling",
        default="exponential",
        choices=["linear", "sigmoidal", "exponential", "reverse_exponential"],
        help="Lambda spacing. Defaults to exponential.",
    )
    optional.add_argument(
        "-l",
        "--start",
        dest="start",
        default="1",
        choices=["1", "0.5", "0"],
        help="Lambda endpoint the simulations start from. Defaults to 1.",
    )
    optional.add_argument(
        "-ts",
        "--timestep",
        dest="timestep",
        default="2fs",
        choices=["1fs", "2fs"],
        help="Simulation timestep. Defaults to 2fs.",
    )
    optional.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        default="298",
        help="Temperature in K; several as 'T1,T2,...'. Defaults to 298.",
    )
    optional.add_argument(
        "-R",
        "--replicates",
        dest="replicates",
        type=int,
        default=10,
        help="Independent repeats per mutation. Defaults to 10.",
    )
    optional.add_argument(
        "-eqs",
        "--eq5-steps",
        dest="eq5_steps",
        type=int,
        default=None,
        help=(
            "Length of the final equilibration (eq5) in steps, overriding the protocol's "
            "own 2.5 ns at the default 2 fs timestep. Use it to match a set that was "
            "equilibrated for longer."
        ),
    )
    optional.add_argument(
        "-ps",
        "--production-steps",
        dest="production_steps",
        type=int,
        default=DEFAULT_PRODUCTION_STEPS,
        help=f"Production length per lambda window in steps. Defaults to {DEFAULT_PRODUCTION_STEPS}.",
    )
    optional.add_argument(
        "-cof",
        "--cofactors",
        dest="cofactors",
        nargs="*",
        help=(
            "Cofactor basenames; each needs a .pdb, .lib and .prm beside the input PDB. They "
            "are added to every prepared sphere."
        ),
    )
    optional.add_argument(
        "-clean",
        "--files-to-clean",
        dest="to_clean",
        nargs="+",
        default=None,
        help="Suffixes to delete once a run finishes, e.g. `-clean dcd`.",
    )
    optional.add_argument(
        "--no-trajectories",
        dest="write_trajectories",
        action="store_false",
        default=True,
        help="Disable DCD output in all generated Q inputs.",
    )
    optional.add_argument(
        "--separate-scaling",
        dest="separate_scaling",
        choices=["on", "off"],
        default="on",
        help="Set Q's separate_scaling option. Defaults to on.",
    )
    optional.add_argument(
        "--manuscript-settings",
        action="store_true",
        help="Force the residue-FEP manuscript protocol settings.",
    )
    seed_options = optional.add_mutually_exclusive_group()
    seed_options.add_argument(
        "-rs",
        "--random_state",
        dest="random_state",
        type=int,
        default=None,
        help="Seed for the per-replicate random seeds. Defaults to None.",
    )
    seed_options.add_argument(
        "--seeds",
        dest="seeds",
        nargs="+",
        type=int,
        default=None,
        help=(
            "Explicit MD random seed for each replicate. The number of values must equal "
            "--replicates; mutually exclusive with --random-state."
        ),
    )
    optional.add_argument(
        "-mp",
        "--mutant-pdbs",
        dest="mutant_pdb_dir",
        default=None,
        help=(
            "Directory of ready-made mutant residue PDBs (<MUTANT><POSITION>.pdb). Defaults "
            "to building them with PyMOL's mutagenesis wizard."
        ),
    )
    optional.add_argument(
        "--legs",
        dest="legs",
        nargs="+",
        default=["protein", "tripeptide"],
        choices=["protein", "tripeptide"],
        help="Legs of the thermodynamic cycle to set up. Defaults to both.",
    )
    optional.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Where the series is written. Defaults to the current directory.",
    )
    optional.add_argument(
        "-n",
        "--validate-only",
        dest="validate_only",
        action="store_true",
        help="Check the structure and every mutation, then stop without writing anything.",
    )
    optional.add_argument(
        "-log",
        "--log-level",
        dest="log",
        default="info",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Log level. Defaults to info.",
    )
    args = parser.parse_args()
    if args.manuscript_settings:
        apply_manuscript_settings(args, include_preparation=True)
    return args


def _qresfep_options(args: argparse.Namespace) -> list[str]:
    """The arguments passed straight through to `qresfep`."""
    options = ["-log", args.log]
    if getattr(args, "manuscript_settings", False):
        options += ["--manuscript-settings"]
    else:
        options += [
            "-t", args.tripeptide_flanks,
            "-w", str(args.windows),
            "-s", args.sampling,
            "-l", args.start,
            "-ts", args.timestep,
            "-T", str(args.temperature),
            "-r", str(args.replicates),
            "-ps", str(args.production_steps),
            "--separate-scaling", args.separate_scaling,
        ]
        if args.eq5_steps is not None:
            options += ["-eqs", str(args.eq5_steps)]
        if args.random_state is not None:
            options += ["-rs", str(args.random_state)]
        if args.seeds is not None:
            options += ["--seeds", *(str(seed) for seed in args.seeds)]

    if args.to_clean:
        options += ["-clean", *args.to_clean]
    if not args.write_trajectories:
        options += ["--no-trajectories"]
    return options


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> list:
    """Set up every mutation of a series.

    Args:
        args: Parsed command line arguments.
        kwargs: Overrides passed straight to
            :class:`~QligFEP.resfep_setup.MutationSeries`.

    Returns:
        One :class:`~QligFEP.resfep_setup.MutationOutcome` per mutation.
    """
    parameters = {}
    validate_only = False
    if args is not None:
        setup_logger(level=args.log.upper())
        validate_only = args.validate_only
        mutations = (
            read_mutations(Path(args.mutations_file), args.chain)
            if args.mutations_file
            else [Mutation.from_string(entry, args.chain) for entry in args.mutations]
        )
        parameters = {
            "structure": Path(args.structure),
            "mutations": mutations,
            "force_field": args.force_field,
            "cluster": args.cluster,
            "radius": args.radius,
            "center": args.center,
            "neutralization_offset": args.neutralization_offset,
            "legs": tuple(args.legs),
            "cofactors": args.cofactors,
            "mutant_pdb_dir": args.mutant_pdb_dir,
            "workdir": Path(args.output_dir) if args.output_dir else None,
            "qresfep_options": _qresfep_options(args),
        }
    parameters.update(kwargs)

    series = MutationSeries(**parameters)
    where = "one sphere per mutated residue"
    if series.center is not None:
        where = f"one sphere for the series, centred on {series.center}"
    logger.info(
        f"{len(series.mutations)} mutation(s) of {series.structure.name}, "
        f"{series.radius:.0f} A radius, {where}"
    )

    if validate_only:
        series.validate()
        logger.info("Every mutation checks out; nothing was written (--validate-only).")
        return []

    outcomes = series.run()

    failed = [outcome for outcome in outcomes if not outcome.ok]
    logger.info(f"Set up {len(outcomes) - len(failed)}/{len(outcomes)} mutation(s)")
    if failed:
        logger.error(
            "These mutations were not set up:\n  "
            + "\n  ".join(f"{o.mutation.name}: {o.error.splitlines()[0]}" for o in failed)
        )
    return outcomes


def main_exe():
    args = parse_arguments()
    try:
        outcomes = main(args)
    except SetupError as error:
        logger.error(str(error))
        sys.exit(1)
    if any(not outcome.ok for outcome in outcomes):
        sys.exit(1)


if __name__ == "__main__":
    main_exe()
