"""Command line interface for setting up amino-acid mutation FEP (QresFEP)."""

import argparse
import datetime
import json
import random
from pathlib import Path
from typing import Optional

from QligFEP import __version__

from ..logger import logger, setup_logger
from ..qresfep import QresFEP
from ..settings.settings import CLUSTER_DICT


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="qresfep",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Generate the input files for an amino-acid mutation FEP, using the dual-topology "
            "protocol. Run `qprep_prot` first: this reads its prep.json, prepared protein PDB "
            "and water.pdb from the working directory.\n\n"
            "The mutant residue must also be present as <MUT><POSITION>.pdb (e.g. ALA39.pdb), "
            "holding that residue alone, positioned on the wild-type backbone.\n\n"
            "A complete calculation needs both legs:\n"
            "  qresfep -m LEU39ALA -mc A -S protein    -FF OPLSAAM -c SNELLIUS\n"
            "  qresfep -m LEU39ALA -mc A -S tripeptide -FF OPLSAAM -c SNELLIUS"
        ),
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-m",
        "--mutation",
        dest="mutation",
        required=True,
        help=(
            "The mutation, as <wild-type><position><mutant> (e.g. LEU39ALA or L39A). "
            "The position is the residue number in the PDB given to `qprep_prot`."
        ),
    )
    required.add_argument(
        "-mc",
        "--mutchain",
        dest="chain",
        required=True,
        help="Chain of the mutated residue, as in the input PDB.",
    )
    required.add_argument(
        "-S",
        "--system",
        dest="system",
        required=True,
        choices=["protein", "tripeptide"],
        help=(
            "Which leg of the thermodynamic cycle to set up: `protein` for the mutation in "
            "the folded protein, `tripeptide` for the capped reference peptide. "
            "ddG_fold is the difference between the two."
        ),
    )
    required.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        required=True,
        choices=list(CLUSTER_DICT.keys()),
        help="Cluster the run scripts are written for; profiles live in settings.py.",
    )

    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "-FF",
        "--forcefield",
        dest="force_field",
        default="OPLSAAM",
        help=(
            "Force field, which must match the one used for `qprep_prot`. Either a name "
            "(OPLSAAM, OPLS2015, OPLS2005, AMBER14sb, CHARMM36) or a path to .lib/.prm files "
            "without the extension. Defaults to OPLSAAM."
        ),
    )
    optional.add_argument(
        "-t",
        "--tripeptide",
        dest="tripeptide_flanks",
        default="A",
        choices=["A", "G", "X", "Z"],
        help=(
            "Residues flanking the mutated one in the reference peptide: A (Ala), G (Gly), "
            "X (none), Z (the native neighbours). Z gives less stable results for "
            "charge-changing mutations. Defaults to A."
        ),
    )
    optional.add_argument(
        "-w",
        "--windows",
        dest="windows",
        type=int,
        default=25,
        help=(
            "Lambda windows *per stage*. A dual-topology mutation runs two stages, so the "
            "total is twice this number. Defaults to 25."
        ),
    )
    optional.add_argument(
        "-s",
        "--sampling",
        dest="sampling",
        default="exponential",
        choices=["linear", "sigmoidal", "exponential", "reverse_exponential"],
        help=(
            "Lambda spacing. `exponential` gives stage 1 exponential and stage 2 "
            "reverse-exponential spacing, concentrating windows where each stage is "
            "steep. Defaults to exponential."
        ),
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
        help="Simulation timestep; SHAKE on hydrogens follows it. Defaults to 2fs.",
    )
    optional.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        default="298",
        help="Temperature in K; several as 'T1,T2,...'. Defaults to 298.",
    )
    optional.add_argument(
        "-r",
        "--replicates",
        dest="replicates",
        type=int,
        default=10,
        help="Independent repeats. Defaults to 10.",
    )
    optional.add_argument(
        "-sh",
        "--shell_rest",
        dest="shell_restraint",
        type=float,
        default=0.0,
        help=(
            "Width of the restrained outer shell, subtracted from the sphere radius. "
            "Recommended for membrane proteins. Defaults to 0.0."
        ),
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
        "-cof",
        "--cofactors",
        dest="cofactors",
        nargs="*",
        help="Cofactor basenames; each needs a .pdb, .lib and .prm in the working directory.",
    )
    optional.add_argument(
        "-clean",
        "--files-to-clean",
        dest="to_clean",
        nargs="+",
        default=None,
        help=(
            "Suffixes to delete once a run finishes, e.g. `-clean dcd` to drop trajectories. "
            "Defaults to keeping everything."
        ),
    )
    optional.add_argument(
        "--no-trajectories",
        dest="write_trajectories",
        action="store_false",
        default=True,
        help="Disable DCD output in every equilibration and production input.",
    )
    optional.add_argument(
        "--coupled-thermostat",
        dest="separate_scaling",
        action="store_false",
        default=True,
        help=(
            "Use the legacy shared solute/solvent heat-bath scaling. This reproduces "
            "the published protocol; corrected runs scale them separately."
        ),
    )
    seed_options = optional.add_mutually_exclusive_group()
    seed_options.add_argument(
        "-rs",
        "--random_state",
        dest="random_state",
        type=int,
        default=None,
        help="Seed for the per-replicate random seeds, for reproducible runs. Defaults to None.",
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
        "-log",
        "--log-level",
        dest="log",
        default="info",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Log level. Defaults to info.",
    )
    return parser.parse_args()


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> Path:
    """Set up one mutation leg.

    Args:
        args: Parsed command line arguments.
        kwargs: Overrides passed straight to :class:`~QligFEP.qresfep.QresFEP`.

    Returns:
        The FEP directory that was created.
    """
    parameters = {}
    if args is not None:
        setup_logger(level=args.log.upper())
        rng = random.Random(args.random_state)
        seeds = (
            list(args.seeds)
            if args.seeds is not None
            else [rng.randint(1, 9999) for _ in range(args.replicates)]
        )
        parameters = {
            "mutation": args.mutation,
            "chain": args.chain,
            "system": args.system,
            "force_field": args.force_field,
            "cluster": args.cluster,
            "windows": args.windows,
            "sampling": args.sampling,
            "start": args.start,
            "timestep": args.timestep,
            "temperature": args.temperature,
            "replicates": args.replicates,
            "tripeptide_flanks": args.tripeptide_flanks,
            "shell_restraint": args.shell_restraint,
            "eq5_steps": args.eq5_steps,
            "cofactors": args.cofactors,
            "to_clean": args.to_clean,
            "write_trajectories": args.write_trajectories,
            "separate_scaling": args.separate_scaling,
            "seeds": seeds,
        }
    parameters.update(kwargs)

    run = QresFEP(**parameters)
    directory = run.run()

    # Record how the directory was produced, mirroring qligfep's fep_config.json.
    # Read off the object rather than the arguments, so the defaults that were
    # actually used are recorded too -- `qresfep_analyze` reads `start` from here to
    # know which direction the stages ran in.
    (directory / "inputfiles" / "resfep_config.json").write_text(
        json.dumps(
            {
                "mutation": run.name,
                "wild_type": run.wild_type,
                "position": run.position,
                "mutant": run.mutant,
                "chain": run.chain,
                "system": run.system,
                "force_field": run.force_field,
                "cluster": run.cluster,
                "windows": run.windows,
                "sampling": run.sampling,
                "start": run.start,
                "timestep": run.timestep,
                "temperature": run.temperature,
                "replicates": run.replicates,
                "tripeptide_flanks": run.tripeptide_flanks,
                "shell_restraint": run.shell_restraint,
                "eq5_steps": run.eq5_steps,
                "cofactors": run.cofactors,
                "to_clean": run.to_clean,
                "write_trajectories": run.write_trajectories,
                "separate_scaling": run.separate_scaling,
                "seeds": run.seeds,
                "q_position": run.q_position,
                "hybrid_residue": run.hybrid_name,
                "sphere_radius": run.radius,
                "sphere_center": run.center,
                "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "QligFEP_version": __version__,
            },
            indent=4,
        )
    )
    logger.info(f"{run.name} ({run.system}) set up in {directory}")
    return directory


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
