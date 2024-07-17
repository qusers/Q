"""Module holding a parse_arguments base function for both qligfep and setupFEP CLI's"""

import argparse

from ..settings.settings import CLUSTER_DICT


def parse_arguments(program: str) -> argparse.Namespace:
    if program == "QligFEP":
        parser = argparse.ArgumentParser(
            prog="QligFEP",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="       == Generate FEP files for dual topology ligand FEP == ",
        )
        parser.add_argument("-l1", "--lig_1", dest="lig1", required=True, help="name of ligand 1", type=str)
        parser.add_argument("-l2", "--lig_2", dest="lig2", required=True, help="name of ligand 2", type=str)
    elif program == "setupFEP":
        parser = argparse.ArgumentParser(
            prog="setupFEP",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(
                "Generate all FEP files for the directory you're working on, according to the "
                "edges input in the json_map file. This includes creating directories for both "
                "water and protein system. Submitting the FEP calculations to the cluster is up to the user. "
                "A minimal example of usage: setupFEP -FF OPLS2015 -c KEBNE -S sigmoidal -r 25 -l 0.5 -w 100"
            ),
        )
    parser.add_argument(
        "-FF",
        "--forcefield",
        dest="FF",
        required=True,
        choices=["OPLS2005", "OPLS2015", "AMBER14sb", "CHARMM36", "CHARMM22", "CHARMM_TEST"],
        help="Forcefield to be used.",
    )
    if program == "QligFEP":
        parser.add_argument(
            "-s",
            "--system",
            dest="system",
            required=True,
            choices=["water", "protein", "vacuum"],
            help="what type of system we are setting up",
        )
    parser.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        required=True,
        help="cluster you want to submit to, cluster specific parameters added to settings.",
        choices=list(CLUSTER_DICT.keys()),
    )
    parser.add_argument(
        "-r",
        "--sphereradius",
        dest="sphereradius",
        required=False,
        default="25",
        help="Size of the simulation sphere. Defaults to 25.",
    )
    parser.add_argument(
        "-b",
        "--cysbond",
        dest="cysbond",
        default="auto",
        help=(
            "Add cystein bonds. Input should be formatted with the atom numbers"
            "(participating in the Cys bond) connected by `_` and with different bonds "
            "separated by `,` as in: `atom1_atom2,atom3_atom4`. Defaults to `auto`, where "
            "cystein bonds will be automatically detected within distance of 1.8 to 2.2 A."
        ),
        type=str,
    )
    parser.add_argument(
        "-l",
        "--start",
        dest="start",
        default="0.5",
        choices=["1", "0.5"],
        help="Starting FEP in the middle or endpoint. Defaults to 0.5.",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        default="298",
        help="Temperature(s), mutliple tempereratures given as 'T1,T2,...,TN'. Defaults to 298K",
    )
    parser.add_argument(
        "-R",
        "--replicates",
        dest="replicates",
        default="10",
        help="How many repeats should be run. Defaults to 10.",
    )
    parser.add_argument(
        "-S",
        "--sampling",
        dest="sampling",
        default="sigmoidal",
        choices=["linear", "sigmoidal", "exponential", "reverse_exponential"],
        help="Lambda spacing type to be used. Defaults to `sigmoidal`.",
    )
    parser.add_argument(
        "-w",
        "--windows",
        dest="windows",
        default="100",
        help="Total number of windows that will be run. Defaults to 100.",
        type=str,
    )
    if program == "QligFEP":
        parser.add_argument(
            "-sc",
            "--softcore",
            dest="softcore",
            default=False,
            action="store_true",
            help="Turn on if you want to use softcore",
        )
    if program == "setupFEP":  # this is only used for setupFEP
        parser.add_argument(
            "-j",
            "--json_map",
            dest="json_map",
            help=(
                "Path for the '.json' QmapFEP file. If not given, the script will "
                "look for a single '.json' file in the current directory and raise "
                "an error if there are more than one."
            ),
            default=None,
        )
    parser.add_argument(
        "-ts",
        "--timestep",
        dest="timestep",
        choices=["1fs", "2fs"],
        default="2fs",
        help="Simulation timestep, default 2fs",
    )
    parser.add_argument(
        "-clean",
        "--files-to-clean",
        dest="to_clean",
        nargs="+",
        default="",
        help=(
            "Files to clean after the simulation. The arguments are given as a list of strings "
            "and the cleaning is done by adding the command `rm -rf *{arg1} *{arg2}` to the job submission. "
            "Usage example: `-clean dcd` will remove all dcd files after the simulation. If left as None, won't clean any files."
        ),
    )
    parser.add_argument(
        "-rest",
        "--restraint_method",
        dest="restraint_method",
        type=str,
        default="chemoverlap_mp",
        choices=[
            "overlap",
            "chemoverlap_mp",
            "chemoverlap_ls",
            "chemoverlap_strict",
        ],
        help=(
            """How to set the restraints to the ligand topologies involved in the perturbation.
              1) `overlap` - old implementation; restraints are applied to whichever atoms
                    within 1 A distance in the perturbation topologies;
              2) `chemoverlap_mp` (more permissive); restraints are applied to atoms which have
                  equivalent ring structures. The restraint is then applied to the respective
                  decorations in the ring, provided same atom types.
              3) `chemoverlap_ls` - (less strict); restraint are applied to atoms which have
                  equivalent ring structures + immediate surround. E.g.: a methyl in only one
                  of the rings would lead to the whole ring left unrestrained.
              4) `chemoverlap_strict` - (strict); restraint are applied to atoms which have
                  equivalent ring structures + immediate surround, but immediate surround
                  will be atom-type sensitive. E.g.: a methyl v.s. a chlorine would lead to
                  the unrestrained ring."""
        ),
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log",
        required=False,
        default="info",
        help="Set the log level for the logger. Defaults to `info`.",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
    )
    return parser.parse_args()
