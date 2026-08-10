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
        parser.add_argument("-l1", "--lig1", dest="lig1", required=True, help="name of ligand 1", type=str)
        parser.add_argument("-l2", "--lig2", dest="lig2", required=True, help="name of ligand 2", type=str)
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
        default="AMBER14sb",
        help=(
            "Protein forcefield to be used. Valid inputs: existing path to a forcefield file without the extensions"
            "(either .lib, .prm, or Path without the extensions will work) or one of the following: "
            "OPLS2005, OPLS2015, AMBER14sb, CHARMM36. Defaults to AMBER14sb."
        ),
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
            "-po",
            "--protein_only",
            action="store_true",
            help="Only generate FEP files for the protein system.",
        )
        parser.add_argument(
            "-wo",
            "--water_only",
            action="store_true",
            help="Only generate FEP files for the water system.",
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
        "--minimize",
        action="store_true",
        help="Run FIRE energy minimization before the first equilibration stage. Disabled by default.",
    )
    parser.add_argument(
        "--production",
        action="store_true",
        help=(
            "Minimize equilibrium FEP output: disable trajectory files and, unless --clean is "
            "specified, remove .en, .re, and dualtop.top files after qfep has produced "
            "energies.csv and qfep.out."
        ),
    )
    parser.add_argument(
        "-clean",
        "--clean",
        "--files-to-clean",
        dest="to_clean",
        nargs="+",
        default=None,
        help=(
            "File suffixes to remove after a successful simulation and qfep run. "
            "Usage example: `--clean en dualtop.top` removes *en and *dualtop.top. "
            "When explicitly supplied with --production, this list replaces the production "
            "defaults, allowing other files (for example .re files) to be retained."
        ),
    )
    parser.add_argument(
        "-rest",
        "--restraint_method",
        dest="restraint_method",
        type=str,
        default="heavyatom_p",
        help=(
            """How to set the restraints to the ligand topologies involved in the perturbation. Defaults to `heavyatom_p`.

            Atom compare method: `heavyatom`, `aromaticity`, `hibridization`, `element`. Setting the first part of the
                string as either of these, will determine how the substituents / ring atoms are treated to be
                defined as equivalent. Heavyatom will only check if both atoms being compared are heavy atoms. Aromaticity
                and hybridization operater similarly, but checking for those properties instead. Element will check
                for atom equivalence based on atomic numbers.
            Surround atom compare: `p` (permissive), `ls` (less strict), `strict`.
                Setting the second part of the string as either of these, will determine if or how the
                direct surrounding atoms to the ring strictures will be taken into account for ring equivalence.
                    1) Permissive: Only the ring atoms are compared.
                    2) Less strict: The ring atoms and their direct surroundings are compared, but element type
                        is ignored.
                    3) Strict: The ring atoms and their direct surroundings are element-wise compared.
           Kartograf atom max distance (optional): int or float to be used by `kartograf` as the maximum distance between
                atoms to be considered for mapping. This is by default set to 0.95 A, but can be changed by passing `_0.95`,
                for example, at the end of the `restraint_method` string."""
        ),
    )
    parser.add_argument(
        "-drf",
        "--distance_restraint_force",
        dest="dr_force",
        type=float,
        default=0.5,
        help=(
            "Force constant applied as distance restrains to ensure space overlap among ligands during FEP. For how/which atoms are mapped "
            "restrained together by this force, check the `restraint_method` parameter description. Defaults to 0.5 kcal/mol/A^2. "
            "restraints forces set through this parameter will be applied for equilibration 5 (eq5) and the md_xxxx_xxxx simulations. "
            "Forces set by this parameter are set to the [distance_restraints] section of `.inp` (input) files."
        ),
    )
    parser.add_argument(
        "-wath",
        "--water-thresh",
        dest="water_thresh",
        type=float,
        default=1.4,
        help=(
            "Threshold (in Angstrom) for removing water molecules that are too close to the ligands "
            "involved in the FEP. Defaults to 1.4 A."
        ),
    )
    parser.add_argument(
        "-wath-ligo",
        "--wath-ligand-only",
        dest="wath_ligand_only",
        action="store_true",
        help=(
            "If set, the water threshold set by the `-wath` parameter will only be applied to the ligand "
            "atoms, not the entire protein structure."
        ),
    )
    parser.add_argument(
        "-rs",
        "--random_state",
        dest="random_state",
        type=int,
        default=None,
        help=(
            "Reproducible random state for the random FEP seed generator. Defaults to None (random FEP seeds)."
        ),
    )
    parser.add_argument(
        "-neq",
        "--neq",
        dest="neq",
        action="store_true",
        help=(
            "Set up a non-equilibrium (NEQ) FEP instead of the windowed equilibrium approach. "
            "Instead of many fixed-lambda windows, NEQ runs the serial `qdyn` engine (fast-switching "
            "mode) to drive lambda "
            "from one endpoint to the other over a single simulation, accumulating the switching "
            "work. Free energies are obtained from BAR over the forward/reverse work distributions "
            "(see `qligfep_neq_analyze`). When set, the windowed parameters `--windows` and "
            "`--sampling` are not used."
        ),
    )
    parser.add_argument(
        "-neqr",
        "--neq-reps",
        dest="neq_reps",
        type=int,
        default=5,
        help=(
            "NEQ only: number of forward/reverse switching pairs run per replicate. Each replicate "
            "keeps a continuous endpoint equilibration and fires a forward and reverse switch from "
            "successive snapshots to decorrelate the work samples. Defaults to 5."
        ),
    )
    parser.add_argument(
        "-neqs",
        "--neq-steps",
        dest="neq_steps",
        type=int,
        default=50000,
        help=(
            "NEQ only: length of each lambda-switching simulation in MD steps. Recommended >16000. "
            "Defaults to 50000."
        ),
    )
    parser.add_argument(
        "-neqes",
        "--neq-eq-steps",
        dest="neq_eq_steps",
        type=int,
        default=1000,
        help=(
            "NEQ only: number of endpoint equilibration steps between successive switches "
            "(spacing). Recommended >250. Defaults to 1000."
        ),
    )
    parser.add_argument(
        "-neqrs",
        "--neq-relax-steps",
        dest="neq_relax_steps",
        type=int,
        default=5000,
        help=(
            "NEQ only: length (MD steps) of the one-time endpoint relaxation run at lambda=0 "
            "and lambda=1 before the first switch, to settle the nearly-decoupled ligand at each "
            "endpoint. Applied to the first switching iteration; later iterations use the shorter "
            "--neq-eq-steps (tEQ) spacing. ~10 ps at 2 fs. Defaults to 5000."
        ),
    )
    parser.add_argument(
        "-L",
        "--neq-steepness",
        dest="neq_L",
        type=float,
        default=8.0,
        help=(
            "NEQ only: steepness `L` of the sigmoidal lambda schedule l(t) = 1/[1+e^(L(t-0.5))]. "
            "Higher L spends more time near lambda=0 and lambda=1; lower L approaches a linear "
            "schedule. Recommended between 4 and 16. Defaults to 8. Ignored when "
            "`--neq-schedule linear` is used."
        ),
    )
    parser.add_argument(
        "-neqsched",
        "--neq-schedule",
        dest="neq_schedule",
        default="sigmoidal",
        choices=["sigmoidal", "linear"],
        help="NEQ only: lambda switching schedule. Defaults to `sigmoidal`.",
    )
    parser.add_argument(
        "-sc",
        "--softcore-method",
        dest="softcore_method",
        type=str,
        default="standard",
        choices=["standard", "beutler_coul", "gapsys"],
        help=(
            "Soft-core method for nonbonded interactions during FEP. "
            "'standard' applies soft-core only to LJ (current default). "
            "'beutler_coul' extends soft-core to Coulomb via a modified effective distance. "
            "'gapsys' uses force-based linearization below a critical radius for both LJ and Coulomb."
        ),
    )
    parser.add_argument(
        "-cm",
        "--charge-method",
        dest="charge_method",
        type=str,
        default="ion_match",
        choices=["none", "ion_match", "coalchemical_water"],
        help=(
            "Strategy for handling FEP edges that change the ligand net charge. "
            "'ion_match' (default) adds Cl-/Na+ counter-ions to the water leg so its total "
            "charge matches the protein leg. "
            "'coalchemical_water' turns a real water in the water leg into a co-alchemical Na+/Cl- "
            "(O <-> ion, H <-> DUM). "
            "'none' disables all neutralization (raw ddG retains the Born artifact)."
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
    if program == "QligFEP":
        parser.add_argument(
            "-pq",
            "--protein_charge",
            dest="protein_charge",
            type=int,
            default=None,
            help=(
                "Total charge from the protein leg's qprep output. Used by setupFEP to pass "
                "the protein system charge so the water leg can match it with counter-ions."
            ),
        )
    return parser.parse_args()
