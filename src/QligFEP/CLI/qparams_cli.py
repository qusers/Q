"""Module containing the QligFEP parameter writing functionalities."""

import argparse

from ..openff2Q import OpenFF2Q

# TODO: implement the other parameter writing functionalities
# see issue: https://github.com/qusers/Q/issues/26


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ligand parameter files for QligFEP according to the chosen forcefield"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input (sdf) file containing one or more ligands.",
    )
    parser.add_argument(
        "-FF",
        "--forcefield",
        dest="forcefield",
        default="OpenFF",
        choices=[
            "OPLS2005",
            "OPLS2015",
            "AMBER14sb",
            "CHARMM36",
            "CHARMM22",
            "CHARMM_TEST",
            "OpenFF",
        ],
        help="Forcefield to be used. Defaults to OpenFF.",
    )
    parser.add_argument(
        "-p",
        "--parallel",
        dest="parallel",
        default=1,
        help=(
            "Number of jobs to parallelize operations. Currently applied to calculate "
            "the charges of the molecules using OpenFF. Defaults to 1."
        ),
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    if args.forcefield == "OpenFF":
        openff2q = OpenFF2Q(args.input, n_jobs=args.parallel)
        openff2q.process_ligands()
    else:
        raise NotImplementedError("Forcefield not supported through this CLI yet")


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
