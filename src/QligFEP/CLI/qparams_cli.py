"""Module containing the QligFEP parameter writing functionalities."""

import argparse

from ..logger import setup_logger
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
        "-lff",
        "--ligand-forcefield",
        dest="lff",
        default="OpenFF",
        choices=[
            "OpenFF",
        ],
        help="Forcefield to be used to parametrize ligands. Defaults to OpenFF.",
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
        type=int,
    )
    parser.add_argument(
        "-am1bcc",
        "--use_am1bcc",
        dest="am1bcc",
        action="store_true",
        help=(
            "Use the AM1-BCC method to calculate partial charges instead of NAGL. "
            "NAGL is the default method and is faster, but AM1-BCC may be needed "
            "for molecules outside NAGL's training domain."
        ),
    )
    parser.add_argument(
        "-ff",
        "--forcefield",
        dest="forcefield",
        default="openff-2.3.0.offxml",
        help=(
            "OpenFF forcefield file to use for ligand parameters. Defaults to openff-2.3.0.offxml. "
            "See https://github.com/openforcefield/openff-forcefields for available forcefields."
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
    parser.add_argument(
        "-pcof",
        "--parametrize-cofactors",
        dest="pcof",
        action="store_true",
        help=(
            "Parametrize cofactors. Use this argument together with `-pff` to create .prm and "
            ".lib files with both protein and cofactor parameters"
        ),
    )
    parser.add_argument(
        "-pff",
        "--protein-forcefield",
        dest="pff",
        default="AMBER14sb",
        help=(
            "To be used with -pcof. Passing a protein forcefield will create both .prm and .lib "
            "files containing the cofactors' and the protein parameters."
        ),
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    setup_logger(level=args.log)
    if args.lff == "OpenFF":
        openff2q = OpenFF2Q(
            args.input, nagl=not args.am1bcc, n_jobs=args.parallel, forcefield=args.forcefield
        )
        openff2q.process_ligands()
        if args.pcof:
            openff2q.write_cofactor_plus_ff_files(args.pff)
        else:
            openff2q.write_ligand_files()
    else:
        raise NotImplementedError("Forcefield not supported through this CLI yet")


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
