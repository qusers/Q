"""Module with the CLI to align ligands to a reference ligand using kcombu."""

import argparse
import sys
from pathlib import Path

import numpy as np
from MolClusterkit.mcs import MCSClustering

from ..lig_aligner import LigandAligner
from ..logger import logger, setup_logger


def parse_arguments() -> argparse.Namespace:
    """Method to parse the arguments."""
    parser = argparse.ArgumentParser(
        description="Lomap wrapper to be used in  QligFEP."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        help=(
            "Input file (sdf) or directory containing the mol files (sdf|mol2) to "
            "be aligned. If sdf, a new dir is created, containing separate sdf files "
            "for each ligand."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help=(
            "Output file (sdf) to store all the aligned ligands. If not provided, defaults to "
            "`{input}_{reference}_aligned.sdf`, stored in the current working directory."
        ),
    )
    parser.add_argument(
        "--reference",
        "-ref",
        type=str,
        help=(
            "Name of the reference ligand to which the ligands will be aligned. If the `input` is "
            "a directory, it should be the name of the file without the `.sdf` suffix. If the "
            "input is a single sdf file with several ligands, the input is the name of that "
            "ligand. If left as None (Default), the program will score the ligands against "
            "each other using a MCS algorithm to find the reference ligand with highest MCS "
            "score with the rest of the ligands."
        ),
        dest="ref",
        default=None,
    )
    parser.add_argument(
        "-rm",
        "--remove_temp_dir",
        dest="rm",
        help=(
            "If set, the temporary directory created to store the separate sdf "
            "files prior to alignment will be removed. Defaults to True."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--parallel",
        dest="parallel",
        default=1,
        help="Number of threads to parallelize operations. Defaults to 1.",
        type=int,
    ),
    parser.add_argument(
        "-log",
        "--log-level",
        dest='log',
        required=False,
        default="info",
        help="Set the log level for the logger. Defaults to `info`.",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
    )
    return parser.parse_args()


def main(args: argparse.Namespace):
    """Main function to align the ligands."""

    # setup the logger with the desired log level
    setup_logger(level=args.log)

    aligner = LigandAligner(
        args.input,
        n_threads=args.parallel,
        tempdir="to_align_ligands",
        delete_tempdir=args.rm,
    )
    if aligner.lig_files != []:
        is_aligned = [f.name.endswith("_aligned.sdf") for f in aligner.lig_files]
        if any(is_aligned):
            logger.warning(
                "Some of the files in the input directory seem to be already aligned. "
                "Please remove them or change the input directory."
            )
            aligned_files = np.where(is_aligned)[0]
            aligned_files = [aligner.lig_files[idx] for idx in aligned_files]
            logger.info(f"Aligned files: {aligned_files}")
            sys.exit(1)


    if args.ref is None:
        smiles_list = [mol.to_smiles() for mol in aligner.molecules]
        mcs_kwargs = {
            "atomCompare": "CompareElements",
            "bondCompare": "CompareAny",
            "ringCompare": "StrictRingFusion",
        }
        mcs_cluster = MCSClustering(smiles_list, **mcs_kwargs)
        logger.info(
            "Reference not provided. Scoring ligands based on Maximum Common "
            "Substructure to find the reference."
        )
        _, simi_matrix = mcs_cluster.compute_similarity_matrix()

        highest_score_idx = simi_matrix.sum(axis=1).argmax()
        reference_name = aligner.lig_names[highest_score_idx]
    else:
        reference_name = args.ref
        if reference_name not in aligner.lig_names:
            logger.error(
                "The reference ligand is not in the input ligands. See the ligands names with the `lig_names` attribute."
            )
            logger.info(f"lig_names; {aligner.lig_names}")
            raise ValueError(
                f"The reference ligand {reference_name} is not in the input ligands."
            )
        logger.info(f"Aligning ligands to reference: {reference_name}")

    # run the alignment using kcombu
    aligner.kcombu_align(reference_name)

    output = Path(args.output)
    if output.exists():
        logger.warning(
            f"The output file {output} already exists. It will be overwritten."
        )

    aligner.output_aligned_ligands(str(output))


def main_exe():
    """Main function to be executed by the CLI."""
    args = parse_arguments()
    main(args)
