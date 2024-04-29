"""Module with the CLI to align ligands to a reference ligand using kcombu."""

import argparse
from pathlib import Path

from MolClusterkit.mcs import MCSClustering

from ..lig_aligner import LigandAligner
from ..logger import logger


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
            "Name of the reference molecule to which the ligands will be aligned. If the `input` is "
            "a directory, it should be the name of the file without the `.sdf` suffix. If the "
            "input is a single sdf file with several molecules, the input is the name of that "
            "molecule. If left as None (Default), the program will score the ligands against "
            "each other using a MCS algorithm to find the reference molecule with highest MCS "
            "score with the rest of the ligands."
        ),
        dest="ref",
        default=None,
    )
    parser.add_argument(
        "-p",
        "--parallel",
        dest="parallel",
        default=4,
        help="Number of jobs to parallelize operations. Defaults to 4.",
        type=int,
    )
    return parser.parse_args()


def main(args: argparse.Namespace):
    """Main function to align the ligands."""
    aligner = LigandAligner(
        args.input, n_threads=4, temp_ligalign_dir="to_align_ligands"
    )

    if args.ref is None:
        smiles_list = [mol.to_smiles() for mol in aligner.molecules]
        mcs_kwargs = {
            "AtomCompare": "CompareAny",
            "BondCompare": "CompareAny",
            "RingCompare": "StrictRingFusion",
        }
        mcs_cluster = MCSClustering(smiles_list, **mcs_kwargs)
        _, simi_matrix = mcs_cluster.compute_similarity_matrix()

        highest_score_idx = simi_matrix.sum(axis=1).argmax()
        reference_name = aligner.lig_names[highest_score_idx]
    else:
        reference_name = args.ref
        if reference_name not in aligner.lig_names:
            logger.error(
                "The reference molecule is not in the input ligands. See the ligands names with the `lig_names` attribute."
            )
            logger.info(f"lig_names; {aligner.lig_names}")
            raise ValueError(
                f"The reference molecule {reference_name} is not in the input ligands."
            )

    # run the alignment using kcombu
    aligner.kcombu_align(reference_name)

    output = Path(args.output)
    if output.exists():
        logger.warning(
            f"The output file {output} already exists. It will be overwritten."
        )

    aligner.output_aligned_molecules(str(output))
