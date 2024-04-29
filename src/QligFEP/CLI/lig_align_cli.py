from ..lig_aligner import LigandAligner
import argparse
from pathlib import Path

def parse_arguments() -> argparse.Namespace:
    """Method to parse the arguments."""
    parser = argparse.ArgumentParser(description='Lomap wrapper to be used in  QligFEP.')
    parser.add_argument('--input',
                        '-i',
                        type=str,
                        help=(
                            'Input file (sdf) or directory containing the mol files (sdf|mol2) to '
                            'be aligned. If sdf, a new dir is created, containing separate sdf files '
                            'for each ligand.'
                            ))
    parser.add_argument('--output',
                        '-o',
                        type=str,
                        help=(
                            'Output file (json) to store the lomap results. If not provided, '
                            'defaults to `lomap.json`, stored in the directory with the separate sdf|mol2 files.'
                            )
                        )
    parser.add_argument('--reference',
                        '-ref',
                        type=str,
                        help=(
                            "Name of the reference molecule to which the ligands will be aligned. If the `input` is "
                            "a directory, it should be the name of the file without the `.sdf` suffix. If the "
                            "input is a single sdf file with several molecules, the input is the name of that "
                            "molecule."
                            ),
                        default=None
                        )
    return parser.parse_args()