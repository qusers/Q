"""Write the PyMOL mutagenesis script for the chain-A T4L examples.

    python make_mutants.py 2LZM_prep.pdb mutations_neutral.txt
    pymol -c mutagenesis.pml
"""

import sys
from pathlib import Path

from QligFEP.resfep_setup import read_mutations, write_mutagenesis_script


CHAIN = "A"


def main(structure: str, mutation_file: str, output: str = "mutagenesis.pml") -> None:
    """Write one standard QresFEP mutagenesis block per listed mutation."""
    mutations = read_mutations(Path(mutation_file), CHAIN)
    write_mutagenesis_script(Path(structure), mutations, Path(output))
    print(f"{output} written for {len(mutations)} mutation(s). Now run: pymol -c {output}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <structure.pdb> <mutations.txt>")
    main(sys.argv[1], sys.argv[2])
