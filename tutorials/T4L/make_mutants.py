"""Write a PyMOL script that builds each mutant side chain as its own PDB file.

QresFEP needs the mutant residue positioned on the wild-type backbone, which is
what PyMOL's mutagenesis wizard produces. One file per mutation, named
``<MUTANT><POSITION>.pdb``.

    python make_mutants.py 2LZM_prep.pdb mutations_neutral.txt
    pymol -c mutagenesis.pml
"""

import re
import sys
from pathlib import Path

_TEMPLATE = """
reinitialize
cmd.load('{structure}')
cmd.wizard('mutagenesis')
cmd.do('refresh_wizard')
cmd.get_wizard().set_mode('{mutant}')
cmd.get_wizard().do_select('resi {position}')
cmd.get_wizard().apply()
save {mutant}{position}.pdb, resi {position}
"""


def main(structure: str, mutation_file: str, output: str = "mutagenesis.pml") -> None:
    """Write the PyMOL script for every mutation listed in `mutation_file`."""
    blocks = []
    for entry in Path(mutation_file).read_text().split():
        parts = re.split(r"(\d+)", entry.strip())
        if len(parts) != 3:
            raise ValueError(f"Could not read {entry!r} as <wild-type><position><mutant>")
        _, position, mutant = parts
        blocks.append(_TEMPLATE.format(structure=structure, mutant=mutant, position=position))

    blocks.append("cmd.quit()\n")
    Path(output).write_text("".join(blocks))
    print(f"{output} written for {len(blocks) - 1} mutation(s). Now run: pymol -c {output}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <structure.pdb> <mutations.txt>")
    main(sys.argv[1], sys.argv[2])
