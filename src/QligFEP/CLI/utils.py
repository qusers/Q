"""Module containing utility functions used by the command line interfaces of QligFEP."""

from itertools import product
from pathlib import Path

from ..logger import logger
from ..pdb_utils import disulfide_search, nest_pdb, unnest_pdb


def get_avail_restraint_methods():
    ring_compare_methods = ["aromaticity", "hybridization", "element"]
    surround_compare_methods = ["p", "ls", "strict"]
    return list(map(lambda x: "_".join(x), product(ring_compare_methods, surround_compare_methods))) + [
        "overlap",
        "kartograf",
    ]


def cysbonds_for_qprep(pdb_file: Path, comment_out: bool = True):
    """Method to search for disulfide bonds in a pdb file and return the bonds in the format
    `!addbond atom1 atom2 y\n` used in the `qprep.inp` file.

    Args:
        pdb_file: path to the pdb file to search for disulfide bonds.
        comment_out: if True, the bonds will be commented out as in `!addbond` so qprep
            will ignore them. If False, the bonds will be uncommented.

    Returns:
        str: string containing the disulfide bonds in the format `!addbond atom1 atom2 y\n`.
    """
    with open(pdb_file) as f:
        pdb_lines = f.readlines()
        npdb = nest_pdb(pdb_lines)
        npdb, cysbonds, renamed = disulfide_search(npdb)
        if renamed:
            logger.warning("Disulfide bonds detected on CYS residues! Renaming to CYX")
            pdbarr = unnest_pdb(npdb)
            with open(pdb_file, "w") as file_out:
                for line in pdbarr:
                    file_out.write(f"{line}")
    addbond = "!addbond" if comment_out else "addbond"
    cysbonds = "".join([f"{addbond} {bond[0]} {bond[1]} y\n" for bond in cysbonds])
    return cysbonds


def handle_cysbonds(input_config: str, pdb_file: Path, comment_out: bool = True):
    """Handle cysbond input. If the input is `auto`, it will search for disulfide bonds
    in the pdb file. If the input is `none`, it will return an empty string so that no
    bonds are added. If the input is a string with the format `atom1_atom2,atom3_atom4`, it
    will return the bonds in the format `!addbond x1 x2 y\n!addbond x3 x4 y\n`.

    Note: atoms in the bond can be defined as either `ResNumber:SG` or `AtomNumber` in Qprep.
    We opt for the former as it's less error-prone.

    Args:
        input_config: configuration for the cysbonds. Can be `auto`, `none`, or a string
            with the format `atom1_atom2,atom3_atom4`.
        pdb_file: path to the pdb file to search for disulfide bonds.
        comment_out: if True, the bonds will be commented out as in `!addbond` so qprep
            will ignore them. If False, the bonds will be uncommented.

    Raises:
        ValueError: if the method fails to parse the input_config

    Returns:
        str: string containing the cysbonds in the format `!addbond at1 at2 y\n`.
    """

    if input_config == "auto":
        qprep_lines = cysbonds_for_qprep(pdb_file, comment_out=comment_out)
    elif input_config.lower() == "none":  # TODO: do this in a smarter way...
        qprep_lines = ""
    elif input_config != "":
        input_config = input_config.split(",")
        qprep_lines = "".join([f"!addbond {b.split('_')[0]} {b.split('_')[1]} y\n" for b in input_config])
    else:
        raise ValueError(f"Invalid cysbond input: {input_config}. Please check the input format.")
    return qprep_lines
