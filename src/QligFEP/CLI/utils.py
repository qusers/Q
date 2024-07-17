"""Module containing utility functions used by the command line interfaces of QligFEP."""

from pathlib import Path

from ..logger import logger
from ..pdb_utils import disulfide_search, nest_pdb, unnest_pdb


def cysbonds_for_qprep(pdb_file: Path, comment_out: bool = True):
    """Method to search for disulfide bonds in a pdb file and return the bonds in the format
    `!addbond at1 at2 y\n` used in the `qprep.inp` file.

    Args:
        pdb_file: path to the pdb file to search for disulfide bonds.
        comment_out: if True, the bonds will be commented out as in `!addbond` so qprep
            will ignore them. If False, the bonds will be uncommented.

    Returns:
        str: string containing the disulfide bonds in the format `!addbond at1 at2 y\n`.
    """
    with open(pdb_file) as f:
        pdb_lines = f.readlines()
        npdb = nest_pdb(pdb_lines)
        npdb, cysbonds, renamed = disulfide_search(npdb)
        if renamed:
            logger.info("Disulfide bonds detected: renaming CYS to CYX")
            pdbarr = unnest_pdb(npdb)
            with open(pdb_file, "w") as file_out:
                for line in pdbarr:
                    file_out.write(f"{line}\n")
    addbond = "!addbond" if comment_out else "addbond"
    cysbonds = "".join([f"{addbond} {atomN[0]} {atomN[1]} y\n" for atomN in cysbonds])
    return cysbonds


def handle_cysbonds(input_config: str, pdb_file: Path, comment_out: bool = True):
    """Handle cysbond input. If the input is `auto`, it will search for disulfide bonds
    in the pdb file. If the input is `none`, it will return an empty string so that no
    bonds are added. If the input is a string with the format `at1_at2,at3_at4`, it
    will return the bonds in the format `!addbond at1 at2 y\n!addbond at3 at4 y\n`.

    Args:
        input_config: configuration for the cysbonds. Can be `auto`, `none`, or a string
            with the format `at1_at2,at3_at4`.
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
    elif input_config != "":
        input_config = input_config.split(",")
        qprep_lines = "".join([f"!addbond {b.split('_')[0]} {b.split('_')[1]} y\n" for b in input_config])
    elif input_config.lower() == "none":  # TODO: do this in a smarter way...
        qprep_lines = ""
    else:
        raise ValueError(f"Invalid cysbond input: {input_config}. Please check the input format.")
    return qprep_lines
