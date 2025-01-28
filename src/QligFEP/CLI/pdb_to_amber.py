"""Module with utility functions (and CLI) to rename .pdb files so that they're compatible with the AMBER forcefield.

!! Note !! The functions:
rename_charged, nest_pdb, unnest_pdb, get_coords, pdb_cleanup, histidine_search,
atom_is_present

are all Python3 adaptations from the original repository / module:
https://github.com/choderalab/mmtools/blob/master/mccetools/rename.py

to use this module, you can use the following code:
`pdb_to_amber -i <pdb_file>`
"""

import argparse
import os
import re
import sys
from pathlib import Path

from ..logger import logger
from ..pdb_utils import (
    nest_pdb,
    read_pdb_to_dataframe,
    unnest_pdb,
    write_dataframe_to_pdb,
)


def reindex_pdb_residues(pdb_path: Path, out_pdb_path: str):
    pdb_df = read_pdb_to_dataframe(pdb_path)
    uniq_indexes = pdb_df.set_index(
        ["residue_seq_number", "residue_name", "chain_id", "insertion_code"]
    ).index
    resn_mapping = {resn: idx for idx, resn in enumerate(uniq_indexes.unique(), 1)}
    pdb_df["residue_seq_number"] = uniq_indexes.map(resn_mapping)
    pdb_df["insertion_code"] = ""
    # pdb_df = pdb_df.assign(residue_seq_number=uniq_indexes.map(resn_mapping))
    write_dataframe_to_pdb(pdb_df, out_pdb_path)


def correct_numbered_atom_names(npdb_i):
    """Corrects atom names that start with numbers by moving the numbers to the end.
    Uses regex to match and extract leading numbers.

    Args:
        npdb_i: nested pdb data structure for a single residue

    Returns:
        Modified npdb_i with corrected atom names
    """

    def process_atom_name(line):
        atom_name = line[12:16].strip()

        # these only exist in AMBER with 2 and 3 for some reason (?)
        sum_after = atom_name in [
            "2HG",
            "1HG",
            "2HB",
            "1HB",
            "1HG1",
            "2HG1",
            "1HA",
            "2HA",
            "1HD",
            "2HD",
            "1HE",
            "2HE",
        ]

        pattern = re.compile(r"^(\d+)([A-Z]+\d*)")
        match = pattern.match(atom_name)

        if not match:
            return line

        # Extract the matched groups
        numbers, letters = match.groups()
        new_atom_name = letters + (str(int(numbers) + 1) if sum_after else numbers)

        # Format according to PDB specifications
        if len(new_atom_name) == 4:
            return line[:12] + new_atom_name + line[16:]
        else:
            return line[:12] + f"{new_atom_name:<4}" + line[16:]

    return [process_atom_name(line) for line in npdb_i]


def correct_amino_acid_atom_names(npdb_i, resname, rename_mapping):
    """corrects the amino acid atom names according to the mapping provided

    Args:
        npdb_i: nested pdb data structure for a single residue
        resname: the residue name
        rename_mapping: a dictionary mapping old names to new names
    """
    if resname in rename_mapping:
        for old_name, new_name in rename_mapping[resname].items():
            npdb_i = [extract_and_replace(x, old_name, new_name) for x in npdb_i]
            # certify that we have the alignment as expected for pdb files
    return npdb_i


def extract_and_replace(line, old_name, new_name):
    """extracts the atom name and replaces it with the new name"""
    atom_name = line[12:16].strip()
    if atom_name != old_name:
        return line
    new_atom_name = atom_name.replace(old_name, new_name).strip()
    if len(new_atom_name) == 4:
        return line[:12] + new_atom_name + line[16:]
    else:
        # return left aligned atom name always with len() == 3 but with a " " in the beginning
        return line[:12] + f" {new_atom_name:<3}" + line[16:]


def fix_pdb(pdb_path: Path, rename_mapping):
    renamed_pdb_path = pdb_path.with_name(pdb_path.stem + "_renamed.pdb")
    with open(pdb_path) as f:
        pdb_lines = f.readlines()

    npdb = nest_pdb(pdb_lines)
    npdb = asp_search(npdb)
    npdb = glu_search(npdb)  # TODO; check if this one is necessary
    npdb = histidine_search(npdb)

    for i, res in enumerate(npdb):
        resname = res[-1][17:21].rstrip()
        if resname == "NMA":  # we use NME in our FF library
            npdb[i] = [x.replace("NMA", "NME") for x in npdb[i]]
            resname = "NME"
        npdb[i] = correct_numbered_atom_names(npdb[i])
        npdb[i] = correct_amino_acid_atom_names(npdb[i], resname, rename_mapping)

    npdb = nc_termini_search(npdb)  # after atom name correction, label N and C termini
    npdb = correct_neutral_arginine(npdb)
    pdb_lines = unnest_pdb(npdb)

    with open(renamed_pdb_path, "w") as f:
        for line in pdb_lines:
            f.write(line)
    return pdb_lines


def correct_neutral_arginine(pdb_arr):
    """
    Updates residue naming for ARN (neutral arginine - AMBER14sb) residues to ensure
    compatibility.

    Will check for HH22 atom in ARN residues. If HH22 is found, it indicates that the
    NH1 and NH2 groups (and their associated hydrogens) should be renamed to match the
    naming convention in the force field (only HH11, HH12, HH21 should exist). This
    renaming ensures that the nitrogen and hydrogen atoms in the side chain of arginine
    are correctly identified and interact as expected according to the simulation parameters.

    Parameters:
    - pdb_lines (list of str): The original PDB content as a list of strings, where each
      string represents a line in the PDB file.

    Returns:
    - list of str: The updated PDB content with corrected atom names for ARN residues,
      ready for use in molecular dynamics simulations.
    """
    # Placeholder for updated PDB lines
    updated_pdb_lines = []

    # Check if any ARN residue contains HH22
    has_hh22 = any("ARN" in line and "HH22" in line for line in pdb_arr)

    # Only proceed with renaming if HH22 is present
    if not has_hh22:
        return pdb_arr  # Return original lines if no HH22 found

    # Temporary mapping for the first pass
    atom_name_replacements = {
        "NH1": "NHT ",  # Temporary name for NH1
        "HH11": "HHT1",
        "NH2": "NH1 ",
        "HH21": "HH11",  # Temporary name for HH11
        "HH22": "HHT2",  # Temporary name for HH22 to avoid direct swap conflict
    }
    atom_name_replacements = {}

    # Reverse mapping for the second pass
    final_name_replacements = {
        "NHT ": "NH2 ",
        "HHT1": "HH21",
        "HHT2": "HH12",  # Assuming HH22 should be renamed to HH12 if present
    }

    for line in pdb_arr:  # First pass: Apply initial renaming
        if line.startswith("ATOM") and "ARN" in line:
            atom_name = line[12:16].strip()  # Extract atom name
            if atom_name in atom_name_replacements:
                new_atom_name = atom_name_replacements[atom_name]
                line = line[:12] + new_atom_name + line[16:]
        updated_pdb_lines.append(line)

    # Second pass: Replace temporary names with their final names
    final_updated_pdb_lines = []
    for line in updated_pdb_lines:
        if line.startswith("ATOM") and "ARN" in line:
            atom_name = line[12:16].strip()  # Extract atom name
            if atom_name in final_name_replacements:
                new_atom_name = final_name_replacements[atom_name].ljust(4)
                line = line[:12] + new_atom_name + line[16:]
        final_updated_pdb_lines.append(line)

    return final_updated_pdb_lines


def rename_charged(npdb):
    """Generate AMBER-specific residue names for charged residues from MCCE residue names. Also fix some problems with atom naming (specifically hydrogens) for charged residues.

    ARGUMENTS
        npdb - "nested" PDB data structure generated by nest_pdb. Modified to reflect AMBER names.

    OPTIONAL ARGUMENTS
        terminology       default 'AMBER', which protonates residues using naming recognized by AMBER and the ffamber ports for GROMACS from the Sorin lab.
                          Optionally, specify 'gAMBER' for the GROMACS AMBER format, which has some slightly different residue naming conventions
                          (for example, LYP -> LYS, LYS -> LYN, and CYN -> CYS, CYS2 -> CYM )

    CHANGE LOG:
    - DLM 7-1-2009: Modified to fix naming of HE1 in GLH to HE2 to conform to ffamber rtp.
    - DFV 20-02-2024: adapted to python3 & added logger.
    """
    for i, res in enumerate(npdb):
        original_resname = res[0][17:21].rstrip()
        resname = res[-1][17:21].rstrip()
        npdb[i] = correct_amino_acid_atom_names(npdb[i], resname)
        new_resname = res[0][17:21].rstrip()  # keep track for the log message
        if original_resname != new_resname:
            logger.info(f"Residue {i+1}: {original_resname} renamed to {new_resname}.")
    return npdb


def histidine_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "HIS":
            HE_present = atom_is_present(npdb[i], "HE2")  # bonded to NE2
            HD_present = atom_is_present(npdb[i], "HD1")  # bonded to ND1

            if HD_present and HE_present:
                npdb[i] = [x.replace("HIS", "HIP") for x in npdb[i]]
            elif HE_present:
                npdb[i] = [x.replace("HIS", "HIE") for x in npdb[i]]
            elif HD_present:
                npdb[i] = [x.replace("HIS", "HID") for x in npdb[i]]
            else:
                raise ValueError("No protons found for histidine.")
    return npdb


def nc_termini_search(npdb):

    NATURAL_AA = (
        "ALA;ARG;ASH;ASN;ASP;CYM;CYS;CYX;GLH;GLN;GLU;GLY;HID;HIE;HIP;"
        "HYP;ILE;LEU;LYN;LYS;MET;PHE;PRO;SER;THR;TRP;TYR;VAL"
    ).split(";")

    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname in NATURAL_AA:
            if resname in ["CYM", "ASH", "GLH", "LYN"]:
                continue  # no parameter for those on C or N terminus
            H3_present = atom_is_present(npdb[i], "H3")  # n-terminus
            OXT_present = atom_is_present(npdb[i], "OXT")  # c-terminus
            if H3_present:
                if resname == "HYP":
                    logger.error(
                        "No parameters available for n-terminal HYP residue!!! Please check your structure"
                    )
                else:
                    npdb[i] = [x.replace(f"{resname} ", f"N{resname}") for x in npdb[i]]
            if OXT_present:
                npdb[i] = [x.replace(f"{resname} ", f"C{resname}") for x in npdb[i]]
            if H3_present and OXT_present:
                raise ValueError(f"residue {npdb[i]} has both H3 and OXT atoms")
    return npdb


def glu_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "GLU":
            HE1_present = atom_is_present(npdb[i], "HE1")
            if HE1_present:
                npdb[i] = [x.replace("GLU", "GLH") for x in npdb[i]]
    return npdb


def asp_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "ASP":
            OD1_present = atom_is_present(npdb[i], "OD1")
            OD2_present = atom_is_present(npdb[i], "OD2")
            if all([OD1_present, OD2_present]):
                npdb[i] = [x.replace("ASP", "ASH") for x in npdb[i]]
    return npdb


def atom_is_present(pdblines, atomname):
    """Returns TRUE if the given atom is present in the given PDB atom lines.

    ARGUMENTS
        pdblines - list of PDB lines
        atomname - the name of the atom to check the existence of

    RETURNS
        is_present - True if the given atom name is present, False otherwise

    """
    atoms = [pdbline[12:16] for pdbline in pdblines]
    return bool(any(atomname in atom for atom in atoms))


def main_exe():
    parser = argparse.ArgumentParser(
        description="Rename amino acids in a PDB file for AMBER forcefield compatibility."
    )
    parser.add_argument("-i", "--input", required=True, help="Input PDB file to process")
    parser.add_argument(
        "-o", "--output", required=False, help="Output name for the processed file", default=None
    )
    args = parser.parse_args()

    input_pdb_path = args.input
    if not os.path.exists(input_pdb_path):
        logger.error(f"Input file {input_pdb_path} does not exist.")
        sys.exit(1)

    logger.info(f"Processing file: {input_pdb_path}")

    try:
        with open(input_pdb_path) as f:
            pdb_lines = f.readlines()

        renamed_pdb_lines = rename_residues(pdb_lines)
        if args.output is not None:
            output_pdb_path = args.output
        else:
            output_pdb_path = input_pdb_path.replace(".pdb", "_renamed.pdb")

        with open(output_pdb_path, "w") as f:
            for line in renamed_pdb_lines:
                f.write(line)

        logger.info(f"Renaming completed. Output file: {output_pdb_path}")

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main_exe()
