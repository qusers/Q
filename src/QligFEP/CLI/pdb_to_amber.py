"""Module with utility functions (and CLI) to rename .pdb files so that they're compatible with the AMBER forcefield.

!! Note !! The functions:
rename_residues, rename_charged, labeledPDB_to_AmberPDB, renumber_atoms,
nest_pdb, unnest_pdb, get_coords, disulfide_search, pdb_cleanup, histidine_search,
atom_is_present

are all Python3 adaptations from the original repository / module:
https://github.com/choderalab/mmtools/blob/master/mccetools/rename.py

to use this module, you can use the following code:
`pdb_to_amber -i <pdb_file>`
"""

import argparse
from typing import List
from ..logger import logger
import os
import sys
import math


def rename_residues(pdbarr):
    logger.debug(f"At start of renaming, pdbarr has {len(pdbarr)} items")
    npdb = nest_pdb(pdbarr)
    logger.debug(f"Nest/unnest leads to {len(unnest_pdb(npdb))}, items")
    npdb = rename_charged(npdb)
    npdb = histidine_search(npdb)
    npdb = disulfide_search(npdb)
    pdbarr = unnest_pdb(npdb)
    pdbarr = correct_neutral_arginine(pdbarr)
    return pdbarr

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
        "NH1": "NHT",  # Temporary name for NH1
        "HH11": "HHT1",  # Temporary name for HH11
        "NH2": "NH1",
        "HH21": "HH11",
        "HH22": "HHT2",  # Temporary name for HH22 to avoid direct swap conflict
    }

    # Reverse mapping for the second pass
    final_name_replacements = {
        "NHT": "NH2",
        "HHT1": "HH21",
        "HHT2": "HH12",  # Assuming HH22 should be renamed to HH12 if present
    }

    for line in pdb_arr: # First pass: Apply initial renaming
        if line.startswith("ATOM") and "ARN" in line:
            atom_name = line[12:16].strip()  # Extract atom name
            if atom_name in atom_name_replacements:
                new_atom_name = atom_name_replacements[atom_name].ljust(4)
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

def correct_amino_acid_atom_names(npdb_i, resname):
    """corrects the amino acid atom names according to the charge state and the force field
    terminology. This function is called within rename_changed() where pdb_i is the
    nested pdb data structure for a single residue. The function modifies the
    residue name in place.

    Args:
        npdb_i: nested pdb data structure for a single residue
        resname_and_state: the residue name and charge state, e.g. 'LYS0' or 'LYS+'
        terminology: which FF terminology will be used. Defaults to 'AMBER'.
    """
    if resname=='GLH':
        npdb_i=[x.replace('HE1  GLH','HE2  GLH') for x in npdb_i]
    elif resname=='ASH':
        npdb_i=[x.replace('HD1  ASH','HD2  ASH') for x in npdb_i]
    elif resname=='LYN':
        npdb_i=[x.replace('HE1  GLH','HE2  GLH') for x in npdb_i]
    elif resname=='ARN':
        npdb_i=[x.replace('HE1  GLH','HE2  GLH') for x in npdb_i]
    elif resname=='LYN':
        # HZ1 and HZ2 should be renamed to HZ2 and HZ3...
        npdb_i=[x.replace('HZ2  LYN','HZ3  LYN') for x in npdb_i]
        npdb_i=[x.replace('HZ1  LYN','HZ2  LYN') for x in npdb_i]
    return npdb_i

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
        original_resname = res[0][17:20]
        resname = res[-1][17:20]
        npdb[i] = correct_amino_acid_atom_names(npdb[i], resname)
        new_resname = res[0][17:20] # keep track for the log message
        if original_resname != new_resname:
            logger.info(f"Residue {i+1}: {original_resname} renamed to {new_resname}.")
    return npdb

def labeledPDB_to_AmberPDB(labeledPDBfile, outPDBfile, renameResidues=True):
    with open(labeledPDBfile, 'r') as fin:
        lines = fin.readlines()
    pdbarr = []
    while lines:
        pdbline = lines.pop(0)
        label = lines.pop(0).strip()
        pdbarr.append(pdbline.strip() + label)
    if renameResidues:
        pdbarr = rename_residues(pdbarr)
    pdbarr = pdb_cleanup(pdbarr)
    with open(outPDBfile, 'w') as fout:
        for line in pdbarr:
            fout.write(f"{line}\n")
            
def renumber_atoms(pdbarr):
    for i, line in enumerate(pdbarr):
        pdbarr[i] = f"{line[:6]}{str(i+1).rjust(5)}{line[11:]}"

def nest_pdb(pdbarr: List[str]) -> List[List[str]]:
    """Organizes a flat list of PDB (Protein Data Bank) file lines into a nested structure
    grouped by residues. This function takes a list of strings, each representing a line
    from a PDB file, and groups these lines by residue. Each residue's lines are collected
    based on continuity of residue identifiers and uniqueness of atom names within the
    residue. This nested structure is useful for operations that require manipulation or
    analysis on a per-residue basis.

    args:
        pdbarr: A list where each element is a string representing a line from a PDB file.

    Returns:
    - nestedpdb : Each inner list contains all the lines from the input 
        corresponding to a single residue. The grouping is based on residue identifiers
        (including residue name, chain identifier, and residue sequence number) and ensures
        that each atom within a residue is unique.

    Notes:
    - The function assumes that the input list is ordered as it would be in a standard PDB file,
        where lines corresponding to atoms of the same residue are consecutive.
    """
    nestedpdb = []
    residue = []
    usedatoms = []
    for line in pdbarr:
        atom = line[12:17].strip()
        if not residue or line[17:27] != residue[-1][17:27] or atom in usedatoms:
            if residue: nestedpdb.append(residue)
            residue = [line]
            usedatoms = [atom]
        else:
            residue.append(line)
            usedatoms.append(atom)
    if residue:
        nestedpdb.append(residue)
    return nestedpdb

def unnest_pdb(npdb):
    return [atm for res in npdb for atm in res]

def get_coords(atomname, residue):
    for line in residue:
        if line[12:16].strip() == atomname.strip():
            return tuple(float(line[i:i+8]) for i in range(30, 54, 8))
    raise ValueError("Atom not found!")

def disulfide_search(npdb, min_dist=1.8, max_dist=2.2):
    residues_to_rename = set()
    for i in range(len(npdb)):
        if npdb[i][0][17:20] not in ['CYS', 'CYD', 'CYX']:
            continue
        iX, iY, iZ = get_coords('SG', npdb[i])
        i_residue_info = npdb[i][0][17:27].strip()  # Extract residue info for logging
        i_atom_number = npdb[i][0][6:11].strip()  # Extract atom number for logging

        for j in range(i + 1, len(npdb)):
            if npdb[j][0][17:20] not in ['CYS', 'CYD', 'CYX']:
                continue
            jX, jY, jZ = get_coords('SG', npdb[j])
            j_residue_info = npdb[j][0][17:27].strip()  # Extract residue info for logging
            j_atom_number = npdb[j][0][6:11].strip()  # Extract atom number for logging

            distance = math.sqrt((iX - jX) ** 2 + (iY - jY) ** 2 + (iZ - jZ) ** 2)
            if min_dist <= distance <= max_dist:
                residues_to_rename.update({i, j})
                # Log the atoms involved in the disulfide bond, including their atom numbers
                logger.info(f"Disulfide bond detected within atoms: {i_atom_number}_{j_atom_number} with distance {distance:.2f} Å.")
                logger.info(f"Bond between residues `{i_residue_info}` and `{j_residue_info}`.")
    
    for i in residues_to_rename:
        npdb[i] = [x.replace('CYS ', 'CYX') if 'CYS' in x or 'CYD' in x else x for x in npdb[i]]
    
    return npdb


def pdb_cleanup(pdbarr):
    updated_pdbarr = []
    for line in pdbarr:
        atom = line[:54]  # Truncate each line to the first 54 characters (up to the end of the coordinates)
        atomsymbol = line[12:16].strip(" 0123456789")[0]  # Extract the atom symbol
        updated_line = f"{atom}  1.00  0.00          {atomsymbol.rjust(2)}"  # Format the line with default occupancy and B-factor, and reposition the atom symbol
        updated_pdbarr.append(updated_line)

    npdb = nest_pdb(updated_pdbarr)  # Nest the PDB for further processing if needed
    
    # Renumber residues and set chain identifiers if required
    for i, residue in enumerate(npdb):
        for j, line in enumerate(residue):
            # Assume chain identifier is blank for simplicity; adjust if handling multi-chain PDBs
            npdb[i][j] = f"{line[:21]} {' ':1}{str(i + 1).rjust(4)}{'    '}{line[30:]}"

    return unnest_pdb(npdb)  # Return the unnested, cleaned-up PDB array

def histidine_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:20]
        if resname == 'HIS':
            HE_present = atom_is_present(npdb[i], 'HE2')  # bonded to NE2
            HD_present = atom_is_present(npdb[i], 'HD1')  # bonded to ND1

            if HD_present and HE_present:
                npdb[i] = [x.replace('HIS', 'HIP') for x in npdb[i]]
            elif HE_present:
                npdb[i] = [x.replace('HIS', 'HIE') for x in npdb[i]]
            elif HD_present:
                npdb[i] = [x.replace('HIS', 'HID') for x in npdb[i]]
            else:
                raise ValueError("No protons found for histidine.")
    return npdb

def atom_is_present(pdblines, atomname):
    """Returns TRUE if the given atom is present in the given PDB atom lines.
    
    ARGUMENTS
        pdblines - list of PDB lines
        atomname - the name of the atom to check the existence of

    RETURNS
        is_present - True if the given atom name is present, False otherwise

    """
    atoms = [pdbline[13:16] for pdbline in pdblines]
    if any(atomname in atom for atom in atoms):
        return True
    else:
        return False


def main_exe():
    parser = argparse.ArgumentParser(description="Rename amino acids in a PDB file for AMBER forcefield compatibility.")
    parser.add_argument('-i', '--input', required=True, help="Input PDB file to process")
    parser.add_argument('-o', '--output', required=False, help="Output name for the processed file", default=None)
    args = parser.parse_args()

    input_pdb_path = args.input
    if not os.path.exists(input_pdb_path):
        logger.error(f"Input file {input_pdb_path} does not exist.")
        sys.exit(1)

    logger.info(f"Processing file: {input_pdb_path}")

    try:
        with open(input_pdb_path, 'r') as f:
            pdb_lines = f.readlines()
        
        renamed_pdb_lines = rename_residues(pdb_lines)
        if args.output is not None:
            output_pdb_path = args.output
        else:
            output_pdb_path = input_pdb_path.replace('.pdb', '_renamed.pdb')

        with open(output_pdb_path, 'w') as f:
            for line in renamed_pdb_lines:
                f.write(line)

        logger.info(f"Renaming completed. Output file: {output_pdb_path}")

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main_exe()
