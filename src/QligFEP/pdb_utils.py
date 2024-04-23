"""Module containing functions for parsing pdb files."""

import math
import numpy as np
from typing import List
from .logger import logger


def pdb_parse_in(line, include=('ATOM','HETATM')):
    """
    Takes a pdb file line and parses it into a list, according to Atomic Coordinate Entry Format 
    v3.3
    """
    at_entry = []
    line = line.strip('\n')
    if line.startswith(include):
        at_entry.append(line[0:6])              #  0 ATOM/HETATM
        at_entry.append(int(line[6:11]))        #  1 ATOM serial number
        at_entry.append(line[12:16].strip())    #  2 ATOM name
        at_entry.append(line[16:17])            #  3 Alternate location indicator
        at_entry.append(line[17:21].strip())    #  4 Residue name
        at_entry.append(line[21:22])            #  5 Chain identifier
        at_entry.append(int(line[22:26]))       #  6 Residue sequence number
        at_entry.append(line[26:27])            #  7 Code for insertion of residue
        at_entry.append(float(line[30:38]))     #  8 Orthogonal coordinates for X
        at_entry.append(float(line[38:46]))     #  9 Orthogonal coordinates for Y
        at_entry.append(float(line[46:54]))     # 10 Orthogonal coordinates for Z
        # These entries can be empty
        try:
            at_entry.append(float(line[54:60])) # 11 Occupancy
            
        except IndexError:
            at_entry.append(0.0)                # 11 Empty Occupancy
            
        try:
            at_entry.append(float(line[60:66])) # 12 Temperature factor
            
        except IndexError:
            at_entry.append(0.0)                # 12 Empty Temperature factor
            
        try:
            at_entry.append(line[76:78])        # 13 Element symbol
            
        except IndexError:
            at_entry.append('  ')               # 13 Empty Element symbol
            
        try:
            at_entry.append(line[78:80])        # 14 Charge on atom
            
        except IndexError:
            at_entry.append('  ')               # 14 Empty charge
        
    else:
        at_entry = line
    
    return at_entry

def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    line = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(*line)
    return line

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

def disulfide_search(npdb, min_dist=1.8, max_dist=2.2):
    residues_to_rename = set()
    for i in range(len(npdb)):
        if npdb[i][0][17:20] not in ['CYS', 'CYD', 'CYX']:
            continue
        iX, iY, iZ = get_coords('SG', npdb[i])
        sg_idx = np.where(np.char.find(npdb[i], 'SG') != -1)[0][0]
        i_residue_info = npdb[i][sg_idx][17:27].strip()  # Extract residue info for logging
        i_atom_number = npdb[i][sg_idx][6:11].strip()  # Extract atom number for logging

        for j in range(i + 1, len(npdb)):
            if npdb[j][0][17:20] not in ['CYS', 'CYD', 'CYX']:
                continue
            jX, jY, jZ = get_coords('SG', npdb[j])
            sg_idx = np.where(np.char.find(npdb[j], 'SG') != -1)[0][0]
            j_residue_info = npdb[j][sg_idx][17:27].strip()  # Extract residue info for logging
            j_atom_number = npdb[j][sg_idx][6:11].strip()  # Extract atom number for logging

            distance = math.sqrt((iX - jX) ** 2 + (iY - jY) ** 2 + (iZ - jZ) ** 2)
            if min_dist <= distance <= max_dist:
                residues_to_rename.update({i, j})
                # Log the atoms involved in the disulfide bond, including their atom numbers
                logger.info(f"Disulfide bond detected within atoms: {i_atom_number}_{j_atom_number} with distance {distance:.2f} Å.")
                logger.info(f"Bond between residues `{i_residue_info}` and `{j_residue_info}`.")
    
    for i in residues_to_rename:
        npdb[i] = [x.replace('CYS ', 'CYX') if 'CYS' in x or 'CYD' in x else x for x in npdb[i]]
    
    return npdb

def get_coords(atomname, residue):
    for line in residue:
        if line[12:16].strip() == atomname.strip():
            return tuple(float(line[i:i+8]) for i in range(30, 54, 8))
    raise ValueError("Atom not found!")