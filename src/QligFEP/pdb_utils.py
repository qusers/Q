"""Module containing functions for parsing pdb files."""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.neighbors import NearestNeighbors

from .logger import logger


def pdb_HOH_nn(pdb_df_query, pdb_df_target, th=2.5, output_file=None):
    """
    Find water oxygen atoms within a distance threshold from protein atoms.

    Args:
        pdb_df_query: DataFrame containing the water molecules.
        pdb_df_target: DataFrame containing the protein atoms.
        th: Distance threshold in Angstroms.
        output_file: Optional path to write the result to a file.

    Returns:
        A DataFrame containing the water oxygen atoms within the distance threshold.
    """
    # Extract coordinates and filter to keep only oxygen atoms from water molecules
    water_oxygen_df = pdb_df_query.query("atom_name == 'O'")
    query_arr = water_oxygen_df[["x", "y", "z"]].values

    # Select only oxygen atoms from the target protein atoms
    protein_oxygen_df = pdb_df_target.query("atom_name == 'O'")
    target_arr = protein_oxygen_df[["x", "y", "z"]].values

    # Use NearestNeighbors with a radius threshold
    knn = NearestNeighbors(radius=th, metric="euclidean", n_jobs=-1)
    knn.fit(query_arr)

    # Find all protein atoms within the radius of water oxygen atoms
    distances, indices = knn.radius_neighbors(target_arr)

    # Flatten indices, ensuring only unique and valid indices are retained
    unique_indices = sorted(set([i for sublist in indices for i in sublist]))
    to_rm_waters = water_oxygen_df.iloc[unique_indices]["residue_seq_number"].tolist()
    logger.info(f"Removing {len(to_rm_waters)} water molecules within {th} Å of protein atoms.")
    final = pdb_df_query[~pdb_df_query["residue_seq_number"].isin(to_rm_waters)].copy()
    # now renumber the atom_serial_number and the residue_seq_number
    startAtom = final["atom_serial_number"].values[0]

    final["atom_serial_number"] = np.arange(startAtom, startAtom + len(final))
    new_residue_seq_numbers = {
        residue: i + 1 for i, residue in enumerate(final["residue_seq_number"].unique())
    }
    final["residue_seq_number"] = final["residue_seq_number"].map(new_residue_seq_numbers)

    # Optionally write results to a PDB file
    if output_file:
        write_dataframe_to_pdb(final, output_file)

    return final


def pdb_parse_in(line, include=("ATOM", "HETATM")):
    """
    Takes a pdb file line and parses it into a list, according to Atomic Coordinate Entry Format
    v3.3
    """
    at_entry = []
    line = line.strip("\n")
    if line.startswith(include):
        at_entry.append(line[0:6])  #  0 ATOM/HETATM
        at_entry.append(int(line[6:11]))  #  1 ATOM serial number
        at_entry.append(line[12:16].strip())  #  2 ATOM name
        at_entry.append(line[16:17])  #  3 Alternate location indicator
        at_entry.append(line[17:21].strip())  #  4 Residue name
        at_entry.append(line[21:22])  #  5 Chain identifier
        at_entry.append(int(line[22:26]))  #  6 Residue sequence number
        at_entry.append(line[26:27])  #  7 Code for insertion of residue
        at_entry.append(float(line[30:38]))  #  8 Orthogonal coordinates for X
        at_entry.append(float(line[38:46]))  #  9 Orthogonal coordinates for Y
        at_entry.append(float(line[46:54]))  # 10 Orthogonal coordinates for Z
        # These entries can be empty
        try:
            at_entry.append(float(line[54:60]))  # 11 Occupancy
        except IndexError:
            at_entry.append(0.0)  # 11 Empty Occupancy
        try:
            at_entry.append(float(line[60:66]))  # 12 Temperature factor
        except IndexError:
            at_entry.append(0.0)  # 12 Empty Temperature factor
        try:
            at_entry.append(line[76:78])  # 13 Element symbol
        except IndexError:
            at_entry.append("  ")  # 13 Empty Element symbol
        try:
            at_entry.append(line[78:80])  # 14 Charge on atom
        except IndexError:
            at_entry.append("  ")  # 14 Empty charge
    else:
        at_entry = line

    return at_entry


def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    line = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        *line
    )
    return line


def nest_pdb(pdbarr: list[str]) -> list[list[str]]:
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
            if residue:
                nestedpdb.append(residue)
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
    cysbonds = []
    for i in range(len(npdb)):
        if npdb[i][0][17:20] not in ["CYS", "CYD", "CYX"]:
            continue
        iX, iY, iZ = get_coords("SG", npdb[i])
        sg_idx = np.where(np.char.find(npdb[i], "SG") != -1)[0][0]
        i_residue_info = npdb[i][sg_idx][17:27].strip()  # Extract residue info for logging
        i_atom_number = npdb[i][sg_idx][6:11].strip()  # Extract atom number for logging

        for j in range(i + 1, len(npdb)):
            if npdb[j][0][17:20] not in ["CYS", "CYD", "CYX"]:
                continue
            jX, jY, jZ = get_coords("SG", npdb[j])
            sg_idx = np.where(np.char.find(npdb[j], "SG") != -1)[0][0]
            j_residue_info = npdb[j][sg_idx][17:27].strip()  # Extract residue info for logging
            j_atom_number = npdb[j][sg_idx][6:11].strip()  # Extract atom number for logging

            distance = math.sqrt((iX - jX) ** 2 + (iY - jY) ** 2 + (iZ - jZ) ** 2)
            if min_dist <= distance <= max_dist:
                residues_to_rename.update({i, j})
                # Log the atoms involved in the disulfide bond, including their atom numbers
                logger.info(
                    f"Disulfide bond detected within atoms: {i_atom_number}_{j_atom_number} with distance {distance:.2f} Å."
                )
                logger.info(f"Bond between residues `{i_residue_info}` and `{j_residue_info}`.")
                cysbonds.append((i_atom_number, j_atom_number))

    for i in residues_to_rename:
        npdb[i] = [x.replace("CYS ", "CYX") if "CYS" in x or "CYD" in x else x for x in npdb[i]]

    return npdb, cysbonds


def get_coords(atomname, residue):
    for line in residue:
        if line[12:16].strip() == atomname.strip():
            return tuple(float(line[i : i + 8]) for i in range(30, 54, 8))
    raise ValueError("Atom not found!")


def _convert_to(value, dtype):
    try:
        return dtype(value)
    except ValueError:
        return value


def _parse_pdb_line(line):
    if line.startswith(("ATOM", "HETATM")):
        parsed_line = [
            line[0:6].strip(),  # record_type
            _convert_to((line[6:11].strip()), int),  # atom_serial_number
            line[12:16].strip(),  # atom_name
            line[16].strip(),  # alt_loc
            line[17:20].strip(),  # residue_name
            line[21].strip(),  # chain_id
            _convert_to((line[22:26].strip()), int),  # residue_seq_number
            line[26].strip(),  # insertion_code
            _convert_to((line[30:38].strip()), float),  # x
            _convert_to((line[38:46].strip()), float),  # y
            _convert_to((line[46:54].strip()), float),  # z
            _convert_to((line[54:60].strip()), float),  # occupancy
            _convert_to((line[60:66].strip()), float),  # temp_factor
            line[72:76].strip(),  # segment_id
            line[76:78].strip(),  # element_symbol
            line[78:80].strip(),  # charge
        ]
        return parsed_line


def read_pdb_to_dataframe(pdb_file):
    columns = [
        "record_type",
        "atom_serial_number",
        "atom_name",
        "alt_loc",
        "residue_name",
        "chain_id",
        "residue_seq_number",
        "insertion_code",
        "x",
        "y",
        "z",
        "occupancy",
        "temp_factor",
        "segment_id",
        "element_symbol",
        "charge",
    ]
    data = []
    if isinstance(pdb_file, list):
        for line in pdb_file:
            result = _parse_pdb_line(line)
            if result is not None:
                data.append(result)
    elif isinstance(pdb_file, (str, Path)):
        assert Path(pdb_file).exists(), f"File {pdb_file} does not exist."
        with open(pdb_file) as file:
            for line in file:
                result = _parse_pdb_line(line)
                if result is not None:
                    data.append(result)

    df = pd.DataFrame(data, columns=columns)
    return df


def write_dataframe_to_pdb(df, output_file):
    with open(output_file, "w") as file:
        for _, row in df.iterrows():
            pdb_line = (
                f"{row['record_type']:<6}{row['atom_serial_number']:>5} "
                f"{row['atom_name']:<4}{row['alt_loc']:<1}{row['residue_name']:>3} "
                f"{row['chain_id']:>1}{row['residue_seq_number']:>4}{row['insertion_code']:>1}   "
                f"{row['x']:>8.3f}{row['y']:>8.3f}{row['z']:>8.3f}{row['occupancy']:>6.2f}"
                f"{row['temp_factor']:>6.2f}          {row['element_symbol']:>2}{row['charge']:>2}\n"
            )
            file.write(pdb_line)


def sdf_to_pdb(in_sdf_file, out_pdb_file):
    """Converts an SDF file to a PDB file.

    Args:
        in_sdf_file: std_in; path with the sdf file
        out_pdb_file: std_out; path for the pdb file
    """
    suppl = Chem.SDMolSupplier(in_sdf_file)
    for mol in suppl:
        if mol is not None:
            # Generate 3D coordinates if not present
            mol_with_h = Chem.AddHs(mol, addCoords=True)
            AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=200)
            with open(out_pdb_file, "w") as f:
                print("overwriting")
                f.write(Chem.MolToPDBBlock(mol_with_h))
            break
