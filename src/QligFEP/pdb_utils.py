"""Module containing functions for parsing pdb files."""

import math
import re
from pathlib import Path
from string import ascii_uppercase
from typing import Optional, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.neighbors import NearestNeighbors

from .logger import logger


def rm_HOH_clash_NN(
    pdb_df_query: pd.DataFrame,
    pdb_df_target: pd.DataFrame,
    th: float = 2.5,
    output_file: Union[str, Path] = None,
    heavy_only: bool = True,
    ligand_only: bool = False,
    header: Optional[str] = None,
    save_removed: bool = False,
):
    """Use a NearestNeighbors approach to find water molecules within a distance threshold
    (Ångström) from an input pdb file (e.g.: a protein-ligand complex), and remove them if
    the atoms are within the threshold distance from the input pdb structure.

    Args:
        pdb_df_query: DataFrame containing the water molecules to be removed.
        pdb_df_target: DataFrame containing the protein atoms.
        th: Distance threshold in Angstroms.
        output_file: Optional path to write the result to a file.
        heavy_only: If True, only consider heavy atoms (i.e., exclude hydrogen atoms).
        ligand_only: If True, only remove waters near ligand atoms (default as LIG and LID).
        header: header to be added to the output pdb file. Defaults to None

    Returns:
        A DataFrame containing the water oxygen atoms within the distance threshold.
    """
    water_query = pdb_df_query.query("atom_name == 'O'") if heavy_only else pdb_df_query
    query_arr = water_query[["x", "y", "z"]].values

    # we ignore ions in the target; since water molecules might be in close proximity to it
    pdb_ions = ["ZN", "SOD", "IOD", "BR", "CL", "CU", "CU1", "NA", "MG", "CA"]  # noqa: F841
    target_df = pdb_df_target.query("~residue_name.isin(@pdb_ions)")
    if heavy_only:
        Hatom_protein_regex = re.compile(r"(?<![NCO])H\d*")  # noqa: F841
        Hatom_ligand_regex = r"^H[A-Z]?\d{0,2}?"  # noqa: F841
        if ligand_only:
            target_arr = target_df.query(
                "(~atom_name.str.match(@Hatom_protein_regex)) & "
                r"(residue_name.isin(['LIG', 'LID']) & (~atom_name.str.contains(@Hatom_ligand_regex)))"
            )[["x", "y", "z"]].values
        else:
            target_arr = target_df.query(
                "(~atom_name.str.match(@Hatom_protein_regex)) | "
                r"(residue_name.isin(['LIG', 'LID']) & (~atom_name.str.contains(@Hatom_ligand_regex)))"
            )[["x", "y", "z"]].values
    else:
        if ligand_only:
            target_arr = target_df.query(r"residue_name.isin(['LIG', 'LID'])")[["x", "y", "z"]].values
        else:
            target_arr = target_df[["x", "y", "z"]].values

    boron_atoms = pdb_df_target.query(  # stricter check for boron-water proximity (crashes in QligFEP)
        r"residue_name.isin(['LIG', 'LID']) & atom_name.str.contains('^B\d{0,2}?')"
    )

    knn = NearestNeighbors(radius=th, metric="euclidean", n_jobs=4)
    knn.fit(query_arr)
    distances, indices = knn.radius_neighbors(target_arr)
    unique_indices = sorted(set([i for sublist in indices for i in sublist]))

    if not boron_atoms.empty:
        boron_th = th + 1
        boron_knn = NearestNeighbors(radius=boron_th, metric="euclidean", n_jobs=4)
        boron_knn.fit(query_arr)
        boron_distances, boron_indices = boron_knn.radius_neighbors(boron_atoms[["x", "y", "z"]].values)
        boron_unique_indices = sorted(set([i for sublist in boron_indices for i in sublist]))
        n_Br_removed = len(water_query.iloc[boron_unique_indices]["residue_seq_number"].tolist())
        if n_Br_removed > 0:
            logger.info(f"Removing {n_Br_removed} waters near Boron atoms - threshold: {boron_th} Å")
            unique_indices = sorted(set(unique_indices + boron_unique_indices))

    to_rm_waters = water_query.iloc[unique_indices]["residue_seq_number"].tolist()
    final = pdb_df_query[~pdb_df_query["residue_seq_number"].isin(to_rm_waters)].copy()
    # now renumber the atom_serial_number and the residue_seq_number
    startAtom = final["atom_serial_number"].values[0]

    n_removed = np.setdiff1d(
        water_query["residue_seq_number"].unique(), final["residue_seq_number"].unique()
    ).shape[0]
    logger.info(
        f"Removed {n_removed} (total) water molecules {'with oxygen atoms' if heavy_only else ''}"
        f" within {th} Å of protein atoms."
    )

    final["atom_serial_number"] = np.arange(startAtom, startAtom + len(final))
    new_residue_seq_numbers = {
        residue: i + 1 for i, residue in enumerate(final["residue_seq_number"].unique())
    }
    final["residue_seq_number"] = final["residue_seq_number"].map(new_residue_seq_numbers)

    # Optionally write results to a PDB file
    if output_file:
        # assign the correct data types
        final["atom_serial_number"] = final["atom_serial_number"].astype(int)
        final["residue_seq_number"] = final["residue_seq_number"].astype(int)
        final["x"] = final["x"].astype(float)
        final["y"] = final["y"].astype(float)
        final["z"] = final["z"].astype(float)
        try:
            final["occupancy"] = final["occupancy"].astype(float)
            final["temp_factor"] = final["temp_factor"].astype(float)
        except ValueError:
            final["occupancy"] = 0
            final["temp_factor"] = 0
        write_dataframe_to_pdb(final, output_file, header=header)
        if save_removed:
            removed_file = Path(output_file).with_name(Path(output_file).stem + "_removed.pdb")
            write_dataframe_to_pdb(
                pdb_df_query[pdb_df_query["residue_seq_number"].isin(to_rm_waters)],
                removed_file,
                header=header,
            )
    return final, n_removed


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
        at_entry.append(line[17:21].strip())  #  4 Residue name - 21 instead of 20 for N- & C- termini
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


def next_chain_id(existing_ids):
    """
    Calculate the next chain ID based on existing IDs.
    Wrap around to 'A' after 'Z', and ensure uniqueness.
    """
    alphabet = list(ascii_uppercase)
    if not existing_ids:
        return "A"
    # Find the highest current chain_id and increment
    highest_id = max([alphabet.index(cid) for cid in existing_ids if cid in alphabet], default=-1)
    next_id_index = (highest_id + 1) % len(alphabet)
    return alphabet[next_id_index]


def append_pdb_to_another(
    main_pdb: Union[pd.DataFrame, str, list[str]],
    to_append_pdb: Union[pd.DataFrame, str, list[str]],
    save_pdb: Optional[str] = None,
    assign_new_chain: bool = False,
    new_ligname: Optional[str] = None,
    ignore_waters: bool = False,
) -> pd.DataFrame:
    """Reads the two pdbs as DataFrames, appends the second to the end of the protein
    file and writes the new pdb file containing both.

    Args:
        main_pdb: main pdb file to receive the appended pdb.
        to_append_pdb: input pdb to be appended to main.
        save_pdb: if desired, the path to save the merged pdb file. Defaults to None.
        assign_new_chain: if True, assigns a new chain ID to the appended part. Defaults to False.
        new_ligname: new residue name for the appended part, if desired. Defaults to None.
        ignore_waters: don't take waters in consideration for atom & residue number assignment.
            Defaults to False.

    Returns:
        DataFrame with the merged structure.
    """
    main_df = read_pdb_to_dataframe(main_pdb) if isinstance(main_pdb, (str, Path)) else main_pdb
    to_append_df = (
        read_pdb_to_dataframe(to_append_pdb) if isinstance(to_append_pdb, (str, Path)) else to_append_pdb
    )

    if assign_new_chain:
        existing_chain_ids = set(main_df["chain_id"].replace({"": np.nan}).dropna().unique())
        if not existing_chain_ids:
            main_df["chain_id"] = "A"  # Default the protein to chain A if no chain_id is present
            new_chain_id = "B"
        else:
            new_chain_id = next_chain_id(existing_chain_ids)

    if ignore_waters:
        last_prot_atom = main_df.query("residue_name != 'HOH'")["atom_serial_number"].astype(int).max()
        last_prot_resn = main_df.query("residue_name != 'HOH'")["residue_seq_number"].astype(int).max()
    else:
        last_prot_atom = main_df["atom_serial_number"].astype(int).max()
        last_prot_resn = main_df["residue_seq_number"].astype(int).max()

    to_append_df = to_append_df.assign(
        atom_serial_number=(to_append_df["atom_serial_number"].astype(int) + last_prot_atom).astype(str),
        residue_seq_number=(to_append_df["residue_seq_number"].astype(int) + last_prot_resn).astype(str),
    )

    if new_ligname is not None:
        to_append_df["residue_name"] = new_ligname
    if assign_new_chain:
        to_append_df["chain_id"] = new_chain_id

    merged_df = pd.concat([main_df, to_append_df], ignore_index=True)
    if save_pdb is not None:
        write_dataframe_to_pdb(merged_df, save_pdb)
    return merged_df


def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line using positional arguments
    """
    line = "{:6s}{:5d} {:<4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
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
                logger.debug(
                    f"Disulfide bond detected within atoms: {i_atom_number}_{j_atom_number} with distance {distance:.2f} Å."
                )
                logger.debug(f"Bond between residues `{i_residue_info}` and `{j_residue_info}`.")
                cysbonds.append((f"{i_residue_info.split()[-1]}:SG", f"{j_residue_info.split()[-1]}:SG"))

    renamed = bool(residues_to_rename)
    for i in residues_to_rename:
        npdb[i] = [x.replace("CYS", "CYX") if "CYS" in x or "CYD" in x else x for x in npdb[i]]

    return npdb, cysbonds, renamed


def get_coords(atomname, residue):
    for line in residue:
        if line[12:16].strip() == atomname.strip():
            return tuple(float(line[i : i + 8]) for i in range(30, 54, 8))
    raise ValueError(f"Atom {atomname} not found in residue {residue}!")


def calculate_distance(atom_coords, center_coords) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom_coords, center_coords)))


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
            line[17:21].strip(),  # residue_name - 21 instead of 20 for N- & C- termini
            line[21:22].strip(),  # chain_id
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


def write_dataframe_to_pdb(df, output_file, header: Optional[str] = None):
    """Save a DataFrame object created from read_pdb_to_dataframe function to a PDB file.

    Args:
        df: DataFrame object containing the parsed PDB file.
        output_file: name of the output file (include .pdb extension).
        header: if desired, a header to be added to the PDB file. Defaults to None.
    """
    with open(output_file, "w") as file:
        if header is not None:
            file.write(f"{header}\n")
        if df.temp_factor.dtype == "float64":
            df.loc[:, "temp_factor"] = df.temp_factor.apply(lambda x: f"{x:.2f}")
        if df.occupancy.dtype == "float64":
            df.loc[:, "occupancy"] = df.occupancy.apply(lambda x: f"{x:.2f}")
        for _, row in df.iterrows():
            pdb_line = (
                f"{row['record_type']:<6}{row['atom_serial_number']:>5} "
                f"{row['atom_name']:<4}{row['alt_loc']:<1}{row['residue_name']:<4}"  # residue_name:>4 for N- & C- termini
                f"{row['chain_id']:>1}{row['residue_seq_number']:>4}{row['insertion_code']:>1}   "
                f"{row['x']:>8.3f}{row['y']:>8.3f}{row['z']:>8.3f}{row['occupancy']:>6}"
                f"{row['temp_factor']:>6}          {row['element_symbol']:>2}{row['charge']:>2}\n"
            )
            # except ValueError as e:
            #     logger.error(f"{row.values}")
            #     raise ValueError from e
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
