"""Module containing functions for parsing pdb files."""

import math
import re
import warnings
from pathlib import Path
from string import ascii_uppercase
from typing import Optional, Union

import MDAnalysis as mda
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


def disulfide_search(npdb, min_dist=1.8, max_dist_cyx=4.0, max_dist_cys=2.5):
    """Detect disulfide bonds using distance-dependent criteria.

    Uses two detection ranges:
    - CYX/CYD residues (already labeled as disulfide): wide range (1.8-4.0Å) to catch
      cases where the disulfide bond was poorly resolved in the structure, leading to
      long SG-SG distances.
    - CYS residues (free thiols): strict range (1.8-2.5Å) to avoid false positives
      from nearby cysteines that happen to be geometrically close.

    Reference for extended range (up to 3.2Å): https://pubs.acs.org/doi/10.1021/jz900214e
    """
    residues_to_rename = set()
    cysbonds = []

    # Collect all cysteine residues with their SG coordinates
    cys_residues = []
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].strip()
        if resname not in ("CYS", "CYD", "CYX"):
            continue
        try:
            x, y, z = get_coords("SG", npdb[i])
        except ValueError:
            continue
        sg_idx = np.where(np.char.find(npdb[i], "SG") != -1)[0][0]
        res_info = npdb[i][sg_idx][17:27].strip()
        atom_number = npdb[i][sg_idx][6:11].strip()
        cys_residues.append(
            {
                "idx": i,
                "resname": resname,
                "coords": (x, y, z),
                "res_info": res_info,
                "atom_number": atom_number,
            }
        )

    # Find disulfide pairs with appropriate distance cutoffs
    for ii, res_i in enumerate(cys_residues):
        for res_j in cys_residues[ii + 1 :]:
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(res_i["coords"], res_j["coords"])))

            # Use wide range only when both residues are CYX/CYD (confirmed disulfide
            # partners). Mixed CYX+CYS pairs use the strict cutoff to avoid false
            # positives from free cysteines that happen to be near a disulfide.
            either_is_cyx = res_i["resname"] in ("CYX", "CYD") or res_j["resname"] in ("CYX", "CYD")
            both_are_cyx = res_i["resname"] in ("CYX", "CYD") and res_j["resname"] in ("CYX", "CYD")
            max_dist = max_dist_cyx if both_are_cyx else max_dist_cys

            if min_dist <= distance <= max_dist:
                residues_to_rename.update({res_i["idx"], res_j["idx"]})
                logger.debug(
                    f"Disulfide bond detected: {res_i['atom_number']}_{res_j['atom_number']} "
                    f"({distance:.2f} Å, {'CYX-labeled' if either_is_cyx else 'distance-based'})"
                )
                logger.debug(f"Bond between residues `{res_i['res_info']}` and `{res_j['res_info']}`.")
                cysbonds.append(
                    (f"{res_i['res_info'].split()[-1]}:SG", f"{res_j['res_info'].split()[-1]}:SG")
                )

    renamed = bool(residues_to_rename)
    for i in residues_to_rename:
        npdb[i] = [x.replace("CYS", "CYX") if "CYS" in x or "CYD" in x else x for x in npdb[i]]
        # Remove HG atom — CYX (disulfide) has no thiol hydrogen
        npdb[i] = [x for x in npdb[i] if " HG " not in x]

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


def write_dataframe_to_pdb(
    df, output_file, header: Optional[str] = None, ter_after_indices: Optional[set[int]] = None
):
    """Save a DataFrame object created from read_pdb_to_dataframe function to a PDB file.

    Automatically inserts TER records at chain_id transitions. Additional TER positions
    (e.g., gaps from residue cropping) can be specified via ter_after_indices.

    Args:
        df: DataFrame object containing the parsed PDB file.
        output_file: name of the output file (include .pdb extension).
        header: if desired, a header to be added to the PDB file. Defaults to None.
        ter_after_indices: extra DataFrame indices after which to write a TER record,
            in addition to auto-detected chain breaks. Defaults to None.
    """
    ter_positions = detect_chain_breaks(df)
    if ter_after_indices is not None:
        ter_positions |= ter_after_indices

    df = df.copy()
    with open(output_file, "w") as file:
        if header is not None:
            file.write(f"{header}\n")
        if df["temp_factor"].dtype == "float64":
            df["temp_factor"] = df["temp_factor"].apply(lambda x: f"{x:.2f}")
        if df["occupancy"].dtype == "float64":
            df["occupancy"] = df["occupancy"].apply(lambda x: f"{x:.2f}")
        for idx, row in df.iterrows():
            pdb_line = (
                f"{row['record_type']:<6}{row['atom_serial_number']:>5} "
                f"{row['atom_name']:<4}{row['alt_loc']:<1}{row['residue_name']:<4}"  # residue_name:>4 for N- & C- termini
                f"{row['chain_id']:>1}{row['residue_seq_number']:>4}{row['insertion_code']:>1}   "
                f"{row['x']:>8.3f}{row['y']:>8.3f}{row['z']:>8.3f}{row['occupancy']:>6}"
                f"{row['temp_factor']:>6}          {row['element_symbol']:>2}{row['charge']:>2}\n"
            )
            file.write(pdb_line)
            if idx in ter_positions:
                file.write(
                    f"TER   {int(row['atom_serial_number']) + 1:>5}      "
                    f"{row['residue_name']:>3} {row['chain_id']:>1}{row['residue_seq_number']:>4}\n"
                )


def detect_chain_breaks(df: pd.DataFrame) -> set[int]:
    """Find DataFrame indices where TER records should be inserted.

    Detects chain boundaries by looking for chain_id changes between
    consecutive ATOM/HETATM records.

    Args:
        df: DataFrame from read_pdb_to_dataframe.

    Returns:
        Set of DataFrame indices corresponding to the last atom before each chain break.
    """
    atom_df = df[df["record_type"].isin(("ATOM", "HETATM"))]
    if atom_df.empty:
        return set()

    chain_ids = atom_df["chain_id"]
    changed = chain_ids != chain_ids.shift(1)
    # Indices where chain_id changes (skip first row — nothing before it)
    change_positions = atom_df.index[changed][1:]

    ter_indices = set()
    for pos in change_positions:
        loc = atom_df.index.get_loc(pos)
        prev_idx = atom_df.index[loc - 1]
        ter_indices.add(prev_idx)

    return ter_indices


def reindex_pdb_residues(pdb_path: Path, out_pdb_path: Path):
    pdb_df = read_pdb_to_dataframe(pdb_path)
    uniq_indexes = pdb_df.set_index(
        ["residue_seq_number", "residue_name", "chain_id", "insertion_code"]
    ).index
    resn_mapping = {resn: idx for idx, resn in enumerate(uniq_indexes.unique(), 1)}
    pdb_df["residue_seq_number"] = uniq_indexes.map(resn_mapping)
    pdb_df["insertion_code"] = ""
    write_dataframe_to_pdb(pdb_df, str(out_pdb_path.resolve().absolute()))


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


def filter_pdb_by_sphere(
    pdb_input: Path,
    pdb_output: Path,
    center: list[float],
    radius: float,
    exclude_residues: set[str] = None,
) -> tuple[int, int]:
    """Filter PDB keeping only fragments with atoms inside the sphere.

    Uses MDAnalysis topology-based fragment detection (bond connectivity)
    rather than relying on chain IDs. Entire molecular fragments are preserved
    if any atom is within the sphere radius.

    Args:
        pdb_input: Input PDB file path
        pdb_output: Output PDB file path
        center: [x, y, z] coordinates of sphere center
        radius: Sphere radius in Angstroms
        exclude_residues: Residue names to exclude from filtering (handled separately).
            Defaults to {"HOH", "LIG", "LID"}.

    Returns:
        Tuple of (original_atom_count, filtered_atom_count)
    """
    if exclude_residues is None:
        exclude_residues = {"HOH", "LIG", "LID"}

    # Custom vdW radii for ions that are common in PDB files but not in MDAnalysis defaults
    # Includes both standard element names and (Q) force field-specific naming conventions
    custom_vdwradii = {
        # Standard element names
        "Na": 2.27,  # Sodium
        "K": 2.75,  # Potassium
        "Cl": 1.75,  # Chloride
        "Ca": 2.31,  # Calcium
        "Mg": 1.73,  # Magnesium
        "Zn": 1.39,  # Zinc
        "Fe": 1.56,  # Iron
        "Cu": 1.40,  # Copper
        "Mn": 1.39,  # Manganese
        # AMBER force field naming conventions
        "SOD": 2.27,  # Sodium (AMBER)
        "MAG": 1.73,  # Magnesium (AMBER)
        "CHL": 1.75,  # Chloride (AMBER)
        "ZIN": 1.39,  # Zinc (AMBER)
        # OPLS force field naming conventions
        "MG2": 1.73,  # Magnesium (OPLS)
        "Na+": 2.27,  # Sodium (OPLS)
        # Non-standard element records with formal charges
        # These occur when formal charges are written in the element column (e.g., "O1-" parsed as "O1")
        "O1": 1.52,  # Oxygen (parsed from O1- element record)
        "N1": 1.55,  # Nitrogen (parsed from N1+ element record)
        # Note - all values come from MDAnalysis itself:
        # https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/guesser/tables.py
    }

    # Suppress expected MDAnalysis warnings (missing CRYST1, missing chain IDs, etc.)
    # These are harmless for our topology-based fragment detection
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unit cell dimensions not found")
        warnings.filterwarnings("ignore", message="Found missing chainIDs")
        warnings.filterwarnings("ignore", message="Found no information for attr")
        warnings.filterwarnings("ignore", message="Element information is missing")

        try:
            # Load Universe without guessing bonds initially
            u = mda.Universe(str(pdb_input))
            # Then guess bonds with custom vdW radii (MDAnalysis 3.0 compatible)
            u.guess_TopologyAttrs(to_guess=["bonds"], vdwradii=custom_vdwradii)
        except Exception as e:
            logger.warning(f"Could not load PDB with MDAnalysis: {e}")
            return 0, 0

        original_count = u.atoms.n_atoms

        if original_count == 0:
            # Empty PDB file
            Path(pdb_output).write_text("END\n")
            return 0, 0

        # Build selection string for the sphere volume
        cx, cy, cz = center
        volume_selection = f"point {cx} {cy} {cz} {radius}"

        # Build exclusion string for specified residues
        exclude_str = " or ".join([f"resname {r}" for r in exclude_residues]) if exclude_residues else None

        # Select: (fragments touching sphere) PLUS (excluded residues that touch sphere)
        # This ensures LIG/LID/HOH in sphere are kept, but their fragments don't pull in distant atoms
        # MDAnalysis's "same fragment as" expands selection to complete connected molecules
        if exclude_str:
            # For non-excluded residues: keep whole fragment if any atom in sphere
            # For excluded residues: only keep them if they're in sphere (no fragment expansion)
            non_excluded_in_sphere = f"(same fragment as (({volume_selection}) and not ({exclude_str})))"
            excluded_in_sphere = f"(({exclude_str}) and ({volume_selection}))"
            selection = f"{non_excluded_in_sphere} or {excluded_in_sphere}"
        else:
            selection = f"same fragment as ({volume_selection})"

        try:
            kept_atoms = u.select_atoms(selection)
        except Exception as e:
            logger.warning(f"MDAnalysis selection failed: {e}. Keeping all atoms.")
            kept_atoms = u.atoms

        filtered_count = kept_atoms.n_atoms

        if filtered_count == 0:
            # All atoms outside sphere
            Path(pdb_output).write_text("END\n")
            return original_count, 0

        # Write filtered structure
        kept_atoms.write(str(pdb_output))

    return original_count, filtered_count


def filter_out_of_sphere_fragments(
    pdb_path: Path,
    center: list[float],
    radius: float,
    exclude_residues: set[str] = None,
) -> tuple[int, int]:
    """Remove molecular fragments completely outside the sphere radius.

    Uses MDAnalysis topology-based fragment detection for robust filtering
    that doesn't depend on chain ID labels. Modifies the PDB file in place.

    Args:
        pdb_path: Path to PDB file (will be modified in place)
        center: [x, y, z] sphere center coordinates
        radius: Sphere radius in Angstroms
        exclude_residues: Residue names to exclude from fragment expansion.
            Defaults to {"HOH", "LIG", "LID"}.

    Returns:
        Tuple of (original_atom_count, filtered_atom_count)
    """
    if exclude_residues is None:
        exclude_residues = {"HOH", "LIG", "LID"}

    # Use a temporary output file, then replace the original
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        orig_count, filt_count = filter_pdb_by_sphere(pdb_path, tmp_path, center, radius, exclude_residues)

        if filt_count < orig_count:
            # Move filtered file to original location
            import shutil

            shutil.move(str(tmp_path), str(pdb_path))
        else:
            # No change needed, remove temp file
            tmp_path.unlink(missing_ok=True)

        return orig_count, filt_count
    except Exception as e:
        # Clean up on error
        tmp_path.unlink(missing_ok=True)
        raise e
