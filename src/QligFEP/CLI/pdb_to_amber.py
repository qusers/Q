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
import math
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
from .aa_atom_rename import rename_mapping

# Monatomic ions: PDB/source residue name -> AMBER14sb library name.
# In AMBER14sb the atom name equals the residue name for every monatomic ion.
ION_RENAME = {
    "MG": "MAG",
    "ZN": "ZIN",
    "NA": "SOD",
    "CL": "CHL",
    "CA": "CAL",
}

# Protein amino acids whose backbone nitrogen carries an amide proton. Proline
# (and hydroxyproline) are excluded because their backbone nitrogen has none.
# Used to scope backbone-proton normalization away from DNA, ions, and caps.
AMINO_ACID_RESIDUES = {
    "ALA", "ARG", "ARN", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX",
    "GLH", "GLN", "GLU", "GLY", "HID", "HIE", "HIP", "HIS", "ILE",
    "LEU", "LYN", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL",
} # fmt: skip


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
            "2HG", "1HG", "2HB", "1HB", "1HG1", "2HG1", "1HA", "2HA", "1HD", "2HD", "1HE", "2HE",
        ] # fmt: skip

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


def correct_amino_acid_atom_names(npdb_i, resname, rename_mapping=rename_mapping):
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


_PROTONATED_CARBOXYLATES = {
    "GLH": {
        "oxygens": ("OE1", "OE2"),
        "hydrogens": ("HE1", "HE2"),
        "target_hydrogen": "HE2",
    },
    "ASH": {
        "oxygens": ("OD1", "OD2"),
        "hydrogens": ("HD1", "HD2"),
        "target_hydrogen": "HD2",
    },
}
_MAX_XH_BOND_DISTANCE = 1.3
_MIN_PARENT_DISTANCE_DIFFERENCE = 0.2


def _atom_coordinates(line):
    """Return the XYZ coordinates from a PDB ATOM/HETATM line."""
    try:
        return tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
    except ValueError as exc:
        raise ValueError(f"Invalid coordinates in PDB line: {line.rstrip()}") from exc


def _rename_atoms(npdb_i, renames):
    """Rename atoms in one pass so exchanges such as OE1 <-> OE2 are safe."""
    renamed = []
    for line in npdb_i:
        atom_name = line[12:16].strip()
        new_name = renames.get(atom_name)
        renamed.append(extract_and_replace(line, atom_name, new_name) if new_name else line)
    return renamed


def normalize_protonated_carboxylate(npdb_i, resname):
    """Normalize protonated Asp/Glu atom identities to the AMBER14sb convention.

    AMBER14sb bonds GLH HE2 to OE2 and ASH HD2 to OD2. Some preparation
    programs use the opposite oxygen numbering, or rename only the acidic
    hydrogen. Geometry identifies the hydroxyl oxygen: when the proton is next
    to oxygen 1, exchange the two oxygen *names* and normalize the hydrogen name
    to the AMBER oxygen-2 convention. Coordinates are never moved.

    Residues without an explicit acidic hydrogen are left alone so qprep can add
    a missing hydrogen when intentionally neutralizing an out-of-sphere residue.
    """
    definition = _PROTONATED_CARBOXYLATES.get(resname)
    if definition is None:
        return npdb_i

    lines_by_name = {}
    for line in npdb_i:
        if line.startswith(("ATOM", "HETATM")):
            lines_by_name.setdefault(line[12:16].strip(), []).append(line)

    hydrogen_names = [name for name in definition["hydrogens"] if name in lines_by_name]
    if not hydrogen_names:
        return npdb_i
    if len(hydrogen_names) != 1 or len(lines_by_name[hydrogen_names[0]]) != 1:
        raise ValueError(f"{resname} must contain exactly one acidic hydrogen, found {hydrogen_names}")

    oxygen1, oxygen2 = definition["oxygens"]
    for oxygen in (oxygen1, oxygen2):
        if len(lines_by_name.get(oxygen, [])) != 1:
            raise ValueError(f"{resname} with an acidic hydrogen must contain exactly one {oxygen} atom")

    hydrogen_name = hydrogen_names[0]
    hydrogen_xyz = _atom_coordinates(lines_by_name[hydrogen_name][0])
    distances = {
        oxygen: math.dist(hydrogen_xyz, _atom_coordinates(lines_by_name[oxygen][0]))
        for oxygen in (oxygen1, oxygen2)
    }
    nearest_oxygen = min(distances, key=distances.get)
    nearest_distance = distances[nearest_oxygen]
    if nearest_distance > _MAX_XH_BOND_DISTANCE:
        raise ValueError(
            f"{resname} acidic hydrogen is not bonded to either side-chain oxygen "
            f"({oxygen1}={distances[oxygen1]:.3f} A, {oxygen2}={distances[oxygen2]:.3f} A)"
        )
    if abs(distances[oxygen1] - distances[oxygen2]) < _MIN_PARENT_DISTANCE_DIFFERENCE:
        raise ValueError(
            f"Ambiguous {resname} acidic-hydrogen assignment "
            f"({oxygen1}={distances[oxygen1]:.3f} A, {oxygen2}={distances[oxygen2]:.3f} A)"
        )

    renames = {hydrogen_name: definition["target_hydrogen"]}
    if nearest_oxygen == oxygen1:
        renames.update({oxygen1: oxygen2, oxygen2: oxygen1})
        residue_id = lines_by_name[hydrogen_name][0][21:27].strip()
        logger.info(f"{resname} {residue_id}: exchanged {oxygen1}/{oxygen2} to match AMBER14sb")
    return _rename_atoms(npdb_i, renames)


def normalize_neutral_arginine_geometry(npdb_i):
    """Exchange ARN NH1/NH2 names when their hydrogens follow the opposite convention."""
    resname = npdb_i[0][17:21].strip()
    if resname != "ARN":
        return npdb_i

    lines_by_name = {line[12:16].strip(): line for line in npdb_i if line.startswith(("ATOM", "HETATM"))}
    required = {"NH1", "NH2", "HH11", "HH12", "HH21"}
    if not required.issubset(lines_by_name):
        return npdb_i

    xyz = {name: _atom_coordinates(lines_by_name[name]) for name in required}
    expected = (("HH11", "NH1"), ("HH12", "NH1"), ("HH21", "NH2"))
    opposite = (("HH11", "NH2"), ("HH12", "NH2"), ("HH21", "NH1"))

    def assignment_is_bonded(assignment):
        return all(
            math.dist(xyz[hydrogen], xyz[nitrogen]) <= _MAX_XH_BOND_DISTANCE
            for hydrogen, nitrogen in assignment
        )

    if assignment_is_bonded(expected):
        return npdb_i
    if assignment_is_bonded(opposite):
        residue_id = lines_by_name["NH1"][21:27].strip()
        logger.info(f"ARN {residue_id}: exchanged NH1/NH2 to match AMBER14sb")
        return _rename_atoms(npdb_i, {"NH1": "NH2", "NH2": "NH1"})

    distances = ", ".join(
        f"{hydrogen}-{nitrogen}={math.dist(xyz[hydrogen], xyz[nitrogen]):.3f} A"
        for hydrogen, nitrogen in expected
    )
    raise ValueError(f"ARN terminal-hydrogen assignment is inconsistent with AMBER14sb ({distances})")


def _normalize_arn_before_nesting(pdb_lines):
    """Normalize ARN terminal nitrogens/hydrogens before duplicate names split a residue.

    The neutral AMBER14sb ARN convention calls the two-hydrogen nitrogen NH1
    (HH11/HH12) and the one-hydrogen nitrogen NH2 (HH21). Source files may use
    the opposite nitrogen numbering or duplicate HH21 names. Geometry and the
    proton count determine the identities without moving coordinates.
    """
    residue_indices = {}
    for index, line in enumerate(pdb_lines):
        if line.startswith(("ATOM", "HETATM")) and line[17:21].strip() == "ARN":
            residue_indices.setdefault(line[17:27], []).append(index)

    for residue_id, indices in residue_indices.items():
        nitrogen_indices = [index for index in indices if pdb_lines[index][12:16].strip() in {"NH1", "NH2"}]
        hydrogen_indices = [
            index for index in indices if pdb_lines[index][12:16].strip() in {"HH11", "HH12", "HH21", "HH22"}
        ]
        # qprep can build atoms omitted from an incomplete PDB. Only normalize geometry when the whole guanidinium
        # group is available; otherwise preserve the supplied names and let qprep fill in the missing atoms.
        if len(nitrogen_indices) != 2 or len(hydrogen_indices) != 3:
            continue

        attached = {index: [] for index in nitrogen_indices}
        for hydrogen_index in hydrogen_indices:
            hydrogen_xyz = _atom_coordinates(pdb_lines[hydrogen_index])
            distances = {
                nitrogen_index: math.dist(hydrogen_xyz, _atom_coordinates(pdb_lines[nitrogen_index]))
                for nitrogen_index in nitrogen_indices
            }
            parent = min(distances, key=distances.get)
            if distances[parent] > _MAX_XH_BOND_DISTANCE:
                raise ValueError(f"ARN {residue_id.strip()} terminal hydrogen is not bonded to NH1 or NH2")
            attached[parent].append(hydrogen_index)

        counts = sorted(len(hydrogens) for hydrogens in attached.values())
        if counts != [1, 2]:
            raise ValueError(
                f"ARN {residue_id.strip()} must have one singly and one doubly protonated terminal nitrogen"
            )

        two_h_nitrogen = next(index for index, hydrogens in attached.items() if len(hydrogens) == 2)
        one_h_nitrogen = next(index for index, hydrogens in attached.items() if len(hydrogens) == 1)
        renames_by_index = {
            two_h_nitrogen: "NH1",
            one_h_nitrogen: "NH2",
            attached[two_h_nitrogen][0]: "HH11",
            attached[two_h_nitrogen][1]: "HH12",
            attached[one_h_nitrogen][0]: "HH21",
        }
        for index, new_name in renames_by_index.items():
            old_name = pdb_lines[index][12:16].strip()
            pdb_lines[index] = extract_and_replace(pdb_lines[index], old_name, new_name)

    return pdb_lines


def _validate_unique_atom_names(pdb_lines):
    """Reject duplicate atom identities before they reach qprep."""
    seen = set()
    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        identity = (line[17:27], line[12:16].strip())
        if identity in seen:
            raise ValueError(
                f"Duplicate atom name {identity[1]} in residue {identity[0].strip()} after normalization"
            )
        seen.add(identity)
    return pdb_lines


def _filter_altloc(pdb_lines):
    """Collapse alternate conformations, keeping the highest-occupancy altloc.

    Structures with alternate locations (column 17 set to A/B/...) carry duplicate
    atoms for the same residue. qprep rejects these with "Too many atoms in residue".
    For each (residue, atom name) we keep the altloc with the highest occupancy
    (ties broken by altloc letter) and clear the altloc indicator on the survivor.
    Atoms with a blank altloc are always kept.
    """
    best = {}  # (residue_id, atom_name) -> (occupancy, altloc_letter, line_index)
    for idx, line in enumerate(pdb_lines):
        if not line.startswith(("ATOM", "HETATM")) or line[16] == " ":
            continue
        key = (line[17:27], line[12:16].strip())
        try:
            occ = float(line[54:60])
        except ValueError:
            occ = 0.0
        alt = line[16]
        current = best.get(key)
        if current is None or occ > current[0] or (occ == current[0] and alt < current[1]):
            best[key] = (occ, alt, idx)

    filtered = []
    for idx, line in enumerate(pdb_lines):
        if not line.startswith(("ATOM", "HETATM")) or line[16] == " ":
            filtered.append(line)
            continue
        key = (line[17:27], line[12:16].strip())
        if idx == best[key][2]:
            filtered.append(line[:16] + " " + line[17:])  # clear altloc indicator
    return filtered


def normalize_backbone_amide_h(npdb_i, resname):
    """Rename a lone, misnamed backbone amide proton to the internal-residue name 'H'.

    A free N-terminus carries three backbone protons (H1, H2, H3 on the NH3+ group)
    and is left untouched so nc_termini_search can relabel it to its N-terminal
    variant. A residue that follows an ACE cap (or is otherwise internal) sometimes
    carries a single backbone proton misnamed 'H1' (or 'H2'); the internal library
    entry expects it to be named 'H', so rename it.
    """
    if resname not in AMINO_ACID_RESIDUES:
        return npdb_i  # DNA bases, ions and caps keep their H1/H2/H3 atoms
    present = {line[12:16].strip() for line in npdb_i}
    nterminal_protons = present & {"H1", "H2", "H3", "H11"}
    if "H" in present or len(nterminal_protons) != 1:
        return npdb_i
    lone = nterminal_protons.pop()
    return [extract_and_replace(line, lone, "H") for line in npdb_i]


def validate_backbone_amide_geometry(npdb_i, resname):
    """Reject explicit backbone amide hydrogens too far from their library parent N."""
    if resname not in AMINO_ACID_RESIDUES:
        return npdb_i
    lines_by_name = {line[12:16].strip(): line for line in npdb_i if line.startswith(("ATOM", "HETATM"))}
    if "H" not in lines_by_name or "N" not in lines_by_name:
        return npdb_i
    distance = math.dist(_atom_coordinates(lines_by_name["H"]), _atom_coordinates(lines_by_name["N"]))
    if distance > _MAX_XH_BOND_DISTANCE:
        residue_id = lines_by_name["H"][21:27].strip()
        raise ValueError(
            f"{resname} {residue_id} backbone H-N distance is {distance:.3f} A; rebuild the hydrogen"
        )
    return npdb_i


def normalize_neutral_lysine_hydrogens(npdb_i, resname):
    """Normalize neutral lysine NZ hydrogens to AMBER14sb HZ2/HZ3 names."""
    if resname != "LYN":
        return npdb_i
    present = {line[12:16].strip() for line in npdb_i}
    hydrogen_names = present & {"HZ1", "HZ2", "HZ3"}
    if hydrogen_names == {"HZ2", "HZ3"} or len(hydrogen_names) < 2:
        return npdb_i
    if hydrogen_names == {"HZ1", "HZ2"}:
        return _rename_atoms(npdb_i, {"HZ1": "HZ2", "HZ2": "HZ3"})
    raise ValueError(f"LYN NZ hydrogen names are inconsistent: {sorted(hydrogen_names)}")


_CTERMINAL_OXYGEN_ALIASES = {
    "O1": "O", "OT1": "O", "OC1": "O",
    "O2": "OXT", "OT2": "OXT", "OC2": "OXT",
} # fmt: skip


def normalize_cterminal_oxygens(npdb_i, resname):
    """Rename C-terminal carboxylate oxygens (O1/O2, OT1/OT2, OC1/OC2) to O/OXT.

    These names appear only on C-terminal residues. The backbone carbonyl oxygen
    becomes 'O' and the extra terminal oxygen becomes 'OXT', but only when the
    standard name is not already present, so an existing O/OXT is never overwritten.
    Running before nc_termini_search lets it detect the C-terminus via OXT.
    """
    if resname not in AMINO_ACID_RESIDUES:
        return npdb_i
    present = {line[12:16].strip() for line in npdb_i}
    renames = {}

    # Some Maestro exports use O for the carbonyl oxygen and O1 for a protonated terminal hydroxyl.
    # AMBER14sb uses the charged C-terminal O/OXT form, so retain the oxygen coordinate as OXT
    # and remove only its directly attached proton.
    if "O" in present and "OXT" not in present:
        extra_oxygens = present & {"O1", "O2", "OT1", "OT2", "OC1", "OC2"}
        if len(extra_oxygens) == 1:
            extra_oxygen = next(iter(extra_oxygens))
            oxygen_line = next(line for line in npdb_i if line[12:16].strip() == extra_oxygen)
            oxygen_xyz = _atom_coordinates(oxygen_line)
            filtered = []
            for line in npdb_i:
                atom_name = line[12:16].strip()
                if (
                    atom_name.startswith("H")
                    and math.dist(_atom_coordinates(line), oxygen_xyz) <= _MAX_XH_BOND_DISTANCE
                ):
                    continue
                filtered.append(line)
            npdb_i = filtered
            present = {line[12:16].strip() for line in npdb_i}
            renames[extra_oxygen] = "OXT"

    for alias, target in _CTERMINAL_OXYGEN_ALIASES.items():
        if alias in present and target not in present and target not in renames.values():
            renames[alias] = target
    for alias, target in renames.items():
        npdb_i = [extract_and_replace(line, alias, target) for line in npdb_i]
    return npdb_i


def _fix_duplicate_backbone_h(pdb_lines):
    """Rename duplicate backbone H atoms to H1, H2 before nest_pdb.

    Maestro exports N-terminal residues with multiple atoms named "H" instead of
    H1, H2, H3. nest_pdb splits on duplicate atom names, so we must rename them
    before nesting to keep the residue intact.
    """
    # Group lines by (chain, resnum, resname) to find duplicates
    residue_h_indices = {}
    for idx, line in enumerate(pdb_lines):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "H":
            key = line[17:27]  # residue identifier (resname + chain + resnum)
            residue_h_indices.setdefault(key, []).append(idx)

    for key, indices in residue_h_indices.items():
        if len(indices) >= 2:
            for count, idx in enumerate(indices, start=1):
                line = pdb_lines[idx]
                pdb_lines[idx] = line[:12] + f" H{count} " + line[16:]
    return pdb_lines


def fix_pdb(pdb_path: Path, rename_mapping=rename_mapping, out_name=None):
    renamed_pdb_path = (
        pdb_path.with_name(pdb_path.stem + "_renamed.pdb") if out_name is None else Path(out_name)
    )
    with open(pdb_path) as f:
        pdb_lines = f.readlines()

    pdb_lines = _filter_altloc(pdb_lines)
    pdb_lines = _fix_duplicate_backbone_h(pdb_lines)
    pdb_lines = _normalize_arn_before_nesting(pdb_lines)
    pdb_lines = _validate_unique_atom_names(pdb_lines)
    npdb = nest_pdb(pdb_lines)
    npdb = asp_search(npdb)
    npdb = glu_search(npdb)  # TODO; check if this one is necessary
    npdb = histidine_search(npdb)

    for i, res in enumerate(npdb):
        resname = res[-1][17:21].strip()
        if resname == "NMA":  # we use NME in our FF library
            npdb[i] = [x.replace("NMA", "NME") for x in npdb[i]]
            resname = "NME"
        elif resname.upper() in ION_RENAME:
            new_name = ION_RENAME[resname.upper()]
            for j, line in enumerate(npdb[i]):
                if line.startswith(("ATOM", "HETATM")):
                    # Atom name (cols 12-15) and residue name (cols 17-20)
                    npdb[i][j] = line[:12] + f"{new_name:<4}" + line[16:17] + f"{new_name:<4}" + line[21:]
            resname = new_name
        # Apply residue-specific mappings before the generic numbered-atom rule.
        # Reversing this order turns 1HB/2HB/3HB into HB2/HB3/HB3 and creates
        # duplicate atom names for methyl groups.
        npdb[i] = correct_amino_acid_atom_names(npdb[i], resname, rename_mapping)
        npdb[i] = correct_numbered_atom_names(npdb[i])
        npdb[i] = normalize_protonated_carboxylate(npdb[i], resname)
        npdb[i] = normalize_neutral_lysine_hydrogens(npdb[i], resname)
        npdb[i] = normalize_cterminal_oxygens(npdb[i], resname)
        npdb[i] = normalize_backbone_amide_h(npdb[i], resname)
        npdb[i] = validate_backbone_amide_geometry(npdb[i], resname)

    npdb = nc_termini_search(npdb)  # after atom name correction, label N and C termini
    npdb = correct_neutral_arginine(npdb)
    npdb = [normalize_neutral_arginine_geometry(residue) for residue in npdb]
    pdb_lines = unnest_pdb(npdb)

    # Renaming can map different source atom names to the same AMBER name.
    # E.g., Input: 1HB, HB1 -> output HB1, HB1
    # Check again to ensure normalization did not create duplicates.
    pdb_lines = _validate_unique_atom_names(pdb_lines)

    with open(renamed_pdb_path, "w") as f:
        for line in pdb_lines:
            f.write(line)
    logger.info(f"Renaming completed. Output file: {renamed_pdb_path}")

    return pdb_lines


def correct_neutral_arginine(npdb):
    """
    Updates residue naming for ARN (neutral arginine - AMBER14sb) residues to ensure
    compatibility.

    Will check for HH22 atom in ARN residues. If HH22 is found, it indicates that the
    NH1 and NH2 groups (and their associated hydrogens) should be renamed to match the
    naming convention in the force field (only HH11, HH12, HH21 should exist). This
    renaming ensures that the nitrogen and hydrogen atoms in the side chain of arginine
    are correctly identified and interact as expected according to the simulation parameters.

    Parameters:
    - npdb (list of list of str): The nested PDB data structure, where each inner list
      represents a residue and each string in the inner list represents a line in the PDB file.

    Returns:
    - list of list of str: The updated nested PDB data structure with corrected atom names
      for ARN residues, ready for use in molecular dynamics simulations.
    """
    atom_name_replacements = {
        "NH1": "NHT",  # temporary
        "HH11": "HHT1",  # temporary
        "NH2": "NH1",
        "HH21": "HH11",
        "HH22": "HH12",
    }
    final_name_replacements = {  # reverse the temporary replacements (T)
        "NHT": "NH2",
        "HHT1": "HH21",
    }

    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "ARN":
            HH22_present = atom_is_present(npdb[i], "HH22")
            if HH22_present:
                for idx, line in enumerate(npdb[i]):
                    atom_name = line[12:16].strip()  # Extract atom name
                    if atom_name in atom_name_replacements:
                        new_atom_name = atom_name_replacements[atom_name]
                        npdb[i][idx] = f"{line[:12]}{new_atom_name:<4}{line[16:]}"

                for idx, line in enumerate(npdb[i]):
                    atom_name = line[12:16].strip()  # Extract atom name
                    if atom_name in final_name_replacements:
                        new_atom_name = final_name_replacements[atom_name]
                        npdb[i][idx] = f"{line[:12]}{new_atom_name:<4}{line[16:]}"
        else:
            continue

    return npdb


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


def _count_atom(npdb_i, atomname):
    """Count how many times an atom name appears in a residue's PDB lines."""
    return sum(1 for line in npdb_i if line[12:16].strip() == atomname)


def _rename_duplicate_h_to_h1_h2(npdb_i):
    """Rename duplicate H atoms to H1, H2 for N-terminal residues from Maestro.

    Maestro sometimes exports N-terminal residues with multiple atoms named "H"
    instead of H1, H2, H3. This renames them sequentially.
    """
    h_count = 0
    for j, line in enumerate(npdb_i):
        if line[12:16].strip() == "H":
            h_count += 1
            npdb_i[j] = line[:12] + f" H{h_count} " + line[16:]
    return npdb_i


def nc_termini_search(npdb):

    NATURAL_AA = [
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE",
        "HIP", "HYP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ] # fmt: skip

    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].strip()
        if resname in NATURAL_AA:
            if resname in ["CYM", "ASH", "GLH", "LYN"]:
                continue  # no parameter for those on C or N terminus

            H3_present = atom_is_present(npdb[i], "H3")  # n-terminus
            H1_present = _count_atom(npdb[i], "H1") >= 1  # exact match (Maestro N-terminal)
            H2_present = _count_atom(npdb[i], "H2") >= 1  # exact match (avoid HH22 false positive)
            duplicate_H = _count_atom(npdb[i], "H") >= 2  # Maestro duplicate H on N-terminal
            OXT_present = atom_is_present(npdb[i], "OXT")  # c-terminus
            HXT_present = _count_atom(npdb[i], "HXT") >= 1  # c-terminus (Maestro convention)

            is_n_terminal = H3_present or (H1_present and H2_present) or duplicate_H
            is_c_terminal = OXT_present

            if is_n_terminal:
                if resname == "HYP":
                    logger.error(
                        "No parameters available for n-terminal HYP residue!!! Please check your structure"
                    )
                else:
                    if duplicate_H:
                        # Re-rename H→H1,H2 (the aa rename mapping may have undone our earlier fix)
                        npdb[i] = _rename_duplicate_h_to_h1_h2(npdb[i])
                    npdb[i] = [x.replace(f"{resname} ", f"N{resname}") for x in npdb[i]]
            if HXT_present:
                # Remove HXT (Maestro C-terminal H). Don't relabel as C-terminal since
                # the PDB lacks OXT (heavy atom) which qprep can't add automatically.
                npdb[i] = [line for line in npdb[i] if line[12:16].strip() != "HXT"]
            if is_c_terminal:
                npdb[i] = [x.replace(f"{resname} ", f"C{resname}") for x in npdb[i]]
            if is_n_terminal and is_c_terminal:
                raise ValueError(f"residue {npdb[i]} has both N-terminal and C-terminal atoms")
    return npdb


def glu_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "GLU":
            # The carboxyl proton of protonated glutamate is HE2; some sources name it HE1.
            HE_present = atom_is_present(npdb[i], "HE1") or atom_is_present(npdb[i], "HE2")
            if HE_present:
                npdb[i] = [x.replace("GLU", "GLH") for x in npdb[i]]
    return npdb


def asp_search(npdb):
    for i in range(len(npdb)):
        resname = npdb[i][0][17:21].rstrip()
        if resname == "ASP":
            # Both names occur in source formats; geometry is normalized later.
            acidic_hydrogen_present = atom_is_present(npdb[i], "HD1") or atom_is_present(npdb[i], "HD2")
            if acidic_hydrogen_present:
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
        fix_pdb(Path(input_pdb_path), out_name=args.output)

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main_exe()
