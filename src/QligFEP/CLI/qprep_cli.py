"""Module containing the command line interface for the qprep fortran program."""

import argparse
import re
from pathlib import Path
from typing import Optional

import numpy as np

from ..IO import get_force_field_paths, run_command
from ..logger import logger, setup_logger
from ..pdb_utils import (
    append_pdb_to_another,
    read_pdb_to_dataframe,
    write_dataframe_to_pdb,
)
from ..settings.settings import CONFIGS, FF_DIR
from .utils import handle_cysbonds

# NOTE: cysbonds will have \n after each bond -> `maketop MKC_p` is in a different line
qprep_inp_content = """rl {ff_lib_path}
rprm {ff_prm_path}
! TO DO Change if protein system is used
rp {pdb_file_path}
! set solute_density 0.05794
! NOTE, this is now large for water system, change for protein system
set solvent_pack {solvent_pack}
boundary 1 {cog} {sphereradius}
solvate {cog} {sphereradius} 1 HOH
{cysbond}maketop MKC_p
writetop dualtop.top
wp top_p.pdb y
rt dualtop.top
mask none
mask not excluded
wp complexnotexcluded.pdb y
q
"""


class QprepError(Exception):
    pass


class QprepAtomLibMissingError(Exception):
    pass


class ProteinNeutralizer:
    """Protein neutralizer for out-of-sphere charged residues in qprep workflow"""

    def __init__(self, center_coords, radius=25.0, boundary_offset=3.0):
        self.center = np.array(center_coords)
        self.radius = float(radius)
        self.boundary_offset = boundary_offset
        self.rest_bound = self.radius - self.boundary_offset  # Neutralize residues OUTSIDE this boundary

        # Define charged residues and their neutral forms + key atoms for distance calculation
        # Format: 'charged': ['neutral', 'key_atom', charge]
        self.charged_residues = {
            "GLU": ["GLH", "CD", -1],  # Glutamic acid -> neutral glutamic acid
            "ASP": ["ASH", "CG", -1],  # Aspartic acid -> neutral aspartic acid
            "ARG": ["ARN", "CZ", +1],  # Arginine -> neutral arginine
            "LYS": ["LYN", "NZ", +1],  # Lysine -> neutral lysine
            "HIP": ["HID", "CG", +1],  # Histidine (protonated) -> histidine (delta protonated)
        }

        # Statistics tracking
        self.stats = {
            "total_charged_residues": 0,
            "residues_outside_boundary": 0,
            "residues_neutralized": 0,
            "salt_bridges_neutralized": 0,
            "original_total_charge": 0,
            "final_total_charge": 0,
            "modifications": {},
            "remaining_outside_charged": [],
        }

    def neutralize_outside_residues(self, pdb_file, salt_bridge_cutoff=4.0):
        """Find and neutralize charged residues outside the sphere boundary"""
        logger.info(f"Neutralizing charged residues outside {self.rest_bound:.1f}Å boundary")

        df = read_pdb_to_dataframe(pdb_file)
        charged_residues_info = self._find_charged_residues(df)

        if not charged_residues_info:
            logger.info("No charged residues found in the PDB file")
            return df, self.stats

        self.stats["total_charged_residues"] = len(charged_residues_info)
        logger.info(f"Found {len(charged_residues_info)} charged residues")

        # Classify residues by distance to center
        outside_residues, inside_residues = self._classify_residues_by_distance(charged_residues_info)
        self.stats["residues_outside_boundary"] = len(outside_residues)

        # Track original charge
        self.stats["original_total_charge"] = sum(res["charge"] for res in charged_residues_info)

        # Find salt bridges between outside and inside residues
        salt_bridge_pairs = self._find_salt_bridges(outside_residues, inside_residues, salt_bridge_cutoff)
        self.stats["salt_bridges_neutralized"] = len(salt_bridge_pairs)

        # Neutralize residues outside boundary and their salt bridge partners
        residues_to_modify = outside_residues + salt_bridge_pairs
        self.stats["residues_neutralized"] = len(residues_to_modify)

        # Modify the dataframe
        modified_df = self._modify_residues(df, residues_to_modify)

        # Track final charge and remaining outside charged residues
        self.stats["final_total_charge"] = sum(
            res["charge"] for res in inside_residues if res not in salt_bridge_pairs
        )

        # Check for remaining charged residues outside boundary
        self._check_remaining_outside_charged(inside_residues, salt_bridge_pairs)

        # Log statistics
        self._log_neutralization_stats()

        return modified_df, self.stats

    def _find_charged_residues(self, df):
        """Find all charged residues in the PDB"""
        charged_residues_info = []

        for res_name in self.charged_residues.keys():
            residues = df[df["residue_name"] == res_name]

            if len(residues) == 0:
                continue

            for (chain, res_num), group in residues.groupby(["chain_id", "residue_seq_number"]):
                key_atom = self.charged_residues[res_name][1]
                charge = self.charged_residues[res_name][2]
                key_atom_row = group[group["atom_name"] == key_atom]  # atom for distance calculation

                if len(key_atom_row) == 0:
                    logger.warning(f"Key atom {key_atom} not found in {res_name} {chain}:{res_num}")
                    continue

                key_atom_coords = key_atom_row[["x", "y", "z"]].values[0]
                distance = np.linalg.norm(key_atom_coords - self.center)

                residues_info = {
                    "residue_name": res_name,
                    "chain_id": chain,
                    "residue_seq_number": res_num,
                    "key_atom": key_atom,
                    "key_atom_coords": key_atom_coords,
                    "distance": distance,
                    "charge": charge,
                    "neutral_form": self.charged_residues[res_name][0],
                }
                charged_residues_info.append(residues_info)

        return charged_residues_info

    def _classify_residues_by_distance(self, charged_residues_info):
        """Classify residues as inside or outside the boundary"""
        outside_residues = []
        inside_residues = []

        for res_info in charged_residues_info:
            if res_info["distance"] > self.rest_bound:
                outside_residues.append(res_info)
                logger.debug(
                    f"Residue {res_info['chain_id']}:{res_info['residue_seq_number']} ({res_info['residue_name']}) "
                    f"outside boundary: {res_info['distance']:.2f}Å > {self.rest_bound:.1f}Å"
                )
            else:
                inside_residues.append(res_info)
                logger.debug(
                    f"Residue {res_info['chain_id']}:{res_info['residue_seq_number']} ({res_info['residue_name']}) "
                    f"inside boundary: {res_info['distance']:.2f}Å <= {self.rest_bound:.1f}Å"
                )

        return outside_residues, inside_residues

    def _find_salt_bridges(self, outside_residues, inside_residues, cutoff):
        """Find salt bridges between outside and inside residues"""
        salt_bridge_partners = []

        for outside_res in outside_residues:
            for inside_res in inside_residues:
                # Only consider oppositely charged residues
                if outside_res["charge"] * inside_res["charge"] >= 0:
                    continue

                distance = np.linalg.norm(outside_res["key_atom_coords"] - inside_res["key_atom_coords"])

                if distance <= cutoff:
                    salt_bridge_partners.append(inside_res)
                    logger.info(
                        f"Salt bridge detected: {outside_res['chain_id']}:{outside_res['residue_seq_number']} "
                        f"({outside_res['residue_name']}) <-> {inside_res['chain_id']}:{inside_res['residue_seq_number']} "
                        f"({inside_res['residue_name']}) at {distance:.2f}Å - neutralizing both"
                    )

        return salt_bridge_partners

    def _modify_residues(self, df, residues_to_modify):
        """Modify residues to their neutral forms and remove appropriate atoms"""
        modified_df = df.copy()

        for res_info in residues_to_modify:
            chain = res_info["chain_id"]
            res_num = res_info["residue_seq_number"]
            old_name = res_info["residue_name"]
            new_name = res_info["neutral_form"]
            distance = res_info["distance"]

            # Change residue nameand remove atoms based on residue type
            mask = (modified_df["chain_id"] == chain) & (modified_df["residue_seq_number"] == res_num)
            modified_df.loc[mask, "residue_name"] = new_name
            atoms_to_remove = self._get_atoms_to_remove(old_name, new_name)

            for atom_name in atoms_to_remove:
                atom_mask = mask & (modified_df["atom_name"] == atom_name)
                indices_to_remove = modified_df[atom_mask].index
                modified_df = modified_df.drop(indices_to_remove)

            mod_key = f"{chain}:{res_num}"
            self.stats["modifications"][mod_key] = {
                "original": old_name,
                "modified": new_name,
                "distance": distance,
                "atoms_removed": atoms_to_remove,
            }

            logger.debug(
                f"Neutralized {old_name} -> {new_name} at {chain}:{res_num} (distance: {distance:.2f}Å)"
            )
            if atoms_to_remove:
                logger.debug(f"Removed atoms: {', '.join(atoms_to_remove)}")

        return modified_df

    def _get_atoms_to_remove(self, old_name, new_name):
        """Get atoms to remove when converting charged to neutral form"""
        atoms_to_remove = []

        if old_name == "ASP" and new_name == "ASH":
            # ASP -> ASH: No atoms removed, HD2 is added (handled by forcefield)
            pass
        elif old_name == "GLU" and new_name == "GLH":
            # GLU -> GLH: No atoms removed, HE2 is added (handled by forcefield)
            pass
        elif old_name == "LYS" and new_name == "LYN":
            # LYS -> LYN: Remove HZ1 hydrogen from NZ
            atoms_to_remove = ["HZ1"]
        elif old_name == "ARG" and new_name == "ARN":
            # ARG -> ARN: Remove HH22 hydrogen from NH2
            atoms_to_remove = ["HH22"]
        elif old_name == "HIP" and new_name == "HID":
            # HIP -> HID: Remove HE2 hydrogen from NE2
            atoms_to_remove = ["HE2"]

        return atoms_to_remove

    def _check_remaining_outside_charged(self, inside_residues, salt_bridge_partners):
        """Check for remaining charged residues outside the boundary"""
        for res in inside_residues:
            if res not in salt_bridge_partners and res["distance"] > self.rest_bound:
                self.stats["remaining_outside_charged"].append(res)

    def _log_neutralization_stats(self):
        """Log neutralization statistics"""
        stats = self.stats

        logger.info(
            f"Neutralization statistics:\n"
            f"  Total charged residues found: {stats['total_charged_residues']}"
            f"  Residues outside boundary ({self.rest_bound:.1f}Å): {stats['residues_outside_boundary']}"
            f"  Salt bridge pairs neutralized: {stats['salt_bridges_neutralized']}"
            f"  Total residues neutralized: {stats['residues_neutralized']}"
            f"  Original total charge: {stats['original_total_charge']:+d}"
            f"  Final total charge: {stats['final_total_charge']:+d}"
        )

        if stats["remaining_outside_charged"]:
            for res in stats["remaining_outside_charged"]:
                logger.warning(
                    f"Charged residue {res['chain_id']}:{res['residue_seq_number']} ({res['residue_name']}) "
                    f"remains outside boundary at {res['distance']:.2f}Å"
                )

        if stats["modifications"]:
            logger.debug("Detailed modifications:")
            for res_id, mod_info in stats["modifications"].items():
                logger.debug(
                    f"  {res_id}: {mod_info['original']} -> {mod_info['modified']} "
                    f"(distance: {mod_info['distance']:.2f}Å)"
                )
                if mod_info["atoms_removed"]:
                    logger.debug(f"    Atoms removed: {', '.join(mod_info['atoms_removed'])}")


def qprep_error_check(qprep_out_path: Path, ff_name) -> None:
    """Check for errors in the qprep.out file and raise an exception if any are found.

    Args:
        qprep_out_path: Path to the qprep.out file.
        ff_name: name of the forcefield to point user to the .lib & .prm files.

    Raises:
        QprepError: ff any errors are found in the qprep.out file.
    """
    error_pat = re.compile(r"ERROR\:\s", re.IGNORECASE)
    missing_lib_pat = re.compile(
        r">>> Atom ...?.? in residue no\.\s+\d+ not found in library entry for [A-Z]+"
    )
    outfile_lines = qprep_out_path.read_text().split("\n")
    error_lines = []
    missing_atomlib_lines = []
    for line in outfile_lines:
        if error_pat.findall(line):
            error_lines.append(line)
            logger.error(
                f"Errors found in qprep output file {qprep_out_path}. Please check if the amino "
                "acids in your pdb file match the residue & atom conventions on the forcefield .lib & .prm files:\n"
                f"{FF_DIR/ ff_name}.prm & {FF_DIR/ ff_name}.lib"
            )
        if missing_lib_pat.findall(line):
            missing_atomlib_lines.append(line)
            logger.error(
                f"Errors found in qprep output file {qprep_out_path}. "
                "Your protein file likely contains atoms that are not present in the forcefield's .lib & .prm files:, \n"
                f"{FF_DIR/ ff_name}.prm & {FF_DIR/ ff_name}.lib"
            )

    if error_lines:
        error_message = {"\n".join(error_lines)}
        raise QprepError(error_message)
    if missing_atomlib_lines:
        error_message = {"\n".join(missing_atomlib_lines)}
        raise QprepAtomLibMissingError(error_message)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Runs qprep to generate the water.pdb, the dualtop.top, and the complexnotexcluded.pdb files for "
            "an input protein.pdb file. In case of failure, please inspect the qprep.log file. "
            "Depending on the forcefield you're using, you might need some additional preparation steps "
            "(e.g.: running `pdb2amber` if you're using that forcefield)."
        )
    )
    parser.add_argument(
        "-i",
        "--input_pdb_file",
        dest="input_pdb_file",
        required=True,
        help="input protein to run qprep",
    )
    parser.add_argument(
        "-FF",
        "--forcefield",
        dest="FF",
        default="AMBER14sb",
        help=(
            "Protein forcefield to be used. Valid inputs: existing path to a forcefield file without the extensions"
            "(either .lib, .prm, or Path without the extensions will work) or one of the following: "
            "OPLS2005, OPLS2015, AMBER14sb, CHARMM36. Defaults to AMBER14sb."
        ),
    )
    parser.add_argument(
        "-cog",
        "--center_of_geometry",
        dest="cog",
        help=(
            "Center of geometry for the protein. The format is 'x y z', where all numbers "
            "contain 3 decimal cases. This center of geometry can be obtained using the `qcog ` "
            "command, but if you include ligands using the `-lig` option, this COG will be "
            "automatically calculated. If you want to calculate the COG manually, you can use "
            "this option. Defaults to None."
        ),
        required=False,
        default=None,
        nargs=3,
        type=str,
    )
    parser.add_argument(
        "-r",
        "--sphereradius",
        dest="sphereradius",
        required=False,
        default=25,
        help="Size of the simulation sphere. If float, only one decimal case will be used. Defaults to 25.",
        type=float,
    )
    parser.add_argument(
        "-b",
        "--cysbond",
        dest="cysbond",
        default="auto",
        help=(
            "Add cystein bonds. Input should be formatted with the atom numbers"
            "(participating in the Cys bond) connected by `_` and with different bonds "
            "separated by `,` as in: `atom1_atom2,atom3_atom4`. Defaults to `auto`, where "
            "cystein bonds will be automatically detected within distance of 1.8 to 2.2 A."
        ),
        type=str,
    )
    parser.add_argument(
        "-sp",
        "--solvent_pack",
        dest="solvent_pack",
        default=3.0,
        help=(
            "Parameter to qprep.inp `set solvent_pack`. According to Q's manual, this value "
            "corresponds to the minimum distance between solute and solvent heavy atoms when "
            "adding solvent (e.g.: HOH) and defaults to 2.4. In QligFEP we use a value of 3.0 "
            "for creating this FEP water sphere. Defaults to 3.0."
        ),
        type=float,
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default="info",
        choices=["info", "debug", "warning", "error", "critical"],
        help=(
            "Set the logging level. Defaults to info. "
            "Choose between: info, debug, warning, error, critical."
        ),
        type=str,
    )
    parser.add_argument(
        "-cof",
        "--cofactors",
        dest="cofactors",
        nargs="*",
        help=(
            "List of cofactors to be added to the system. Inputs should be one or more "
            "pdb files containing the cofactors to be added."
        ),
    )
    parser.add_argument(
        "-skip-n",
        "--skip-neutralization",
        dest="skip_neutralization",
        action="store_true",
        help=(
            "Skip neutralizing charged residues outside the spherical boundary. Highly advised against."
            "This neutralization is crucial for maintaining the stability of the system and prevents "
            "electrostatic interference from distant charges. Defaults to False."
        ),
    )
    parser.add_argument(
        "-nbo",
        "--neutralize_boundary_offset",
        dest="neutralize_boundary_offset",
        type=float,
        default=3.0,
        help=(
            "Distance offset from sphere radius to define neutralization boundary. "
            "Residues outside (radius - offset) will be neutralized. Defaults to 3.0Å."
        ),
    )
    parser.add_argument(
        "-sbc",
        "--salt_bridge_cutoff",
        dest="salt_bridge_cutoff",
        type=float,
        default=4.0,
        help=(
            "Distance cutoff for detecting salt bridges between inside and outside "
            "residues. Salt bridge partners will also be neutralized. Defaults to 4.0Å."
        ),
    )
    return parser.parse_args()


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> None:
    """Either runs the qprep program with the given arguments via **kwargs or parses
    the arguments and runs the program.

    Args:
        args: argparse Namespace containing the arguments for the qprep program. Defaults to None
    """
    cwd = Path.cwd()
    if args.log_level != "info":
        setup_logger(args.log_level.upper())
    qprep_path = CONFIGS["QPREP"]
    logger.debug(f"Running qprep from path: {qprep_path}")
    sphereradius = f"{args.sphereradius:.1f}"
    formatted_solvent_pack = f"{args.solvent_pack:.1f}"

    cog = " ".join(args.cog)
    logger.debug(f"COG is {cog}")

    pdb_file = str(cwd / args.input_pdb_file)
    pdb_path = cwd / args.input_pdb_file

    ff_lib_path, ff_prm_path = get_force_field_paths(args.FF)

    if args.cofactors:  # append cofactors to the protein if any are passed...
        pdb_data = read_pdb_to_dataframe(pdb_file)
        for cofactor in args.cofactors:
            pdb_data = append_pdb_to_another(pdb_data, cwd / cofactor, ignore_waters=True)
        cofactor_path = pdb_path.with_name(f"{pdb_path.stem}_plus_cofactors.pdb")
        write_dataframe_to_pdb(pdb_data, cofactor_path)
        pdb_path = cofactor_path
        pdb_file = str(pdb_path)

    qprep_inp_path = cwd / "qprep.inp"
    qprep_out_path = cwd / "qprep.out"

    cysbonds = handle_cysbonds(args.cysbond, pdb_file, comment_out=True)

    # write out without crystal waters - (will be in sphere)
    protein_df = read_pdb_to_dataframe(pdb_file)
    crystal_waters_df = protein_df.query("residue_name == 'HOH'")
    if not crystal_waters_df.empty:
        fname = args.input_pdb_file.split(".")[0]
        write_dataframe_to_pdb(
            protein_df.query("residue_name != 'HOH'"), Path(pdb_file).with_stem(f"{fname}_noHOH")
        )

    neutralization_stats = None
    if not args.skip_neutralization:
        logger.info("Neutralizing charged residues outside spherical boundary")
        center_coords = [float(coord) for coord in args.cog]
        neutralizer = ProteinNeutralizer(center_coords, args.sphereradius, args.neutralize_boundary_offset)
        # Neutralize the protein and get statistics
        neutralized_df, neutralization_stats = neutralizer.neutralize_outside_residues(
            pdb_file, args.salt_bridge_cutoff
        )
        # Write the neutralized protein file
        neutralized_pdb_path = pdb_path.with_name(f"{pdb_path.stem}_neutralized.pdb")
        write_dataframe_to_pdb(neutralized_df, neutralized_pdb_path)
        # Update the pdb_file path to use the neutralized version
        pdb_file = str(neutralized_pdb_path)
        pdb_path = neutralized_pdb_path
        logger.info(f"Neutralized protein saved as: {neutralized_pdb_path}")

    if args is not None:
        param_dict = {
            "pdb_file_path": pdb_file,
            "cog": cog,
            "ff_lib_path": ff_lib_path,
            "ff_prm_path": ff_prm_path,
            "sphereradius": sphereradius,
            "cysbond": cysbonds,
            "solvent_pack": formatted_solvent_pack,
        }
    else:
        param_dict = {**kwargs}

    if qprep_inp_path.exists():
        logger.warning("qprep.inp already exists!! Overwriting...")
    with qprep_inp_path.open("w") as qprep_inp_f:
        qprep_inp_f.write(qprep_inp_content.format(**param_dict))

    options = " < qprep.inp > qprep.out"
    logger.debug(f"Running command {qprep_path} {options}")
    run_command(qprep_path, options, string=True)
    qprep_error_check(qprep_out_path, args.FF)
    logger.info("qprep run finished. Check the output `qprep.out` for more information.")

    # Log neutralization summary if performed
    if neutralization_stats and not args.skip_neutralization:
        logger.info(
            "NEUTRALIZATION SUMMARY\n"
            f"Total charged residues processed: {neutralization_stats['total_charged_residues']}\n"
            f"Residues outside boundary: {neutralization_stats['residues_outside_boundary']}\n"
            f"Salt bridge pairs neutralized: {neutralization_stats['salt_bridges_neutralized']}\n"
            f"Total residues neutralized: {neutralization_stats['residues_neutralized']}\n"
            f"Charge change: {neutralization_stats['original_total_charge']:+d} -> {neutralization_stats['final_total_charge']:+d}\n"
        )
        if neutralization_stats["remaining_outside_charged"]:
            logger.warning(
                f"Warning: {len(neutralization_stats['remaining_outside_charged'])} charged residues remain outside boundary"
            )

    waterfile = Path(cwd / "water.pdb")
    # Write water file and deal with possible errors
    if not Path("complexnotexcluded.pdb").exists():
        logger.error(
            "`complexnotexcluded.pdb` file not found. This is as sign qprep didn't "
            "run correctly. Check the outoput in your console and try again..."
        )
        logger.info(
            "If your console contains something like `libgfortran.so.5: cannot "
            "open shared object file: No such file or directory`, you might need to load "
            "some module in your HPC system that you used to compile Q."
        )
        raise FileNotFoundError("complexnotexcluded.pdb file not found. Something went wrong")
    with open("complexnotexcluded.pdb") as f:
        lines = f.readlines()
    with open("water.pdb", "w") as f:
        water_header = f"TITLE      Water Sphere Generated with Qprep: COG {cog}"
        f.write(f"{water_header}\n")
        for line in lines:
            if line.startswith("ATOM") and line[17:20] == "HOH":
                f.write(line)
        logger.info("water.pdb file created.")

    # Now that the water file is created, we remove water molecules outside the sphere radius
    cog = [float(i) for i in cog.split()]
    water_df = read_pdb_to_dataframe(waterfile)
    oxygen_subset = water_df.query('atom_name == "O"')
    euclidean_distances = oxygen_subset[["x", "y", "z"]].sub(cog).pow(2).sum(1).apply(np.sqrt)
    outside = np.where(euclidean_distances > args.sphereradius * 1.05)[0]  # we add a tolerance of 5%
    outside_HOH_residues = oxygen_subset.iloc[outside].residue_seq_number.unique()  # noqa: F841
    if outside.shape[0] > 0:
        logger.warning(f"Found {outside.shape[0]} water molecules outside the sphere radius.")
        logger.warning("Removing these water molecules from the water.pdb file.")
        todrop_idxs = water_df.query("residue_seq_number in @outside_HOH_residues").index
        water_df.drop(index=todrop_idxs, inplace=True)
        new_distances = (
            water_df.query('atom_name == "O"')[["x", "y", "z"]].sub(cog).pow(2).sum(1).apply(np.sqrt)
        )
        logger.debug(f"Final highest distance is {new_distances.max():.2f} A")
        write_dataframe_to_pdb(water_df, waterfile, header=water_header)
    else:
        logger.info("All water molecules are inside the sphere radius.")
        logger.debug(f"Final highest distance to COG is {euclidean_distances.max():.2f} A")


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
