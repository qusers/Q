"""Module containing the command line interface for the qprep fortran program."""

import argparse
from pathlib import Path
from typing import Optional

import numpy as np

from .. import sphere_prep
from ..IO import get_force_field_paths, parse_lib, parse_qprep_total_charge, run_qprep
from ..logger import logger, setup_logger
from ..pdb_utils import (
    append_pdb_to_another,
    filter_out_of_sphere_fragments,
    read_pdb_to_dataframe,
    reindex_pdb_residues,
    write_dataframe_to_pdb,
)
from ..settings.settings import CONFIGS
from ..templates import QprepProteinParameters, render_qprep_protein_input
from ._popc_utils import (
    convert_pymemdyn_to_unified_dataframe,
    df_to_pdb_corrected_element,
    has_pop_residues,
)
from .utils import handle_cysbonds


class Neutralizer:
    """Neutralizes out-of-sphere charged residues (protein and DNA) in the qprep workflow."""

    def __init__(self, center_coords, radius=25.0, boundary_offset=3.0, force_field="AMBER14sb"):
        self.center = np.array(center_coords)
        self.radius = float(radius)
        self.boundary_offset = boundary_offset
        self.rest_bound = self.radius - self.boundary_offset  # Neutralize residues OUTSIDE this boundary
        self.force_field = force_field
        self._library = None

        # Define charged residues and their neutral forms + key atoms for distance calculation
        # Format: 'charged': ['neutral', 'key_atom', charge]
        self.charged_residues = {
            # Protein residues
            "GLU": ["GLH", "CD", -1],  # Glutamic acid -> neutral glutamic acid
            "ASP": ["ASH", "CG", -1],  # Aspartic acid -> neutral aspartic acid
            "ARG": ["ARN", "CZ", +1],  # Arginine -> neutral arginine
            "LYS": ["LYN", "NZ", +1],  # Lysine -> neutral lysine
            "HIP": ["HID", "CG", +1],  # Histidine (protonated) -> histidine (delta protonated)
            # DNA nucleotides (phosphate backbone carries -1 charge each)
            "DA": ["DAN", "C1'", -1],
            "DC": ["DCN", "C1'", -1],
            "DG": ["DGN", "C1'", -1],
            "DT": ["DTN", "C1'", -1],
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
        df = read_pdb_to_dataframe(pdb_file)
        return self.neutralize_outside_residues_dataframe(df, salt_bridge_cutoff)

    def neutralize_outside_residues_dataframe(self, df, salt_bridge_cutoff=4.0):
        """Find and neutralize charged residues outside the sphere boundary for a DataFrame"""
        logger.info(f"Neutralizing charged residues outside {self.rest_bound:.1f}Å boundary")

        # Terminal neutralization runs unconditionally (independent of side chain charges)
        modified_df = self._neutralize_nterminals(df)
        modified_df = self._neutralize_cterminals(modified_df)

        charged_residues_info = self._find_charged_residues(modified_df)

        if not charged_residues_info:
            logger.info("No charged residues found in the PDB data")
            return modified_df, self.stats

        # Separate protein and DNA charged residues for different handling
        protein_charged = [r for r in charged_residues_info if r["residue_name"] not in self._DNA_RESIDUES]
        dna_charged = [r for r in charged_residues_info if r["residue_name"] in self._DNA_RESIDUES]

        self.stats["total_charged_residues"] = len(charged_residues_info)
        logger.info(
            f"Found {len(charged_residues_info)} charged residues ({len(protein_charged)} protein, {len(dna_charged)} DNA)"
        )

        # Classify PROTEIN residues by distance to center (uses rest_bound)
        outside_residues, inside_residues = self._classify_residues_by_distance(protein_charged)
        self.stats["residues_outside_boundary"] = len(outside_residues)

        # Track original charge
        self.stats["original_total_charge"] = sum(res["charge"] for res in charged_residues_info)

        # Find salt bridges between outside and inside protein residues
        salt_bridge_pairs = self._find_salt_bridges(outside_residues, inside_residues, salt_bridge_cutoff)
        self.stats["salt_bridges_neutralized"] = len(salt_bridge_pairs)

        protein_to_modify = outside_residues + salt_bridge_pairs

        # Find DNA nucleotides outside the neutralization boundary (C1' > rest_bound)
        dna_to_neutralize = self._find_outside_dna(df)

        self.stats["residues_neutralized"] = len(protein_to_modify) + len(dna_to_neutralize)

        # Modify the dataframe: protein side chains neutralized at rest_bound,
        # DNA outside boundary removed with 5' terminal capping
        modified_df = self._modify_residues(modified_df, protein_to_modify)
        if dna_to_neutralize:
            modified_df = self._remove_and_cap_outside_dna(modified_df, dna_to_neutralize)

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

        for res_name in self.charged_residues:
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
                    "insertion_code": group["insertion_code"].iloc[0],
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

            # Change residue name and remove atoms based on residue type. Match on
            # insertion_code too: a cap (e.g. NME 167A) can share chain+seq_number
            # with the residue being neutralized and must not be rewritten.
            mask = (
                (modified_df["chain_id"] == chain)
                & (modified_df["residue_seq_number"] == res_num)
                & (modified_df["insertion_code"] == res_info["insertion_code"])
            )
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

    _DNA_RESIDUES = {"DA", "DC", "DG", "DT"}

    def _neutralize_nterminals(self, df):
        """Convert N-terminal residues outside the neutralization boundary to internal forms.

        In AMBER14sb, N-terminal residues are named "N" + internal name (e.g. NILE → ILE)
        and carry +1 charge from the protonated NH3+ group. Converting to internal form
        removes this charge by renaming H1→H and removing H2, H3.
        """
        modified_df = df.copy()
        nterm_mask = modified_df["residue_name"].str.startswith("N") & (
            modified_df["residue_name"].str.len() > 3
        )
        if not nterm_mask.any():
            return modified_df

        nterm_residues = modified_df[nterm_mask].groupby(["chain_id", "residue_seq_number", "residue_name"])
        n_converted = 0

        for (chain, res_num, nterm_name), group in nterm_residues:
            n_atom = group[group["atom_name"] == "N"]
            if len(n_atom) == 0:
                continue
            dist = np.linalg.norm(n_atom[["x", "y", "z"]].values[0] - self.center)
            if dist <= self.rest_bound:
                continue

            internal_name = nterm_name[1:]  # NILE → ILE, NGLY → GLY, etc.
            ins_code = group["insertion_code"].iloc[0]
            mask = (
                (modified_df["chain_id"] == chain)
                & (modified_df["residue_seq_number"] == res_num)
                & (modified_df["insertion_code"] == ins_code)
            )
            modified_df.loc[mask, "residue_name"] = internal_name
            h1_mask = mask & (modified_df["atom_name"] == "H1")
            if internal_name == "PRO":
                # Proline's N is tertiary (bonded to CA, CD, prev-C) and has no backbone amide H.
                modified_df = modified_df.drop(modified_df[h1_mask].index)
            else:
                modified_df.loc[h1_mask, "atom_name"] = "H"
            for hatom in ["H2", "H3"]:
                h_mask = mask & (modified_df["atom_name"] == hatom)
                modified_df = modified_df.drop(modified_df[h_mask].index)
            logger.debug(
                f"Neutralized N-terminal {nterm_name} -> {internal_name} "
                f"at {chain}:{res_num} (N: {dist:.2f}Å)"
            )
            n_converted += 1

        if n_converted > 0:
            logger.info(f"Neutralized {n_converted} N-terminal residues outside boundary")
        return modified_df

    def _neutralize_cterminals(self, df):
        """Convert C-terminal residues outside the neutralization boundary to internal forms.

        In AMBER14sb, C-terminal residues are named "C" + internal name (e.g. CALA → ALA)
        and carry -1 charge from the deprotonated COO- group. Converting to internal form
        removes this charge by stripping the "C" prefix and removing OXT.
        """
        modified_df = df.copy()
        cterm_mask = modified_df["residue_name"].str.startswith("C") & (
            modified_df["residue_name"].str.len() > 3
        )
        if not cterm_mask.any():
            return modified_df

        cterm_residues = modified_df[cterm_mask].groupby(["chain_id", "residue_seq_number", "residue_name"])
        n_converted = 0

        for (chain, res_num, cterm_name), group in cterm_residues:
            n_atom = group[group["atom_name"] == "N"]
            if len(n_atom) == 0:
                continue
            dist = np.linalg.norm(n_atom[["x", "y", "z"]].values[0] - self.center)
            if dist <= self.rest_bound:
                continue

            internal_name = cterm_name[1:]  # CALA → ALA, CGLY → GLY, etc.
            ins_code = group["insertion_code"].iloc[0]
            mask = (
                (modified_df["chain_id"] == chain)
                & (modified_df["residue_seq_number"] == res_num)
                & (modified_df["insertion_code"] == ins_code)
            )
            modified_df.loc[mask, "residue_name"] = internal_name
            oxt_mask = mask & (modified_df["atom_name"] == "OXT")
            modified_df = modified_df.drop(modified_df[oxt_mask].index)
            logger.debug(
                f"Neutralized C-terminal {cterm_name} -> {internal_name} "
                f"at {chain}:{res_num} (N: {dist:.2f}Å)"
            )
            n_converted += 1

        if n_converted > 0:
            logger.info(f"Neutralized {n_converted} C-terminal residues outside boundary")
        return modified_df

    def _find_outside_dna(self, df):
        """Find DNA nucleotides whose C1' atom is outside the sphere radius.

        DNA in the restrained shell (between rest_bound and radius) is kept
        with its native charges, consistent with how protein residues in the
        shell are handled. Only DNA beyond the sphere radius is removed, as
        Q would exclude these via charge-group-based exclusion, leaving their
        phosphate charges as static artifacts.
        """
        outside = []
        for res_name in self._DNA_RESIDUES:
            residues = df[df["residue_name"] == res_name]
            if len(residues) == 0:
                continue
            for (chain, res_num), group in residues.groupby(["chain_id", "residue_seq_number"]):
                c1_row = group[group["atom_name"] == "C1'"]
                if len(c1_row) == 0:
                    logger.warning(f"C1' not found in {res_name} {chain}:{res_num}")
                    continue
                c1_coords = c1_row[["x", "y", "z"]].values[0]
                c1_dist = np.linalg.norm(c1_coords - self.center)
                if c1_dist > self.radius:
                    outside.append(
                        {
                            "residue_name": res_name,
                            "chain_id": chain,
                            "residue_seq_number": res_num,
                            "charge": self.charged_residues[res_name][2],
                            "neutral_form": self.charged_residues[res_name][0],
                            "distance": c1_dist,
                        }
                    )
                    logger.debug(
                        f"DNA {res_name} {chain}:{res_num} outside sphere "
                        f"(C1': {c1_dist:.2f}Å > {self.radius:.1f}Å)"
                    )
        if outside:
            logger.info(f"Found {len(outside)} DNA nucleotides outside sphere radius")
        return outside

    _DNA_5PRIME = {"DA": "DA5", "DC": "DC5", "DG": "DG5", "DT": "DT5"}

    def _remove_and_cap_outside_dna(self, df, fully_outside_residues):
        """Remove fully-outside DNA and cap the inside boundary residues.

        All DNA nucleotides fully outside the sphere are removed. The last
        inside-sphere residue on each chain boundary is converted to a
        5' terminal form (removing the phosphodiester linkage to the now-absent
        upstream residue) when the removed DNA was on its 5' side.
        """
        modified_df = df.copy()
        if not fully_outside_residues:
            return modified_df

        # Group outside residues by chain
        outside_by_chain = {}
        for r in fully_outside_residues:
            outside_by_chain.setdefault(r["chain_id"], []).append(r)

        # Find all DNA residue numbers per chain (inside + outside)
        all_dna_by_chain = {}
        for res_name in self._DNA_RESIDUES:
            for (chain, res_num), _ in df[df["residue_name"] == res_name].groupby(
                ["chain_id", "residue_seq_number"]
            ):
                all_dna_by_chain.setdefault(chain, set()).add(res_num)

        # Remove all fully-outside DNA
        for r in fully_outside_residues:
            mask = (modified_df["chain_id"] == r["chain_id"]) & (
                modified_df["residue_seq_number"] == r["residue_seq_number"]
            )
            modified_df = modified_df.drop(modified_df[mask].index)
            logger.debug(f"Removed out-of-sphere DNA {r['chain_id']}:{r['residue_seq_number']}")

        # Cap inside boundary residues where removed DNA was upstream (5' side)
        n_caps = 0
        for chain, outside_list in outside_by_chain.items():
            outside_resnums = {r["residue_seq_number"] for r in outside_list}
            inside_resnums = sorted(all_dna_by_chain.get(chain, set()) - outside_resnums)
            if not inside_resnums:
                continue

            first_inside = min(inside_resnums)
            has_removed_upstream = any(rn < first_inside for rn in outside_resnums)
            if not has_removed_upstream:
                continue

            mask = (modified_df["chain_id"] == chain) & (modified_df["residue_seq_number"] == first_inside)
            old_name = modified_df.loc[mask, "residue_name"].iloc[0]
            if old_name in self._DNA_5PRIME:
                new_name = self._DNA_5PRIME[old_name]
                modified_df.loc[mask, "residue_name"] = new_name
                for atom in ["P", "OP1", "OP2", "HP"]:
                    atom_mask = mask & (modified_df["atom_name"] == atom)
                    modified_df = modified_df.drop(modified_df[atom_mask].index)
                logger.debug(f"5' terminal cap: {old_name} -> {new_name} at {chain}:{first_inside}")
                n_caps += 1

        logger.info(
            f"Removed {len(fully_outside_residues)} out-of-sphere DNA nucleotides, "
            f"applied {n_caps} 5' terminal caps"
        )
        return modified_df

    # Which proton a neutral form drops is force-field specific: AMBER14sb's LYN keeps
    # HZ2/HZ3 while OPLS and CHARMM keep HZ1/HZ2, and ARN differs the same way. These are
    # the last-resort answers for residues absent from the library (the DNA forms, mainly).
    _FALLBACK_ATOMS_TO_REMOVE = {
        ("LYS", "LYN"): ["HZ1"],
        ("ARG", "ARN"): ["HH22"],
        ("HIP", "HID"): ["HE2"],
    }

    def _get_atoms_to_remove(self, old_name, new_name):
        """Return the atoms the neutral form does not have.

        Read off the force field library rather than hardcoded, so the right proton is
        dropped whichever naming the library uses. Deprotonation only ever removes atoms;
        forms that gain one (ASP -> ASH, GLU -> GLH) let qprep add it.
        """
        if self._library is None:
            try:
                self._library = parse_lib(self.force_field)
            except Exception as exc:  # unknown/custom force field: fall back to the table
                logger.debug(f"Could not read the {self.force_field} library ({exc})")
                self._library = {}

        charged = self._library.get(old_name)
        neutral = self._library.get(new_name)
        if charged and neutral:
            neutral_atoms = {atom["name"] for atom in neutral["atoms"]}
            return [atom["name"] for atom in charged["atoms"] if atom["name"] not in neutral_atoms]

        return list(self._FALLBACK_ATOMS_TO_REMOVE.get((old_name, new_name), []))

    def _check_remaining_outside_charged(self, inside_residues, salt_bridge_partners):
        """Check for remaining charged residues outside the boundary"""
        for res in inside_residues:
            if res not in salt_bridge_partners and res["distance"] > self.rest_bound:
                self.stats["remaining_outside_charged"].append(res)

    def _log_neutralization_stats(self):
        """Log neutralization statistics"""
        stats = self.stats

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
            "OPLS2005, OPLS2015, OPLSAAM, AMBER14sb, CHARMM36. Defaults to AMBER14sb."
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
    parser.add_argument(
        "-skip-ff",
        "--skip-fragment-filter",
        dest="skip_fragment_filter",
        action="store_true",
        help=(
            "Skip filtering molecular fragments outside the simulation sphere. "
            "Fragment filtering removes entire chains, cofactors, or lipids that are "
            "completely outside the sphere radius, reducing system size for faster "
            "downstream FEP setup."
        ),
    )
    parser.add_argument(
        "--strip-crystal-waters",
        dest="strip_crystal_waters",
        action="store_true",
        help=(
            "Remove crystallographic HOH residues before solvating. By default they are "
            "preserved and qprep removes added solvent that overlaps them."
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

    # Track processing steps for file naming
    processing_steps = []
    original_stem = pdb_path.stem

    ff_lib_path, ff_prm_path = get_force_field_paths(args.FF)

    # Step 1: Keep crystallographic waters by default. qprep recognises input HOH
    # residues as solvent and its normal solvation pass rejects added waters that
    # overlap them. This preserves experimentally resolved hydrogen-bond networks
    # instead of replacing them with nearby grid waters.
    pdb_data = read_pdb_to_dataframe(pdb_file)
    crystal_waters_df = pdb_data.query("residue_name == 'HOH'")
    if not crystal_waters_df.empty:
        water_keys = ["chain_id", "residue_seq_number", "insertion_code"]
        crystal_water_count = len(crystal_waters_df[water_keys].drop_duplicates())
        if getattr(args, "strip_crystal_waters", False):
            logger.info(f"Removing {crystal_water_count} crystallographic water molecules")
            pdb_data = pdb_data.query("residue_name != 'HOH'")
            processing_steps.append("removed crystallographic waters")
        else:
            logger.info(
                f"Preserving {crystal_water_count} crystallographic water molecules; "
                "qprep will resolve overlaps with added solvent"
            )

    # Step 2: Add cofactors if provided
    if args.cofactors:
        logger.info(f"Adding {len(args.cofactors)} cofactor(s)")
        for cofactor in args.cofactors:
            pdb_data = append_pdb_to_another(pdb_data, cwd / cofactor, ignore_waters=True)
        processing_steps.append("added cofactors")

    # Step 3: Neutralization
    neutralization_stats = None
    if not args.skip_neutralization:
        logger.info("Neutralizing charged residues outside spherical boundary")
        center_coords = [float(coord) for coord in args.cog]
        neutralizer = Neutralizer(
            center_coords, args.sphereradius, args.neutralize_boundary_offset, force_field=args.FF
        )
        # Neutralize the protein and get statistics
        pdb_data, neutralization_stats = neutralizer.neutralize_outside_residues_dataframe(
            pdb_data, args.salt_bridge_cutoff
        )
        processing_steps.append("neutralized out-of-sphere residues")
        logger.info("Charged residues neutralized")

    # Step 4: POPC lipid renaming (only for AMBER14sb with POP residues)
    lipid_conversion_performed = False
    pop_count = 0
    if args.FF == "AMBER14sb":
        has_pop, pop_count = has_pop_residues(pdb_data)
        if has_pop:
            logger.info(f"Detected {pop_count} POP lipid residue(s) with AMBER14sb forcefield")
            # Check if POP atoms are already in Q/Lipid21 convention (e.g., from memprep -oqc)
            # Pymemdyn uses N4/P8; Q/Lipid21 uses N31/P31 — these never overlap
            pop_atoms = set(pdb_data[pdb_data["residue_name"] == "POP"]["atom_name"])
            pymemdyn_markers = {"N4", "P8", "C5", "C6", "O7"}
            unified_markers = {"N31", "P31", "O31", "O32"}
            if unified_markers & pop_atoms:
                logger.info(
                    "POP atoms already in Q/Lipid21 convention (e.g., from memprep -oqc) — skipping conversion"
                )
            elif pymemdyn_markers & pop_atoms:
                logger.info("Converting pymemdyn POPC atom names to unified format compatible with AMBER14sb")
                pdb_data = convert_pymemdyn_to_unified_dataframe(pdb_data)
                processing_steps.append("renamed lipid atoms to match Amber14sb convention")
                lipid_conversion_performed = True
            else:
                logger.warning(
                    "POP atoms detected but naming convention is unclear — skipping conversion to be safe"
                )

    # Write the processed protein to a predictable filename
    if processing_steps:
        processed_pdb_path = pdb_path.with_name("protein_processed.pdb")
        logger.info(f"Processing steps applied: {', '.join(processing_steps)}")
    else:
        processed_pdb_path = pdb_path

    # Use element-corrected writing if lipids were converted
    if lipid_conversion_performed:
        df_to_pdb_corrected_element(pdb_data, processed_pdb_path)
    else:
        write_dataframe_to_pdb(pdb_data, processed_pdb_path)
    pdb_file = str(processed_pdb_path)
    pdb_path = processed_pdb_path
    logger.info(f"Final processed protein saved as: {processed_pdb_path}")

    # Step 5: Reindex residues so that every residue has a unique number.
    # Q/qprep ignores chain IDs, so multi-chain systems (e.g. protein + DNA)
    # that reuse residue numbers across chains will collide without this step.
    original_numbering = reindex_pdb_residues(processed_pdb_path, processed_pdb_path)

    # Step 6: Filter out-of-sphere molecular fragments (chains, cofactors, lipids)
    if not args.skip_fragment_filter:
        center_coords = [float(coord) for coord in args.cog]
        orig_count, filt_count = filter_out_of_sphere_fragments(
            processed_pdb_path, center_coords, float(args.sphereradius)
        )
        if filt_count < orig_count:
            logger.info(
                f"Fragment filtering: {orig_count} → {filt_count} atoms "
                f"({orig_count - filt_count} removed)"
            )
        else:
            logger.info("Fragment filtering: no fragments outside sphere to remove")

    qprep_inp_path = cwd / "qprep.inp"
    qprep_out_path = cwd / "qprep.out"

    cysbonds = handle_cysbonds(args.cysbond, pdb_file, comment_out=True)

    params = QprepProteinParameters(
        ff_lib_path=ff_lib_path,
        ff_prm_path=ff_prm_path,
        pdb_file_path=pdb_file,
        cog=cog,
        sphere_radius=sphereradius,
        solvent_pack=formatted_solvent_pack,
        cysbonds=cysbonds,
    )

    if qprep_inp_path.exists():
        logger.warning("qprep.inp already exists!! Overwriting...")
    with qprep_inp_path.open("w") as qprep_inp_f:
        qprep_inp_f.write(render_qprep_protein_input(params))

    logger.debug(f"Running qprep from {qprep_path}")
    run_qprep(qprep_path, "qprep.inp", "qprep.out", args.FF)
    logger.info("qprep run finished. Check the output `qprep.out` for more information.")

    # Log lipid conversion summary if performed
    if lipid_conversion_performed:
        logger.info(
            "LIPID CONVERSION SUMMARY\n"
            f"Converted {pop_count} POP lipid residue(s) from pymemdyn to Amber-unified format. "
            "With the correct heavy atom naming, qprep will add the missing hydrogen atoms. "
            "For this, check the output top_p.pdb from this program."
        )

    # Log neutralization summary if performed and charged residues were found
    if (
        neutralization_stats
        and not args.skip_neutralization
        and neutralization_stats["total_charged_residues"] > 0
    ):
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
    if outside.shape[0] > 0:
        outside_HOH_residues = oxygen_subset.iloc[outside].residue_seq_number.unique()  # noqa: F841
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

    # Record what this run established about the sphere. Setups that follow (QresFEP in
    # particular) need the centre, the enclosed charge, the disulfides and the PDB -> Q
    # residue numbering, none of which survive in the PDB files alone.
    prep = sphere_prep.collect(
        input_pdb=cwd / args.input_pdb_file,
        prepared_pdb=processed_pdb_path,
        force_field=args.FF,
        center=cog,
        radius=args.sphereradius,
        total_charge=parse_qprep_total_charge(qprep_out_path),
        cysbond_lines=cysbonds,
        topology_pdb=cwd / "top_p.pdb",
        original_numbering=original_numbering,
        neutralization_offset=args.neutralize_boundary_offset,
    )
    prep_path = prep.write(cwd)
    logger.info(
        f"{prep_path.name} written: sphere charge {prep.total_charge:+d}, "
        f"{len(prep.residues)} solute residues, {len(prep.disulfides)} disulfide(s)."
    )


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
