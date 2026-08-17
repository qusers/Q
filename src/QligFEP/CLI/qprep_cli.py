"""Module containing the command line interface for the qprep fortran program."""

import argparse
from pathlib import Path
from typing import Optional

import numpy as np

from ..IO import get_force_field_paths, parse_lib, parse_prm_options, run_qprep
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
    """Neutralize ionizable groups in Q's outer SCAAS boundary layer.

    Boundary references follow the selected force field: Q uses the first atom
    of each explicit charge group when ``switch_atoms`` is on (OPLS), and the
    geometric charge-group center when it is off (AMBER14sb).
    """

    def __init__(self, center_coords, radius=25.0, boundary_offset=3.0, force_field="AMBER14sb"):
        self.center = np.asarray(center_coords, dtype=float)
        self.radius = float(radius)
        self.boundary_offset = float(boundary_offset)
        if not 0.0 <= self.boundary_offset <= self.radius:
            raise ValueError("boundary_offset must be between zero and the sphere radius")
        self.rest_bound = self.radius - self.boundary_offset
        self.force_field = force_field
        if "CHARMM" in Path(str(force_field)).stem.upper():
            raise NotImplementedError(
                "Force-field-sensitive boundary neutralization is not implemented for CHARMM; "
                "use AMBER14sb or OPLS, or explicitly skip neutralization."
            )
        self.residue_library = parse_lib(force_field)
        switch_atoms = parse_prm_options(force_field).get("switch_atoms")
        if switch_atoms not in {"on", "off"}:
            raise ValueError(f"Force field {force_field!r} does not define 'switch_atoms on|off'")
        self.switch_atoms = switch_atoms == "on"

        # Charged/neutral residue names are chemical states. Their atom sets,
        # charges, charge groups, and switching atoms come from the force field.
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

        self._charged_heavy_atoms = {
            "GLU": {"OE1", "OE2"},
            "ASP": {"OD1", "OD2"},
            "ARG": {"NH1", "NH2"},
            "LYS": {"NZ"},
            "HIP": {"ND1", "NE2"},
        }

        # Statistics tracking
        self.stats = {
            "forcefield": str(force_field),
            "switch_atoms": self.switch_atoms,
            "boundary_offset": self.boundary_offset,
            "boundary_radius": self.rest_bound,
            "total_charged_residues": 0,
            "residues_outside_boundary": 0,
            "residues_neutralized": 0,
            "terminals_neutralized": 0,
            "salt_bridges_neutralized": 0,
            "boundary_salt_bridges": [],
            "original_total_charge": 0,
            "final_total_charge": 0,
            "modifications": {},
            "remaining_outside_charged": [],
        }

    def _library_atom_names(self, residue_name):
        entry = self.residue_library.get(residue_name)
        if entry is None:
            raise ValueError(f"Residue {residue_name!r} is missing from the {self.force_field} library")
        return [atom["name"] for atom in entry["atoms"]]

    def _charge_groups(self, residue_name):
        entry = self.residue_library.get(residue_name)
        if entry is None:
            raise ValueError(f"Residue {residue_name!r} is missing from the {self.force_field} library")
        return entry["charge_groups"] or [self._library_atom_names(residue_name)]

    def _formal_charge_group(self, residue_name, formal_charge):
        entry = self.residue_library[residue_name]
        charges = {atom["name"]: atom["charge"] for atom in entry["atoms"]}
        groups = self._charge_groups(residue_name)
        matching = [
            atom_names
            for atom_names in groups
            if formal_charge * sum(charges.get(name, 0.0) for name in atom_names) > 0.5
        ]
        if not matching:
            raise ValueError(
                f"Could not identify the formal-charge group for {residue_name} in {self.force_field}"
            )
        return max(matching, key=lambda names: abs(sum(charges.get(name, 0.0) for name in names)))

    def _boundary_reference(self, pdb_group, charge_group_atoms, residue_label):
        available = pdb_group.set_index("atom_name", drop=False)
        if self.switch_atoms:
            switch_atom = charge_group_atoms[0]
            if switch_atom not in available.index:
                raise ValueError(
                    f"Q switch atom {switch_atom} is missing from {residue_label}; "
                    "cannot classify it at the spherical boundary"
                )
            row = available.loc[switch_atom]
            if getattr(row, "ndim", 1) > 1:
                row = row.iloc[0]
            return row[["x", "y", "z"]].to_numpy(dtype=float), f"switch atom {switch_atom}"

        present = pdb_group[pdb_group["atom_name"].isin(charge_group_atoms)]
        if present.empty:
            raise ValueError(f"No charge-group atoms for {residue_label} are present in the PDB")
        missing = set(charge_group_atoms) - set(present["atom_name"])
        if missing:
            logger.debug(
                f"Computing {residue_label} charge-group center without missing atoms: "
                f"{', '.join(sorted(missing))}"
            )
        return present[["x", "y", "z"]].to_numpy(dtype=float).mean(axis=0), "charge-group center"

    def _template_charge(self, residue_name):
        entry = self.residue_library.get(residue_name)
        if entry is None:
            return 0.0
        return sum(atom["charge"] for atom in entry["atoms"])

    def _total_library_charge(self, df):
        total = 0.0
        keys = ["chain_id", "residue_seq_number", "insertion_code"]
        for _, group in df.groupby(keys, dropna=False):
            total += self._template_charge(group["residue_name"].iloc[0])
        return int(round(total))

    def neutralize_outside_residues(self, pdb_file, salt_bridge_cutoff=4.0):
        """Find and neutralize charged residues outside the sphere boundary"""
        df = read_pdb_to_dataframe(pdb_file)
        return self.neutralize_outside_residues_dataframe(df, salt_bridge_cutoff)

    def neutralize_outside_residues_dataframe(self, df, salt_bridge_cutoff=4.0):
        """Neutralize excluded charge groups and direct included salt-bridge partners."""
        reference = "switch atom" if self.switch_atoms else "charge-group center"
        logger.info(
            f"Neutralizing charged residues whose Q {reference} is at or beyond "
            f"{self.rest_bound:.1f}Å (outer {self.boundary_offset:.1f}Å SCAAS layer)"
        )
        self.stats["original_total_charge"] = self._total_library_charge(df)

        # Terminal neutralization runs unconditionally (independent of side chain charges)
        modified_df = self._neutralize_nterminals(df)
        modified_df = self._neutralize_cterminals(modified_df)

        charged_residues_info = self._find_charged_residues(modified_df)

        if not charged_residues_info:
            self.stats["residues_neutralized"] = self.stats["terminals_neutralized"]
            self.stats["final_total_charge"] = self._total_library_charge(modified_df)
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

        # Find salt bridges between outside and inside protein residues
        salt_bridge_pairs = self._find_salt_bridges(outside_residues, inside_residues, salt_bridge_cutoff)
        self.stats["salt_bridges_neutralized"] = len(salt_bridge_pairs)

        protein_to_modify = outside_residues + salt_bridge_pairs

        # Find DNA nucleotides outside the neutralization boundary (C1' > rest_bound)
        dna_to_neutralize = self._find_outside_dna(df)

        self.stats["residues_neutralized"] = (
            len(protein_to_modify) + len(dna_to_neutralize) + self.stats["terminals_neutralized"]
        )

        # Modify the dataframe: protein side chains neutralized at rest_bound,
        # DNA outside boundary removed with 5' terminal capping
        modified_df = self._modify_residues(modified_df, protein_to_modify)
        if dna_to_neutralize:
            modified_df = self._remove_and_cap_outside_dna(modified_df, dna_to_neutralize)

        # Track the full force-field template charge after all transformations.
        self.stats["final_total_charge"] = self._total_library_charge(modified_df)

        # Check for remaining charged residues outside boundary
        self._check_remaining_outside_charged(inside_residues, salt_bridge_pairs)

        # Log statistics
        self._log_neutralization_stats()

        return modified_df, self.stats

    def _find_charged_residues(self, df):
        """Find charged residues and calculate their force-field boundary references."""
        charged_residues_info = []
        group_keys = ["chain_id", "residue_seq_number", "insertion_code"]

        for res_name, (neutral_form, dna_key_atom, charge) in self.charged_residues.items():
            residues = df[df["residue_name"] == res_name]
            for (chain, res_num, insertion_code), group in residues.groupby(group_keys, dropna=False):
                if res_name in self._DNA_RESIDUES:
                    reference_atoms = [dna_key_atom]
                else:
                    reference_atoms = self._formal_charge_group(res_name, charge)

                reference_coords, reference_kind = self._boundary_reference(
                    group, reference_atoms, f"{res_name} {chain}:{res_num}{insertion_code}"
                )
                distance = float(np.linalg.norm(reference_coords - self.center))

                charged_atom_names = self._charged_heavy_atoms.get(res_name, {dna_key_atom})
                charged_atoms = group[group["atom_name"].isin(charged_atom_names)]
                charged_atom_coords = charged_atoms[["x", "y", "z"]].to_numpy(dtype=float)
                if len(charged_atom_coords) == 0:
                    charged_atom_coords = np.array([reference_coords])

                charged_residues_info.append(
                    {
                        "residue_name": res_name,
                        "chain_id": chain,
                        "residue_seq_number": res_num,
                        "insertion_code": insertion_code,
                        "boundary_reference": reference_kind,
                        "boundary_reference_coords": reference_coords,
                        "charged_atom_coords": charged_atom_coords,
                        "distance": distance,
                        "charge": charge,
                        "neutral_form": neutral_form,
                    }
                )

        return charged_residues_info

    def _classify_residues_by_distance(self, charged_residues_info):
        """Classify residues as inside or outside the boundary"""
        outside_residues = []
        inside_residues = []

        for res_info in charged_residues_info:
            if res_info["distance"] >= self.rest_bound:
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
        """Find direct salt bridges crossing the included/excluded boundary."""
        salt_bridge_partners = {}

        for outside_res in outside_residues:
            for inside_res in inside_residues:
                if outside_res["charge"] * inside_res["charge"] >= 0:
                    continue

                deltas = (
                    outside_res["charged_atom_coords"][:, None, :]
                    - inside_res["charged_atom_coords"][None, :, :]
                )
                distance = float(np.linalg.norm(deltas, axis=2).min())

                if distance <= cutoff:
                    identity = (
                        inside_res["chain_id"],
                        inside_res["residue_seq_number"],
                        inside_res["insertion_code"],
                    )
                    salt_bridge_partners[identity] = inside_res
                    bridge = {
                        "excluded": {
                            "chain_id": outside_res["chain_id"],
                            "residue_seq_number": int(outside_res["residue_seq_number"]),
                            "residue_name": outside_res["residue_name"],
                        },
                        "included": {
                            "chain_id": inside_res["chain_id"],
                            "residue_seq_number": int(inside_res["residue_seq_number"]),
                            "residue_name": inside_res["residue_name"],
                        },
                        "distance": distance,
                    }
                    if bridge not in self.stats["boundary_salt_bridges"]:
                        self.stats["boundary_salt_bridges"].append(bridge)
                    logger.info(
                        f"Boundary-crossing salt bridge: "
                        f"{outside_res['chain_id']}:{outside_res['residue_seq_number']} "
                        f"({outside_res['residue_name']}) <-> "
                        f"{inside_res['chain_id']}:{inside_res['residue_seq_number']} "
                        f"({inside_res['residue_name']}) at {distance:.2f}Å; "
                        "neutralizing the included partner"
                    )

        return list(salt_bridge_partners.values())

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
                indices_to_remove = atom_mask[atom_mask].index.intersection(modified_df.index)
                modified_df = modified_df.drop(indices_to_remove)

            mod_key = f"{chain}:{res_num}{res_info['insertion_code']}"
            self.stats["modifications"][mod_key] = {
                "original": old_name,
                "modified": new_name,
                "distance": distance,
                "boundary_reference": res_info["boundary_reference"],
                "atoms_removed": atoms_to_remove,
            }

            logger.debug(
                f"Neutralized {old_name} -> {new_name} at {chain}:{res_num} (distance: {distance:.2f}Å)"
            )
            if atoms_to_remove:
                logger.debug(f"Removed atoms: {', '.join(atoms_to_remove)}")

        return modified_df

    _DNA_RESIDUES = {"DA", "DC", "DG", "DT"}

    _INTERNAL_PROTEIN_RESIDUES = {
        "ALA", "ARG", "ARN", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLU", "GLH",
        "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    } #fmt: skip

    def _terminal_neutral_form(self, terminal_name, terminal_charge):
        """Find the internal template after removing an N- or C-terminal charge.

        Selection uses template charge and atom-set similarity, which also covers
        the historical OPLS2005 names such as NAR+, NLYP, CAR+, and CLYP.
        """
        source_atoms = set(self._library_atom_names(terminal_name))
        target_charge = self._template_charge(terminal_name) - terminal_charge
        candidates = []
        for name in self._INTERNAL_PROTEIN_RESIDUES & self.residue_library.keys():
            charge_error = abs(self._template_charge(name) - target_charge)
            if charge_error > 0.25:
                continue
            target_atoms = set(self._library_atom_names(name))
            similarity = len(source_atoms & target_atoms) / len(source_atoms | target_atoms)
            candidates.append((similarity, -charge_error, name))
        if not candidates:
            raise ValueError(
                f"Could not find an internal template for terminal residue {terminal_name} "
                f"in {self.force_field}"
            )
        return max(candidates)[2]

    def _terminal_charge_group(self, residue_name, marker_atom):
        matching = [group for group in self._charge_groups(residue_name) if marker_atom in group]
        if len(matching) != 1:
            raise ValueError(
                f"Expected one {residue_name} charge group containing {marker_atom}, found {len(matching)}"
            )
        return matching[0]

    def _convert_residue_template(self, df, mask, old_name, new_name, preserve_one_hydrogen=False):
        """Rename a residue and remove atoms absent from its new force-field template."""
        new_atoms = self._library_atom_names(new_name)
        present_atoms = list(df.loc[mask, "atom_name"])
        atoms_to_remove = [name for name in present_atoms if name not in set(new_atoms)]

        # N-terminal force fields use H1/H2/H3 or HT1/HT2/HT3 whereas the
        # internal residue uses H. Preserve one equivalent N-H coordinate.
        if preserve_one_hydrogen and "H" in new_atoms and "H" not in present_atoms:
            old_n_hydrogens = [name for name in ("H1", "HT1", "H2", "HT2") if name in atoms_to_remove]
            if old_n_hydrogens:
                source = old_n_hydrogens[0]
                df.loc[mask & (df["atom_name"] == source), "atom_name"] = "H"
                atoms_to_remove.remove(source)

        df.loc[mask, "residue_name"] = new_name
        for atom_name in atoms_to_remove:
            current_mask = mask.reindex(df.index, fill_value=False)
            df = df.drop(df[current_mask & (df["atom_name"] == atom_name)].index)
        return df, atoms_to_remove

    def _neutralize_nterminals(self, df):
        """Convert charged N-terminal templates outside the boundary to internal forms."""
        modified_df = df.copy()
        candidate_names = {
            name
            for name in modified_df["residue_name"].unique()
            if name.startswith("N") and len(name) > 3 and name in self.residue_library
        }
        nterm_mask = modified_df["residue_name"].isin(candidate_names)
        group_keys = ["chain_id", "residue_seq_number", "insertion_code", "residue_name"]
        n_converted = 0

        for (chain, res_num, insertion_code, nterm_name), group in modified_df[nterm_mask].groupby(
            group_keys, dropna=False
        ):
            reference_atoms = self._terminal_charge_group(nterm_name, "N")
            reference_coords, reference_kind = self._boundary_reference(
                group, reference_atoms, f"{nterm_name} {chain}:{res_num}{insertion_code}"
            )
            dist = float(np.linalg.norm(reference_coords - self.center))
            if dist < self.rest_bound:
                continue

            internal_name = self._terminal_neutral_form(nterm_name, terminal_charge=+1.0)
            mask = (
                (modified_df["chain_id"] == chain)
                & (modified_df["residue_seq_number"] == res_num)
                & (modified_df["insertion_code"] == insertion_code)
            )
            modified_df, atoms_removed = self._convert_residue_template(
                modified_df, mask, nterm_name, internal_name, preserve_one_hydrogen=True
            )
            self.stats["modifications"][f"{chain}:{res_num}{insertion_code}"] = {
                "original": nterm_name,
                "modified": internal_name,
                "distance": dist,
                "boundary_reference": reference_kind,
                "atoms_removed": atoms_removed,
            }
            logger.debug(
                f"Neutralized N-terminal {nterm_name} -> {internal_name} at {chain}:{res_num} "
                f"({reference_kind}: {dist:.2f}Å)"
            )
            n_converted += 1

        self.stats["terminals_neutralized"] += n_converted
        if n_converted:
            logger.info(f"Neutralized {n_converted} N-terminal residues outside boundary")
        return modified_df

    def _neutralize_cterminals(self, df):
        """Convert charged C-terminal templates outside the boundary to internal forms."""
        modified_df = df.copy()
        candidate_names = {
            name
            for name in modified_df["residue_name"].unique()
            if name.startswith("C") and len(name) > 3 and name in self.residue_library
        }
        cterm_mask = modified_df["residue_name"].isin(candidate_names)
        group_keys = ["chain_id", "residue_seq_number", "insertion_code", "residue_name"]
        n_converted = 0

        for (chain, res_num, insertion_code, cterm_name), group in modified_df[cterm_mask].groupby(
            group_keys, dropna=False
        ):
            reference_atoms = self._terminal_charge_group(cterm_name, "OXT")
            reference_coords, reference_kind = self._boundary_reference(
                group, reference_atoms, f"{cterm_name} {chain}:{res_num}{insertion_code}"
            )
            dist = float(np.linalg.norm(reference_coords - self.center))
            if dist < self.rest_bound:
                continue

            internal_name = self._terminal_neutral_form(cterm_name, terminal_charge=-1.0)
            mask = (
                (modified_df["chain_id"] == chain)
                & (modified_df["residue_seq_number"] == res_num)
                & (modified_df["insertion_code"] == insertion_code)
            )
            modified_df, atoms_removed = self._convert_residue_template(
                modified_df, mask, cterm_name, internal_name
            )
            self.stats["modifications"][f"{chain}:{res_num}{insertion_code}"] = {
                "original": cterm_name,
                "modified": internal_name,
                "distance": dist,
                "boundary_reference": reference_kind,
                "atoms_removed": atoms_removed,
            }
            logger.debug(
                f"Neutralized C-terminal {cterm_name} -> {internal_name} at {chain}:{res_num} "
                f"({reference_kind}: {dist:.2f}Å)"
            )
            n_converted += 1

        self.stats["terminals_neutralized"] += n_converted
        if n_converted:
            logger.info(f"Neutralized {n_converted} C-terminal residues outside boundary")
        return modified_df

    def _find_outside_dna(self, df):
        """Find DNA nucleotides whose C1' atom is outside the sphere radius.

        DNA in the outer SCAAS layer is kept because changing a nucleotide to
        a neutral protein-style template would break its backbone chemistry.
        Only nucleotides beyond the simulation sphere are removed and capped.
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
                    modified_df = modified_df.drop(atom_mask[atom_mask].index.intersection(modified_df.index))
                logger.debug(f"5' terminal cap: {old_name} -> {new_name} at {chain}:{first_inside}")
                n_caps += 1

        logger.info(
            f"Removed {len(fully_outside_residues)} out-of-sphere DNA nucleotides, "
            f"applied {n_caps} 5' terminal caps"
        )
        return modified_df

    def _get_atoms_to_remove(self, old_name, new_name):
        """Return atoms present only in the charged force-field template."""
        old_atoms = self._library_atom_names(old_name)
        new_atoms = set(self._library_atom_names(new_name))
        atoms_to_remove = [name for name in old_atoms if name not in new_atoms]
        unexpected = [name for name in atoms_to_remove if not name.startswith("H")]
        if unexpected:
            raise ValueError(
                f"Refusing {old_name}->{new_name}: neutralization would remove heavy atoms "
                f"{', '.join(unexpected)}"
            )
        return atoms_to_remove

    def _check_remaining_outside_charged(self, inside_residues, salt_bridge_partners):
        """Check for remaining charged residues outside the boundary"""
        for res in inside_residues:
            if res not in salt_bridge_partners and res["distance"] >= self.rest_bound:
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
            "Width of the outer SCAAS layer in which ionizable protein groups are neutralized. "
            "Q's force-field switching atom or charge-group center is compared with "
            "(radius - offset). Defaults to 3.0 Å, matching Q's water-polarization layer "
            "and the QresFEP preparation protocol. Use 0 for qprep's nominal exclusion boundary."
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
    return parser.parse_args()


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> Optional[dict]:
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

    # Step 1: Remove crystal waters first (they will be replaced by sphere waters)
    pdb_data = read_pdb_to_dataframe(pdb_file)
    crystal_waters_df = pdb_data.query("residue_name == 'HOH'")
    if not crystal_waters_df.empty:
        logger.info(f"Removing {len(crystal_waters_df)} crystal water molecules")
        pdb_data = pdb_data.query("residue_name != 'HOH'")
        processing_steps.append("extracted waters to water.pdb")

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
            center_coords,
            args.sphereradius,
            args.neutralize_boundary_offset,
            force_field=args.FF,
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
    reindex_pdb_residues(processed_pdb_path, processed_pdb_path)

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

    return neutralization_stats


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
