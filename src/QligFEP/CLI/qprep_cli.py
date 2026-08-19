"""Module containing the command line interface for the qprep fortran program."""

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

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

Coordinates = Sequence[float]
ChargedResidueSpec = tuple[str, str, int]
ResidueInfo = dict[str, Any]
NeutralizationStats = dict[str, Any]
ForceFieldEntry = dict[str, Any]


class Neutralizer:
    """Neutralize ionizable groups in Q's outer SCAAS boundary layer.

    Boundary references follow the selected force field: Q uses the first atom
    of each explicit charge group when ``switch_atoms`` is on (OPLS), and the
    geometric charge-group center when it is off (AMBER14sb).
    """

    def __init__(
        self,
        center_coords: Coordinates,
        radius: float = 25.0,
        boundary_offset: float = 3.0,
        force_field: str = "AMBER14sb",
    ) -> None:
        self.center = np.asarray(center_coords, dtype=float)
        self.radius = float(radius)
        self.boundary_offset = float(boundary_offset)
        if not 0.0 <= self.boundary_offset <= self.radius:
            raise ValueError("boundary_offset must be between zero and the sphere radius")
        self.rest_bound = self.radius - self.boundary_offset
        self.force_field = force_field
        self.residue_library: dict[str, ForceFieldEntry] = parse_lib(force_field)
        # Q defaults to switching atoms when this option is not specified.
        # For explicit charge groups, the first listed atom determines whether
        # the entire group is included or excluded at the spherical boundary.
        switch_atoms = parse_prm_options(force_field).get("switch_atoms", "on")
        if switch_atoms not in {"on", "off"}:
            raise ValueError(f"Force field {force_field!r} must define 'switch_atoms' as on or off")
        self.switch_atoms = switch_atoms == "on"

        # Charged/neutral residue names are chemical states. Their atom sets,
        # charges, charge groups, and switching atoms come from the force field.
        self.charged_residues: dict[str, ChargedResidueSpec] = {
            # Protein residues
            "GLU": ("GLH", "CD", -1),  # Glutamic acid -> neutral glutamic acid
            "ASP": ("ASH", "CG", -1),  # Aspartic acid -> neutral aspartic acid
            "ARG": ("ARN", "CZ", +1),  # Arginine -> neutral arginine
            "LYS": ("LYN", "NZ", +1),  # Lysine -> neutral lysine
            "HIP": ("HID", "CG", +1),  # Histidine (protonated) -> histidine (delta protonated)
            # DNA nucleotides (phosphate backbone carries -1 charge each)
            "DA": ("DAN", "C1'", -1),
            "DC": ("DCN", "C1'", -1),
            "DG": ("DGN", "C1'", -1),
            "DT": ("DTN", "C1'", -1),
        }

        self._charged_heavy_atoms: dict[str, set[str]] = {
            "GLU": {"OE1", "OE2"},
            "ASP": {"OD1", "OD2"},
            "ARG": {"NH1", "NH2"},
            "LYS": {"NZ"},
            "HIP": {"ND1", "NE2"},
        }

        # Statistics tracking
        self.stats: NeutralizationStats = {
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

    _DNA_RESIDUES = {"DA", "DC", "DG", "DT"}
    _TERMINAL_ALIASES = {
        "N": {"H": ("H1", "HT1", "H2", "HT2")},
        "C": {"O": ("OT1",)},
    } # fmt: skip

    _INTERNAL_PROTEIN_RESIDUES = {
        "ALA", "ARG", "ARN", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLU", "GLH",
        "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    } #fmt: skip

    def _library_entry(self, residue_name: str) -> ForceFieldEntry:
        entry = self.residue_library.get(residue_name)
        if entry is None:
            raise ValueError(f"Residue {residue_name!r} is missing from the {self.force_field} library")
        return entry

    def _library_atom_names(self, residue_name: str) -> list[str]:
        return [atom["name"] for atom in self._library_entry(residue_name)["atoms"]]

    def _charge_groups(self, residue_name: str) -> list[list[str]]:
        entry = self._library_entry(residue_name)
        return entry["charge_groups"] or [[atom["name"] for atom in entry["atoms"]]]

    def _formal_charge_group(
        self,
        residue_name: str,
        formal_charge: int,
        marker_atoms: set[str] | None = None,
    ) -> list[str]:
        entry = self._library_entry(residue_name)
        charges = {atom["name"]: atom["charge"] for atom in entry["atoms"]}
        matching = []
        for atom_names in self._charge_groups(residue_name):
            group_charge = sum(charges.get(name, 0.0) for name in atom_names)
            contains_marker = marker_atoms is None or bool(marker_atoms & set(atom_names))
            if abs(group_charge - formal_charge) <= 1e-6 and contains_marker:
                matching.append(atom_names)
        if len(matching) != 1:
            raise ValueError(
                f"Expected one {formal_charge:+d} formal-charge group for {residue_name} "
                f"in {self.force_field}, found {len(matching)}"
            )
        return matching[0]

    def _boundary_reference(
        self,
        pdb_group: pd.DataFrame,
        charge_group_atoms: Sequence[str],
        residue_label: str,
    ) -> tuple[np.ndarray, str]:
        if self.switch_atoms:
            available = pdb_group.set_index("atom_name", drop=False)
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

    def _template_charge(self, residue_name: str) -> float:
        entry = self.residue_library.get(residue_name)
        if entry is None:
            return 0.0
        return sum(atom["charge"] for atom in entry["atoms"])

    def _has_integer_charge_groups(self, residue_name: str) -> bool:
        entry = self._library_entry(residue_name)
        charges = {atom["name"]: atom["charge"] for atom in entry["atoms"]}
        for group in self._charge_groups(residue_name):
            group_charge = sum(charges.get(atom, 0.0) for atom in group)
            if abs(group_charge - round(group_charge)) > 1e-6:
                return False
        return True

    def _terminal_sidechain_neutral_form(
        self,
        terminal_name: str,
        neutral_form: str,
        charge: int,
    ) -> str:
        """Choose a valid terminal neutral template, falling back to the internal form."""
        prefixed_neutral = f"{terminal_name[0]}{neutral_form}"
        expected_charge = self._template_charge(terminal_name) - charge
        if (
            prefixed_neutral in self.residue_library
            and self._has_integer_charge_groups(prefixed_neutral)
            and abs(self._template_charge(prefixed_neutral) - expected_charge) <= 1e-6
        ):
            return prefixed_neutral
        return neutral_form

    def _charged_reference_atoms(
        self,
        residue_name: str,
        key_atom: str,
        charge: int,
    ) -> tuple[list[str], bool]:
        """Return the sidechain boundary atoms and whether they bypass charge groups."""
        if residue_name in self._DNA_RESIDUES:
            return [key_atom], False

        base_name = residue_name[1:] if len(residue_name) > 3 else residue_name
        entry = self._library_entry(residue_name)
        if entry["charge_groups"] or base_name == residue_name:
            markers = self._charged_heavy_atoms.get(base_name) if base_name != residue_name else None
            return self._formal_charge_group(residue_name, charge, markers), False

        # AMBER has no explicit groups. A terminal template's total charge can
        # mask an oppositely charged sidechain, so use the charged-site atoms.
        atoms = self._charged_heavy_atoms.get(base_name)
        return sorted(atoms) if atoms else [key_atom], True

    def _total_library_charge(self, df: pd.DataFrame) -> int:
        keys = ["chain_id", "residue_seq_number", "insertion_code"]
        residue_names = df.drop_duplicates(keys)["residue_name"]
        return int(round(sum(self._template_charge(name) for name in residue_names)))

    def neutralize_outside_residues(
        self,
        pdb_file: str | Path,
        salt_bridge_cutoff: float = 4.0,
    ) -> tuple[pd.DataFrame, NeutralizationStats]:
        """Find and neutralize charged residues outside the sphere boundary"""
        df = read_pdb_to_dataframe(pdb_file)
        return self.neutralize_outside_residues_dataframe(df, salt_bridge_cutoff)

    def neutralize_outside_residues_dataframe(
        self,
        df: pd.DataFrame,
        salt_bridge_cutoff: float = 4.0,
    ) -> tuple[pd.DataFrame, NeutralizationStats]:
        """Neutralize excluded charge groups and direct included salt-bridge partners."""
        reference = "switch atom" if self.switch_atoms else "charge-group center"
        logger.info(
            f"Neutralizing charged residues whose Q {reference} is at or beyond "
            f"{self.rest_bound:.1f}Å (outer {self.boundary_offset:.1f}Å SCAAS layer)"
        )
        self.stats["original_total_charge"] = self._total_library_charge(df)

        # Classify sidechain charges before terminal charges. A terminal template
        # can contain both sites (for example, CHARMM NGLU has +1 and -1 groups),
        # and removing the terminal charge first can otherwise hide its sidechain.
        charged_residues_info = self._find_charged_residues(df)
        protein_charged = [
            residue for residue in charged_residues_info if residue["residue_name"] not in self._DNA_RESIDUES
        ]
        dna_count = len(charged_residues_info) - len(protein_charged)
        self.stats["total_charged_residues"] = len(charged_residues_info)
        logger.info(
            f"Found {len(charged_residues_info)} charged residues "
            f"({len(protein_charged)} protein, {dna_count} DNA)"
        )

        outside_residues, inside_residues = self._classify_residues_by_distance(protein_charged)
        self.stats["residues_outside_boundary"] = len(outside_residues)
        salt_bridge_pairs = self._find_salt_bridges(outside_residues, inside_residues, salt_bridge_cutoff)
        self.stats["salt_bridges_neutralized"] = len(salt_bridge_pairs)

        modified_df = self._modify_residues(df, outside_residues + salt_bridge_pairs)
        modified_df = self._neutralize_nterminals(modified_df)
        modified_df = self._neutralize_cterminals(modified_df)

        # DNA uses its separate full-sphere truncation policy.
        dna_to_neutralize = self._find_outside_dna(df)
        if dna_to_neutralize:
            modified_df = self._remove_and_cap_outside_dna(modified_df, dna_to_neutralize)

        # A multiply charged terminal residue can undergo both a sidechain and a
        # terminal conversion; count its residue identity only once.
        self.stats["residues_neutralized"] = len(self.stats["modifications"]) + len(dna_to_neutralize)
        self.stats["final_total_charge"] = self._total_library_charge(modified_df)

        # Log statistics
        self._log_neutralization_stats()

        return modified_df, self.stats

    def _find_charged_residues(self, df: pd.DataFrame) -> list[ResidueInfo]:
        """Find charged residues and calculate their force-field boundary references."""
        charged_residues_info = []
        group_keys = ["chain_id", "residue_seq_number", "insertion_code"]

        residue_specs = list(self.charged_residues.items())
        for charged_name, (neutral_form, key_atom, charge) in self.charged_residues.items():
            if charged_name in self._DNA_RESIDUES:
                continue
            for prefix in ("N", "C", "n", "c"):
                terminal_name = f"{prefix}{charged_name}"
                if terminal_name in self.residue_library:
                    terminal_neutral = self._terminal_sidechain_neutral_form(
                        terminal_name, neutral_form, charge
                    )
                    residue_specs.append((terminal_name, (terminal_neutral, key_atom, charge)))

        for res_name, (neutral_form, key_atom, charge) in residue_specs:
            residues = df[df["residue_name"] == res_name]
            for (chain, res_num, insertion_code), group in residues.groupby(group_keys, dropna=False):
                reference_atoms, uses_charged_site = self._charged_reference_atoms(res_name, key_atom, charge)
                reference_coords, reference_kind = self._boundary_reference(
                    group,
                    reference_atoms,
                    f"{res_name} {chain}:{res_num}{insertion_code}",
                )
                if uses_charged_site:
                    reference_kind = "charged-site center"
                distance = float(np.linalg.norm(reference_coords - self.center))

                base_name = res_name[1:] if len(res_name) > 3 else res_name
                charged_atom_names = self._charged_heavy_atoms.get(base_name, {key_atom})
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
                        "charged_atom_coords": charged_atom_coords,
                        "distance": distance,
                        "charge": charge,
                        "neutral_form": neutral_form,
                    }
                )

        return charged_residues_info

    def _classify_residues_by_distance(
        self,
        charged_residues_info: Sequence[ResidueInfo],
    ) -> tuple[list[ResidueInfo], list[ResidueInfo]]:
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

    def _find_salt_bridges(
        self,
        outside_residues: Sequence[ResidueInfo],
        inside_residues: Sequence[ResidueInfo],
        cutoff: float,
    ) -> list[ResidueInfo]:
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

    @staticmethod
    def _residue_mask(
        df: pd.DataFrame,
        chain: str,
        residue_number: int,
        insertion_code: str,
    ) -> pd.Series:
        return (
            (df["chain_id"] == chain)
            & (df["residue_seq_number"] == residue_number)
            & (df["insertion_code"] == insertion_code)
        )

    def _record_modification(
        self,
        chain: str,
        residue_number: int,
        insertion_code: str,
        old_name: str,
        new_name: str,
        distance: float,
        boundary_reference: str,
        atoms_removed: list[str],
    ) -> None:
        identity = f"{chain}:{residue_number}{insertion_code}"
        existing = self.stats["modifications"].get(identity)
        if existing is None:
            self.stats["modifications"][identity] = {
                "original": old_name,
                "modified": new_name,
                "distance": distance,
                "boundary_reference": boundary_reference,
                "atoms_removed": atoms_removed,
            }
            return

        existing["modified"] = new_name
        existing["distance"] = max(existing["distance"], distance)
        if boundary_reference not in existing["boundary_reference"]:
            existing["boundary_reference"] += f"; {boundary_reference}"
        existing["atoms_removed"] = list(dict.fromkeys(existing["atoms_removed"] + atoms_removed))

    def _modify_residues(
        self,
        df: pd.DataFrame,
        residues_to_modify: Sequence[ResidueInfo],
    ) -> pd.DataFrame:
        """Modify residues to their neutral forms and remove appropriate atoms."""
        modified_df = df.copy()

        for res_info in residues_to_modify:
            chain = res_info["chain_id"]
            residue_number = res_info["residue_seq_number"]
            insertion_code = res_info["insertion_code"]
            old_name = res_info["residue_name"]
            new_name = res_info["neutral_form"]
            distance = res_info["distance"]

            # Include the insertion code because a cap such as NME 167A can
            # share its chain and sequence number with residue 167.
            mask = self._residue_mask(modified_df, chain, residue_number, insertion_code)
            aliases = self._TERMINAL_ALIASES.get(old_name[0].upper()) if len(old_name) > 3 else None
            modified_df, atoms_to_remove = self._convert_residue_template(
                modified_df,
                mask,
                new_name,
                preserved_aliases=aliases,
            )

            self._record_modification(
                chain,
                residue_number,
                insertion_code,
                old_name,
                new_name,
                distance,
                res_info["boundary_reference"],
                atoms_to_remove,
            )

            logger.debug(
                f"Neutralized {old_name} -> {new_name} at {chain}:{residue_number} "
                f"(distance: {distance:.2f}Å)"
            )
            if atoms_to_remove:
                logger.debug(f"Removed atoms: {', '.join(atoms_to_remove)}")

        return modified_df

    def _terminal_neutral_form(self, terminal_name: str, terminal_charge: float) -> str:
        """Find the internal template after removing an N- or C-terminal charge.

        Selection uses template charge and atom-set similarity, which also covers
        the historical OPLS2005 names such as NAR+, NLYP, CAR+, and CLYP.
        """
        source_atoms = set(self._library_atom_names(terminal_name))
        target_charge = self._template_charge(terminal_name) - terminal_charge
        direct_name = terminal_name[1:]
        if (
            direct_name in self._INTERNAL_PROTEIN_RESIDUES
            and direct_name in self.residue_library
            and abs(self._template_charge(direct_name) - target_charge) <= 0.25
        ):
            return direct_name

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

    def _terminal_charge_group(
        self,
        residue_name: str,
        marker_atoms: Sequence[str],
    ) -> list[str]:
        matching = [
            group
            for group in self._charge_groups(residue_name)
            if any(marker in group for marker in marker_atoms)
        ]
        if len(matching) != 1:
            markers = ", ".join(marker_atoms)
            raise ValueError(
                f"Expected one {residue_name} charge group containing one of "
                f"{markers}, found {len(matching)}"
            )
        return matching[0]

    def _convert_residue_template(
        self,
        df: pd.DataFrame,
        mask: pd.Series,
        new_name: str,
        preserved_aliases: dict[str, Sequence[str]] | None = None,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Rename a residue and remove atoms absent from its new force-field template."""
        new_atoms = set(self._library_atom_names(new_name))
        present_atoms = list(df.loc[mask, "atom_name"])
        atoms_to_remove = [name for name in present_atoms if name not in new_atoms]

        for target, sources in (preserved_aliases or {}).items():
            if target not in new_atoms or target in present_atoms:
                continue
            source = next((name for name in sources if name in atoms_to_remove), None)
            if source is not None:
                df.loc[mask & (df["atom_name"] == source), "atom_name"] = target
                atoms_to_remove.remove(source)

        df.loc[mask, "residue_name"] = new_name
        if atoms_to_remove:
            remove_mask = mask & df["atom_name"].isin(atoms_to_remove)
            df = df.drop(df.index[remove_mask])
        return df, atoms_to_remove

    def _neutralize_terminals(
        self,
        df: pd.DataFrame,
        *,
        prefix: str,
        marker_atoms: Sequence[str],
        terminal_charge: float,
        preserved_aliases: dict[str, Sequence[str]] | None = None,
    ) -> pd.DataFrame:
        """Convert one kind of charged terminal template to its internal form."""
        modified_df = df.copy()
        candidate_names = {
            name
            for name in modified_df["residue_name"].unique()
            if name.startswith(prefix) and len(name) > 3 and name in self.residue_library
        }
        terminal_mask = modified_df["residue_name"].isin(candidate_names)
        group_keys = ["chain_id", "residue_seq_number", "insertion_code", "residue_name"]
        converted = 0

        for (chain, residue_number, insertion_code, terminal_name), group in modified_df[
            terminal_mask
        ].groupby(group_keys, dropna=False):
            reference_atoms = self._terminal_charge_group(terminal_name, marker_atoms)
            reference_coords, reference_kind = self._boundary_reference(
                group,
                reference_atoms,
                f"{terminal_name} {chain}:{residue_number}{insertion_code}",
            )
            distance = float(np.linalg.norm(reference_coords - self.center))
            if distance < self.rest_bound:
                continue

            internal_name = self._terminal_neutral_form(terminal_name, terminal_charge)
            mask = self._residue_mask(modified_df, chain, residue_number, insertion_code)
            modified_df, atoms_removed = self._convert_residue_template(
                modified_df,
                mask,
                internal_name,
                preserved_aliases=preserved_aliases,
            )
            self._record_modification(
                chain,
                residue_number,
                insertion_code,
                terminal_name,
                internal_name,
                distance,
                reference_kind,
                atoms_removed,
            )
            logger.debug(
                f"Neutralized {prefix}-terminal {terminal_name} -> {internal_name} "
                f"at {chain}:{residue_number} ({reference_kind}: {distance:.2f}Å)"
            )
            converted += 1

        self.stats["terminals_neutralized"] += converted
        if converted:
            logger.info(f"Neutralized {converted} {prefix}-terminal residues outside boundary")
        return modified_df

    def _neutralize_nterminals(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert charged N-terminal templates outside the boundary to internal forms."""
        return self._neutralize_terminals(
            df,
            prefix="N",
            marker_atoms=("N",),
            terminal_charge=+1.0,
            preserved_aliases=self._TERMINAL_ALIASES["N"],
        )

    def _neutralize_cterminals(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert charged C-terminal templates outside the boundary to internal forms."""
        return self._neutralize_terminals(
            df,
            prefix="C",
            marker_atoms=("OXT", "OT2"),
            terminal_charge=-1.0,
            preserved_aliases=self._TERMINAL_ALIASES["C"],
        )

    def _find_outside_dna(self, df: pd.DataFrame) -> list[ResidueInfo]:
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
    _DNA_3PRIME = {"DA": "DA3", "DC": "DC3", "DG": "DG3", "DT": "DT3"}

    def _remove_and_cap_outside_dna(
        self,
        df: pd.DataFrame,
        fully_outside_residues: Sequence[ResidueInfo],
    ) -> pd.DataFrame:
        """Remove fully-outside DNA and cap each retained strand fragment.

        A retained residue next to removed upstream DNA gets a 5' terminal
        template, while one next to removed downstream DNA gets a 3' template.
        The complementary terminal charges keep each cut fragment integral.
        """
        modified_df = df.copy()
        if not fully_outside_residues:
            return modified_df

        # Group outside residues by chain
        outside_by_chain = {}
        for r in fully_outside_residues:
            outside_by_chain.setdefault(r["chain_id"], []).append(r)

        # Find all DNA residue numbers per chain (inside + outside).
        all_dna_by_chain = {}
        dna = df[df["residue_name"].isin(self._DNA_RESIDUES)]
        for (chain, residue_number), _ in dna.groupby(["chain_id", "residue_seq_number"]):
            all_dna_by_chain.setdefault(chain, set()).add(residue_number)

        # Remove all fully-outside DNA
        for r in fully_outside_residues:
            mask = (modified_df["chain_id"] == r["chain_id"]) & (
                modified_df["residue_seq_number"] == r["residue_seq_number"]
            )
            modified_df = modified_df.drop(modified_df[mask].index)
            logger.debug(f"Removed out-of-sphere DNA {r['chain_id']}:{r['residue_seq_number']}")

        n_5prime_caps = 0
        n_3prime_caps = 0
        for chain, outside_list in outside_by_chain.items():
            outside_resnums = {r["residue_seq_number"] for r in outside_list}
            ordered_resnums = sorted(all_dna_by_chain.get(chain, set()))

            for index, residue_number in enumerate(ordered_resnums):
                if residue_number in outside_resnums:
                    continue
                removed_upstream = index > 0 and ordered_resnums[index - 1] in outside_resnums
                removed_downstream = (
                    index < len(ordered_resnums) - 1 and ordered_resnums[index + 1] in outside_resnums
                )
                if not removed_upstream and not removed_downstream:
                    continue

                mask = (modified_df["chain_id"] == chain) & (
                    modified_df["residue_seq_number"] == residue_number
                )
                old_name = modified_df.loc[mask, "residue_name"].iloc[0]
                if removed_upstream and removed_downstream:
                    new_name = self.charged_residues[old_name][0]
                elif removed_upstream:
                    new_name = self._DNA_5PRIME[old_name]
                else:
                    new_name = self._DNA_3PRIME[old_name]

                modified_df.loc[mask, "residue_name"] = new_name
                if removed_upstream:
                    remove_mask = mask & modified_df["atom_name"].isin({"P", "OP1", "OP2", "HP"})
                    modified_df = modified_df.drop(modified_df.index[remove_mask])
                    n_5prime_caps += 1
                if removed_downstream:
                    n_3prime_caps += 1
                logger.debug(f"DNA boundary cap: {old_name} -> {new_name} at {chain}:{residue_number}")

        logger.info(
            f"Removed {len(fully_outside_residues)} out-of-sphere DNA nucleotides; "
            f"applied {n_5prime_caps} 5' and {n_3prime_caps} 3' terminal caps"
        )
        return modified_df

    def _get_atoms_to_remove(self, old_name: str, new_name: str) -> list[str]:
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

    def _log_neutralization_stats(self) -> None:
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


def main(
    args: argparse.Namespace | None = None,
    **kwargs: Any,
) -> NeutralizationStats | None:
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


def main_exe() -> None:
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
