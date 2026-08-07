import contextlib
import glob
import os
import re
import shutil
import stat
from pathlib import Path
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from .CLI.utils import handle_cysbonds, parse_restraint_method
from .counter_ions import minimize_coulomb_on_sphere
from .functions import COG, overlapping_pairs, sigmoid
from .IO import get_force_field_paths, qprep_error_check, replace, run_qprep
from .logger import logger
from .pdb_utils import (
    calculate_distance,
    filter_out_of_sphere_fragments,
    pdb_parse_in,
    pdb_parse_out,
    read_pdb_to_dataframe,
    reindex_pdb_residues,
    rm_HOH_clash_NN,
    shift_from_ligand_collision,
)
from .settings.settings import CLUSTER_DICT, CONFIGS
from .templates import (
    QprepFEPParameters,
    format_energy_files,
    get_equilibration_configs,
    get_neq_endpoint_config,
    get_production_config,
    render_md_input,
    render_qfep_input,
    render_qprep_fep_input,
)
from .templates.run_local import LocalRunConfig, render_local_run
from .templates.sections import (
    format_distance_restraints,
    format_sequence_restraint,
    format_water_restraint,
)

WATER_ATOM_TYPES = {
    "AMBER14sb": ("OW", "HW"),
    "OPLS2015": ("OT", "HT"),
    "CHARMM36": ("OT", "HT"),
}

TIP3P_CHARGE_O = -0.834
TIP3P_CHARGE_H = 0.417

CHLORIDE_NAME = {
    "AMBER14sb": "CHL",
    "OPLS2015": "CLA",
    "CHARMM36": "CLA",
}

# Residue name for the co-alchemical counter-water. Must not be HOH, since
# qprep treats HOH as solvent and the counter-water needs to be a solute.
COUNTER_WATER_RESNAME = "CWT"


def lrf_required_for_edge(same_charge: "bool | None") -> bool:
    """Whether this edge must run with LRF on.

    A charge-changing edge (``same_charge`` is False) changes the in-sphere net charge, whose
    long-range (~1/r) Coulomb a bare cutoff truncates — adding a large, state-dependent artifact
    that the post-hoc finite-sphere Born correction cannot remove (it is not the continuum
    ``-C*Q**2`` monopole term; it is a cutoff-truncation error). LRF (reaction field) restores the
    long-range interaction. Neutral edges (``same_charge`` True) cancel between states and are fine
    with a plain cutoff, so the configured default is kept.
    """
    return same_charge is False


class QligFEP:
    """
    Create dual topology FEP files based on two ligands
    """

    def __init__(
        self,
        lig1: str,
        lig2: str,
        FF: str,
        system: str,
        cluster: str = "TETRA",
        sphereradius: str = "25",
        cysbond: str = "none",
        start: Literal["0.0", "0.5"] = "0.0",
        temperature: str = "298",
        replicates: str = "10",
        sampling: Literal["sigmoidal", "linear", "exponential", "reverse_exponential"] = "sigmoidal",
        timestep: Literal["1fs", "2fs"] = "2fs",
        to_clean: Optional[list[str]] = None,
        water_thresh: Union[float, int] = 1.4,
        dr_force: float = 0.5,
        random_state: Optional[int] = 42,
        wath_ligand_only: bool = False,
        neq: bool = False,
        neq_reps: int = 5,
        neq_steps: int = 50000,
        neq_eq_steps: int = 1000,
        neq_relax_steps: int = 5000,
        neq_L: float = 8.0,
        neq_schedule: Literal["sigmoidal", "linear"] = "sigmoidal",
        protein_charge: Optional[int] = None,
        charge_method: str = "ion_match",
        minimize: bool = False,
    ):
        self.timestep = timestep
        self.minimize = minimize
        self.lig1 = lig1
        self.lig2 = lig2
        self.FF = FF
        self.lib_file, self.prm_file = get_force_field_paths(FF)
        self.system = system
        self.rootdir = os.getcwd()
        self.cluster = cluster
        self.sphereradius = sphereradius
        self.cysbond = cysbond
        self.start = start
        self.include = ["ATOM", "HETATM"]
        self.temperature = temperature
        self.replicates = replicates
        self.sampling = sampling
        self.to_clean = to_clean
        self.water_thresh = water_thresh
        self.dr_force = dr_force  # dr for distance restraint
        self.wath_ligand_only = wath_ligand_only
        self.neq = neq
        self.neq_reps = neq_reps
        self.neq_steps = neq_steps
        self.neq_eq_steps = neq_eq_steps
        self.neq_relax_steps = neq_relax_steps
        self.neq_L = neq_L
        self.neq_schedule = neq_schedule
        if self.neq and self.cluster in ("LOCAL", "LOCALP"):
            # NEQ setups are generated from the SLURM run_neq.sh template, which needs the
            # qdynp/qdyn engine paths that only the cluster profiles carry; there is no local NEQ runner.
            raise ValueError(
                f"--neq is not supported on cluster '{self.cluster}': the {self.cluster} profile "
                "does not define QDYN_NEQ. Use a SLURM cluster profile for non-equilibrium setups."
            )
        # Temporary until flag is here
        self.ABS = False  # True
        self.ABS_waters = []
        self.write_dir = None
        self.pdb_fname = f"{self.lig1}_{self.lig2}.pdb"
        self.seeds = self.set_seeds(random_state)
        self.n_counter_ions = 0
        # counter ion type to neutralized charge-changing perturbations (e.g., SOD or CLA)
        self.ion_type = None
        self.protein_charge = protein_charge
        valid_charge_methods = ("none", "ion_match", "coalchemical_water")
        if charge_method not in valid_charge_methods:
            raise ValueError(f"charge_method={charge_method!r} not in {valid_charge_methods}")
        self.charge_method = charge_method
        # Populated by read_files() once formal charges are known.
        self.same_charge: Optional[bool] = None
        # Co-alchemical water state, populated by place_counter_water() when
        # charge_method == "coalchemical_water". Each entry is a dict with
        # keys "topology_indices" (3 ints) and "qatoms" (3 Q-atom descriptors).
        self.counter_water_atoms: list = []
        # VdW lines for ion + water atom types, written into the FEP file.
        self.counter_ion_vdw: list = []

        if self.system == "protein":
            # Get last atom and residue from complexfile!
            txt_lines = Path("protein.pdb").read_text().splitlines()
            for line in reversed(txt_lines):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        resnr = int(line[22:26])
                        atnr = int(line[6:11])
                    except (IndexError, ValueError):
                        continue
                    break
            self.residueoffset = resnr
            self.atomoffset = atnr
            # NOTE: atomoffset is updated in `write_FEP_file` as Q might protonate
            # amino acids outside the sphere, making the residue offset incorrect.
        else:
            self.atomoffset = 0
            self.residueoffset = 0

    def set_seeds(self, random_state):
        """Set the seeds for reproduciblity"""
        if random_state is None:
            return np.random.default_rng().integers(0, 32767, size=int(self.replicates))
        rng = np.random.default_rng(random_state)
        return rng.integers(0, 32767, size=int(self.replicates))

    def makedir(self):
        lignames = f"{self.lig1}_{self.lig2}"
        directory = self.rootdir + "/FEP_" + lignames
        if not os.path.exists(directory):
            os.makedirs(directory)

        if not os.path.exists(directory + "/inputfiles"):
            os.makedirs(directory + "/inputfiles")

        self.write_dir = directory
        return directory

    def read_files(self):
        changes_1 = {}
        changes_2 = {}
        charges = []
        atomtypes = []
        merged_molsize = 0

        with open(self.lig1 + ".lib") as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 0:
                    if line[0] == "[atoms]":
                        block = 1
                        continue
                    if line[0] == "[bonds]":
                        block = 2

                if block == 1 and len(line) > 0:
                    # construct for FEP file
                    merged_molsize = merged_molsize + 1
                    charges.append([merged_molsize, line[3], "0.000"])
                    atomtypes.append([merged_molsize, line[2], "DUM"])

                if block == 2:
                    break

            molsize_lig1 = len(atomtypes)

        with open(self.lig2 + ".lib") as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 0:
                    if line[0] == "[atoms]":
                        block = 1
                        continue
                    if line[0] == "[bonds]":
                        block = 2

                if block == 1 and len(line) > 0:
                    # construct for FEP file
                    merged_molsize = merged_molsize + 1
                    charges.append([merged_molsize, "0.000", line[3]])

                    # adjustments to be made for lib and prm files
                    cnt = 0
                    for i in [line[1], line[2]]:
                        cnt = cnt + 1
                        if "AMBER14sb" in self.FF or "CHARMM36" in self.FF:
                            j = "X" + i
                        else:
                            match = re.match(r"([a-z]+)([0-9]+)", i, re.I)
                            if match:
                                items = match.groups()
                                j = str(items[0]) + str(int(items[1]) + int(molsize_lig1))

                        if cnt == 1:
                            changes_1[i] = j
                        if cnt == 2:
                            changes_2[i] = j
                            atomtypes.append([merged_molsize, "DUM", j])

        molsize_lig2 = merged_molsize - molsize_lig1

        # Compute formal charges from partial charge sums
        self.charge_lig1 = round(sum(float(c[1]) for c in charges[:molsize_lig1]))
        self.charge_lig2 = round(sum(float(c[2]) for c in charges[molsize_lig1:]))
        self.same_charge = self.charge_lig1 == self.charge_lig2

        return ([changes_1, changes_2], [charges, atomtypes], [molsize_lig1, molsize_lig2])

    def change_lib(self, replacements, writedir):
        replacements["LIG"] = "LID"
        pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")

        with open(self.lig2 + ".lib") as infile:
            file_replaced = []
            for line in infile:
                line2 = pattern.sub(lambda x: replacements[x.group()], line)
                file_replaced.append(line2)

        with open(writedir + "/" + self.lig2 + "_renumber.lib", "w") as outfile:
            for line in file_replaced:
                outfile.write(line)

        shutil.copy(self.lig1 + ".lib", writedir + "/" + self.lig1 + ".lib")

    def change_prm(self, replacements, writedir):
        pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
        file1 = glob.glob(self.lig1 + ".prm")[0]
        file2 = glob.glob(self.lig2 + ".prm")[0]
        prm_file = self.prm_file
        prm_merged = {"vdw": [], "bonds": [], "angle": [], "torsion": [], "improper": []}

        for file in [file1, file2]:
            with open(file) as infile:
                block = 0
                for line in infile:
                    if file == file2:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                    if line == "[atom_types]\n":
                        block = 1
                        continue
                    elif line == "[bonds]\n":
                        block = 2
                        continue
                    elif line == "[angles]\n":
                        block = 3
                        continue
                    elif line == "[torsions]\n":
                        block = 4
                        continue
                    if line == "[impropers]\n":
                        block = 5
                        continue
                    if block == 1:
                        prm_merged["vdw"].append(line)

                    elif block == 2:
                        prm_merged["bonds"].append(line)

                    elif block == 3:
                        prm_merged["angle"].append(line)

                    elif block == 4:
                        prm_merged["torsion"].append(line)

                    elif block == 5:
                        prm_merged["improper"].append(line)

        prm_fname = f"{writedir}/{self.FF}_{self.lig1}_{self.lig2}_merged.prm"
        with open(prm_file) as infile, open(prm_fname, "w") as outfile:
            for line in infile:
                block = 0
                outfile.write(line)
                if len(line) > 1:
                    if line == "! Ligand vdW parameters\n":
                        block = 1
                    elif line == "! Ligand bond parameters\n":
                        block = 2
                    elif line == "! Ligand angle parameters\n":
                        block = 3
                    elif line == "! Ligand torsion parameters\n":
                        block = 4
                    elif line == "! Ligand improper parameters\n":
                        block = 5
                # Read the parameters in from file and store them
                if block == 1:
                    for line in prm_merged["vdw"]:
                        outfile.write(line)

                elif block == 2:
                    for line in prm_merged["bonds"]:
                        outfile.write(line)
                elif block == 3:
                    for line in prm_merged["angle"]:
                        outfile.write(line)
                elif block == 4:
                    for line in prm_merged["torsion"]:
                        outfile.write(line)

                elif block == 5:
                    for line in prm_merged["improper"]:
                        outfile.write(line)

        # AND return the vdW list for the FEP file
        FEP_vdw = []
        for line in prm_merged["vdw"]:
            if len(line) > 1 and line[0] != "!" and line[0:1]:
                line = line.split()
                line2 = f"{line[0]:10}{line[1]:10}{line[3]:10}{str(0):10}{str(0):10}{line[4]:10}{line[5]:10}{line[6]:10}"
                FEP_vdw.append(line2)
        return FEP_vdw

    def write_FEP_file(
        self, change_charges, change_vdw, FEP_vdw, writedir, lig_size1, lig_size2, softcore_method="standard"
    ):
        lig_size1 = int(lig_size1)
        lig_size2 = int(lig_size2)
        lig_tot = lig_size1 + lig_size2
        exclude_residues = ["HOH", "LIG", "LID", COUNTER_WATER_RESNAME]
        if self.system == "water" and self.n_counter_ions > 0 and self.ion_type:
            # Method 1 (real ion in water leg) excludes the ion residue too.
            exclude_residues.append(self.ion_type)
        self.atomoffset = (
            read_pdb_to_dataframe(Path(writedir) / "top_p.pdb")
            .query("~residue_name.isin(@exclude_residues)")
            .shape[0]
        )

        # Counter-water Q-atom topology indices (3 per counter-water). Empty
        # unless place_counter_water() populated them (Method 3 only).
        cw_topo_indices = []
        if self.counter_water_atoms:
            for cw in self.counter_water_atoms:
                cw_topo_indices.extend(cw["topology_indices"])

        with open(writedir + "/FEP1.fep", "w") as outfile:
            total_atoms = len(change_charges)
            outfile.write("!info: " + self.lig1 + " --> " + self.lig2 + "\n")
            outfile.write("[FEP]\n")
            outfile.write("states 2\n")
            outfile.write("softcore_use_max_potential on\n")
            if softcore_method != "standard":
                outfile.write(f"softcore_method {softcore_method}\n")
            outfile.write("\n")

            # defining the atom order taken user given offset into account
            outfile.write("[atoms]\n")
            for i in range(1, total_atoms + 1):
                outfile.write(f"{str(i):5}{str(i + self.atomoffset):5}\n")
            for idx, topo_idx in enumerate(cw_topo_indices):
                q_idx = total_atoms + 1 + idx
                outfile.write(f"{str(q_idx):5}{str(topo_idx):5}\n")
            outfile.write("\n\n")

            # changing charges
            outfile.write("[change_charges]\n")

            for line in change_charges:
                outfile.write(f"{line[0]:<5}{line[1]:>10}{line[2]:>10}\n")
            if self.counter_water_atoms:
                q_idx = total_atoms
                for cw in self.counter_water_atoms:
                    for qatom in cw["qatoms"]:
                        q_idx += 1
                        outfile.write(f"{q_idx:<5}{qatom['charge_s1']:>10}{qatom['charge_s2']:>10}\n")
            outfile.write("\n\n")

            # add the Q atomtypes
            outfile.write("[atom_types]\n")
            for line in FEP_vdw:
                outfile.write(line + "\n")
            if self.counter_ion_vdw:
                for vdw_line in self.counter_ion_vdw:
                    outfile.write(vdw_line + "\n")

            outfile.write("DUM       0.0000    0.0000    0         0         0.0000    0.0000    1.0080")
            outfile.write("\n\n")

            outfile.write("[softcore]\n")
            # ADD softcore
            for i in range(1, lig_size1 + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i), "0", "20"))

            for i in range(1 + lig_size1, lig_tot + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i), "20", "0"))
            if self.counter_water_atoms:
                q_idx = total_atoms
                for cw in self.counter_water_atoms:
                    for qatom in cw["qatoms"]:
                        q_idx += 1
                        outfile.write(
                            "{:<5}{:>10}{:>10}\n".format(
                                str(q_idx), qatom["softcore_s1"], qatom["softcore_s2"]
                            )
                        )

            outfile.write("\n\n")

            # changing atom types
            outfile.write("[change_atoms]\n")
            for line in change_vdw:
                outfile.write(f"{line[0]:<5}{line[1]:>10}{line[2]:>10}\n")
            if self.counter_water_atoms:
                q_idx = total_atoms
                for cw in self.counter_water_atoms:
                    for qatom in cw["qatoms"]:
                        q_idx += 1
                        outfile.write(f"{q_idx:<5}{qatom['type_s1']:>10}{qatom['type_s2']:>10}\n")

    def merge_pdbs(self, writedir):
        replacements = {}
        replacements["LIG"] = "LID"
        file_replaced = []
        occupied = []
        atnr = self.atomoffset
        with open(self.lig2 + ".pdb") as infile:
            for line in infile:
                if line.split()[0].strip() in self.include:
                    atom1 = pdb_parse_in(line)
                    atom1[4] = "LID"
                    occupied.append((atom1[8], atom1[9], atom1[10]))
                    file_replaced.append(pdb_parse_out(atom1) + "\n")

        protein_contents = ""
        if self.system == "protein":
            with open("protein.pdb") as protfile:
                protein_contents = protfile.read()
            occupied += [
                (a[8], a[9], a[10])
                for a in map(pdb_parse_in, protein_contents.splitlines())
                if isinstance(a, list)
            ]

        with open(f"{self.lig1}.pdb") as ligfile:
            lig1_atoms = [a for a in map(pdb_parse_in, ligfile) if isinstance(a, list)]

        # Lig1 (LIG) is translated as a rigid body so no atom lands on a lig2 (LID) or
        # on other alchemical atom. Q crashes on atoms occupying the exact same space.
        # The previous solution was to apply a fixed 0.001 offset, but this can cause
        # ligand clashes in a few rare cases. This solution is the same, but increments
        # are checked in a while loop so that clashes are fully avoided.
        shift = shift_from_ligand_collision([(a[8], a[9], a[10]) for a in lig1_atoms], occupied)

        with open(f"{writedir}/{self.pdb_fname}", "w") as outfile:
            if protein_contents:
                outfile.write(protein_contents)
                if not protein_contents.endswith("\n"):
                    outfile.write("\n")
            for atom1 in lig1_atoms:
                atnr += 1  # The atoms are not allowed to overlap in Q
                atom1[1] = atom1[1] + self.atomoffset
                atom1[6] = atom1[6] + self.residueoffset
                atom1[8] = atom1[8] + shift
                atom1[9] = atom1[9] + shift
                atom1[10] = atom1[10] + shift
                outfile.write(pdb_parse_out(atom1) + "\n")

            self.residueoffset = self.residueoffset + 2
            resnr = f"{self.residueoffset:4}"
            for line in file_replaced:
                atnr = atnr + 1
                atchange = f"{atnr:5}"
                line = line[0:6] + atchange + line[11:22] + resnr + line[26:]
                outfile.write(line)

    def place_counter_ions(self, writedir: str) -> int:
        """Place counter-ions in the combined PDB to neutralize charge differences.

        Only applicable for water system when ligands have different formal charges.
        Ions are placed on a sphere (radius = sphereradius - 11) centered on the
        ligand COG to maximize separation via Coulomb minimization.

        When protein_charge is provided (via setupFEP), the water leg is set up to
        match the protein leg's total system charge so that electrostatic artifacts
        from the charged sphere cancel in the ddG calculation.

        Args:
            writedir: Directory containing the combined PDB (inputfiles subdirectory).

        Returns:
            Number of counter-ions placed.
        """
        delta_q = self.charge_lig1 - self.charge_lig2
        if delta_q == 0 or self.system != "water":
            return 0
        if self.charge_method != "ion_match":
            return 0

        # Temporary: in the future let's have standard 2-letter code for ions and not
        # these made up names for Q because it demands for 3-letter residue names...
        CHLORIDE_NAME = {
            "AMBER14sb": "CHL",
            "OPLS2015": "CLA",
            "CHARMM36": "CLA",
        }
        chloride = CHLORIDE_NAME.get(self.FF, "CLA")

        if self.protein_charge is not None:
            water_charge = self.charge_lig1 + self.charge_lig2
            ions_charge = self.protein_charge - water_charge
            n_ions = abs(ions_charge)
            self.ion_type = chloride if ions_charge < 0 else "SOD"
            logger.info(
                f"Matching protein leg charge {self.protein_charge}: "
                f"placing {n_ions} {self.ion_type} ion(s)"
            )
        else:
            n_ions = abs(delta_q)
            self.ion_type = chloride if delta_q > 0 else "SOD"

        self.n_counter_ions = n_ions

        cog = np.array(COG(self.lig1 + ".pdb"))
        ion_radius = int(self.sphereradius) - 11
        ion_xyz = minimize_coulomb_on_sphere(n_ions, ion_radius, cog, seed=42)

        # Read existing PDB to get last atom serial and residue number
        pdb_path = Path(writedir) / self.pdb_fname
        pdb_text = pdb_path.read_text()
        last_atnr = 0
        last_resnr = 0
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    last_atnr = max(last_atnr, int(line[6:11]))
                    last_resnr = max(last_resnr, int(line[22:26]))
                except (ValueError, IndexError):
                    continue

        # Append ion ATOM lines to the combined PDB
        with open(pdb_path, "a") as outfile:
            for i in range(n_ions):
                atnr = last_atnr + 1 + i
                resnr = last_resnr + 1 + i
                x, y, z = ion_xyz[i]
                ion_name = self.ion_type
                line = (
                    f"ATOM  {atnr:5d}  {ion_name:<3s} {ion_name:<4s} {resnr:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00\n"
                )
                outfile.write(line)

        return n_ions

    def place_counter_water(self, writedir: str) -> int:
        """Co-alchemical counter-water placement (Method 3).

        Appends one or more water residues (CWT) to the merged PDB as solute
        atoms. In the FEP, each counter-water oxygen perturbs between a real
        water oxygen and a Na+/Cl- ion, with its two hydrogens perturbing
        between real H and DUM. The net effect is a charge-conserving swap
        across the two FEP endpoints, with no softcore needed (both ends have
        finite VdW for O <-> ion, and the H <-> DUM pair is zero VdW on both
        sides). Must be called BEFORE qprep so the waters become solute atoms.

        Args:
            writedir: Inputfiles directory containing the merged PDB.

        Returns:
            Number of counter-waters placed.
        """
        delta_q = self.charge_lig1 - self.charge_lig2
        if delta_q == 0:
            return 0
        if self.charge_method != "coalchemical_water":
            return 0

        chloride = CHLORIDE_NAME.get(self.FF, "CLA")
        if self.FF not in WATER_ATOM_TYPES:
            raise ValueError(
                f"coalchemical_water mode does not have water atom types defined for FF={self.FF}; "
                f"supported FFs: {sorted(WATER_ATOM_TYPES)}"
            )
        water_o, water_h = WATER_ATOM_TYPES[self.FF]
        n_ions = abs(delta_q)

        # Choose state assignment that minimizes absolute system charge at each
        # endpoint. If lig2 has smaller |q|, the ion is "real" in state 1 (lig1
        # endpoint); otherwise it's real in state 2.
        ion_real_in_state1 = abs(self.charge_lig2) <= abs(self.charge_lig1)
        if ion_real_in_state1:
            ion_real_charge = (self.charge_lig2 - self.charge_lig1) / n_ions
        else:
            ion_real_charge = (self.charge_lig1 - self.charge_lig2) / n_ions
        self.ion_type = chloride if ion_real_charge < 0 else "SOD"
        ion_real_charge = round(ion_real_charge)

        if ion_real_in_state1:
            o_qatom = {
                "charge_s1": float(ion_real_charge),
                "charge_s2": TIP3P_CHARGE_O,
                "type_s1": self.ion_type,
                "type_s2": water_o,
                "softcore_s1": "0",
                "softcore_s2": "0",
            }
            h_qatom = {
                "charge_s1": 0.0,
                "charge_s2": TIP3P_CHARGE_H,
                "type_s1": "DUM",
                "type_s2": water_h,
                "softcore_s1": "0",
                "softcore_s2": "0",
            }
        else:
            o_qatom = {
                "charge_s1": TIP3P_CHARGE_O,
                "charge_s2": float(ion_real_charge),
                "type_s1": water_o,
                "type_s2": self.ion_type,
                "softcore_s1": "0",
                "softcore_s2": "0",
            }
            h_qatom = {
                "charge_s1": TIP3P_CHARGE_H,
                "charge_s2": 0.0,
                "type_s1": water_h,
                "type_s2": "DUM",
                "softcore_s1": "0",
                "softcore_s2": "0",
            }

        vdw_lines = [self._read_vdw_from_ff(self.ion_type), self._read_vdw_from_ff(water_o)]

        with contextlib.suppress(ValueError):
            # Some FFs omit HW from [atom_types] because its VdW is zero.
            vdw_lines.append(self._read_vdw_from_ff(water_h))

        pdb_path = Path(writedir) / self.pdb_fname
        pdb_df = read_pdb_to_dataframe(pdb_path)
        lig_df = pdb_df[pdb_df["residue_name"] == "LIG"]
        lig_cog = np.array([lig_df["x"].mean(), lig_df["y"].mean(), lig_df["z"].mean()])

        sphere_r = int(self.sphereradius)
        target_dist = max(sphere_r - 11.0, 10.0)
        if sphere_r - 5.0 - 10.0 < 3.0:
            logger.warning(
                f"Sphere radius {sphere_r}A is narrow for counter-water placement: "
                f"valid distance band is only {sphere_r - 5.0 - 10.0:.1f}A wide. "
                f"Consider a larger sphere for charge-changing perturbations."
            )

        placement_dir = self._counter_water_direction(pdb_df, lig_cog)

        positions = []
        for i in range(n_ions):
            if i == 0:
                positions.append(lig_cog + placement_dir * target_dist)
                continue
            angle = 2.0 * np.pi * i / n_ions
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            perp = np.cross(placement_dir, [0, 0, 1])
            if np.linalg.norm(perp) < 1e-6:
                perp = np.cross(placement_dir, [0, 1, 0])
            perp = perp / np.linalg.norm(perp)
            rotated = placement_dir * cos_a + perp * sin_a
            rotated = rotated / np.linalg.norm(rotated)
            positions.append(lig_cog + rotated * target_dist)

        last_serial = int(pdb_df["atom_serial_number"].max())
        last_resseq = int(pdb_df["residue_seq_number"].max())
        resname = COUNTER_WATER_RESNAME

        def _pdb_atom(serial, name, resseq, pos, element):
            return pdb_parse_out(
                [
                    "ATOM  ",
                    serial,
                    name,
                    " ",
                    resname,
                    " ",
                    resseq,
                    " ",
                    *pos,
                    1.00,
                    0.00,
                    element,
                    "  ",
                ]
            )

        self.counter_water_atoms = []
        with open(pdb_path, "a") as fh:
            for i, o_pos in enumerate(positions):
                resseq = last_resseq + 1 + i
                o_serial = last_serial + 1 + i * 3
                h1_serial = o_serial + 1
                h2_serial = o_serial + 2
                h1_pos, h2_pos = self._tip3p_hydrogen_positions(o_pos, placement_dir)
                fh.write(_pdb_atom(o_serial, "O", resseq, o_pos, " O") + "\n")
                fh.write(_pdb_atom(h1_serial, "H1", resseq, h1_pos, " H") + "\n")
                fh.write(_pdb_atom(h2_serial, "H2", resseq, h2_pos, " H") + "\n")
                self.counter_water_atoms.append(
                    {
                        "topology_indices": [o_serial, h1_serial, h2_serial],
                        "qatoms": [o_qatom, h_qatom.copy(), h_qatom.copy()],
                    }
                )

        self.n_counter_ions = n_ions
        self.counter_ion_vdw = vdw_lines
        logger.info(
            f"Placed {n_ions} counter-water(s) for {self.ion_type} swap "
            f"({o_qatom['type_s1']} -> {o_qatom['type_s2']})"
        )
        return n_ions

    def update_counter_water_indices(self, writedir: str) -> None:
        """Refresh counter-water topology indices after qprep renumbers atoms.

        qprep emits a fresh top_p.pdb whose serial numbers may differ from the
        ones we wrote into the merged PDB. We locate the CWT residues by name
        and overwrite each entry's topology_indices.
        """
        if not self.counter_water_atoms:
            return
        pdb_df = read_pdb_to_dataframe(Path(writedir) / "top_p.pdb")
        cwt_df = pdb_df[pdb_df["residue_name"] == COUNTER_WATER_RESNAME]
        if len(cwt_df) == 0:
            raise ValueError(
                f"No {COUNTER_WATER_RESNAME} residues found in top_p.pdb after qprep. "
                "Check that cwt.lib is loaded in qprep.inp."
            )
        cwt_resseqs = sorted(cwt_df["residue_seq_number"].unique())
        if len(cwt_resseqs) != len(self.counter_water_atoms):
            raise ValueError(
                f"Expected {len(self.counter_water_atoms)} CWT residues in top_p.pdb, "
                f"found {len(cwt_resseqs)}"
            )
        for i, resseq in enumerate(cwt_resseqs):
            res_atoms = cwt_df[cwt_df["residue_seq_number"] == resseq].sort_values("atom_serial_number")
            topo_indices = res_atoms["atom_serial_number"].astype(int).tolist()
            self.counter_water_atoms[i]["topology_indices"] = topo_indices

    @staticmethod
    def _tip3p_hydrogen_positions(o_pos, direction):
        """Compute the two H positions for a TIP3P water given the O and a
        reference direction. TIP3P: O-H = 0.9572 A, H-O-H = 104.52 deg."""
        bond_length = 0.9572
        half_angle = np.radians(104.52 / 2.0)
        z = direction / np.linalg.norm(direction)
        perp = np.cross(z, [1, 0, 0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(z, [0, 1, 0])
        x = perp / np.linalg.norm(perp)
        h1_dir = z * np.cos(half_angle) + x * np.sin(half_angle)
        h2_dir = z * np.cos(half_angle) - x * np.sin(half_angle)
        return o_pos + bond_length * h1_dir, o_pos + bond_length * h2_dir

    def _counter_water_direction(self, pdb_df, lig_cog):
        """Pick a placement direction for the counter-water that points away
        from the protein bulk (or +x in pure water systems)."""
        non_lig_residues = {"LIG", "LID", "HOH"}
        protein_df = pdb_df[~pdb_df["residue_name"].isin(non_lig_residues)]
        if len(protein_df) > 0:
            protein_com = np.array([protein_df["x"].mean(), protein_df["y"].mean(), protein_df["z"].mean()])
            direction = lig_cog - protein_com
            norm = np.linalg.norm(direction)
            if norm > 1e-6:
                return direction / norm
        return np.array([1.0, 0.0, 0.0])

    def _read_vdw_from_ff(self, atom_type: str) -> str:
        """Pull a VdW parameter line for `atom_type` from the FF .prm file and
        format it for the FEP [atom_types] block."""
        in_atom_types = False
        with open(self.prm_file) as fh:
            for line in fh:
                if line.strip() == "[atom_types]":
                    in_atom_types = True
                    continue
                if in_atom_types and line.startswith("["):
                    break
                if in_atom_types:
                    parts = line.split()
                    if parts and parts[0] == atom_type:
                        return (
                            f"{parts[0]:10}{parts[1]:10}{parts[3]:10}"
                            f"{str(0):10}{str(0):10}"
                            f"{parts[4]:10}{parts[5]:10}{parts[6]:10}"
                        )
        raise ValueError(f"Atom type {atom_type!r} not found in {self.prm_file}")

    def _format_counter_water_restraint(self) -> str:
        """Sequence-restraint lines pinning each counter-water oxygen, applied
        across all eq and production stages. Returns "" if no counter-waters."""
        if not self.counter_water_atoms:
            return ""
        lines = []
        for cw in self.counter_water_atoms:
            o_idx = cw["topology_indices"][0]
            lines.append(format_sequence_restraint(o_idx, o_idx, force=1.0))
        return "\n".join(lines)

    def write_water_pdb(self, writedir):
        header = self.sphereradius + ".0 SPHERE\n"
        with open("water.pdb") as infile, open(writedir + "/water.pdb", "w") as outfile:
            outfile.write(header)
            for line in infile:
                if line.startswith("TITLE"):  # qprep doesn't accept titles
                    continue
                outfile.write(line)

    def get_lambdas(self, windows, sampling):
        # Constructing the lambda partition scheme
        windows = int(windows)
        step = int(windows / 2)
        lambdas = []
        lmbda_1 = []
        lmbda_2 = []
        k_dic = {"sigmoidal": -1.1, "linear": 1000, "exponential": -1.1, "reverse_exponential": 1.1}
        k = k_dic[sampling]

        if sampling == "sigmoidal":
            for i in range(0, step + 1):
                # lambda 2 + lambda 1 should be 1.0
                lmbda1 = f"{0.5 * (sigmoid(float(i) / float(step), k) + 1):.3f}"
                lmbda2 = f"{1.0 - float(lmbda1):.3f}"
                lmbda_1.append(lmbda1)
                lmbda_2.append(lmbda2)

            lmbda_2 = lmbda_2[1:]

            for i in reversed(lmbda_2):
                lambdas.append(i)

            for i in lmbda_1:
                lambdas.append(i)

        else:
            for i in range(0, windows + 1):
                lmbda = f"{sigmoid(float(i) / float(windows), k):.3f}"
                lambdas.append(lmbda)

        lambdas = lambdas[::-1]
        return lambdas

    def set_restraints(self, writedir, restraint_method, strict_check: bool = True) -> list[list[int]]:
        """Function to set the restraints for FEP. Originally, this was performed on
        overlapping atoms, but based on our observations this was changed to a more
        chemistry-aware method, implemented under `QligFEP.restraints.restraint_setter`.

        The configuration on how these restraints will be applied depend on three strings, passed into
        `method` as `{ring_compare_method}_{surround_compare_method}_{atom_max_distance}`. Alternatively, the user
        can opt for `overlap` which simply restrains atoms within 1 A from each other, or `kartograf` to use
        the package's functionality without further post-processing of the mappings.

        Explanation:
            Ring atom compare: `aromaticity`, `hibridization`, `element`. Setting the first part of the
                string as either of these, will determine how the substituents / ring atoms are treated to be
                defined as equivalent.

            Surround atom compare: `p` (permissive), `ls` (less strict), `strict`.
                Setting the second part of the string as either of these, will determine if or how the
                direct surrounding atoms to the ring strictures will be taken into account for ring equivalence.

            - Permissive: Only the ring atoms are compared.
            - Less strict: The ring atoms and their direct surroundings are compared, but element type
                is ignored.
            - Strict: The ring atoms and their direct surroundings are element-wise compared.

            Kartograf atom max distance (optional): int or float to be used by `kartograf` as the maximum distance between
                atoms to be considered for mapping. This is by default set to 0.95 A, but can be changed by passing `_0.95`,
                for example, at the end of the `restraint_method` string.

        Args:
            writedir: directory to get the input files from, e.g.: FEP_lig1_lig2/inputfiles.
            strict_check: whether to assert the atom indexes are correctly assigned.

        Returns:
            list: list of overlapping atoms.
        """
        pdbfile = writedir + f"/inputfiles/{self.lig1}_{self.lig2}.pdb"
        if self.system == "protein":
            pdb_df = read_pdb_to_dataframe(Path(pdbfile).parent / "top_p.pdb")
        elif self.system == "water":
            pdb_df = read_pdb_to_dataframe(pdbfile)
        else:
            raise ValueError(
                f"System {self.system} not supported by this distance "
                "restraint method. Please use 'protein' or 'water'."
            )
        subset_lig1 = pdb_df.query("residue_name == 'LIG'").reset_index()
        subset_lig2 = pdb_df.query("residue_name == 'LID'").reset_index()
        logger.debug(f"Subset LIG shape: {subset_lig1.shape}")
        logger.debug(f"Subset LID shape: {subset_lig2.shape}")
        if restraint_method == "overlap":
            reslist = ["LIG", "LID"]
            torestraint_list = overlapping_pairs(pdbfile, reslist)
        else:
            # These preparation dependencies are distributed through conda-forge rather
            # than PyPI. Import them only for restraint methods that need them so the
            # analysis and input-rendering modules remain usable in a pip/uv environment.
            from .restraints.restraint_setter import RestraintSetter

            parsed = parse_restraint_method(restraint_method)
            atom_max_distance = parsed.pop("kartograf_max_atom_distance", 0.95)
            parent_write_dir = Path(writedir).parent
            lig1_path = parent_write_dir / f"{self.lig1}.sdf"
            lig2_path = parent_write_dir / f"{self.lig2}.sdf"
            if not lig1_path.exists() or not lig2_path.exists():
                logger.error(
                    "Using restraint methods other than `overlap` requires the sdf of the ligands to also be in the perturbation directory."
                )
                raise FileNotFoundError(
                    f"Could not find the sdf files for the ligands in the perturbation directory: {lig1_path}, {lig2_path}"
                )
            else:
                logger.debug(f'Loading sdf for restraint calculation:\nlig1:"{lig1_path}"\nlig2"{lig2_path}"')
                rsetter = RestraintSetter(lig1_path, lig2_path, kartograf_max_atom_distance=atom_max_distance)
                restraints = rsetter.set_restraints(**parsed)
                logger.debug(f"Restraints set using {restraint_method} method. Parameters: {parsed}")
                if strict_check:  # Good to check in case sdf in directory doesn't belong to the structure
                    rdLig1 = rsetter.molA.to_rdkit()
                    rdLig2 = rsetter.molB.to_rdkit()
                    for AtomIdx_Lig1, AtomIdx_Lig2 in restraints.items():
                        try:
                            rowLig1 = subset_lig1.iloc[AtomIdx_Lig1]
                            rowLig2 = subset_lig2.iloc[AtomIdx_Lig2]
                        except IndexError:
                            qprep_out_path = Path(writedir) / "qprep.out"
                            logger.error(
                                f"Index error for atom {AtomIdx_Lig1} in LIG or {AtomIdx_Lig2} in LID. "
                                "It's likely coming from a qprep failure.\nPlease check the qprep output at "
                                f"{qprep_out_path}."
                            )
                            qprep_error_check(qprep_out_path, self.FF)
                        atom1_in_pdb = rowLig1["atom_name"].strip("1234567890")
                        atom1_in_rdkit = rdLig1.GetAtomWithIdx(AtomIdx_Lig1).GetSymbol()
                        atom2_in_pdb = rowLig2["atom_name"].strip("1234567890")
                        atom2_in_rdkit = rdLig2.GetAtomWithIdx(AtomIdx_Lig2).GetSymbol()
                        assert atom1_in_pdb == atom1_in_rdkit
                        assert atom2_in_pdb == atom2_in_rdkit
                # convert the numbers accordingly
                pdb_atoms_lig1 = subset_lig1["atom_serial_number"].values
                pdb_atoms_lig2 = subset_lig2["atom_serial_number"].values
                torestraint_list = [[pdb_atoms_lig1[k], pdb_atoms_lig2[v]] for k, v in restraints.items()]
        rest_atom_count = len(torestraint_list)
        rest_pct_lig1 = rest_atom_count / subset_lig1.shape[0]
        rest_pct_lig2 = rest_atom_count / subset_lig2.shape[0]
        if rest_pct_lig1 > 0.3 or rest_pct_lig2 > 0.3 or rest_atom_count < 6:
            logger.warning(
                f"{rest_atom_count} restraints set with `{restraint_method}` account for {rest_pct_lig1:.2%} "
                f"of lig1 and {rest_pct_lig2:.2%} of lig2. Make sure this is intendend. Too few restraints might "
                "lead to crashes or unstable perturbations."
            )
        return torestraint_list

    def _write_equilibration_files(
        self,
        writedir,
        lig_size1,
        lig_size2,
        overlapping_atoms,
        eq_lambda1,
        eq_lambda2,
        cw_restraint="",
        wall_restraints="",
    ):
        """Write the eq1-eq5 equilibration input files from the templates module.

        Shared by the windowed and non-equilibrium setups. Returns the list of written
        file names.
        """
        file_list = []
        for i, eq_config in enumerate(
            get_equilibration_configs(
                self.timestep,
                int(self.sphereradius),
                minimize=self.minimize,
            )
        ):
            dr_str, seq_str = self._format_restraints_for_eq(
                overlapping_atoms, lig_size1, lig_size2, eq_config, self.dr_force
            )
            if cw_restraint:
                seq_str = f"{seq_str}\n{cw_restraint}" if seq_str else cw_restraint

            content = render_md_input(
                params=eq_config.params,
                lambda1=eq_lambda1,
                lambda2=eq_lambda2,
                trajectory_file=f"{eq_config.name}.dcd",
                final_file=f"{eq_config.name}.re",
                restart_file=f"eq{i}.re" if i > 0 else None,
                distance_restraints=dr_str,
                sequence_restraints=seq_str,
                wall_restraints=wall_restraints,
                is_eq1=(i == 0),
            )
            (Path(writedir) / f"{eq_config.name}.inp").write_text(content)
            file_list.append(f"{eq_config.name}.inp")
            logger.debug(f"Writing {eq_config.name}.inp")
        return file_list

    def _format_restraints_for_eq(
        self,
        overlapping_atoms: list,
        lig_size1: int,
        lig_size2: int,
        eq_config,
        dr_force: float,
    ) -> tuple[str, str]:
        """Format distance and sequence restraints for an equilibration stage.

        Args:
            overlapping_atoms: List of [atom1, atom2] pairs for distance restraints
            lig_size1: Number of atoms in ligand 1
            lig_size2: Number of atoms in ligand 2
            eq_config: EquilibrationConfig for this stage
            dr_force: Distance restraint force to use for eq5 (production dr_force)

        Returns:
            Tuple of (distance_restraints_str, sequence_restraints_str)
        """
        # Distance restraints: use config force, or dr_force for eq5
        force = (
            eq_config.distance_restraint_force if eq_config.distance_restraint_force is not None else dr_force
        )
        distance_restraints = format_distance_restraints([(a, b) for a, b in overlapping_atoms], force=force)

        # Sequence restraints depend on system type and stage
        if eq_config.use_water_restraint:
            # eq5/production: water systems use water restraint, protein systems have none
            if self.system in ("water", "vacuum"):
                sequence_restraints = format_water_restraint(
                    self.atomoffset + 1,
                    self.atomoffset + lig_size1 + lig_size2,
                    force=1.0,
                )
            else:  # protein
                sequence_restraints = ""
        else:
            # eq1-eq4: water/vacuum systems use sequence restraints, protein systems also do
            if self.system in ("water", "vacuum", "protein"):
                sequence_restraints = format_sequence_restraint(
                    self.atomoffset + 1,
                    self.atomoffset + lig_size1 + lig_size2,
                    force=eq_config.sequence_restraint_force,
                )
            else:
                sequence_restraints = ""

        return distance_restraints, sequence_restraints

    def write_md_files(
        self,
        lambdas: list[str],
        writedir: str,
        lig_size1: int,
        lig_size2: int,
        overlapping_atoms: list,
        wall_restraints: str = "",
    ) -> list[list[str]]:
        """Write equilibration and production MD input files using templates module.

        Args:
            lambdas: List of lambda values for the FEP windows
            writedir: Directory to write files to (inputfiles subdirectory)
            lig_size1: Number of atoms in ligand 1
            lig_size2: Number of atoms in ligand 2
            overlapping_atoms: List of [atom1, atom2] pairs for distance restraints
            wall_restraints: Pre-formatted wall restraints for counter-ions

        Returns:
            List of three file lists: [eq_files + md_start, forward_lambda_files, reverse_lambda_files]
        """
        file_list2 = []  # forward lambda files
        file_list3 = []  # reverse lambda files

        # Sequence restraint pinning the counter-water oxygen(s), if any.
        # Appended to every eq and production stage's sequence_restraints.
        cw_restraint = self._format_counter_water_restraint()

        # Determine equilibration lambdas based on start mode
        if self.start == "0.5":
            eq_lambda1, eq_lambda2 = "0.500", "0.500"
        else:  # start == "0.0" or "1"
            eq_lambda1, eq_lambda2 = "1.000", "0.000"

        prod_config = get_production_config(self.timestep, int(self.sphereradius), self.dr_force)

        # eq files + initial md file
        file_list1 = self._write_equilibration_files(
            writedir,
            lig_size1,
            lig_size2,
            overlapping_atoms,
            eq_lambda1,
            eq_lambda2,
            cw_restraint,
            wall_restraints,
        )

        # Format restraints for production files
        prod_dr_str = format_distance_restraints([(a, b) for a, b in overlapping_atoms], force=self.dr_force)
        if self.system in ("water", "vacuum"):
            prod_seq_str = format_water_restraint(
                self.atomoffset + 1,
                self.atomoffset + lig_size1 + lig_size2,
                force=1.0,
            )
        else:  # protein
            prod_seq_str = ""
        if cw_restraint:
            prod_seq_str = f"{prod_seq_str}\n{cw_restraint}" if prod_seq_str else cw_restraint

        # Write production files based on start mode
        if self.start == "0.5":
            file_list1, file_list2, file_list3 = self._write_production_05(
                lambdas, writedir, prod_config, prod_dr_str, prod_seq_str, file_list1, wall_restraints
            )
        else:
            file_list1, file_list2, file_list3 = self._write_production_1(
                lambdas, writedir, prod_config, prod_dr_str, prod_seq_str, file_list1, wall_restraints
            )

        return [file_list1, file_list2, file_list3]

    def _write_production_05(
        self,
        lambdas: list[str],
        writedir: str,
        prod_config,
        dr_str: str,
        seq_str: str,
        file_list1: list[str],
        wall_str: str = "",
    ) -> tuple[list[str], list[str], list[str]]:
        """Write production MD files for start=0.5 mode.

        Starts from 0.500/0.500, then branches in both directions.
        """
        file_list2 = []
        file_list3 = []

        # Split lambdas into forward and reverse from center
        lambda_1 = []
        lambda_2 = []
        block = 0
        for lmbda in lambdas:
            if lmbda == "0.500":
                block = 1
            if block == 0:
                lambda_1.append(lmbda)
            if block == 1:
                lambda_2.append(lmbda)

        lambda_1 = lambda_1[::-1]
        lambda_2 = lambda_2[1:]

        # Write initial md_0500_0500.inp
        content = render_md_input(
            params=prod_config.params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="md_0500_0500.dcd",
            final_file="md_0500_0500.re",
            restart_file="eq5.re",
            energy_file="md_0500_0500.en",
            distance_restraints=dr_str,
            sequence_restraints=seq_str,
            wall_restraints=wall_str,
        )
        filepath = Path(writedir) / "md_0500_0500.inp"
        filepath.write_text(content)
        file_list1.append("md_0500_0500.inp")

        # Write lambda window files
        for index, lambda_list in enumerate([lambda_1, lambda_2], start=1):
            filename_N = "md_0500_0500"

            for filenr, _ in enumerate(lambda_list):
                if index == 1:
                    l1 = lambda_1[filenr]
                    l2 = lambda_2[filenr]
                else:
                    l1 = lambda_2[filenr]
                    l2 = lambda_1[filenr]

                filename = f"md_{l1.replace('.', '')}_{l2.replace('.', '')}"

                content = render_md_input(
                    params=prod_config.params,
                    lambda1=l1,
                    lambda2=l2,
                    trajectory_file=f"{filename}.dcd",
                    final_file=f"{filename}.re",
                    restart_file=f"{filename_N}.re",
                    energy_file=f"{filename}.en",
                    distance_restraints=dr_str,
                    sequence_restraints=seq_str,
                    wall_restraints=wall_str,
                )
                filepath = Path(writedir) / f"{filename}.inp"
                filepath.write_text(content)

                filename_N = filename

                if index == 1:
                    file_list2.append(f"{filename}.inp")
                else:
                    file_list3.append(f"{filename}.inp")

        return file_list1, file_list2, file_list3

    def _write_production_1(
        self,
        lambdas: list[str],
        writedir: str,
        prod_config,
        dr_str: str,
        seq_str: str,
        file_list1: list[str],
        wall_str: str = "",
    ) -> tuple[list[str], list[str], list[str]]:
        """Write production MD files for start=0.0/1 mode.

        Starts from 1.000/0.000 and goes through all windows sequentially.
        """
        file_list2 = []
        file_list3 = []
        total_lambda = len(lambdas)

        # Write initial md_1000_0000.inp
        content = render_md_input(
            params=prod_config.params,
            lambda1="1.000",
            lambda2="0.000",
            trajectory_file="md_1000_0000.dcd",
            final_file="md_1000_0000.re",
            restart_file="eq5.re",
            energy_file="md_1000_0000.en",
            distance_restraints=dr_str,
            sequence_restraints=seq_str,
            wall_restraints=wall_str,
        )
        filepath = Path(writedir) / "md_1000_0000.inp"
        filepath.write_text(content)
        file_list1.append("md_1000_0000.inp")

        filename_N = "md_1000_0000"
        filenr = 0

        for lmbda in lambdas:
            if lmbda == "1.000":
                continue

            step_n = total_lambda - filenr - 2
            l1 = lmbda
            l2 = lambdas[step_n]

            filename = f"md_{l1.replace('.', '')}_{l2.replace('.', '')}"

            content = render_md_input(
                params=prod_config.params,
                lambda1=l1,
                lambda2=l2,
                trajectory_file=f"{filename}.dcd",
                final_file=f"{filename}.re",
                restart_file=f"{filename_N}.re",
                energy_file=f"{filename}.en",
                distance_restraints=dr_str,
                sequence_restraints=seq_str,
                wall_restraints=wall_str,
            )
            filepath = Path(writedir) / f"{filename}.inp"
            filepath.write_text(content)

            filename_N = filename
            filenr += 1
            file_list2.append(f"{filename}.inp")

        return file_list1, file_list2, file_list3

    def _write_endpoint_file(
        self, file_out, lambda1, lambda2, steps, distance_restraints, sequence_restraints, lambda_scaling=None
    ):
        """Write a single non-equilibrium endpoint MD input file via render_md_input.

        Used for both the endpoint equilibration files (eq6_*) and the lambda-switching
        files (neq_*). The per-replicate restart/final file names and the temperature stay
        as RESTARTFILE/FINALFILE/T_VAR placeholders that the run script fills in. When
        lambda_scaling is given it becomes the [lambda_scaling] section that drives lambda
        over the run (a switching run); without it the file is a plain equilibrium MD at
        the endpoint.
        """
        params = get_neq_endpoint_config(self.timestep, int(self.sphereradius), steps)
        content = render_md_input(
            params=params,
            lambda1=lambda1,
            lambda2=lambda2,
            trajectory_file="neq.dcd",
            final_file="FINALFILE",
            restart_file="RESTARTFILE",
            distance_restraints=distance_restraints,
            sequence_restraints=sequence_restraints,
            lambda_scaling=lambda_scaling,
        )
        Path(file_out).write_text(content)

    def write_MD_neq(self, writedir, lig_size1, lig_size2, overlapping_atoms):
        """Write the non-equilibrium input files: eq1-5 (equilibration), relax_{0,1}
        (one-time endpoint relaxation), eq6_{0,1} (endpoint equilibration spacing) and
        neq_{0,1} (lambda switching). State 0 starts at lambda 0.0->1.0 and state 1 at
        1.0->0.0; the run script relaxes each endpoint once, then chains the restarts and
        repeats the switches per replicate. Returns the list of written file names.
        """
        file_list = self._write_equilibration_files(
            writedir, lig_size1, lig_size2, overlapping_atoms, "0.500", "0.500"
        )

        dr_str = format_distance_restraints([(a, b) for a, b in overlapping_atoms], force=self.dr_force)
        if self.system in ("water", "vacuum"):
            seq_str = format_water_restraint(
                self.atomoffset + 1, self.atomoffset + lig_size1 + lig_size2, force=1.0
            )
        else:  # protein
            seq_str = ""

        lambda_scaling = f"scaling_parameter          {self.neq_schedule}\nL_sigmoid        {self.neq_L}"
        endpoint_lambdas = {"0": ("0.000", "1.000"), "1": ("1.000", "0.000")}
        for state, (lambda1, lambda2) in endpoint_lambdas.items():
            self._write_endpoint_file(
                writedir + f"/relax_{state}.inp", lambda1, lambda2, self.neq_relax_steps, dr_str, seq_str
            )
            self._write_endpoint_file(
                writedir + f"/eq6_{state}.inp", lambda1, lambda2, self.neq_eq_steps, dr_str, seq_str
            )
            self._write_endpoint_file(
                writedir + f"/neq_{state}.inp",
                lambda1,
                lambda2,
                self.neq_steps,
                dr_str,
                seq_str,
                lambda_scaling,
            )
            file_list += [f"relax_{state}.inp", f"eq6_{state}.inp", f"neq_{state}.inp"]
        return file_list

    def write_neq_runfile(self, writedir, file_list):
        """Write the SLURM run script for a non-equilibrium FEP. Each array task runs
        eq1-5, then loops `neq_reps` iterations of a forward/reverse lambda switch pair with the
        serial qdyn engine (NEQ mode selected by the [lambda_scaling] section in neq_*.inp).

        Each iteration produces a forward and a reverse switch, so `neq_reps` iterations give
        ``2*neq_reps`` independent switches. They run concurrently, so the allocation is sized to
        hold them all in a single wave: ``max(2*neq_reps, cluster_default)``. Flooring at the
        cluster default keeps at least the billed core count (some clusters bill a fixed minimum,
        and the equilibration runs the parallel qdynp engine across all of them); growing to the
        switch count when the protocol needs more avoids a starved tail wave that would double the
        switch-phase wall time.
        """
        if self.neq_reps < 1:
            raise ValueError(f"neq_reps must be >= 1 to size the switch allocation, got {self.neq_reps}")

        src = CONFIGS["INPUT_DIR"] + "/run_neq.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"

        replacements = dict(CLUSTER_DICT[self.cluster])
        replacements["FEPS"] = "FEP1.fep"
        replacements["NEQ_REPS"] = str(self.neq_reps)
        n_switches = 2 * self.neq_reps
        ncores = max(n_switches, int(replacements["NTASKS"]))
        replacements["NTASKS"] = str(ncores)
        logger.info(f"NEQ run: {n_switches} switches per replicate on {ncores} cores (one concurrent wave).")
        with open(src) as infile, open(tgt, "w") as outfile:
            for line in infile:
                if line.strip() == "#SBATCH --array=1-TOTAL_JOBS":
                    replacements["TOTAL_JOBS"] = str(self.replicates)
                if line.strip() == "temperatures=(TEMP_VAR)":
                    replacements["TEMP_VAR"] = str(self.temperature)
                if line.strip() == "seeds=(RANDOM_SEEDS)":
                    replacements["RANDOM_SEEDS"] = " ".join([str(s) for s in self.seeds])
                if line.strip() == "#SBATCH -A ACCOUNT":
                    try:  # Try to take account info - not for all clusters!
                        replacements["ACCOUNT"]
                    except KeyError:
                        line = ""
                if line.strip() == "#SBATCH -J JOBNAME":
                    if self.cluster == "DARDEL":  # TODO: refactor this...
                        outfile.write("#SBATCH -p shared\n")
                    elif self.cluster == "SNELLIUS":
                        outfile.write("#SBATCH -p rome\n")
                    try:
                        if self.system == "water":
                            jobname = "w_"
                        elif self.system == "protein":
                            jobname = "p_"
                        elif self.system == "vacuum":
                            jobname = "v_"
                        jobname += self.lig1 + "_" + self.lig2
                        replacements["JOBNAME"] = jobname
                    except Exception as e:
                        logger.error(f"Something went wrong while defining the jobname:\n{e}")
                        line = ""
                outline = replace(line, replacements)
                if outline.strip().startswith("#SBATCH --mem-per-cpu=512") and self.cluster == "CSB":
                    continue
                outfile.write(outline)
                if line.strip() == "#CLEANUP" and self.to_clean is not None:
                    rm_line = "rm -r " + " ".join(["*" + x for x in self.to_clean]) + "\n"
                    outfile.write(rm_line)

        try:
            st = os.stat(tgt)
            os.chmod(tgt, st.st_mode | stat.S_IEXEC)
        except OSError:
            logger.warning(f"Could not change permission for {tgt}")

    def write_submitfile(self, writedir):
        if self.cluster in ("LOCAL", "LOCALP"):
            return
        replacements = {}
        replacements["RUNFILE"] = "run" + self.cluster + ".sh"
        submit_in = CONFIGS["ROOT_DIR"] + "/INPUTS/FEP_submit.sh"
        submit_out = writedir + ("/FEP_submit.sh")
        with open(submit_in) as infile, open(submit_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)

        try:
            st = os.stat(submit_out)
            os.chmod(submit_out, st.st_mode | stat.S_IEXEC)
        except OSError:
            print("WARNING: Could not change permission for " + submit_out)

    def _write_local_runfile(self, writedir, file_list):
        """Generate a local run script using the Python template."""
        cluster_config = CLUSTER_DICT[self.cluster]
        eq_basenames = [Path(f).stem for f in sorted(glob.glob(writedir + "/eq*.inp"))]

        # Collect MD basenames in execution order
        if self.start == "1":
            md_basenames = [Path(f).stem for f in reversed(sorted(glob.glob(writedir + "/md*.inp")))]
        elif self.start == "0.5":
            md_1 = file_list[1]
            md_2 = file_list[2]
            md_basenames = ["md_0500_0500"]
            for i in range(len(md_1)):
                md_basenames.append(Path(md_1[i]).stem)
                md_basenames.append(Path(md_2[i]).stem)

        config = LocalRunConfig(
            qdyn_path=cluster_config["QDYN"],
            qfep_path=cluster_config["QFEP"],
            use_mpi=cluster_config["USE_MPI"],
            ntasks=int(cluster_config["NTASKS"]),
            temperatures=self.temperature.split(","),
            seeds=self.seeds,
            eq_files=eq_basenames,
            md_files=md_basenames,
            fep_files=["FEP1.fep"],
            final_md_restart="md_0000_1000.re",
            cleanup_patterns=self.to_clean if self.to_clean else None,
        )

        script_content = render_local_run(config)
        tgt = writedir + "/run" + self.cluster + ".sh"
        with open(tgt, "w") as outfile:
            outfile.write(script_content)

        try:
            st = os.stat(tgt)
            os.chmod(tgt, st.st_mode | stat.S_IEXEC)
        except Exception:
            logger.warning(f"Could not change permission for {tgt}")

    def write_runfile(self, writedir, file_list):
        if self.cluster in ("LOCAL", "LOCALP"):
            self._write_local_runfile(writedir, file_list)
            return

        src = CONFIGS["INPUT_DIR"] + "/run.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"
        EQ_files = sorted(glob.glob(writedir + "/eq*.inp"))

        if self.start == "1":
            MD_files = reversed(sorted(glob.glob(writedir + "/md*.inp")))

        elif self.start == "0.5":
            md_1 = file_list[1]
            md_2 = file_list[2]

        # --map-by core pins one rank per core in the job's cpuset. Above two ranks, mpirun defaults to
        # spreading ranks evenly across sockets, ignoring how many cores SLURM placed on each; on an
        # uneven split (e.g. 6+2 of 8) that maps more ranks to the small socket than it has cores, and
        # --bind-to core then refuses to launch ("binding more processes than cpus on a resource").
        mpirun = "time mpirun -n $SLURM_NTASKS --map-by core --bind-to core $qdyn"

        replacements = dict(CLUSTER_DICT[self.cluster])
        replacements["FEPS"] = "FEP1.fep"
        with open(src) as infile, open(tgt, "w") as outfile:
            for line in infile:
                if line.strip() == "#SBATCH --array=1-TOTAL_JOBS":
                    replacements["TOTAL_JOBS"] = str(self.replicates)
                if line.strip() == "temperatures=(TEMP_VAR)":
                    replacements["TEMP_VAR"] = str(self.temperature)
                if line.strip() == "seeds=(RANDOM_SEEDS)":
                    replacements["RANDOM_SEEDS"] = " ".join([str(s) for s in self.seeds])
                if line.strip() == "#SBATCH -A ACCOUNT":
                    try:  # Try to take account info - not for all clusters!
                        replacements["ACCOUNT"]
                    except KeyError:
                        line = ""
                if line.strip() == "#SBATCH -J JOBNAME":
                    if self.cluster == "DARDEL":  # TODO: refactor this...
                        outfile.write("#SBATCH -p shared\n")
                    elif self.cluster == "SNELLIUS":
                        outfile.write("#SBATCH -p rome\n")
                    try:
                        if self.system == "water":
                            jobname = "w_"
                        elif self.system == "protein":
                            jobname = "p_"
                        elif self.system == "vacuum":
                            jobname = "v_"
                        jobname += self.lig1 + "_" + self.lig2
                        replacements["JOBNAME"] = jobname
                    except Exception as e:
                        logger.error(f"Something went wrong while defining the jobname:\n{e}")
                        line = ""
                outline = replace(line, replacements)
                # This configuration is not available for CSB
                if outline.strip().startswith("#SBATCH --mem-per-cpu=512") and self.cluster == "CSB":
                    continue
                outfile.write(outline)
                if line.strip() == "#EQ_FILES":
                    for line in EQ_files:
                        file_base = Path(line).stem
                        outfile.write(f"{mpirun} {file_base}.inp > {file_base}.log\n")

                if line.strip() == "#RUN_FILES":
                    if self.start == "1":
                        for line in MD_files:
                            file_base = line.split("/")[-1][:-4]
                            outfile.write(f"{mpirun} {file_base}.inp > {file_base}.log\n")

                    elif self.start == "0.5":
                        outfile.write(f"{mpirun} md_0500_0500.inp > md_0500_0500.log\n\n")
                        for md1, md2 in zip(md_1, md_2):
                            outfile.write(f"{mpirun} {md1[:-4]}.inp > {md1[:-4]}.log\n")
                            outfile.write(f"{mpirun} {md2[:-4]}.inp > {md2[:-4]}.log\n")
                            outfile.write("\n")
                if line.strip() == "#CLEANUP" and self.to_clean is not None:
                    replacements["CLEANUP"] = "#Cleaned {} files\n".format(" ".join(self.to_clean))
                    outline = replace(line, replacements)
                    rm_line = "rm -r " + " ".join(["*" + x for x in self.to_clean]) + "\n"
                    outfile.write(rm_line)
                    outfile.write(outline[1:])

    def write_qfep(self, windows, lambdas):
        qfep_out = self.write_dir + "/inputfiles/qfep.inp"
        # TO DO: multiple files will need to be written out for temperature range
        energy_files = format_energy_files(lambdas)
        content = render_qfep_input(
            total_lambdas=len(lambdas),
            temperature=float(self.temperature),
            windows=windows,
            energy_files=energy_files,
        )
        with open(qfep_out, "w") as outfile:
            outfile.write(content)

    def avoid_water_protein_clashes(self, writedir, header: Optional[str] = None, save_removed: bool = False):
        """Function to remove water molecules too close to protein & ligands | ligands (water leg).
        Thresholds are the distances in Ångström from the protein & ligands | ligands atoms
        to the nearest heavy atom in the water molecule (HOH).

        The threshold used for the removal of these clashing water molecules are defined in the
        `self.water_thresh` attribute.

        Args:
            writedir: directory in which QligFEP will write the input files.
            header: header to be added to the water.pdb file. This detail is needed in Qprep's
                current version (2024/07/21), as it doesn't recognize the radius when merging an
                existing water.pdb file to the topology.
            save_removed: whether to save the removed water molecules to a file. This is automatically
                set to True in the `CLI` in case `--log-level` is passed as either `debug` or `trace`.
        """
        waterfile = Path(writedir) / "water.pdb"
        protfile = Path(writedir) / self.pdb_fname
        threshold = self.water_thresh
        system_to_log = "protein & ligands" if self.system == "protein" else "ligands"

        logger.info(f"Removing water molecules too close to {system_to_log} - threshold: {threshold} A.")
        _, n_removed = rm_HOH_clash_NN(
            pdb_df_query=read_pdb_to_dataframe(waterfile),
            pdb_df_target=read_pdb_to_dataframe(protfile),
            th=threshold,
            output_file=waterfile,
            heavy_only=True,
            ligand_only=self.wath_ligand_only,
            header=header,
            save_removed=save_removed,
        )

    def _get_cog_from_water(self, writedir: str) -> list[float]:
        """Extract COG from water.pdb TITLE line.

        Args:
            writedir: directory where inputfiles are written (e.g., FEP_lig1_lig2/inputfiles)

        Returns:
            list[float]: [x, y, z] coordinates of center of geometry
        """
        cog_regex = re.compile(r"([-]?\d{1,3}\.\d{3})\s+([-]?\d{1,3}\.\d{3})\s+([-]?\d{1,3}\.\d{3})")
        water_pdb_path = Path(writedir).parents[1] / "water.pdb"

        if water_pdb_path.exists():
            first_lines = water_pdb_path.read_text().split("\n")[:3]
            for line in first_lines:
                if line.startswith("TITLE"):
                    matches = cog_regex.search(line)
                    if matches:
                        return [float(x) for x in matches.groups()]

        # Fallback: calculate from ligand
        logger.debug("COG not found in water.pdb TITLE, calculating from ligand.")
        return list(COG(self.lig1 + ".pdb"))

    def _update_offsets_from_pdb(self, pdb_path: Path) -> None:
        """Update atomoffset and residueoffset from filtered PDB.

        After filtering fragments from the protein PDB, the atom and residue
        numbering may have changed. This method recalculates the offsets.

        Args:
            pdb_path: Path to the (filtered) merged PDB file
        """
        pdb_df = read_pdb_to_dataframe(pdb_path)
        protein_atoms = pdb_df[~pdb_df["residue_name"].isin(["HOH", "LIG", "LID"])]
        if not protein_atoms.empty:
            self.atomoffset = int(protein_atoms["atom_serial_number"].max())
            self.residueoffset = int(protein_atoms["residue_seq_number"].max())
            logger.debug(f"Updated offsets: atomoffset={self.atomoffset}, residueoffset={self.residueoffset}")

    def filter_protein_fragments(self, writedir: str) -> tuple[int, int]:
        """Filter merged PDB to remove molecular fragments completely outside sphere.

        Uses MDAnalysis topology-based fragment detection. Entire chains, cofactors,
        lipids, or other molecular units are removed if ALL their atoms are outside
        the simulation sphere.

        Args:
            writedir: directory where inputfiles are written (e.g., FEP_lig1_lig2/inputfiles)

        Returns:
            Tuple of (original_atom_count, filtered_atom_count)
        """
        cog = self._get_cog_from_water(writedir)
        merged_pdb = Path(writedir) / self.pdb_fname

        orig_count, filt_count = filter_out_of_sphere_fragments(merged_pdb, cog, float(self.sphereradius))

        # Update offsets if structure was filtered
        if filt_count < orig_count:
            self._update_offsets_from_pdb(merged_pdb)
            logger.info(
                f"Filtered structure: {orig_count} → {filt_count} atoms ({orig_count - filt_count} removed)"
            )

        return orig_count, filt_count

    def write_qprep(self, writedir):
        """Write the qprep.inp file for Q. If the system is water, the center of geometry
        will be calculated from the lig1's atoms coordinates. If the system is protein, it
        will try to extract the center of geometry from the TITLE line in the water.pdb file,
        added by the `qprep_prot` program.

        Args:
            writedir: directory in which QligFEP will write the input files.
        """
        self.cog = self._get_cog_from_water(writedir)
        center = f"{self.cog[0]} {self.cog[1]} {self.cog[2]}"

        qprep_out = writedir + "/qprep.inp"

        # Determine solvent specification
        if self.system == "vacuum":
            solvent = ""
            solvate = False
        elif self.system == "water":
            solvent = "1 HOH"
            solvate = True
        else:  # protein
            solvent = "4 water.pdb"
            solvate = True

        pdb_df = read_pdb_to_dataframe(Path(writedir) / self.pdb_fname)
        density = self.get_sphere_density(pdb_df) if self.system == "protein" else 0.05794

        # We reindex the residues prior to defining the cysbonds because Q considers
        # the first residue to be always 1, regardless of the numbering in the PDB file.
        reindex_pdb_residues(Path(writedir) / self.pdb_fname, Path(writedir) / self.pdb_fname)
        cysbond_str = handle_cysbonds(
            self.cysbond, Path(writedir) / self.pdb_fname, comment_out=(self.system != "protein")
        )
        if self.system == "protein" and self.cysbond == "auto" and cysbond_str != "":
            # cysbond shouldn't be there if the AA is out of the sphere radius
            new_cysbond_str = ""
            for line in cysbond_str.strip().split("\n"):
                parts = line.split()
                resn1, at1 = parts[1].split(":")
                resn1 = int(resn1)
                resn2, at2 = parts[2].split(":")
                resn2 = int(resn2)

                atom1 = pdb_df.query("residue_seq_number == @resn1 & atom_name == @at1")
                atom2 = pdb_df.query("residue_seq_number == @resn2 & atom_name == @at2")
                if not atom1.empty and not atom2.empty:
                    atom1_coords = atom1[["x", "y", "z"]].values[0]
                    atom2_coords = atom2[["x", "y", "z"]].values[0]

                    dist1 = calculate_distance(atom1_coords, self.cog)
                    dist2 = calculate_distance(atom2_coords, self.cog)
                    logger.debug(f"{resn1}:{at1} and {resn2}:{at2} within {dist1} and {dist2} of the COG.")
                    if dist1 <= int(self.sphereradius) and dist2 <= int(self.sphereradius):
                        new_cysbond_str += line + "\n"
                    else:
                        logger.info(
                            f"Excluding cysbond {line}; one or both atoms are outside the sphere radius."
                        )
                else:
                    logger.warning(f"Atom information not found for bond {line}.")
            if new_cysbond_str:
                cysbond_str = new_cysbond_str

        params = QprepFEPParameters(
            ff_lib=self.lib_file,
            lig1_lib=self.lig1 + ".lib",
            lig2_lib=self.lig2 + "_renumber.lib",
            ligand_prm=self.FF + "_" + self.lig1 + "_" + self.lig2 + "_merged.prm",
            ligand_pdb=self.lig1 + "_" + self.lig2 + ".pdb",
            center=center,
            sphere_radius=self.sphereradius,
            solute_density=f"{density:.5f}",
            solvent=solvent,
            cysbonds=cysbond_str,
            solvate=solvate,
        )

        content = render_qprep_fep_input(params)
        with open(qprep_out, "w") as outfile:
            outfile.write(content)

    def get_sphere_density(self, pdb_df: pd.DataFrame) -> float:
        """Calculate the solute density for the FEP system taking into consideration the different
        densities for protein and lipids. The density is calculated as the weighted average
        of the densities of the protein and lipids, where the weights are the number of atoms
        in each.

        Args:
            pdb_df: DataFrame containing the PDB file information.

        Returns:
            float: The density of the system.
        """
        # regex with negative lookbehind to avoid matching NH, CH, OH (heavy atoms)
        H_atom_regex = re.compile(r"(?<![NCO])H\d*")  # noqa: F841
        protein_vol = 0.05794  # A**-3
        lipid_vol = 0.03431  # A**-3 from octane

        knn = NearestNeighbors(radius=float(self.sphereradius), metric="euclidean", n_jobs=4)
        atom_coord_arr = pdb_df[["x", "y", "z"]].values
        knn.fit(np.array(atom_coord_arr))
        _, indices = knn.radius_neighbors(np.array(self.cog).reshape(1, -1))
        heavy_atoms_df = (
            pdb_df.iloc[indices[0]]
            .query("residue_name != 'HOH'")
            .query("~atom_name.str.match(@H_atom_regex)")
        )
        n_lipid_at = heavy_atoms_df.query("residue_name == 'POP'").shape[0]
        if n_lipid_at == 0:
            return 0.05794
        n_protein_at = heavy_atoms_df.shape[0] - n_lipid_at

        density = (n_protein_at * protein_vol + n_lipid_at * lipid_vol) / heavy_atoms_df.shape[0]
        if density != 0.05794:
            logger.info(f"Calculated solute density: {density:.5f} g/cm^3")

        return density

    def qprep(self, writedir):
        os.chdir(writedir)
        cluster_options = CLUSTER_DICT[self.cluster]
        qprep_path = cluster_options["QPREP"]
        logger.info(f"Running QPREP from path {qprep_path}")
        run_qprep(qprep_path, "qprep.inp", "qprep.out", self.FF)
        os.chdir("../../")
