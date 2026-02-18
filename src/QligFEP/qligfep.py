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

from .CLI.utils import handle_cysbonds
from .counter_ions import minimize_coulomb_on_sphere
from .functions import COG, kT, overlapping_pairs, sigmoid
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
)
from .restraints.restraint_setter import RestraintSetter
from .scoring import parse_restraint_method
from .settings.settings import CLUSTER_DICT, CONFIGS
from .templates import (
    QprepFEPParameters,
    format_energy_files,
    get_equilibration_configs,
    get_production_config,
    render_md_input,
    render_qfep_input,
    render_qprep_fep_input,
)
from .templates.sections import (
    format_distance_restraints,
    format_sequence_restraint,
    format_wall_restraints,
    format_water_restraint,
)


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
    ):
        self.timestep = timestep
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

    def replace(self, string, replacements):
        pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
        replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
        return replaced_string

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
        exclude_residues = ["HOH", "LIG", "LID"]
        if self.system == "water" and self.n_counter_ions > 0:
            exclude_residues.append(self.ion_type)
        self.atomoffset = (
            read_pdb_to_dataframe(Path(writedir) / "top_p.pdb")
            .query("~residue_name.isin(@exclude_residues)")
            .shape[0]
        )

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
            outfile.write("\n\n")

            # changing charges
            outfile.write("[change_charges]\n")

            for line in change_charges:
                outfile.write(f"{line[0]:<5}{line[1]:>10}{line[2]:>10}\n")
            outfile.write("\n\n")

            # add the Q atomtypes
            outfile.write("[atom_types]\n")
            for line in FEP_vdw:
                outfile.write(line + "\n")

            outfile.write("DUM       0.0000    0.0000    0         0         0.0000    0.0000    1.0080")
            outfile.write("\n\n")

            outfile.write("[softcore]\n")
            # ADD softcore
            for i in range(1, lig_size1 + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i), "0", "20"))

            for i in range(1 + lig_size1, lig_tot + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i), "20", "0"))

            outfile.write("\n\n")

            # changing atom types
            outfile.write("[change_atoms]\n")
            for line in change_vdw:
                outfile.write(f"{line[0]:<5}{line[1]:>10}{line[2]:>10}\n")

    def merge_pdbs(self, writedir):
        replacements = {}
        replacements["LIG"] = "LID"
        file_replaced = []
        atnr = self.atomoffset
        with open(self.lig2 + ".pdb") as infile:
            for line in infile:
                if line.split()[0].strip() in self.include:
                    atom1 = pdb_parse_in(line)
                    atom1[4] = "LID"
                    line = pdb_parse_out(atom1) + "\n"
                    file_replaced.append(line)

        with open(f"{self.lig1}.pdb") as infile, open(f"{writedir}/{self.pdb_fname}", "w") as outfile:
            if self.system == "protein":
                with open("protein.pdb") as protfile:
                    contents = protfile.read()
                    outfile.write(contents)
                    if contents and not contents.endswith("\n"):
                        outfile.write("\n")
            for line in infile:
                if line.split()[0].strip() in self.include:
                    resnr = int(line[22:26])
                    atnr += 1  # The atoms are not allowed to overlap in Q
                    atom1 = pdb_parse_in(line)
                    atom1[1] = atom1[1] + self.atomoffset
                    atom1[6] = atom1[6] + self.residueoffset
                    atom1[8] = float(atom1[8]) + 0.001
                    atom1[9] = float(atom1[9]) + 0.001
                    atom1[10] = float(atom1[10]) + 0.001
                    line = pdb_parse_out(atom1) + "\n"
                    outfile.write(line)

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
                lmbda1 = f"{0.5 * (sigmoid(float(i) / float(step), k) + 1):.3f}"
                lmbda2 = f"{0.5 * (-sigmoid(float(i) / float(step), k) + 1):.3f}"
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
        if restraint_method == "overlap":
            reslist = ["LIG", "LID"]
            torestraint_list = overlapping_pairs(pdbfile, reslist)
        else:
            parsed = parse_restraint_method(restraint_method)
            atom_max_distance = parsed.pop("kartograf_max_atom_distance", 0.95)
            parent_write_dir = Path(writedir).parent

            if self.system == "protein":  # In this case order of elements in PDB file is: prot, LIG, LID, HOH
                pdb_df = read_pdb_to_dataframe(Path(pdbfile).parent / "top_p.pdb")
            elif self.system == "water":  # Here both top_p.pdb
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

    def _set_common_md_replacements(self, lig_size1, lig_size2, eq_lambda):
        """Populate self.replacements with the atom ranges, sphere size, equilibration
        lambdas and water restraint shared by every MD/equilibration input file.

        Returns the extra [sequence_restraints] lines (only populated for the ABS water
        case; an empty list otherwise).
        """
        replacements = self.replacements
        lig_total = lig_size1 + lig_size2
        cnt = -1
        restlist = []
        replacements["ATOM_START_LIG1"] = f"{self.atomoffset + 1:<6}"
        replacements["ATOM_END_LIG1"] = f"{self.atomoffset + lig_size1:<7}"
        replacements["ATOM_START_LIG2"] = f"{self.atomoffset + lig_size1 + 1:<6}"
        replacements["ATOM_END_LIG2"] = f"{self.atomoffset + lig_size1 + lig_size2:<7}"
        replacements["SPHERE"] = self.sphereradius
        replacements["ATOM_END"] = f"{self.atomoffset + lig_total:<6}"
        replacements["EQ_LAMBDA"] = eq_lambda

        if self.system == "water" or self.system == "vacuum":
            if self.ABS is False:
                replacements["WATER_RESTRAINT"] = "{:<7}{:<7} 1.0 0 1   \n".format(
                    self.atomoffset + 1, self.atomoffset + lig_size1 + lig_size2
                )

            else:
                replacements["WATER_RESTRAINT"] = "{:<7}{:<7} 1.0 0 1   \n".format(
                    self.atomoffset + 1, self.atomoffset + lig_size1
                )

                for i in range(
                    self.atomoffset + 1 + lig_size1,
                    self.atomoffset + 2 + lig_size1 + lig_size2,
                ):
                    cnt += 1
                    if cnt == 0:
                        rest = f"{i:<7}{i:<7} 1.0 0 1   \n"
                        restlist.append(rest)

                    if cnt == 2:
                        cnt = -1

        elif self.system == "protein":
            replacements["WATER_RESTRAINT"] = ""
        return restlist

    def write_eq_files(self, writedir, overlapping_atoms, restlist):
        """Write the equilibration input files eq1-5.inp from the templates, injecting
        the distance and sequence restraints. Requires self.replacements to be populated
        (see _set_common_md_replacements). Returns the list of written file names.
        """
        file_list = []
        for eq_file_in in sorted(glob.glob(CONFIGS["ROOT_DIR"] + "/INPUTS/eq[1-5].inp")):
            eq_file = os.path.basename(eq_file_in)
            rest_force = 1.5 if eq_file != "eq5.inp" else self.dr_force  # 1.5 for eq1-4
            logger.debug(f"Writing {eq_file}")
            eq_file_out = writedir + "/" + eq_file

            with open(eq_file_in) as infile, open(eq_file_out, "w") as outfile:
                for line in infile:
                    line = replace(line, self.replacements)
                    outfile.write(line)
                    if line == "[distance_restraints]\n":
                        for atompair in overlapping_atoms:
                            outfile.write(f"{atompair[0]:d} {atompair[1]:d} 0.0 0.1 {rest_force:.1f} 0\n")

                    if line == "[sequence_restraints]\n":
                        for restline in restlist:
                            outfile.write(restline)
            file_list.append(eq_file)
        return file_list

    def write_MD_05(self, lambdas, writedir, lig_size1, lig_size2, overlapping_atoms):
        replacements = self.replacements
        file_list2 = []
        file_list3 = []
        lambda_1 = []
        lambda_2 = []
        block = 0
        index = 0

        for line in lambdas:
            if line == "0.500":
                block = 1

            if block == 0:
                lambda_1.append(line)

            if block == 1:
                lambda_2.append(line)

        lambda_1 = lambda_1[::-1]
        lambda_2 = lambda_2[1:]

        restlist = self._set_common_md_replacements(lig_size1, lig_size2, "0.500 0.500")
        file_list1 = self.write_eq_files(writedir, overlapping_atoms, restlist)

        # WRITING THE FEP MOLECULAR DYNAMICS INPUT FILES (e.g.: md_0500_0500.inp)
        file_in = CONFIGS["INPUT_DIR"] + "/md_0500_0500.inp"
        file_out = writedir + "/md_0500_0500.inp"
        with open(file_in) as infile, open(file_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)
                if line == "[distance_restraints]\n":
                    for line in overlapping_atoms:
                        outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 {self.dr_force:.1f} 0\n")

                if line == "[sequence_restraints]\n":
                    for line in restlist:
                        outfile.write(line)
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
        file_list1 = []  # eq files + initial md file
        file_list2 = []  # forward lambda files
        file_list3 = []  # reverse lambda files

        # Determine equilibration lambdas based on start mode
        if self.start == "0.5":
            eq_lambda1, eq_lambda2 = "0.500", "0.500"
        else:  # start == "0.0" or "1"
            eq_lambda1, eq_lambda2 = "1.000", "0.000"

        # Get configurations
        eq_configs = get_equilibration_configs(self.timestep, int(self.sphereradius))
        prod_config = get_production_config(self.timestep, int(self.sphereradius), self.dr_force)

        # Write equilibration files (eq1-eq5)
        for i, eq_config in enumerate(eq_configs):
            dr_str, seq_str = self._format_restraints_for_eq(
                overlapping_atoms, lig_size1, lig_size2, eq_config, self.dr_force
            )

            restart_file = f"eq{i}.re" if i > 0 else None

            content = render_md_input(
                params=eq_config.params,
                lambda1=eq_lambda1,
                lambda2=eq_lambda2,
                trajectory_file=f"{eq_config.name}.dcd",
                final_file=f"{eq_config.name}.re",
                restart_file=restart_file,
                distance_restraints=dr_str,
                sequence_restraints=seq_str,
                wall_restraints=wall_restraints,
                is_eq1=(i == 0),
            )

            filepath = Path(writedir) / f"{eq_config.name}.inp"
            filepath.write_text(content)
            file_list1.append(f"{eq_config.name}.inp")
            logger.debug(f"Writing {eq_config.name}.inp")

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
        self, file_out, eq_lambda, steps, overlapping_atoms, restlist, lambda_scaling=None
    ):
        """Write a single endpoint MD input file from the neq_endpoint.inp template.

        Used for both the endpoint equilibration files (eq6_*) and the lambda-switching
        files (neq_*). The per-replicate restart/final file names and the temperature stay
        as RESTART_VAR/FINAL_VAR/T_VAR placeholders that the run script fills in. When
        lambda_scaling is given, its lines are appended as the [lambda_scaling] section
        that drives lambda over the course of the simulation (turning the file into a
        switching run); without it the file is a plain equilibrium MD at the endpoint.
        """
        replacements = dict(self.replacements)
        replacements["EQ_LAMBDA"] = eq_lambda
        replacements["STEPS_VAR"] = str(steps)
        replacements["OUTPUT_VAR"] = "10"
        template = CONFIGS["INPUT_DIR"] + "/neq_endpoint.inp"
        with open(template) as infile, open(file_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)
                if line == "[distance_restraints]\n":
                    for atompair in overlapping_atoms:
                        outfile.write(f"{atompair[0]:d} {atompair[1]:d} 0.0 0.1 {self.dr_force:.1f} 0\n")
                if line == "[sequence_restraints]\n":
                    for restline in restlist:
                        outfile.write(restline)
            if lambda_scaling is not None:
                outfile.write("\n".join(lambda_scaling) + "\n")

    def write_MD_neq(self, writedir, lig_size1, lig_size2, overlapping_atoms):
        """Write the non-equilibrium input files: eq1-5 (equilibration), relax_{0,1}
        (one-time endpoint relaxation), eq6_{0,1} (endpoint equilibration spacing) and
        neq_{0,1} (lambda switching). State 0 starts at lambda 0.0->1.0 and state 1 at
        1.0->0.0; the run script relaxes each endpoint once, then chains the restarts and
        repeats the switches per replicate. Returns the list of written file names.
        """
        restlist = self._set_common_md_replacements(lig_size1, lig_size2, "0.500 0.500")
        file_list = self.write_eq_files(writedir, overlapping_atoms, restlist)

        endpoint_lambdas = {"0": "0.000 1.000", "1": "1.000 0.000"}
        lambda_scaling = [
            "",
            "[lambda_scaling]",
            f"scaling_parameter          {self.neq_schedule}",
            f"L_sigmoid        {self.neq_L}",
        ]
        for state, eq_lambda in endpoint_lambdas.items():
            relax_out = writedir + f"/relax_{state}.inp"
            self._write_endpoint_file(relax_out, eq_lambda, self.neq_relax_steps, overlapping_atoms, restlist)
            eq6_out = writedir + f"/eq6_{state}.inp"
            self._write_endpoint_file(eq6_out, eq_lambda, self.neq_eq_steps, overlapping_atoms, restlist)
            neq_out = writedir + f"/neq_{state}.inp"
            self._write_endpoint_file(
                neq_out, eq_lambda, self.neq_steps, overlapping_atoms, restlist, lambda_scaling
            )
            file_list += [f"relax_{state}.inp", f"eq6_{state}.inp", f"neq_{state}.inp"]
        return file_list

    def write_neq_runfile(self, writedir, file_list):
        """Write the SLURM run script for a non-equilibrium FEP. Each array task runs
        eq1-5, then loops `neq_reps` forward/reverse lambda switches with qdyn_neq.
        """
        src = CONFIGS["INPUT_DIR"] + "/run_neq.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"

        replacements = CLUSTER_DICT[self.cluster]
        replacements["FEPS"] = "FEP1.fep"
        replacements["NEQ_REPS"] = str(self.neq_reps)
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
        except:
            print("WARNING: Could not change permission for " + submit_out)

    def write_runfile(self, writedir, file_list):

        src = CONFIGS["INPUT_DIR"] + "/run.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"
        EQ_files = sorted(glob.glob(writedir + "/eq*.inp"))

        if self.start == "1":
            MD_files = reversed(sorted(glob.glob(writedir + "/md*.inp")))

        elif self.start == "0.5":
            md_1 = file_list[1]
            md_2 = file_list[2]

        replacements = CLUSTER_DICT[self.cluster]
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
                        outline = f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {file_base}.inp > {file_base}.log\n"
                        outfile.write(outline)

                if line.strip() == "#RUN_FILES":
                    if self.start == "1":
                        for line in MD_files:
                            file_base = line.split("/")[-1][:-4]
                            outline = (
                                f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {file_base}.inp"
                                f" > {file_base}.log\n"
                            )
                        outfile.write(outline)

                    elif self.start == "0.5":
                        outline = "time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0500_0500.inp > md_0500_0500.log\n\n"
                        outfile.write(outline)
                        for i, md in enumerate(md_1):
                            outline1 = f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {md_1[i][:-4]}.inp > {md_1[i][:-4]}.log\n"
                            outline2 = f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {md_2[i][:-4]}.inp > {md_2[i][:-4]}.log\n"

                            outfile.write(outline1)
                            outfile.write(outline2)
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
