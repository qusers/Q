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

from .CLI.qprep_cli import QprepError, qprep_error_check
from .CLI.utils import get_avail_restraint_methods, handle_cysbonds
from .functions import COG, kT, overlapping_pairs, sigmoid
from .IO import get_force_field_paths, replace, run_command
from .logger import logger
from .pdb_utils import (
    calculate_distance,
    pdb_parse_in,
    pdb_parse_out,
    read_pdb_to_dataframe,
    rm_HOH_clash_NN,
)
from .restraints.restraint_setter import RestraintSetter
from .settings.settings import CLUSTER_DICT, CONFIGS


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
    ):
        self.replacements = {}  # TODO: make this explicit in the future
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
        # Temporary until flag is here
        self.ABS = False  # True
        self.ABS_waters = []
        self.write_dir = None
        self.pdb_fname = f"{self.lig1}_{self.lig2}.pdb"
        self.seeds = self.set_seeds(random_state)

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

    def set_timestep(self):
        if self.timestep == "1fs":
            logger.debug("Using 1fs timestep")
            self.replacements["NSTEPS1"] = "100000"
            self.replacements["NSTEPS2"] = "10000"
            self.replacements["STEPSIZE"] = "1.0"
            self.replacements["STEPTOGGLE"] = "off"

        elif self.timestep == "2fs":
            logger.debug("Using 2fs timestep")
            self.replacements["NSTEPS1"] = "50000"
            self.replacements["NSTEPS2"] = "5000"
            self.replacements["STEPSIZE"] = "2.0"
            self.replacements["STEPTOGGLE"] = "on"

        else:
            raise ValueError("Timestep not recognized")

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

    def write_FEP_file(self, change_charges, change_vdw, FEP_vdw, writedir, lig_size1, lig_size2):
        lig_size1 = int(lig_size1)
        lig_size2 = int(lig_size2)
        lig_tot = lig_size1 + lig_size2
        self.atomoffset = (
            read_pdb_to_dataframe(Path(writedir) / "top_p.pdb")
            .query("~residue_name.isin(['HOH', 'LIG', 'LID'])")
            .shape[0]
        )

        with open(writedir + "/FEP1.fep", "w") as outfile:
            total_atoms = len(change_charges)
            outfile.write("!info: " + self.lig1 + " --> " + self.lig2 + "\n")
            outfile.write("[FEP]\n")
            outfile.write("states 2\n")
            outfile.write("softcore_use_max_potential on\n\n")

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
        pattern = r"_(\d+\.?\d*)"  # check for the optional atom max distance
        match = re.search(pattern, restraint_method)
        if match:
            atom_max_distance = float(match.group(1))
            restraint_method = re.sub(pattern, "", restraint_method)
        else:
            atom_max_distance = 0.95

        avail_methods = get_avail_restraint_methods()
        if restraint_method not in avail_methods:
            raise ValueError(f"Method {restraint_method} not recognized. Please use one of {avail_methods}")

        pdbfile = writedir + f"/inputfiles/{self.lig1}_{self.lig2}.pdb"
        if restraint_method == "overlap":
            reslist = ["LIG", "LID"]
            torestraint_list = overlapping_pairs(pdbfile, reslist)

            if self.ABS:
                with open(pdbfile) as infile:
                    for line in infile:
                        if line[13].strip() == "O":
                            line = pdb_parse_in(line)
                            self.ABS_waters.append(int(line[1]) + self.atomoffset)

        else:
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
            subset_lig1 = pdb_df.query("residue_name == 'LIG'")
            subset_lig2 = pdb_df.query("residue_name == 'LID'")
            logger.debug(f"lig1.shape: {subset_lig1.shape}")
            logger.debug(f"lig2.shape: {subset_lig2.shape}")
            if any([subset_lig1.shape[0] == 0, subset_lig2.shape[0] == 0]):
                raise QprepError(
                    "Something went wrong while preparing the protein edge. Please have a look at the "
                    f"input files and at the qprep.out within {Path(writedir / 'qprep.inp')}"
                )
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
                if restraint_method == "kartograf":
                    restraints = rsetter.set_restraints(kartograf_native=True)
                else:
                    atom_compare_method, permissiveness_lvl = restraint_method.split("_")
                    if permissiveness_lvl == "p":
                        params = {"strict_surround": False}
                    elif permissiveness_lvl == "ls":
                        params = {"strict_surround": True, "ignore_surround_atom_type": True}
                    elif permissiveness_lvl == "strict":
                        params = {"strict_surround": True, "ignore_surround_atom_type": False}
                    restraints = rsetter.set_restraints(atom_compare_method=atom_compare_method, **params)
                    logger.debug(f"Restraints set using {restraint_method} method. Parameters: {params}")
                if strict_check:  # Good to check in case sdf in directory doesn't belong to the structure
                    rdLig1 = rsetter.molA.to_rdkit()
                    rdLig2 = rsetter.molB.to_rdkit()
                    for AtomIdx_Lig1, AtomIdx_Lig2 in restraints.items():
                        rowLig1 = subset_lig1.iloc[AtomIdx_Lig1]
                        rowLig2 = subset_lig2.iloc[AtomIdx_Lig2]
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

    def write_MD_05(self, lambdas, writedir, lig_size1, lig_size2, overlapping_atoms):
        replacements = self.replacements
        file_list1 = []
        file_list2 = []
        file_list3 = []
        lig_total = lig_size1 + lig_size2
        lambda_1 = []
        lambda_2 = []
        block = 0
        index = 0
        cnt = -1
        restlist = []

        for line in lambdas:
            if line == "0.500":
                block = 1

            if block == 0:
                lambda_1.append(line)

            if block == 1:
                lambda_2.append(line)

        lambda_1 = lambda_1[::-1]
        lambda_2 = lambda_2[1:]
        replacements["ATOM_START_LIG1"] = f"{self.atomoffset + 1:<6}"
        replacements["ATOM_END_LIG1"] = f"{self.atomoffset + lig_size1:<7}"
        replacements["ATOM_START_LIG2"] = f"{self.atomoffset + lig_size1 + 1:<6}"
        replacements["ATOM_END_LIG2"] = f"{self.atomoffset + lig_size1 + lig_size2:<7}"
        replacements["SPHERE"] = self.sphereradius
        replacements["ATOM_END"] = f"{self.atomoffset + lig_total:<6}"
        replacements["EQ_LAMBDA"] = "0.500 0.500"

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

        # WRITING THE EQUILIBRATION INPUT FILES (eq1-5.inp), NOT PART OF THE FEP YET
        for eq_file_in in sorted(glob.glob(CONFIGS["ROOT_DIR"] + "/INPUTS/eq*.inp")):
            eq_file = eq_file_in.split("/")[-1:][0]
            rest_force = 1.5 if eq_file != "eq5.inp" else self.dr_force  # 1.5 for eq1-4
            logger.debug(f"Writing {eq_file}")
            eq_file_out = writedir + "/" + eq_file

            with open(eq_file_in) as infile, open(eq_file_out, "w") as outfile:
                for line in infile:
                    line = replace(line, replacements)
                    outfile.write(line)
                    if line == "[distance_restraints]\n":
                        for line in overlapping_atoms:
                            outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 {rest_force:.1f} 0\n")

                    if line == "[sequence_restraints]\n":
                        for line in restlist:
                            outfile.write(line)
            file_list1.append(eq_file)

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
        file_list1.append("md_0500_0500.inp")

        for lambdas in [lambda_1, lambda_2]:
            index += 1
            filename_N = "md_0500_0500"
            filenr = -1

            for line in lambdas:
                filenr += 1
                if index == 1:
                    lambda1 = lambda_1[filenr]
                    lambda2 = lambda_2[filenr]

                elif index == 2:
                    lambda1 = lambda_2[filenr]
                    lambda2 = lambda_1[filenr]

                filename = "md_" + lambda1.replace(".", "") + "_" + lambda2.replace(".", "")
                replacements["FLOAT_LAMBDA1"] = lambda1
                replacements["FLOAT_LAMBDA2"] = lambda2
                replacements["FILE"] = filename
                replacements["FILE_N"] = filename_N

                # Consider putting this in a function seeing as it is called multiple times
                pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
                file_in = CONFIGS["INPUT_DIR"] + "/md_XXXX_XXXX.inp"
                file_out = writedir + "/" + filename + ".inp"

                with open(file_in) as infile, open(file_out, "w") as outfile:
                    for line in infile:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                        outfile.write(line)
                        if line == "[distance_restraints]\n":
                            for line in overlapping_atoms:
                                outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 {self.dr_force:.1f} 0\n")

                    if line == "[sequence_restraints]\n":
                        for line in restlist:
                            outfile.write(line)
                filename_N = filename

                if index == 1:
                    file_list2.append(filename + ".inp")

                elif index == 2:
                    file_list3.append(filename + ".inp")
        return [file_list1, file_list2, file_list3]

    def write_MD_1(self, lambdas, writedir, lig_size1, lig_size2, overlapping_atoms):
        replacements = self.replacements
        totallambda = len(lambdas)
        file_list_1 = []
        file_list_2 = []
        file_list_3 = []
        replacements = {}
        lig_total = lig_size1 + lig_size2

        replacements["ATOM_START_LIG1"] = f"{self.atomoffset + 1:<6}"
        replacements["ATOM_END_LIG1"] = f"{self.atomoffset + lig_size1:<7}"
        replacements["ATOM_START_LIG2"] = f"{self.atomoffset + lig_size1 + 1:<6}"
        replacements["ATOM_END_LIG2"] = f"{self.atomoffset + lig_size1 + lig_size2:<7}"
        replacements["SPHERE"] = self.sphereradius
        replacements["ATOM_END"] = f"{self.atomoffset + lig_total:<6}"
        replacements["EQ_LAMBDA"] = "1.000 0.000"

        if self.system == "water" or self.system == "vacuum":
            replacements["WATER_RESTRAINT"] = "{:<7}{:<7} 1.0 0 1   ".format(
                self.atomoffset + 1, self.atomoffset + lig_size1 + lig_size2
            )
        elif self.system == "protein":
            replacements["WATER_RESTRAINT"] = ""

        for eq_file_in in sorted(glob.glob(CONFIGS["ROOT_DIR"] + "/INPUTS/eq*.inp")):
            eq_file = eq_file_in.split("/")[-1:][0]
            eq_file_out = writedir + "/" + eq_file
            with open(eq_file_in) as infile:
                with open(eq_file_out, "w") as outfile:
                    for line in infile:
                        line = replace(line, replacements)
                        outfile.write(line)
                        if line == "[distance_restraints]\n":
                            for line in overlapping_atoms:
                                outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.2 0.5 0\n")
                file_list_1.append(eq_file)

        file_in = CONFIGS["INPUT_DIR"] + "/md_1000_0000.inp"
        file_out = writedir + "/md_1000_0000.inp"
        with open(file_in) as infile, open(file_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)
                if line == "[distance_restraints]\n":
                    for line in overlapping_atoms:
                        outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.2 0.5 0\n")

        file_list_1.append("md_1000_0000.inp")
        filenr = 0

        for lambd in lambdas:
            if lambd == "1.000":
                filename_N = "md_1000_0000"
                continue
            else:
                step_n = totallambda - filenr - 2

                lambda1 = lambd
                lambda2 = lambdas[step_n]
                filename = "md_" + lambda1.replace(".", "") + "_" + lambda2.replace(".", "")
                replacements["FLOAT_LAMBDA1"] = lambda1
                replacements["FLOAT_LAMBDA2"] = lambda2
                replacements["FILE"] = filename
                replacements["FILE_N"] = filename_N

                # Move to functio
                pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
                file_in = CONFIGS["INPUT_DIR"] + "/md_XXXX_XXXX.inp"
                file_out = writedir + "/" + filename + ".inp"

                with open(file_in) as infile, open(file_out, "w") as outfile:
                    for line in infile:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                        outfile.write(line)
                        if line == "[distance_restraints]\n":
                            for line in overlapping_atoms:
                                outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 0.5 0\n")

                filename_N = filename
                filenr += 1
                file_list_2.append(filename + ".inp")

        return [file_list_1, file_list_2, file_list_3]

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
        qfep_in = CONFIGS["ROOT_DIR"] + "/INPUTS/qfep.inp"
        qfep_out = self.write_dir + "/inputfiles/qfep.inp"
        i = 0
        total_l = len(lambdas)

        # TO DO: multiple files will need to be written out for temperature range
        k_T = kT(float(self.temperature))
        replacements = {}
        replacements["kT"] = str(k_T)
        replacements["WINDOWS"] = str(windows)
        replacements["TOTAL_L"] = str(total_l)
        with open(qfep_in) as infile, open(qfep_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)

            if line == "!ENERGY_FILES\n":
                for i in range(0, total_l):
                    j = -(i + 1)
                    lambda1 = lambdas[i]
                    lambda2 = lambdas[j]
                    filename = "md_" + lambda1.replace(".", "") + "_" + lambda2.replace(".", "") + ".en\n"
                    outfile.write(filename)

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

    def write_qprep(self, writedir):
        """Write the qprep.inp file for Q. If the system is water, the center of geometry
        will be calculated from the lig1's atoms coordinates. If the system is protein, it
        will try to extract the center of geometry from the TITLE line in the water.pdb file,
        added by the `qprep_prot` program.

        Args:
            writedir: directory in which QligFEP will write the input files.
        """
        replacements = {}
        cog = None
        cog_regex = re.compile(r"([-]?\d{1,3}\.\d{3})\s+([-]?\d{1,3}\.\d{3})\s+([-]?\d{1,3}\.\d{3})")
        # see if the water.pdb file has a COG from qprep
        first_lines = (Path(writedir).parents[1] / "water.pdb").read_text().split("\n")[:3]
        for line in first_lines:
            if line.startswith("TITLE"):
                _matches = cog_regex.search(line)
                center = " ".join(_matches.groups()) if _matches else None
                if _matches is None:
                    logger.warning(f"Failed to extract COG from line:\n{line}")
                    logger.info("Will calculate the COG from the ligand atoms.")
        center = COG(self.lig1 + ".pdb") if cog is None or self.system == "water" else cog
        center = f"{center[0]} {center[1]} {center[2]}"
        self.cog = [float(coord) for coord in center.split()]
        qprep_in = CONFIGS["ROOT_DIR"] + "/INPUTS/qprep.inp"
        qprep_out = writedir + "/qprep.inp"
        replacements["FF_LIB"] = self.lib_file
        replacements["LIG1"] = self.lig1 + ".lib"
        replacements["LIG2"] = self.lig2 + "_renumber.lib"
        replacements["LIGPRM"] = self.FF + "_" + self.lig1 + "_" + self.lig2 + "_merged.prm"
        replacements["LIGPDB"] = self.lig1 + "_" + self.lig2 + ".pdb"
        replacements["CENTER"] = center
        replacements["SPHERE"] = self.sphereradius
        if self.system == "vacuum":
            replacements["solvate"] = "!solvate"
        if self.system == "water":
            replacements["SOLVENT"] = "1 HOH"
        if self.system == "protein":
            replacements["SOLVENT"] = "4 water.pdb"

        pdb_df = read_pdb_to_dataframe(Path(writedir) / self.pdb_fname)
        density = self.get_sphere_density(pdb_df) if self.system == "protein" else 0.05794
        replacements["SOLUTEDENS"] = f"{density:.5f}"

        with open(qprep_in) as infile, open(qprep_out, "w") as outfile:
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
                        logger.debug(
                            f"{resn1}:{at1} and {resn2}:{at2} within {dist1} and {dist2} of the COG."
                        )
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
            for line in infile:
                line = replace(line, replacements)
                if line == "!addbond at1 at2 y\n" and cysbond_str != "":
                    outfile.write(cysbond_str)
                    continue
                outfile.write(line)

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
        qprep = cluster_options["QPREP"]
        logger.info(f"Running QPREP from path {qprep}")
        options = " < qprep.inp > qprep.out"
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        run_command(qprep, options, string=True)
        qprep_error_check(Path("qprep.out"), self.FF)
        os.chdir("../../")
