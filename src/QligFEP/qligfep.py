import glob
import os
import re
import shutil
import stat
from pathlib import Path
from typing import Literal, Optional

import numpy as np

from .CLI.qprep_cli import qprep_error_check
from .CLI.utils import get_avail_restraint_methods, handle_cysbonds
from .functions import COG, kT, overlapping_pairs, sigmoid
from .IO import replace, run_command
from .logger import logger
from .pdb_utils import (
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
        softcore: bool = False,  # Not implemented yet
        to_clean: Optional[list[str]] = None,
        random_state: Optional[int] = 42,
        **kwargs,
    ):
        self.softcore = softcore
        self.replacements = {}  # TODO: make this explicit in the future
        self.timestep = timestep
        self.lig1 = lig1
        self.lig2 = lig2
        self.FF = FF
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
        # Temporary until flag is here
        self.ABS = False  # True
        self.ABS_waters = []
        self.write_dir = None
        self.pdb_fname = f"{self.lig1}_{self.lig2}.pdb"
        self.seeds = self.set_seeds(random_state)

        if self.system == "protein":
            # Get last atom and residue from complexfile!
            with open("protein.pdb") as infile:
                for line in infile:
                    try:
                        resnr = int(line[22:26])
                        atnr = int(line[6:11])
                    except IndexError:
                        continue
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
                        if self.FF == "AMBER14sb" or self.FF == "CHARMM36":
                            j = "X" + i
                        else:
                            match = re.match(r"([a-z]+)([0-9]+)", i, re.I)
                            if match:
                                items = match.groups()
                                j = str(items[0]) + str(int(items[1]) + int(molsize_lig1))

                        if self.FF == "CHARMM_TEST":
                            j = "X_" + i

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
            block = 0
            for line in file_replaced:
                outfile.write(line)

        shutil.copy(self.lig1 + ".lib", writedir + "/" + self.lig1 + ".lib")

    def change_prm(self, replacements, writedir):
        pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
        file1 = glob.glob(self.lig1 + ".prm")[0]
        file2 = glob.glob(self.lig2 + ".prm")[0]
        prm_file = CONFIGS["FF_DIR"] + "/" + self.FF + ".prm"
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
                line2 = "{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}".format(
                    line[0], line[1], line[3], str(0), str(0), line[4], line[5], line[6]
                )
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
                    for line in protfile:
                        outfile.write(line)
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

    def set_restraints(self, writedir, restraint_method, strict_check: bool = True):
        """Function to set the restraints for FEP. Originally, this was performed on
        overlapping atoms, but based on our observations this was changed to a more
        chemistry-aware method, implemented under `QligFEP.restraints.restraint_setter`.

        The configuration on how these restraints will be applied depend on two strings, passed into
        `method` as `{ring_compare_method}_{surround_compare_method}`. Alternatively, the user can
        opt for `overlap` which simply restrains atoms within 1 A from each other.

        Explanation:
            Ring atom compare: `aromaticity`, `hibridization`, `element`. Setting the first part of the
                string as either of these, will determine how the ring atoms are treated to be defined as
                equivalent.
            Surround atom compare: `p` (permissive), `ls` (less strict), `strict`.
                Setting the second part of the string as either of these, will determine if or how the
                direct surrounding atoms to the ring strictures will be taken into account for ring equivalence.
                    - Permissive: Only the ring atoms are compared.
                    - Less strict: The ring atoms and their direct surroundings are compared, but element type
                        is ignored.
                    - Strict: The ring atoms and their direct surroundings are element-wise compared.

        Args:
            writedir: directory to get the input files from, e.g.: FEP_lig1_lig2/inputfiles.
            strict_check: whether to assert the atom indexes are correctly assigned.

        Returns:
            list: list of overlapping atoms.
        """
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
            lig1_path = parent_write_dir / f"{self.lig1}.sdf"
            lig2_path = parent_write_dir / f"{self.lig2}.sdf"
            if not lig1_path.exists() or not lig2_path.exists():
                logger.error(
                    "Loading ligands from PDB files for the `chemoverlap` method is not supported yet."
                )
                raise FileNotFoundError(
                    "If you're using the `chemoverlap` method, you need to have the `sdf` ligands in the "
                    f"same directory you're using for the FEP calculations: {parent_write_dir}"
                )
            else:
                logger.debug(f'Loading sdf for restraint calculation:\nlig1:"{lig1_path}"\nlig2"{lig2_path}"')
                rsetter = RestraintSetter(lig1_path, lig2_path)
                ring_atom_compare = restraint_method.split("_")[0]
                surround_atom_compare = restraint_method.split("_")[1]
                if surround_atom_compare == "p":
                    strict_surround = False
                    ignore_surround_atom_type = False
                else:
                    strict_surround = True
                    ignore_surround_atom_type = surround_atom_compare == "ls"
                restraints = rsetter.set_restraints(
                    ring_compare_method=ring_atom_compare,
                    strict_surround=strict_surround,
                    ignore_surround_atom_type=ignore_surround_atom_type,
                )
                if strict_check:  # TODO: This should be moved to the tests in the future...
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

        for eq_file_in in sorted(glob.glob(CONFIGS["ROOT_DIR"] + "/INPUTS/eq*.inp")):
            eq_file = eq_file_in.split("/")[-1:][0]
            logger.debug(f"Writing {eq_file}")
            eq_file_out = writedir + "/" + eq_file

            with open(eq_file_in) as infile:
                with open(eq_file_out, "w") as outfile:
                    for line in infile:
                        line = replace(line, replacements)
                        outfile.write(line)
                        if line == "[distance_restraints]\n":
                            for line in overlapping_atoms:
                                outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 3.0 0\n")

                        if line == "[sequence_restraints]\n":
                            for line in restlist:
                                outfile.write(line)

                file_list1.append(eq_file)

        file_in = CONFIGS["INPUT_DIR"] + "/md_0500_0500.inp"
        file_out = writedir + "/md_0500_0500.inp"
        with open(file_in) as infile, open(file_out, "w") as outfile:
            for line in infile:
                line = replace(line, replacements)
                outfile.write(line)
                if line == "[distance_restraints]\n":
                    for line in overlapping_atoms:
                        outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 0.5 0\n")

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
                                outfile.write(f"{line[0]:d} {line[1]:d} 0.0 0.1 0.5 0\n")

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
        lambda_1 = []
        lambda_2 = []

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
        replacements["TEMP_VAR"] = str(self.temperature)
        replacements["RUN_VAR"] = str(self.replicates)
        replacements["RUNFILE"] = "run" + self.cluster + ".sh"
        replacements["SEED_VAR"] = " ".join([str(s) for s in self.seeds])
        if self.softcore:
            submit_in = CONFIGS["ROOT_DIR"] + "/INPUTS/FEP_submit_sc.sh"
        else:
            submit_in = CONFIGS["ROOT_DIR"] + "/INPUTS/FEP_submit.sh"
        submit_out = writedir + ("/FEP_submit.sh")
        with open(submit_in) as infile, open(submit_out, "w") as outfile:
            for line in infile:
                try:
                    line = replace(line, replacements)
                except TypeError:
                    print("heyyy")
                outfile.write(line)

        try:
            st = os.stat(submit_out)
            os.chmod(submit_out, st.st_mode | stat.S_IEXEC)

        except:
            print("WARNING: Could not change permission for " + submit_out)

    def write_runfile(self, writedir, file_list):

        ntasks = CLUSTER_DICT[self.cluster]["NTASKS"]
        src = CONFIGS["INPUT_DIR"] + "/run.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"
        EQ_files = sorted(glob.glob(writedir + "/eq*.inp"))

        if self.start == "1":
            MD_files = reversed(sorted(glob.glob(writedir + "/md*.inp")))

        elif self.start == "0.5":
            md_1 = file_list[1]
            md_2 = file_list[2]

        replacements = CLUSTER_DICT[self.cluster]
        if self.start == "1" and self.softcore is True:  # not implemented yet
            replacements["FEPS"] = "FEP1.fep FEP2.fep FEP3.fep"
        else:
            replacements["FEPS"] = "FEP1.fep"
        with open(src) as infile, open(tgt, "w") as outfile:
            for line in infile:
                if line.strip() == "#SBATCH -A ACCOUNT":
                    try:  # Try to take account info - not for all clusters!
                        replacements["ACCOUNT"]
                    except KeyError:
                        line = ""
                if line.strip() == "#SBATCH -J JOBNAME":
                    if self.cluster == 'DARDEL': # TODO: refactor this...
                        outfile.write('#SBATCH -p shared\n')
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
                outfile.write(outline)
                if line.strip() == "#EQ_FILES":
                    for line in EQ_files:
                        file_base = Path(line).stem
                        outline = f"time srun -n $SLURM_NTASKS $qdyn {file_base}.inp > {file_base}.log\n"
                        outfile.write(outline)

                if line.strip() == "#RUN_FILES":
                    if self.start == "1":
                        for line in MD_files:
                            file_base = line.split("/")[-1][:-4]
                            outline = "time srun -n $SLURM_NTASKS $qdyn {}.inp" " > {}.log\n".format(
                                file_base, file_base
                            )
                        outfile.write(outline)

                    elif self.start == "0.5":
                        outline = "time srun -n $SLURM_NTASKS $qdyn {}.inp" " > {}.log\n\n".format(
                            "md_0500_0500", "md_0500_0500"
                        )
                        outfile.write(outline)
                        for i, md in enumerate(md_1):
                            outline1 = "time srun -n $SLURM_NTASKS $qdyn {}.inp" " > {}.log\n".format(
                                md_1[i][:-4], md_1[i][:-4]
                            )

                            outline2 = "time srun -n $SLURM_NTASKS $qdyn {}.inp" " > {}.log\n".format(
                                md_2[i][:-4], md_2[i][:-4]
                            )

                            outfile.write(outline1)
                            outfile.write(outline2)
                            outfile.write("\n")
                if line.strip() == "#CLEANUP" and self.to_clean is not None:
                    replacements["CLEANUP"] = "#Cleaned {} files\n".format(" ".join(self.to_clean))
                    outline = replace(line, replacements)
                    for suffix in self.to_clean:
                        if suffix is None:
                            break
                        outfile.write(f"rm -f *{suffix}\n")
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

    def avoid_water_protein_clashes(
        self, writedir, prot_th=1.4, water_th=1.4, header: Optional[str] = None, save_removed: bool = False
    ):
        """Function to remove water molecules too close to protein & ligands | ligands.
        Thresholds are the distances in Ångström from the protein & ligands | ligands atoms
        to the nearest atom in the water molecule (HOH).

        Args:
            writedir: directory in which QligFEP will write the input files.
            prot_th: threshold (Å) to remove water molecules clashing with protein & ligand
                on the protein leg. Defaults to 1.4.
            water_th: threshold (Å) to remove water molecules clashing with the ligands
                in the water leg. Defaults to 1.4.
            header: header to be added to the water.pdb file. This detail is needed in Qprep's
                current version (2024/07/21), as it doesn't recognize the radius when merging an
                existing water.pdb file to the topology.
        """
        waterfile = Path(writedir) / "water.pdb"
        protfile = Path(writedir) / self.pdb_fname
        threshold = prot_th if self.system == "protein" else water_th
        system_to_log = "protein & ligands" if self.system == "protein" else "ligands"

        logger.info(f"Removing water molecules too close to {system_to_log} - threshold: {threshold} A.")
        _, n_removed = rm_HOH_clash_NN(
            pdb_df_query=read_pdb_to_dataframe(waterfile),
            pdb_df_target=read_pdb_to_dataframe(protfile),
            th=threshold,
            output_file=waterfile,
            heavy_only=True,
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
        qprep_in = CONFIGS["ROOT_DIR"] + "/INPUTS/qprep.inp"
        qprep_out = writedir + "/qprep.inp"
        replacements["FF_LIB"] = CONFIGS["ROOT_DIR"] + "/FF/" + self.FF + ".lib"
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

        with open(qprep_in) as infile, open(qprep_out, "w") as outfile:
            # FIXME: it might be that the cystein bond is out of the sphere (excluded from TOP). Need to check that
            cysbond_str = handle_cysbonds(
                self.cysbond, Path(writedir) / self.pdb_fname, comment_out=(self.system != "protein")
            )
            for line in infile:
                line = replace(line, replacements)
                if line == "!addbond at1 at2 y\n" and cysbond_str != "":
                    outfile.write(cysbond_str)
                    continue
                outfile.write(line)

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
