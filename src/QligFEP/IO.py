import os
import re
import shlex
import stat
from subprocess import check_output

import numpy as np
import pandas as pd

from .functions import sigmoid
from .logger import logger
from .settings.settings import CONFIGS

qfep_error_regex = re.compile(r"ERROR:")

## Some useful objects TO DO add GLH etc.
charged_res = {"HIS": {"HD1": "HID", "HE2": "HIE"}, "GLU": {"HE2": "GLH"}, "ASP": {"HD2": "ASH"}}


def replace(string, replacements):
    pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string


def run_command(executable, options, string=False):
    """
    Takes three variables, the executable location and its options as strings and a tag if the
    options need to be split or not (e.g. Q runs with one string), and runs the program.
    Returns the output of that program as an unformatted string.
    """
    if string is False:
        args = shlex.split(executable + options)
        out = check_output(args)
        print(" ".join(args))
    else:
        os.system(executable + options)
        out = None

    return out


def AA(AA):
    """
    Handy dictionary to convert 3 letter AA code to one and vice versa
    """
    threeAA = {
        "CYS": "C",
        "CYX": "C",
        "ASH": "D",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYN": "K",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HID": "H",
        "HIP": "H",
        "HIE": "H",
        "HIS": "H",
        "LEU": "L",
        "ARN": "R",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLH": "E",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
    }

    fourAA = {
        "CCYS": "C",
        "CASP": "D",
        "CASH": "H",
        "CSER": "S",
        "CGLN": "Q",
        "CLYN": "K",
        "CLYS": "K",
        "CILE": "I",
        "CPRO": "P",
        "CTHR": "T",
        "CPHE": "F",
        "CASN": "N",
        "CGLY": "G",
        "CHIE": "H",
        "CHID": "H",
        "CHIP": "H",
        "CLEU": "L",
        "CARG": "R",
        "CARN": "R",
        "CTRP": "W",
        "CALA": "A",
        "CVAL": "V",
        "CGLU": "E",
        "CGLH": "E",
        "CTYR": "Y",
        "CMET": "M",
    }

    oneAA = {
        "C": "CYS",
        "D": "ASP",
        "S": "SER",
        "Q": "GLN",
        "K": "LYS",
        "I": "ILE",
        "P": "PRO",
        "T": "THR",
        "F": "PHE",
        "N": "ASN",
        "G": "GLY",
        "H": "HID",
        "L": "LEU",
        "R": "ARG",
        "W": "TRP",
        "A": "ALA",
        "V": "VAL",
        "E": "GLU",
        "Y": "TYR",
        "M": "MET",
    }

    if len(AA) == 4:
        AA = fourAA[AA]
        return AA

    if len(AA) == 3:
        AA = threeAA[AA]
        return AA

    if len(AA) == 1:
        AA = oneAA[AA]
        return AA


def read_prm(prmfiles):
    """
    Takes a list of Q .prm files and merges them, first file is the referene .prm file
    Returns a dicitonary of the merged .prm files
    """
    block = 0
    prm = {
        "[options]": [],
        "[atom_types]": [],
        "[bonds]": [],
        "[angles]": [],
        "[torsions]": [],
        "[impropers]": [],
    }

    for filename in prmfiles:
        with open(filename) as infile:
            for line in infile:
                if line == "[options]\n":
                    block = 1
                    continue
                elif line == "[atom_types]\n":
                    block = 2
                    continue
                elif line == "[bonds]\n":
                    block = 3
                    continue
                elif line == "[angles]\n":
                    block = 4
                    continue
                elif line == "[torsions]\n":
                    block = 5
                    continue
                if line == "[impropers]\n":
                    block = 6
                    continue

                if block == 1:
                    prm["[options]"].append(line)

                if block == 2:
                    prm["[atom_types]"].append(line)

                elif block == 3:
                    prm["[bonds]"].append(line)

                elif block == 4:
                    prm["[angles]"].append(line)

                elif block == 5:
                    prm["[torsions]"].append(line)

                elif block == 6:
                    prm["[impropers]"].append(line)

    return prm


def get_lambdas(windows, sampling):
    windows = int(windows)
    step = int(windows / 2)
    lambdas = []
    lmbda_1 = []
    lmbda_2 = []
    k_dic = {"sigmoidal": -1.1, "linear": 1000, "exponential": -1.1, "reverse_exponential": 1.1}
    k = k_dic[sampling]

    if sampling == "sigmoidal":
        for i in range(0, step + 1):
            lmbda1 = f"{0.5 * (sigmoid(float(i)/float(step), k) + 1):.3f}"
            lmbda2 = f"{0.5 * (-sigmoid(float(i)/float(step), k) + 1):.3f}"
            lmbda_1.append(lmbda1)
            lmbda_2.append(lmbda2)

        lmbda_2 = lmbda_2[1:]

        for i in reversed(lmbda_2):
            lambdas.append(i)

        for i in lmbda_1:
            lambdas.append(i)

    else:
        for i in range(0, windows + 1):
            lmbda = f"{sigmoid(float(i)/float(windows), k):.3f}"
            lambdas.append(lmbda)

    lambdas = lambdas[::-1]
    return lambdas


def write_submitfile(writedir, replacements):
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


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def read_qfep(qfep):
    """
    Reads a given qfep.out file.

    returns [Zwanzig, dGfr, dGr, TI, OS, BAR]
    """
    with open(qfep) as infile:
        block = 0
        for line in infile:
            try:
                if qfep_error_regex.findall(line):
                    error_main = qfep_error_regex.findall(line)[0]
                    error_body = "".join([next(infile) for _ in range(2)])
                    raise OSError(f"QFEP ERROR !! {error_main}\n{error_body}")
            except StopIteration as e:
                logger.info("Reached the end of the file before capturing the full error body.")
                raise OSError(f"QFEP ERROR !! {error_main}") from e

            line = line.split()
            if len(line) > 3:

                if line[3] == "Free":
                    block = 1

                if line[3] == "Termodynamic":
                    # continue
                    block = 2

                if line[3] == "Overlap":
                    block = 3

                if line[3] == "BAR":
                    block = 4

                if line[3] == "Reaction":
                    block = 0

            if len(line) > 1:
                if block == 1:
                    if line[0] == "1.000000":
                        Zwanzig_r = float(line[4])

                    elif line[0] == "0.000000":
                        Zwanzig_f = float(line[2])

                        Zwanzig = np.nan if line[5] == "-Infinity" else float(line[5])

                if block == 2 and line[0] == "0.000000":
                    try:
                        thermo_integration = line[2]
                        if line[2] == "-Infinity":
                            thermo_integration = np.nan
                    except IndexError:
                        thermo_integration = np.nan  # TODO: this line is never reached... # noqa: F841

                if block == 3 and line[0] == "0.000000":
                    overlap_sampling = np.nan if line[2] == "-Infinity" else float(line[2])

                if block == 4 and line[0] == "0.000000":
                    bar = np.nan if line[2] == "-Infinity" else float(line[2])

    return [Zwanzig, Zwanzig_f, Zwanzig_r, overlap_sampling, bar]


def read_qfep_verbose(file_path):
    """
    Reads a qfep output file and extracts data from Part 0 and Part 6 into pandas DataFrames.

    Args:
        file_path (str): Path to the qfep output file.

    Returns:
        tuple: A tuple containing two pandas DataFrames:
               - part0_df: DataFrame for Part 0 data, or None if Part 0 is not found.
               - part6_df: DataFrame for Part 6 data, or None if Part 6 is not found.
    """
    part0_start = ">>>>> Enter a new file name (or nothing to cancel): "
    part0_data = []
    part6_data = []
    in_part0 = False
    in_part6 = False

    try:
        with open(file_path) as f:
            for line in f:
                line = line.strip()

                if line.startswith("# Part 0: Average energies for all states in all files"):
                    in_part0 = True
                    in_part6 = False  # Ensure other parts are not active
                    continue
                elif line.startswith("# Part 1:"):  # Stop reading Part 0 when Part 1 starts
                    in_part0 = False
                elif line.startswith("# Part 6: BAR Bennet:"):
                    in_part6 = True
                    in_part0 = False  # Ensure other parts are not active
                    continue
                elif line.startswith("# Part"):  # Stop reading Part 6 when next Part starts
                    in_part6 = False

                if in_part0:
                    if line.startswith("--> Name of file number"):
                        continue  # Skip header lines for file numbers
                    if line.startswith(">>>>> Failed to open"):  # Handle error lines
                        continue  # Skip error lines
                    if line.startswith("#") or not line:
                        continue  # Skip comment lines and empty lines
                    if line.startswith(part0_start):
                        line = line.replace(part0_start, "")  # Remove the start of the line

                    parts = line.split()
                    if len(parts) >= 17:  # Check if line has enough parts to parse
                        try:
                            file_name = parts[0]
                            state = int(parts[1])
                            pts = int(parts[2])
                            lambda_val = float(parts[3])
                            eqtot = float(parts[4]) if parts[4] != "********" else float("nan")
                            eqbond = float(parts[5]) if parts[5] != "********" else float("nan")
                            eqang = float(parts[6]) if parts[6] != "********" else float("nan")
                            eqtor = float(parts[7]) if parts[7] != "********" else float("nan")
                            eqimp = float(parts[8]) if parts[8] != "********" else float("nan")
                            eqel = float(parts[9]) if parts[9] != "********" else float("nan")
                            eqvdw = float(parts[10]) if parts[10] != "********" else float("nan")
                            eel_qq = float(parts[11]) if parts[11] != "********" else float("nan")
                            evdw_qq = float(parts[12]) if parts[12] != "********" else float("nan")
                            eel_qp = float(parts[13]) if parts[13] != "********" else float("nan")
                            evdw_qp = float(parts[14]) if parts[14] != "********" else float("nan")
                            eel_qw = float(parts[15]) if parts[15] != "********" else float("nan")
                            evdw_qw = float(parts[16]) if parts[16] != "********" else float("nan")
                            eqrstr = (
                                float(parts[17])
                                if len(parts) > 17 and parts[17] != "********"
                                else float("nan")
                            )

                            part0_data.append(
                                [
                                    file_name,
                                    state,
                                    pts,
                                    lambda_val,
                                    eqtot,
                                    eqbond,
                                    eqang,
                                    eqtor,
                                    eqimp,
                                    eqel,
                                    eqvdw,
                                    eel_qq,
                                    evdw_qq,
                                    eel_qp,
                                    evdw_qp,
                                    eel_qw,
                                    evdw_qw,
                                    eqrstr,
                                ]
                            )
                        except ValueError:
                            print(f"Warning: Skipping line due to parsing error in Part 0: {line}")
                            continue  # Skip line if parsing fails

                elif in_part6:
                    if line.startswith("#") or not line:
                        continue  # Skip comment lines and empty lines

                    parts = line.split()
                    if len(parts) == 3:  # Check if line has enough parts to parse for Part 6
                        try:
                            lambda_val_p6 = float(parts[0])
                            dg = float(parts[1])
                            sum_dg = float(parts[2])
                            part6_data.append([lambda_val_p6, dg, sum_dg])
                        except ValueError:
                            print(f"Warning: Skipping line due to parsing error in Part 6: {line}")
                            continue  # Skip line if parsing fails

    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        return None, None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None, None

    part0_df = None
    if part0_data:
        part0_columns = [
            "file_name",
            "state",
            "pts",
            "lambda_val",
            "EQtot",
            "EQbond",
            "EQang",
            "EQtor",
            "EQimp",
            "EQel",
            "EQvdW",
            "Eel_qq",
            "EvdW_qq",
            "Eel_qp",
            "EvdW_qp",
            "Eel_qw",
            "EvdW_qw",
            "Eqrstr",
        ]
        part0_df = pd.DataFrame(part0_data, columns=part0_columns)

    part6_df = None
    if part6_data:
        part6_columns = ["lambda_val", "dG", "sum_dG"]
        part6_df = pd.DataFrame(part6_data, columns=part6_columns)

    return part0_df, part6_df
