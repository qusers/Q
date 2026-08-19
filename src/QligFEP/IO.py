import os
import re
import shlex
import stat
import subprocess
from pathlib import Path
from typing import NamedTuple, Optional

import numpy as np
import pandas as pd

from .functions import sigmoid
from .logger import logger
from .settings.settings import CONFIGS, FF_DIR

qfep_error_regex = re.compile(r"ERROR:")

# Failure markers -> the status label reported for a run, shared by both FEP analyzers
# (analyze_FEP and analyze_neq) so equilibrium and NEQ runs are classified identically.
#
# Dict order decides the label when a log carries several markers.
RUN_FAILURE_KEYWORDS = {
    "A request was made to bind to": "MPI_LAUNCH_FAILED",
    "There are not enough slots available": "MPI_LAUNCH_FAILED",
    "mpirun was unable to": "MPI_LAUNCH_FAILED",
    "DUE TO TIME LIMIT": "TIMEOUT",
    "CANCELLED": "CANCELLED",
    "Out Of Memory": "OOM",
    "abnormally": "CRASHED",
}


class SlurmRunInfo(NamedTuple):
    """Per-replicate run metadata parsed from one ``slurm*.out`` footer."""

    runtime: str
    seed: Optional[str]
    replicate: Optional[str]
    status: str


def read_slurm_diagnostics(slurm_out) -> SlurmRunInfo:
    """Parse one run.sh / run_neq.sh ``slurm*.out`` footer into a ``SlurmRunInfo``.

    Runtime, seed and replicate come from the standardized ``#    Runtime:`` /
    ``#    Random seed:`` / ``#    Replicate Number:`` footer the run scripts write. Reading the
    footer -- rather than the ``Parameters`` line or the slurm filename -- means the true
    ``run_num`` is used, so multi-temperature array jobs (where the filename's ``%a`` differs
    from the replicate) are labeled correctly. ``status`` is ``"SUCCESS"`` unless a known failure
    marker (``RUN_FAILURE_KEYWORDS``) appears anywhere in the log -- the first marker in dict order
    wins -- or, failing that, the footer reports a non-zero ``#    Exit status:``. Logs written
    before that footer line existed have none, and fall back to keyword matching alone.
    """
    with open(slurm_out, errors="ignore") as handle:
        text = handle.read()
    runtime, seed, replicate, status = "", None, None, "SUCCESS"
    exit_status = None
    for line in text.splitlines():
        if line.startswith("#    Runtime:"):
            runtime = line.split()[-1].strip()
        elif line.startswith("#    Random seed:"):
            seed = line.split()[-1].strip()
        elif line.startswith("#    Replicate Number:"):
            replicate = line.split()[-1].strip()
        elif line.startswith("#    Exit status:"):
            exit_status = line.split()[-1].strip()
    for keyword, label in RUN_FAILURE_KEYWORDS.items():
        if keyword in text:
            status = label
            break
    # A keyword names a specific cause, so it is preferred; the reported exit code is the catch-all
    # for failures we have no marker for.
    if status == "SUCCESS" and exit_status not in (None, "0"):
        status = "FAILED"
    return SlurmRunInfo(runtime, seed, replicate, status)


class QprepError(Exception):
    """Raised when qprep encounters a general error."""

    pass


class QprepAtomLibMissingError(Exception):
    """Raised when qprep encounters atom/residue not found in library."""

    pass


def ddG_json_path(mapping_json: str) -> str:
    """Return the ``<stem>_ddG.json`` sibling of a mapping JSON path.

    Derived from the path components rather than by string substitution, so a mapping file
    that is not named ``*.json`` still gets a distinct output instead of being overwritten,
    and a ``.json`` earlier in the path is left alone.
    """
    path = Path(mapping_json)
    return str(path.with_name(f"{path.stem}_ddG.json"))


def parse_qprep_total_charge(qprep_out_path: Path) -> int:
    """Parse the total system charge from a qprep.out file.

    Reads the "total charge of not excluded: X.00" line and returns the integer charge.

    Args:
        qprep_out_path: Path to the qprep.out file.

    Returns:
        Integer total charge of the system within the simulation sphere.

    Raises:
        ValueError: If the charge line is not found in the file.
    """
    charge_pat = re.compile(r"total charge of not excluded:\s*([-]?\d+\.\d+)")
    for line in qprep_out_path.read_text().splitlines():
        match = charge_pat.search(line)
        if match:
            return round(float(match.group(1)))
    raise ValueError(f"Could not find 'total charge of not excluded' in {qprep_out_path}")


def qprep_error_check(qprep_out_path: Path, ff_name: str) -> None:
    """Check for errors in the qprep.out file and raise an exception if any are found.

    Args:
        qprep_out_path: Path to the qprep.out file.
        ff_name: name of the forcefield to point user to the .lib & .prm files.

    Raises:
        QprepError: If any errors are found in the qprep.out file.
        QprepAtomLibMissingError: If atom/residue library mismatches are found.
    """
    error_pat = re.compile(r"ERROR\:\s", re.IGNORECASE)
    missing_lib_pat = re.compile(
        r">>> Atom ...?.? in residue no\.\s+\d+ not found in library entry for [A-Z]+"
        r"|>>> Heavy atom ...?.? missing in residue\s+ [0-9]+"
    )
    outfile_lines = qprep_out_path.read_text().split("\n")
    error_lines = []
    missing_atomlib_lines = []
    for line in outfile_lines:
        if error_pat.findall(line):
            error_lines.append(line)
        if missing_lib_pat.findall(line):
            missing_atomlib_lines.append(line)

    if error_lines:
        error_message = "\n".join(error_lines)
        logger.error(
            f"Found {len(error_lines)} error(s) in {qprep_out_path}. Check residue & atom names "
            f"against {FF_DIR / ff_name}.lib:\n{error_message}"
        )
        raise QprepError(error_message)
    if missing_atomlib_lines:
        error_message = "\n".join(missing_atomlib_lines)
        logger.error(
            f"Found {len(missing_atomlib_lines)} atom/residue mismatch(es) in {qprep_out_path}. "
            f"Check against {FF_DIR / ff_name}.lib:\n{error_message}"
        )
        raise QprepAtomLibMissingError(error_message)


def run_qprep(
    qprep_path: str,
    input_file: str = "qprep.inp",
    output_file: str = "qprep.out",
    ff_name: str = "AMBER14sb",
) -> subprocess.CompletedProcess:
    """Run qprep and check for errors.

    Args:
        qprep_path: Path to qprep executable.
        input_file: Input file name (default: qprep.inp).
        output_file: Output file name for inspection (default: qprep.out).
        ff_name: Force field name for error messages.

    Returns:
        CompletedProcess result.

    Raises:
        QprepError: On general qprep errors.
        QprepAtomLibMissingError: On atom/residue library mismatches.
    """
    cmd = f"{qprep_path} {input_file} > {output_file}"
    result = subprocess.run(cmd, shell=True, text=True)

    # Always check for errors in output file
    qprep_error_check(Path(output_file), ff_name)

    return result


## Some useful objects TO DO add GLH etc.
charged_res = {"HIS": {"HD1": "HID", "HE2": "HIE"}, "GLU": {"HE2": "GLH"}, "ASP": {"HD2": "ASH"}}


def get_force_field_paths(force_field: str):
    """Return the paths to the .lib and .prm files for the given force field. Inputs can
    either be Path-like objects or strings, as for QligFEP-implemented forcefields (AMBER14sb,
    CHARMM36, OPLS2005, OPLS2015, OPLSAAM).

    Args:
        force_field: Either a string with the name of the force field or a Path-like object
            pointing to the `lib`, the `prm` or both files (no extension).

    Raises:
        FileNotFoundError: In case the .lib or .prm file is not found in the given directory.

    Returns:
        tuple: A tuple containing the paths to the .lib and .prm files.
    """
    _available = ["AMBER14sb", "CHARMM36", "OPLS2005", "OPLS2015", "OPLSAAM"]
    if force_field in _available:
        lib_file = CONFIGS["FF_DIR"] + "/" + force_field + ".lib"
        prm_file = CONFIGS["FF_DIR"] + "/" + force_field + ".prm"
    else:
        ff_path_obj = Path(force_field)
        extensions = [".lib", ".prm"]
        for ext in extensions:
            if not ff_path_obj.with_suffix(ext).exists():
                logger.error(
                    "If passing a Path-like object, both .lib and .prm files must be present in the same directory."
                )
                raise FileNotFoundError(f"Could not find the {ext} file for the force field: {ff_path_obj}")
        lib_file = str(ff_path_obj.with_suffix(".lib").resolve().absolute())
        prm_file = str(ff_path_obj.with_suffix(".prm").resolve().absolute())
    return lib_file, prm_file


def parse_lib(force_field: str = "AMBER14sb") -> dict:
    """Parse a Q force field .lib file into a dict of residue entries.

    Args:
        force_field: Force field name (e.g., 'AMBER14sb') or path to .lib file.

    Returns:
        Dict mapping residue names to their entries. Each entry has:
        - 'atoms': list of dicts with 'name', 'type', 'charge'
        - 'charge_groups': lists of atom names in Q charge-group order
        - 'comment': the text after the residue name on the header line
    """
    lib_path, _ = get_force_field_paths(force_field)
    residues = {}
    current_res = None
    section = None

    with open(lib_path) as f:
        for line in f:
            line = line.rstrip("\n")
            # Residue header: {RESNAME}  ! comment
            m = re.match(r"^\{([^}]+)\}\s*(.*)", line)
            if m:
                current_res = m.group(1)
                comment = m.group(2).lstrip("! ").strip()
                residues[current_res] = {"atoms": [], "charge_groups": [], "comment": comment}
                section = None
                continue
            if current_res is None:
                continue
            stripped = line.strip()
            section_match = re.match(r"^\[([^]]+)\]", stripped)
            if section_match:
                section = section_match.group(1).lower()
                continue
            if not stripped or stripped.startswith("*"):
                continue
            content = stripped.split("!", 1)[0].strip()
            if not content:
                continue
            parts = content.split()
            if section == "atoms" and len(parts) >= 4:
                try:
                    charge = float(parts[3])
                except ValueError:
                    continue
                residues[current_res]["atoms"].append({"name": parts[1], "type": parts[2], "charge": charge})
            elif section == "charge_groups":
                residues[current_res]["charge_groups"].append(parts)
    return residues


def parse_prm_options(force_field: str = "AMBER14sb") -> dict[str, str]:
    """Parse the ``[options]`` section of a Q force-field parameter file."""
    _, prm_path = get_force_field_paths(force_field)
    options = {}
    section = None
    with open(prm_path) as f:
        for line in f:
            stripped = line.strip()
            section_match = re.match(r"^\[([^]]+)\]", stripped)
            if section_match:
                section = section_match.group(1).lower()
                continue
            if section != "options" or not stripped or stripped.startswith(("!", "*", "#")):
                continue
            parts = stripped.split("!", 1)[0].split()
            if len(parts) >= 2:
                options[parts[0].lower()] = parts[1].lower()
    return options


def lookup_residue(query: str, force_field: str = "AMBER14sb"):
    """Look up residue(s) in a Q force field .lib file. Prints matching entries.

    Usage from CLI:
        python -c "from QligFEP.IO import lookup_residue; lookup_residue('SOD')"
        python -c "from QligFEP.IO import lookup_residue; lookup_residue('HI')"  # partial match

    Args:
        query: Residue name or partial name to search for (case-insensitive).
        force_field: Force field name. Default: AMBER14sb.
    """
    residues = parse_lib(force_field)
    query_upper = query.upper()
    matches = {k: v for k, v in residues.items() if query_upper in k}
    if not matches:
        print(f"No residues matching '{query}' in {force_field}.lib")
        return

    for resname, entry in sorted(matches.items()):
        total_charge = sum(a["charge"] for a in entry["atoms"])
        print(f"\n{resname}  ({entry['comment']})  total_charge={total_charge:.2f}")
        print(f"  {'atom_name':<10s} {'atom_type':<10s} {'charge':>8s}")
        for atom in entry["atoms"]:
            print(f"  {atom['name']:<10s} {atom['type']:<10s} {atom['charge']:>8.4f}")


def replace(string, replacements):
    pattern = re.compile(r"\b(" + "|".join(replacements.keys()) + r")\b")
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string


def run_command(executable, options, string=False):
    """
    Runs a command and returns the result.

    Args:
        executable: The executable path
        options: Command options as a string
        string: If True, run as a shell command (for programs like Q that need one string).
                If False, split the command into args.

    Returns:
        For string=False: stdout as bytes (for backward compatibility)
        For string=True: CompletedProcess result
    """
    if string:
        cmd = f"{executable}{options}"
        result = subprocess.run(cmd, shell=True, text=True)
        return result
    else:
        args = shlex.split(f"{executable}{options}")
        print(" ".join(args))
        result = subprocess.run(args, capture_output=True, check=True)
        return result.stdout


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


def parse_prm(prm_contents: list[str]) -> dict[str, list[str]]:
    sections = {}
    current_section = None
    current_content = []

    for line in prm_contents:
        if line.startswith("[") and line.endswith("]"):
            if current_section is not None:  # If we were processing a section, save it
                sections[current_section] = "\n".join(current_content)
                current_content = []

            current_section = line.strip("[]")
            continue

        # If we're in a section, add the line to current content
        if current_section is not None:
            current_content.append(line)

    # save the last section
    if current_section is not None and current_content:
        sections[current_section] = "\n".join(current_content)

    return sections


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
            # lambda 2 + lambda 1 should be 1.0
            lmbda1 = f"{0.5 * (sigmoid(float(i)/float(step), k) + 1):.3f}"
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
