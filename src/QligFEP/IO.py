import os
import re
import shlex
import stat
from subprocess import check_output

import numpy as np

from .functions import sigmoid
from .logger import logger
from .settings.settings import CONFIGS

qfep_error_regex = re.compile(r'ERROR:')

## Some useful objects TO DO add GLH etc.
charged_res = {'HIS': {'HD1' : 'HID',
                       'HE2' : 'HIE'},
               'GLU': {'HE2' : 'GLH'},
               'ASP': {'HD2' : 'ASH'}
              }

def replace(string, replacements):
    pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string

def run_command(executable, options, string = False):
    """
    Takes three variables, the executable location and its options as strings and a tag if the
    options need to be split or not (e.g. Q runs with one string), and runs the program.
    Returns the output of that program as an unformatted string.
    """
    if string is False:
        args = shlex.split(executable + options)
        out = check_output(args)
        print(' '.join(args))
    else:
        os.system(executable + options)
        out = None

    return out

def AA(AA):
    """
    Handy dictionary to convert 3 letter AA code to one and vice versa
    """
    threeAA = {'CYS': 'C', 'CYX': 'C', 'ASH': 'D', 'ASP': 'D', 'SER': 'S',
               'GLN': 'Q', 'LYN': 'K', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
               'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HID': 'H',
               'HIP': 'H', 'HIE': 'H', 'HIS': 'H', 'LEU': 'L', 'ARN': 'R',
               'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLH': 'E',
               'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
              }

    fourAA = { 'CCYS': 'C', 'CASP': 'D', 'CASH': 'H', 'CSER': 'S',
               'CGLN': 'Q', 'CLYN': 'K', 'CLYS': 'K', 'CILE': 'I',
               'CPRO': 'P', 'CTHR': 'T', 'CPHE': 'F', 'CASN': 'N',
               'CGLY': 'G', 'CHIE': 'H', 'CHID': 'H', 'CHIP': 'H',
               'CLEU': 'L', 'CARG': 'R', 'CARN': 'R', 'CTRP': 'W',
               'CALA': 'A', 'CVAL': 'V', 'CGLU': 'E', 'CGLH': 'E',
               'CTYR': 'Y', 'CMET': 'M'
             }

    oneAA = {  'C' : 'CYS', 'D' : 'ASP', 'S' : 'SER', 'Q' : 'GLN',
               'K' : 'LYS', 'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR',
               'F' : 'PHE', 'N' : 'ASN', 'G' : 'GLY', 'H' : 'HID',
               'L' : 'LEU', 'R' : 'ARG', 'W' : 'TRP', 'A' : 'ALA',
               'V' : 'VAL', 'E' : 'GLU', 'Y' : 'TYR', 'M' : 'MET'
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
    prm = {'[options]':[],
            '[atom_types]':[],
            '[bonds]':[],
            '[angles]':[],
            '[torsions]':[],
            '[impropers]':[]}

    for filename in prmfiles:
        with open(filename) as infile:
            for line in infile:
                if line == '[options]\n':
                    block = 1
                    continue
                elif line == '[atom_types]\n':
                    block = 2
                    continue
                elif line == '[bonds]\n':
                    block = 3
                    continue
                elif line == '[angles]\n':
                    block = 4
                    continue
                elif line == '[torsions]\n':
                    block = 5
                    continue
                if line == '[impropers]\n':
                    block = 6
                    continue

                if block == 1:
                    prm['[options]'].append(line)

                if block == 2:
                    prm['[atom_types]'].append(line)

                elif block == 3:
                    prm['[bonds]'].append(line)

                elif block == 4:
                    prm['[angles]'].append(line)

                elif block == 5:
                    prm['[torsions]'].append(line)

                elif block == 6:
                    prm['[impropers]'].append(line)

    return prm

def get_lambdas(windows, sampling):
    windows = int(windows)
    step = int(windows/2)
    lambdas = []
    lmbda_1 = []
    lmbda_2 = []
    k_dic = {'sigmoidal':-1.1,
             'linear':1000,
             'exponential':-1.1,
             'reverse_exponential':1.1
            }
    k = k_dic[sampling]

    if sampling == 'sigmoidal':
        for i in range(0, step + 1):
            lmbda1 = f'{0.5 * (sigmoid(float(i)/float(step), k) + 1):.3f}'
            lmbda2 = f'{0.5 * (-sigmoid(float(i)/float(step), k) + 1):.3f}'
            lmbda_1.append(lmbda1)
            lmbda_2.append(lmbda2)

        lmbda_2 = lmbda_2[1:]

        for i in reversed(lmbda_2):
            lambdas.append(i)

        for i in lmbda_1:
            lambdas.append(i)

    else:
        for i in range(0, windows + 1):
            lmbda = f'{sigmoid(float(i)/float(windows), k):.3f}'
            lambdas.append(lmbda)

    lambdas = lambdas[::-1]
    return lambdas

def write_submitfile(writedir, replacements):
    submit_in = CONFIGS['ROOT_DIR'] + '/INPUTS/FEP_submit.sh'
    submit_out = writedir + ('/FEP_submit.sh')
    with open(submit_in) as infile, open (submit_out, 'w') as outfile:
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
                    error_body = ''.join([next(infile) for _ in range(2)])
                    raise OSError(f'QFEP ERROR !! {error_main}\n{error_body}')
            except StopIteration as e:
                logger.info('Reached the end of the file before capturing the full error body.')
                raise OSError(f'QFEP ERROR !! {error_main}') from e

            line = line.split()
            if len(line) > 3:

                if line[3] == 'Free':
                    block = 1

                if line[3] == 'Termodynamic':
                    #continue
                    block = 2

                if line[3] == 'Overlap':
                    block = 3

                if line[3] == 'BAR':
                    block = 4

                if line[3] == 'Reaction':
                    block = 0

            if len(line) > 1:
                if block == 1:
                    if line[0] == '1.000000':
                        Zwanzig_r = float(line[4])

                    elif line[0] == '0.000000':
                        Zwanzig_f = float(line[2])

                        Zwanzig = np.nan if line[5] == '-Infinity' else float(line[5])

                if block == 2 and line[0] == '0.000000':
                    try:
                        thermo_integration = line[2]
                        if line[2] == '-Infinity':
                            thermo_integration = np.nan
                    except IndexError:
                        thermo_integration = np.nan # TODO: this line is never reached... # noqa: F841

                if block == 3 and line[0] == '0.000000':
                    overlap_sampling = np.nan if line[2] == '-Infinity' else float(line[2])

                if block == 4 and line[0] == '0.000000':
                    bar = np.nan if line[2] == '-Infinity' else float(line[2])

    return [Zwanzig, Zwanzig_f, Zwanzig_r, overlap_sampling, bar]

def read_qfep_verbose(qfep):
    """Reads a given qfep.out file and outputs a two-dimensional numpy array with the following structure:

    Args:
        qfep: qpfe.out file to be parsed.

    Raises:
        IOError: if the qfep.out file contains an error message and no content to be parsed.

    Returns:
        A two-dimensional numpy array with the following structure:
    >>> [[Zwanzig, dGfr, dGr, OS, BAR]   lambda 1
    >>>               ....               lambda ..
    >>>  [Zwanzig, dGfr, dGr, OS, BAR]]  lambda n
    """
    array = [[],[],[],[],[]]
    with open(qfep) as infile:
        block = 0
        for line in infile:
            if qfep_error_regex.findall(line):
                error_main = qfep_error_regex.findall(line)[0]
                # error_body = ''.join([next(infile) for i in range(1)])
                raise OSError(f'QFEP ERROR !! {error_main}\n')

            line = line.split()
            if len(line) > 3:
                if line[3] == 'Free':
                    block = 1

                if line[3] == 'Termodynamic':
                    #continue
                    block = 2

                if line[3] == 'Overlap':
                    block = 3

                if line[3] == 'BAR':
                    block = 4

                if line[3] == 'Reaction':
                    block = 0

            if len(line) > 1:
                if line[0] == '#':
                    continue

                if block == 1:
                    try:
                        array[0].append(float(line[5]))
                    except IndexError:
                        array[0].append(np.nan)
                    try:
                        array[1].append(float(line[4]))
                    except IndexError:
                        array[1].append(np.nan)
                    try:
                        array[2].append(float(line[2]))
                    except IndexError:
                        array[2].append(np.nan)

                if block == 3:
                    try:
                        array[3].append(float(line[2]))
                    except IndexError:
                        array[3].append(np.nan)

                if block == 4:
                    try:
                        array[4].append(float(line[2]))
                    except IndexError:
                        array[4].append(np.nan)
    arr = np.array(array)
    if np.isnan(arr).any():
        logger.warning('This function ran into exceptions; could use this file for debugging!')
    return arr.T
