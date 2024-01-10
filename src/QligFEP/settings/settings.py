from QligFEP import SRC
from multiprocessing import cpu_count

ROOT_DIR = SRC / "QligFEP"

# The directories to the input FF and run related input files are given here
FF_DIR = ROOT_DIR / "FF"
INPUT_DIR = ROOT_DIR / "INPUTS"
# Dicionary of locations of Q executables
Q_DIR = SRC / "q6"
BIN = Q_DIR / "bin/q6"

if not BIN.exists():
    raise FileNotFoundError(
        f"Could not find {BIN}.\n"
        "Please make sure you have compiled q6 by running "
        f"`make all && make mpi` in the {Q_DIR} directory."
    )

# some example schrodinger directory
SCHROD_DIR = '/mnt/c/Program\ Files/Schrodinger2022-2/'

# quick fix to run .exe files on wsl
EXE = '.exe'

Q_PATHS = {
    'QDYN'       : 'qdyn=' + str((BIN / 'qdynp').absolute()),
    'QPREP'      : str(BIN / 'qprep'),
    'QFEP'       : str(BIN / 'qfep'),
    'QCALC'      : str(BIN / 'qcalc'),
}

CONFIGS = {
    "ROOT_DIR" : str(ROOT_DIR),
    "FF_DIR" : str(FF_DIR),
    "INPUT_DIR" : str(INPUT_DIR),
    "Q_DIR" : str(Q_DIR),
    "BIN" : str(BIN),
    "SCHROD_DIR" : str(SCHROD_DIR),
    **Q_PATHS
}
# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '16',
       'TIME'         : '0-12:00:00',  # d-hh:mm:ss
       'MODULES'      : 'module purge\nmodule load Q/6.0.1_ME\n',
       **Q_PATHS
      }

# CLUSTER INPUTS. To add your own cluster, use the same input as below
SNELLIUS = {'NODES'        : '1',
            'NTASKS'       : '16',
            'TIME'         : '0-12:00:00',  # d-hh:mm:ss
            'MODULES'      : 'module load 2021\n module load gompi/2021a',
            **Q_PATHS
      }

ALICE = {'MAINDIR'      : Q_DIR,
         'NODES'        : '1',
         'NTASKS'       : '24',
         'TIME'         : '0-3:00:00',  # d-hh:mm:ss
         'MODULES'      : 'module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1',
         **Q_PATHS
        }


HEBBE = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-02:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load GCC/5.4.0-2.26\nmodule load OpenMPI/1.10.3\n', # Add a \n for every added module
         'ACCOUNT'    : 'SNIC2018-2-3',
         **Q_PATHS
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '28',
         'TIME'       : '0-04:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load gompi/2017b\n', # Add a \n for every added module
         'ACCOUNT'    : 'SNIC2018-2-3',
         **Q_PATHS
        }

STALLO = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-12:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28\n', # Add a \n for every added module
         'ACCOUNT'    : 'nn4654K',
         **Q_PATHS
        }

UPPMAX = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : 'gcc/9.2.0\nopenmpi/4.0.2\n', # Add a \n for every added module
         'ACCOUNT'    : 'snic2018-2-3',
         **Q_PATHS
        }

TETRA  = {'NODES'      : '1',
          'NTASKS'     : '32',
          'TIME'       : '0-24:00:00',  # d-hh:mm:ss
          'MODULES'    : '\n', # Add a \n for every added module
          'ACCOUNT'    : 'snic2022-3-2',
          **Q_PATHS
        }

LOCAL = {'NODES'      : 1,
         'NTASKS'     : cpu_count(),
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         **Q_PATHS
        }

CLUSTER_DICT = {
    "CSB" : CSB,
    "SNELLIUS" : SNELLIUS,
    "ALICE" : ALICE,
    "HEBBE" : HEBBE,
    "KEBNE" : KEBNE,
    "STALLO" : STALLO,
    "UPPMAX" : UPPMAX,
    "TETRA" : TETRA,
    "LOCAL" : LOCAL
}