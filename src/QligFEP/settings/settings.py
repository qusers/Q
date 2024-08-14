from multiprocessing import cpu_count

from QligFEP import SRC


def nljoin(list_strings):  # nl for new line
    return "\n".join(list_strings) + "\n"


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
SCHROD_DIR = r"/mnt/c/Program\sFiles/Schrodinger2022-2/"

# quick fix to run .exe files on wsl
EXE = ".exe"

Q_PATHS = {
    "QDYN": "qdyn=" + str((BIN / "qdynp").absolute()),
    "QPREP": str(BIN / "qprep"),
    "QFEP": str(BIN / "qfep"),
    "QCALC": str(BIN / "qcalc"),
}

CONFIGS = {
    "ROOT_DIR": str(ROOT_DIR),
    "FF_DIR": str(FF_DIR),
    "INPUT_DIR": str(INPUT_DIR),
    "Q_DIR": str(Q_DIR),
    "BIN": str(BIN),
    "SCHROD_DIR": str(SCHROD_DIR),
    **Q_PATHS,
}
# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {
    "NODES": "1",
    "NTASKS": "16",
    "TIME": "0-12:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module purge", "module load Q/6.0.1_ME"]),
    **Q_PATHS,
}

# CLUSTER INPUTS. To add your own cluster, use the same input as below
SNELLIUS = {
    "NODES": "1",
    "NTASKS": "16",
    "TIME": "0-12:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module load 2021", "module load gompi/2021a"]),
    **Q_PATHS,
}

ALICE = {
    "MAINDIR": Q_DIR,
    "NODES": "1",
    "NTASKS": "24",
    "TIME": "0-3:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1"]),
    **Q_PATHS,
}


HEBBE = {
    "NODES": "1",
    "NTASKS": "20",
    "TIME": "0-02:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module load GCC/5.4.0-2.26", "module load OpenMPI/1.10.3"]),
    "ACCOUNT": "SNIC2018-2-3",
    **Q_PATHS,
}

KEBNE = {
    "NODES": "1",
    "NTASKS": "28",
    "TIME": "0-04:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module load gompi/2017b"]),
    "ACCOUNT": "SNIC2018-2-3",
    **Q_PATHS,
}

STALLO = {
    "NODES": "1",
    "NTASKS": "20",
    "TIME": "0-12:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28"]),
    "ACCOUNT": "nn4654K",
    **Q_PATHS,
}

UPPMAX = {
    "NODES": "1",
    "NTASKS": "20",
    "TIME": "0-24:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["gcc/9.2.0", "openmpi/4.0.2"]),
    "ACCOUNT": "snic2018-2-3",
    **Q_PATHS,
}

TETRA = {
    "NODES": "1",
    "NTASKS": "8",
    "TIME": "0-06:00:00",  # d-hh:mm:ss
    # Commented out because we're not compiling Q during the job runtime
    "MODULES": "# module load buildenv-gcc/2022a-eb\n",
    "ACCOUNT": "naiss2023-3-5",
    **Q_PATHS,
}

DARDEL = {
    "NODES": "1",
    "NTASKS": "8",
    "TIME": "0-06:00:00",  # d-hh:mm:ss
    "MODULES": nljoin(["# module purge", "# module load cpe/23.12", "# module load PrgEnv-gnu/8.5.0"]),
    "ACCOUNT": "naiss2023-3-5",
    **Q_PATHS,
}

LOCAL = {
    "NODES": 1,
    "NTASKS": str(cpu_count()),
    "TIME": "0-24:00:00",  # d-hh:mm:ss
    "MODULES": "\n",
    **Q_PATHS,
}

CLUSTER_DICT = {
    "CSB": CSB,
    "SNELLIUS": SNELLIUS,
    "ALICE": ALICE,
    "HEBBE": HEBBE,
    "KEBNE": KEBNE,
    "STALLO": STALLO,
    "UPPMAX": UPPMAX,
    "TETRA": TETRA,
    "LOCAL": LOCAL,
}
