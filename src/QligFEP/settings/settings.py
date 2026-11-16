import os
from collections.abc import Mapping
from pathlib import Path

from QligFEP import SRC


def nljoin(list_strings):  # nl for new line
    return "\n".join(list_strings) + "\n"


ROOT_DIR = SRC / "QligFEP"

# The directories to the input FF and run related input files are given here
FF_DIR = ROOT_DIR / "FF"
INPUT_DIR = ROOT_DIR / "INPUTS"
# Dicionary of locations of Q executables
Q_DIR = SRC / "q6"
# $QLIGFEP_Q_BIN points this checkout at binaries compiled elsewhere (another clone, a
# shared install); by default each checkout uses the engine built inside it.
BIN = Path(os.environ["QLIGFEP_Q_BIN"]).expanduser() if os.environ.get("QLIGFEP_Q_BIN") else Q_DIR / "bin/q6"

if not BIN.exists():
    raise FileNotFoundError(
        f"Could not find {BIN}.\n"
        "The qprep binary was not compiled during installation.\n"
        "This can happen if gfortran was not available during pip install.\n\n"
        "To fix this, please:\n"
        "  1. Install gfortran:\n"
        "     micromamba install gfortran -c conda-forge\n"
        "  2. Recompile qprep manually:\n"
        f"     cd {Q_DIR}\n"
        "     make qprep COMP=gcc\n"
        "     mkdir -p bin/q6 && mv qprep bin/q6/\n\n"
        "For HPC users who need MPI support, also run:\n"
        f"     cd {Q_DIR}\n"
        "     make mpi COMP=gcc\n"
    )

# some example schrodinger directory
SCHROD_DIR = r"/mnt/c/Program\sFiles/Schrodinger2022-2/"

# quick fix to run .exe files on wsl
EXE = ".exe"

Q_PATHS = {
    # Windowed FEP run templates use $qdyn as the profile's primary engine (the MPI
    # qdynp on clusters, the serial qdyn on LOCAL).
    "QDYN": "qdyn=" + str((BIN / "qdynp").absolute()),
    # Non-equilibrium runs drive both engines by name: $qdynp (parallel) for the
    # fixed-lambda equilibration and $qdyn (serial) for the lambda switches, which select
    # NEQ mode from the input's [lambda_scaling] section. The standalone qdyn_neq binary
    # no longer exists -- the switch is the consolidated serial qdyn.
    "QDYNP": "qdynp=" + str((BIN / "qdynp").absolute()),
    "QDYN_NEQ": "qdyn=" + str((BIN / "qdyn").absolute()),
    "QPREP": str(BIN / "qprep"),
    "QFEP": str(BIN / "qfep"),
    "QCALC": str(BIN / "qcalc"),
}

# Bare paths for local execution (Python template handles formatting, not shell assignment)
Q_PATHS_LOCAL = {
    "QDYN": str((BIN / "qdyn").absolute()),
    "QPREP": str(BIN / "qprep"),
    "QFEP": str(BIN / "qfep"),
    "QCALC": str(BIN / "qcalc"),
}

Q_PATHS_LOCALP = {
    "QDYN": str((BIN / "qdynp").absolute()),
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


class _ProfileClusterDict(Mapping):
    """Read-only ``CLUSTER_DICT`` view over the YAML execution profiles.

    Cluster settings now live in ``QligFEP/config/defaults.yaml`` and the user/project
    configuration files (see ``qligfep_config``). This mapping keeps the historical
    ``CLUSTER_DICT[name]["NTASKS"]`` access working for code not yet ported to
    :class:`QligFEP.config.Profile`. It is resolved on first use so that importing this
    module never fails on a broken user configuration file.
    """

    def _config(self):
        from QligFEP.config import load_config

        return load_config()

    def __getitem__(self, name):
        from QligFEP.config import ConfigError

        try:
            return self._config().get_profile(name).as_legacy_dict()
        except ConfigError as error:
            if name not in self:
                raise KeyError(name) from error
            raise

    def __iter__(self):
        return iter(self._config().profile_names())

    def __len__(self):
        return len(self._config().profile_names())


CLUSTER_DICT = _ProfileClusterDict()
