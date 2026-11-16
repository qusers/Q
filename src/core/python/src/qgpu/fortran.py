"""The qgpu Engine API on top of Fortran Q (qdyn / qdynp), optionally under MPI.

    engine = FortranEngine(config, qdyn="qdynp", mpi=("mpirun", "-np", "8"))
    engine.run()
    engine.update_config(lambdas=(0.5, 0.5), energy=10, energy_file="md_0500_0500.en")
    engine.run()                          # continues from the previous run's restart file

Fortran Q is a separate (MPI) program, so it cannot live inside the Python process the
way the CUDA engine does. Every run() writes one .inp file in the work folder and starts
qdyn on it; the state is carried from run to run in Q restart files. The .inp files are
an implementation detail: the workflow is written with Config and update_config() only.

This module is pure Python: it needs no native extension, CUDA or GPU.
"""

from collections.abc import Sequence
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import tempfile

from .config import _INP_KEYS, _RESTRAINTS, Config, _abs
from .engine import Engine
from .errors import ConstraintError, InputError, QgpuError

# Config field -> the key Fortran Q reads, where it differs from the Config name's .inp key.
_FORTRAN_KEYS = {"polarisation": "polarization", "polarisation_force": "polarization_force"}
# Config field -> the only value Fortran Q supports (None = unset is always fine).
_UNSUPPORTED = {"thermostat": "berendsen", "hmr": False}
# Not Fortran keys: validated with _UNSUPPORTED, never written.
_SKIPPED = frozenset({"thermostat", "hmr", "hmr_target_mass"})
# Only the first run starts from Config.restart / new velocities; later runs continue.
_FIRST_RUN_ONLY = frozenset({"seed", "initial_temperature"})
_SECTION_ORDER = ("md", "cut-offs", "sphere", "solvent", "intervals")
_RESTRAINT_SECTIONS = {
    "seq_restraints": "sequence_restraints",
    "dist_restraints": "distance_restraints",
    "ang_restraints": "angle_restraints",
    "wall_restraints": "wall_restraints",
    "pos_restraints": "atom_restraints",
}


class FortranEngine(Engine):
    """Engine with the same interface, run by Fortran Q.

    qdyn:    the qdyn or qdynp executable (name on PATH, or a path).
    mpi:     launcher prefix, e.g. ("mpirun", "-np", "8") or "srun"; empty = run qdyn directly.
    workdir: where the .inp, log and restart files of every run are kept. Default: a
             temporary folder that is removed with the engine.
    """

    def __init__(self, config: Config, qdyn="qdyn", mpi: Sequence[str] | str = (), workdir=None):
        if not isinstance(config, Config):
            raise TypeError(f"config must be a qgpu.Config, got {type(config).__name__}")
        executable = shutil.which(str(Path(qdyn).expanduser()))
        if executable is None:
            raise InputError(f"qdyn executable not found: {qdyn}")
        self._command = [*(shlex.split(mpi) if isinstance(mpi, str) else map(str, mpi)),
                         str(Path(executable).resolve())]
        if workdir is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="qgpu_fortran_")
            workdir = self._tmp.name
        self._workdir = Path(workdir).expanduser().resolve()
        self._workdir.mkdir(parents=True, exist_ok=True)

        self._topology = self._link(config.topology, "topology.top")
        self._fep = self._link(config.fep, "fep.fep") if config.fep else ""
        self._restart = self._link(config.restart, "start.re") if config.restart else ""
        self._n_atoms = _topology_atoms(Path(_abs(config.topology)))
        self._n_states = _fep_states(Path(_abs(config.fep))) if config.fep else 0
        self._n_runs = 0
        self.log_file = None  # log of the latest run()
        self._configure(config)
        self._config = config

    @property
    def precision(self):
        return "double"

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def workdir(self):
        return self._workdir

    def run(self):
        """Run config.steps steps from the current state, in one qdyn process.

        After an error the engine is still at the state before this run().
        """
        config = self._config
        name = f"run_{self._n_runs + 1:04d}"
        outputs = {"final": f"{name}.re"}
        if (config.energy or 0) > 0:
            outputs["energy"] = f"{name}.en"
        if (config.trajectory or 0) > 0:
            outputs["trajectory"] = f"{name}.dcd"
        (self._workdir / f"{name}.inp").write_text(self._inp(config, outputs))

        log = self._workdir / f"{name}.log"
        with open(log, "w") as out:
            try:
                status = subprocess.run([*self._command, f"{name}.inp"], cwd=self._workdir,
                                        stdout=out, stderr=subprocess.STDOUT).returncode
            except OSError as e:
                raise InputError(f"Could not start {shlex.join(self._command)}: {e}") from e
        self.log_file = log
        text = log.read_text(errors="replace")
        final = self._workdir / outputs["final"]
        if status != 0 or "ABNORMAL TERMINATION" in text or not final.is_file():
            # qdyn also writes a restart file when it dies: it is not used, the next
            # run() starts from the same restart file as this one
            raise _error(text, status, log)

        self._n_runs += 1
        self._restart = outputs["final"]
        for key, target in (("energy", config.energy_file), ("trajectory", config.trajectory_file)):
            if key in outputs:
                _place(shutil.move, self._workdir / outputs[key], target)
        if config.final_file:
            _place(shutil.copyfile, final, config.final_file)

    def get_state(self):
        """Coordinates, velocities and water-polarisation state, as the CUDA engine returns them."""
        if not self._restart:
            raise InputError("No state yet: call run() first, or start from Config(restart=...).")
        return read_restart(self._workdir / self._restart)

    def set_state(self, state):
        if len(state["coords"]) != self._n_atoms or len(state["velocities"]) != self._n_atoms:
            raise InputError(f"State must have {self._n_atoms} atoms.")
        self._n_runs += 1  # a new name: never overwrite a restart file that a state came from
        name = f"state_{self._n_runs:04d}.re"
        write_restart(self._workdir / name, state)
        self._restart = name

    def _configure(self, config: Config):
        for name, supported in _UNSUPPORTED.items():
            if getattr(config, name) not in (None, supported):
                raise InputError(f"Fortran Q does not support {name}={getattr(config, name)!r}.")
        # qdyn itself runs with the wrong number of lambdas without complaining
        if config.lambdas is not None and len(config.lambdas) != self._n_states:
            raise InputError(f"Expected {self._n_states} lambdas, got {len(config.lambdas)}")
        if not config.restart and (config.initial_temperature is None or not config.seed):
            raise InputError("A fresh start (no restart) needs initial_temperature and a seed > 0.")

    def _link(self, source, name):
        """Give an input file a short name in the work folder (Fortran file names are limited)."""
        source = Path(_abs(source))
        if not source.is_file():
            raise InputError(f"File does not exist: {source}")
        link = self._workdir / name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(source)
        return name

    def _inp(self, config: Config, outputs):
        first = self._n_runs == 0 and self._restart in ("", "start.re")
        sections = {}
        for name, (section, key, kind) in _INP_KEYS.items():
            value = getattr(config, name)
            if value is None or name in _SKIPPED or (name in _FIRST_RUN_ONLY and not first):
                continue
            sections.setdefault(section, {})[_FORTRAN_KEYS.get(name, key)] = _value(value, kind)
        if config.solute_constraint_algorithm is not None:
            sections["md"]["constraint_algorithm"] = \
                f"{config.solute_constraint_algorithm} {config.solvent_constraint_algorithm}"

        files = {"topology": self._topology, "fep": self._fep, "restart": self._restart, **outputs}
        sections["files"] = {key: path for key, path in files.items() if path}

        lines = []
        for section in (*_SECTION_ORDER, "files"):
            if section in sections:
                lines.append(f"[{section}]")
                lines += [f"{key:<26}{value}" for key, value in sections[section].items()]
                lines.append("")
        if config.lambdas:
            lines += ["[lambdas]", " ".join(repr(x) for x in config.lambdas), ""]
        if config.trajectory_atoms and "trajectory" in outputs:
            lines += ["[trajectory_atoms]", config.trajectory_atoms, ""]
        for name, section in _RESTRAINT_SECTIONS.items():
            rows = [_restraint_row(r) for r in getattr(config, name)]
            if rows:
                lines += [f"[{section}]", *rows, ""]
        return "\n".join(lines)


def _value(value, kind):
    if isinstance(value, str):
        return value
    if kind is bool:
        return "on" if value else "off"
    return repr(kind(value))


def _restraint_row(r):
    """One restraint as an .inp row, in the column order Fortran Q reads (md.f90, initialize)."""
    kind = type(r).__name__
    if kind == "SeqRestraint":
        values = (r.ai, r.aj, r.k, int(r.ih), r.to_center)
    elif kind == "PosRestraint":
        values = (r.a, r.x.x, r.x.y, r.x.z, r.k.x, r.k.y, r.k.z, r.ipsi)
    elif kind == "DistRestraint":
        values = (r.ai, r.aj, r.d1, r.d2, r.k, r.ipsi)
    elif kind == "AngRestraint":
        values = (r.ai, r.aj, r.ak, r.ang, r.k, r.ipsi)
    else:
        values = (r.ai, r.aj, r.d, r.k, r.a_morse, r.d_morse, int(r.ih))
    return " ".join(repr(v) for v in values)


def _place(operation, source, target):
    target = Path(_abs(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    operation(source, target)


def _error(text, status, log):
    lines = [line.strip() for line in text.splitlines() if line.strip(" !=\t")]
    reasons = [line for line in lines if line.startswith(">>>") and "WARNING" not in line.upper()]
    for i, line in enumerate(lines):
        if "ABNORMAL TERMINATION" in line:
            reasons += lines[i + 1:i + 4]
            break
    message = "; ".join(reasons) or "\n".join(lines[-5:]) or "no output"
    message = f"qdyn failed (exit status {status}): {message} [log: {log}]"
    lowered = message.lower()
    if "invalid data in input" in lowered or "error opening input" in lowered:
        return InputError(message)
    if any(word in lowered for word in ("shake", "lincs", "settle", "constraint")):
        return ConstraintError(message)
    return QgpuError(message)


def _topology_atoms(path):
    with open(path, errors="replace") as f:
        for line in f:
            if "Total no. of atoms" in line:
                return int(line.split()[0])
    raise InputError(f"Could not read the number of atoms from topology {path}")


def _fep_states(path):
    """Number of states in a FEP file: `states` in its [FEP] section (default 1, as in Fortran Q)."""
    for row in _sections(path).get("fep", []):
        if row[0].lower() == "states" and len(row) > 1:
            return int(row[1])
    return 1


# ---- Q restart files: Fortran unformatted records (4-byte length before and after)
#   1. nat3 (int32), coordinates (nat3 float64)
#   2. nat3 (int32), velocities (nat3 float64)
#   3. n_shells (int32), theta_corr (n_shells float32)   only with water polarisation restraints
# (md.f90, write_xfin). Periodic-box records are not supported: QGPU is sphere-only too.

def _records(data):
    position = 0
    while position < len(data):
        (size,) = struct.unpack_from("<i", data, position)
        end = position + 4 + size
        if size <= 0 or end + 4 > len(data) or struct.unpack_from("<i", data, end)[0] != size:
            raise InputError("Invalid or truncated Q restart file.")
        yield data[position + 4:end]
        position = end + 4


def read_restart(path):
    """A Q restart file as {"coords": (n, 3), "velocities": (n, 3), "theta_corr": (n_shells,)}."""
    import numpy as np

    records = list(_records(Path(path).read_bytes()))
    if len(records) not in (2, 3):
        raise InputError(f"Unsupported Q restart file (periodic box?): {path}")
    vectors = []
    for record in records[:2]:
        (nat3,) = struct.unpack_from("<i", record)
        if nat3 <= 0 or nat3 % 3 or len(record) != 4 + 8 * nat3:
            raise InputError(f"Invalid Q restart file: {path}")
        vectors.append(np.frombuffer(record, "<f8", offset=4).reshape(-1, 3).copy())
    theta = np.frombuffer(records[2], "<f4", offset=4).astype(float) if len(records) == 3 else np.zeros(0)
    return {"coords": vectors[0], "velocities": vectors[1], "theta_corr": theta}


def write_restart(path, state):
    import numpy as np

    def record(count, values):
        payload = struct.pack("<i", count) + values.tobytes()
        size = struct.pack("<i", len(payload))
        return size + payload + size

    data = b""
    for key in ("coords", "velocities"):
        values = np.ascontiguousarray(state[key], "<f8").reshape(-1)
        data += record(len(values), values)
    theta = np.ascontiguousarray(state["theta_corr"], "<f4").reshape(-1)
    if len(theta):
        data += record(len(theta), theta)
    Path(path).write_bytes(data)


# ---- reading .inp files without the C++ parser (Config.from_inp on Fortran-only installs)

_BOOL = {"on": True, "off": False, "true": True, "false": False, ".true.": True, ".false.": False}
# Older spellings that Fortran Q still accepts.
_KEY_ALIASES = {"shake_solvent": "constrain_solvent", "shake_solute": "constrain_solute",
                "shake_hydrogens": "constrain_hydrogens", "polarization": "polarisation",
                "polarization_force": "polarisation_force"}
_SECTION_ALIASES = {"cut_offs": "cut-offs", "water": "solvent"}
_FILE_FIELDS = {"topology": "topology", "fep": "fep", "restart": "restart", "final": "final_file",
                "energy": "energy_file", "trajectory": "trajectory_file"}
# Restraint row columns as they are in the file -> constructor order of the restraint class.
_ROW_TYPES = {
    "seq_restraints": (int, int, float, int, int),
    "dist_restraints": (int, int, float, float, float, int),
    "ang_restraints": (int, int, int, float, float, int),
    "wall_restraints": (int, int, float, float, float, float, int),
    "pos_restraints": (int, float, float, float, float, float, float, int),
}


def _sections(path):
    sections, rows = {}, None
    for line in Path(path).read_text().splitlines():
        for mark in "!#*":
            line = line.split(mark)[0]
        line = line.strip()
        if line.startswith("["):
            name = line.strip("[]").strip().lower()
            rows = sections.setdefault(_SECTION_ALIASES.get(name, name), [])
        elif line and rows is not None:
            rows.append(line.split())
    return sections


def read_inp(path):
    """The Config values in a Q .inp file. Settings that the file leaves out are not returned."""
    path = Path(path).expanduser().resolve()
    sections = _sections(path)
    keys = {(section, key): (name, kind) for name, (section, key, kind) in _INP_KEYS.items()}
    values = {}
    for section, rows in sections.items():
        for row in rows:
            key = _KEY_ALIASES.get(row[0].lower(), row[0].lower())
            if (section, key) in keys and len(row) > 1:
                name, kind = keys[section, key]
                values[name] = _BOOL[row[1].lower()] if kind is bool else kind(row[1])
            elif section == "md" and key == "constraint_algorithm" and len(row) > 2:
                values["solute_constraint_algorithm"] = row[1].lower()
                values["solvent_constraint_algorithm"] = row[2].lower()
            elif section == "files" and key in _FILE_FIELDS and len(row) > 1:
                values[_FILE_FIELDS[key]] = str((path.parent / row[1]).resolve())
    if sections.get("lambdas"):
        values["lambdas"] = tuple(float(x) for row in sections["lambdas"] for x in row)
    if sections.get("trajectory_atoms"):
        values["trajectory_atoms"] = " ".join(sections["trajectory_atoms"][0])
    for name in _RESTRAINTS:
        values[name] = tuple(_read_restraint(name, row) for row in sections.get(_RESTRAINT_SECTIONS[name], []))
    return values


def _read_restraint(name, row):
    kinds = _ROW_TYPES[name]
    v = [kind(x) for kind, x in zip(kinds, row)]
    if name == "seq_restraints":
        return (*v[:3], bool(v[3]) if len(v) > 3 else False, v[4] if len(v) > 4 else 0)
    if len(v) != len(kinds):
        raise InputError(f"Invalid {_RESTRAINT_SECTIONS[name]} row: {' '.join(row)}")
    if name == "ang_restraints":      # file: i j k angle fc state -> AngRestraint(ai, aj, ak, ipsi, ang, k)
        return (*v[:3], v[5], v[3], v[4])
    if name == "pos_restraints":      # file: i x y z fx fy fz state -> PosRestraint(a, ipsi, x, k)
        return (v[0], v[7], tuple(v[1:4]), tuple(v[4:7]))
    if name == "wall_restraints":
        return (*v[:6], bool(v[6]))
    return tuple(v)
