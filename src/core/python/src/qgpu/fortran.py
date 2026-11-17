"""The qgpu Engine API on top of Fortran Q (qdyn / qdynp), optionally under MPI.

    with FortranEngine(config, qdyn="qdynp", mpi=("mpirun", "-np", "8")) as engine:
        engine.run()
        engine.update_config(lambdas=(0.5, 0.5), energy=10, energy_file="md_0500_0500.en")
        engine.run()                      # continues from the previous run's restart file
        engine.update_config(lambdas=(0.6, 0.4), energy_file="")
        engine.run()
        energies = engine.get_energies()  # the run's energies in memory; no file is kept

Fortran Q is a separate (MPI) program, so it cannot live inside the Python process the
way the CUDA engine does. Every run() starts qdyn on a short-lived .inp file, and the state
is carried from run to run in one Q restart file. Both are an implementation detail and are
removed again: the workflow is written with Config and update_config() only, and what
remains on disk is what the Config asked for (energy_file, trajectory_file, final_file).
Energies need no file: get_energies() returns those of the latest run().

This module is pure Python: it needs no native extension, CUDA or GPU.
"""

from collections.abc import Sequence
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys
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
    workdir: where the engine keeps its own files. Default: a temporary folder.
    keep_files: keep the .inp, log and restart file of every run in workdir, for debugging.
             By default a run leaves only the restart file the next run starts from; the
             files of a failed run are kept.
    persistent: True: one `qdyn --serve` process serves every run(); the state stays in its
             memory, and no restart file or process start comes between runs. False: one qdyn
             process per run(), the state carried in a restart file. "auto" (default): True
             when the qdyn binary supports --serve, else False. Results are the same.
             After a failed run() the process is gone; the next run() starts a new one from
             the state the failed run started from.

    qdyn's output is printed to stdout, as the CUDA engine's is.
    Use the engine in a `with` block, or call close(), to remove its files at the end.
    """

    def __init__(self, config: Config, qdyn="qdyn", mpi: Sequence[str] | str = (), workdir=None,
                 keep_files=False, persistent="auto"):
        if persistent not in (True, False, "auto"):
            raise InputError(f"persistent must be True, False or 'auto', not {persistent!r}")
        self._persistent = persistent
        self._server = None
        self._served_state = False  # the server's memory holds the current state
        self._banner = []
        self._n_processes = 0
        if not isinstance(config, Config):
            raise TypeError(f"config must be a qgpu.Config, got {type(config).__name__}")
        executable = shutil.which(str(Path(qdyn).expanduser()))
        if executable is None:
            raise InputError(f"qdyn executable not found: {qdyn}")
        self._command = [*(shlex.split(mpi) if isinstance(mpi, str) else map(str, mpi)),
                         str(Path(executable).resolve())]
        self._keep = keep_files
        self._tmp = None
        if workdir is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="qgpu_fortran_")
            workdir = self._tmp.name
        self._workdir = Path(workdir).expanduser().resolve()
        self._created = not self._workdir.exists()
        self._workdir.mkdir(parents=True, exist_ok=True)
        self._links = []

        # under its own name: qdyn writes the topology file name into trajectory headers
        self._topology = self._link(config.topology, Path(config.topology).name)
        self._fep = self._link(config.fep, "fep.fep") if config.fep else ""
        self._restart = self._link(config.restart, "start.re") if config.restart else ""
        self._n_atoms = _topology_atoms(Path(_abs(config.topology)))
        self._n_states = _fep_states(Path(_abs(config.fep))) if config.fep else 0
        self._n_runs = 0
        self._energies = None
        self.log_file = None  # log of the latest run(): only of a failed one, unless keep_files
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
        """Run config.steps steps from the current state.

        After an error the engine is still at the state before this run().
        """
        if self._persistent is not False:
            if self._server is None and not self._start_server():
                return self._run_process()  # "auto" with a qdyn that cannot serve
            return self._run_served()
        return self._run_process()

    def _run_process(self):
        """One run in a qdyn process of its own, carrying the state in a restart file."""
        config = self._config
        name = f"run_{self._n_runs + 1:04d}"
        outputs = {"final": f"{name}.re"}
        if (config.energy or 0) > 0:
            outputs["energy"] = f"{name}.en"
        if (config.trajectory or 0) > 0:
            outputs["trajectory"] = f"{name}.dcd"
        (self._workdir / f"{name}.inp").write_text(self._inp(config, outputs))

        log = self._workdir / f"{name}.log"
        try:
            process = subprocess.Popen([*self._command, f"{name}.inp"], cwd=self._workdir, text=True,
                                       errors="replace", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except OSError as e:
            raise InputError(f"Could not start {shlex.join(self._command)}: {e}") from e
        self._n_processes += 1
        # qdyn's output goes to this process's stdout as it is produced, like the CUDA
        # engine's, and to the log file of this run.
        lines = []
        with process, open(log, "w", errors="replace") as out:
            for line in process.stdout:
                lines.append(line)
                out.write(line)
                sys.stdout.write(line)
            sys.stdout.flush()
        status = process.returncode
        self.log_file = log
        text = "".join(lines)
        final = self._workdir / outputs["final"]
        if status != 0 or "ABNORMAL TERMINATION" in text or not final.is_file():
            # qdyn also writes a restart file when it dies: it is not used, the next
            # run() starts from the same restart file as this one
            raise _error(text, status, log)

        self._n_runs += 1
        self._advance(outputs["final"])
        self._energies = read_energy_file(self._workdir / outputs["energy"]) if "energy" in outputs else None
        if not self._keep:
            (self._workdir / f"{name}.inp").unlink()
            log.unlink()
            self.log_file = None
        for key, target in (("energy", config.energy_file), ("trajectory", config.trajectory_file)):
            if key in outputs:
                if target:
                    _place(shutil.move, self._workdir / outputs[key], target)
                else:
                    (self._workdir / outputs[key]).unlink()
        if config.final_file:
            _place(shutil.copyfile, final, config.final_file)

    def get_energies(self):
        """The energies the latest run() recorded (every `energy` steps), as an :class:`Energies`.

        Available whether or not the Config names an energy_file.
        """
        if self._energies is None:
            raise InputError("No energies: the latest run() did not record any (energy interval 0), "
                             "or there was no run() yet.")
        return self._energies

    def get_state(self):
        """Coordinates, velocities and water-polarisation state, as the CUDA engine returns them."""
        if self._server is not None and self._served_state:
            path = self._workdir / "get_state.re"
            self._command_server(f"save {path.name}")
            state = read_restart(path)
            path.unlink()
            return state
        if not self._restart:
            raise InputError("No state yet: call run() first, or start from Config(restart=...).")
        return read_restart(self._workdir / self._restart)

    def set_state(self, state):
        if len(state["coords"]) != self._n_atoms or len(state["velocities"]) != self._n_atoms:
            raise InputError(f"State must have {self._n_atoms} atoms.")
        self._n_runs += 1  # a new name: never overwrite a restart file that a state came from
        name = f"state_{self._n_runs:04d}.re"
        write_restart(self._workdir / name, state)
        if self._server is not None:
            self._command_server(f"load {name}")
            self._served_state = True
            (self._workdir / name).unlink()
            self._advance("")
        else:
            self._advance(name)

    @property
    def n_processes(self):
        """How many qdyn processes this engine has started so far."""
        return self._n_processes

    def close(self):
        """Remove the engine's own files. Outputs named in the Config are not touched."""
        self._stop_server()
        if self._tmp is not None:
            self._tmp.cleanup()
        elif not self._keep:
            # stage.inp: the input `qdyn --serve` reads, overwritten by every run
            for name in (self._restart, *self._links, "stage.inp"):
                if name:
                    (self._workdir / name).unlink(missing_ok=True)
            if self._created and not any(self._workdir.iterdir()):
                self._workdir.rmdir()
        self._restart = ""

    def __enter__(self):
        return self

    def __exit__(self, *error):
        # after an error the files stay, to see what happened
        if error[0] is None or self._tmp is not None:
            self.close()

    # ---- persistent mode: one `qdyn --serve` process for every run()

    def _start_server(self):
        """Start `qdyn --serve`. False (persistent="auto" only) if this qdyn cannot serve."""
        try:
            process = subprocess.Popen([*self._command, "--serve", "."], cwd=self._workdir, text=True,
                                       errors="replace", stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, bufsize=1)
        except OSError as e:
            raise InputError(f"Could not start {shlex.join(self._command)}: {e}") from e
        self._n_processes += 1
        banner = []
        for line in process.stdout:
            if line.startswith("@qdyn ready"):
                self._server, self._banner, self._served_state = process, banner, False
                return True
            banner.append(line)
        process.wait()
        if self._persistent is True:
            raise InputError(f"{self._command[-1]} cannot run as a persistent process (qdyn --serve): "
                             f"{''.join(banner[-5:]).strip() or 'no output'}. Rebuild Q, or use persistent=False.")
        self._persistent = False  # "auto": this qdyn predates --serve
        return False

    def _command_server(self, command, body=None, echo=None):
        """Send one command and wait for its answer. Returns qdyn's output in between."""
        process = self._server
        try:
            process.stdin.write(command + "\n")
            if body is not None:
                process.stdin.write(body.rstrip("\n") + "\n@end\n")
            process.stdin.flush()
        except BrokenPipeError:
            pass  # it died; its output below says why
        lines = []
        for line in process.stdout:
            if line.startswith("@qdyn done"):
                return lines
            if line.startswith("@qdyn error"):
                raise InputError(f"qdyn --serve: {line[len('@qdyn error'):].strip()}")
            lines.append(line)
            if echo is not None:
                echo(line)
        # the process ended: a stage died, or the launcher failed
        status = process.wait()
        self._server = None
        raise _ServerEnded("".join(lines), status)

    def _stop_server(self):
        process, self._server = self._server, None
        if process is None:
            return
        try:
            process.stdin.write("quit\n")
            process.stdin.flush()
            process.stdin.close()
        except (BrokenPipeError, ValueError):
            pass
        for _ in process.stdout:
            pass
        process.wait()

    def _run_served(self):
        config = self._config
        name = f"run_{self._n_runs + 1:04d}"
        outputs = {}
        if (config.energy or 0) > 0:
            outputs["energy"] = f"{name}.en"
        if (config.trajectory or 0) > 0:
            outputs["trajectory"] = f"{name}.dcd"
        # the process starts from memory once it holds a state, from the restart file before
        restart = "" if self._served_state else self._restart
        text = self._inp(config, outputs, restart=restart)
        if self._keep:
            (self._workdir / f"{name}.inp").write_text(text)

        log = self._workdir / f"{name}.log"
        with open(log, "w", errors="replace") as out:
            def echo(line):
                out.write(line)
                sys.stdout.write(line)

            for line in self._banner:  # the process's startup output belongs to its first run
                echo(line)
            self._banner = []
            try:
                lines = self._command_server("run", text, echo)
            except _ServerEnded as ended:
                sys.stdout.flush()
                self.log_file = log
                self._recover()
                raise _error(ended.text, ended.status, log) from None
            finally:
                sys.stdout.flush()
        text_out = "".join(lines)
        if "terminated normally" not in text_out:
            self.log_file = log
            raise _error(text_out, 0, log)

        self._n_runs += 1
        self._served_state = True
        self._advance("")  # the state is in the process now; a start file is no longer needed
        self._energies = read_energy_file(self._workdir / outputs["energy"]) if "energy" in outputs else None
        self.log_file = log
        if not self._keep:
            log.unlink()
            self.log_file = None
        for key, target in (("energy", config.energy_file), ("trajectory", config.trajectory_file)):
            if key in outputs:
                if target:
                    _place(shutil.move, self._workdir / outputs[key], target)
                else:
                    (self._workdir / outputs[key]).unlink()
        if config.final_file:
            target = Path(_abs(config.final_file))
            target.parent.mkdir(parents=True, exist_ok=True)
            self._command_server(f"save {target}")
        if self._keep:  # as in process mode: the restart file of every run
            self._command_server(f"save {name}.re")

    def _recover(self):
        """After a served run died: the state it started from, for the next process."""
        recovery = self._workdir / "recovery.re"
        if self._served_state:
            if not recovery.is_file():
                self._restart = ""
                self._served_state = False
                raise QgpuError("qdyn --serve ended without leaving the state of the last good run "
                                "(recovery.re); this engine cannot continue.")
            self._n_runs += 1
            name = f"state_{self._n_runs:04d}.re"
            recovery.rename(self._workdir / name)
            self._advance(name)
        self._served_state = False

    def _advance(self, restart):
        """`restart` is now the state; the file it replaces is no longer needed."""
        previous, self._restart = self._restart, restart
        if previous and previous != restart and not self._keep:
            (self._workdir / previous).unlink(missing_ok=True)

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
        self._links.append(name)
        return name

    def _inp(self, config: Config, outputs, restart=None):
        restart = self._restart if restart is None else restart
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

        files = {"topology": self._topology, "fep": self._fep, "restart": restart, **outputs}
        sections["files"] = {key: path for key, path in files.items() if path}
        for section, items in config.extra.items():
            sections.setdefault(section, {}).update({key: _value(v, None) for key, v in items.items()})

        lines = []
        order = (*_SECTION_ORDER, "files")
        for section in (*order, *(s for s in config.extra if s not in order)):
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


class _ServerEnded(Exception):
    """The `qdyn --serve` process ended while running a command."""

    def __init__(self, text, status):
        super().__init__(status)
        self.text, self.status = text, status


def _value(value, kind):
    if isinstance(value, str):
        return value
    if kind is bool or isinstance(value, bool):
        return "on" if value else "off"
    return repr(kind(value)) if kind else str(value)


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


# ---- Q energy files (.en): Fortran unformatted records (4-byte length before and after).
# No header; one frame every `energy` steps (md.f90 calls put_ene of nrgy.f90):
#   per state:  istate (int32), q_energies (15 float64)
#   then:       offdiag_save * n_offdiag (int32 i, int32 j, float64 Hij, float64 rkl), possibly empty
# qfep reads the same layout back with get_ene.

#: q_energies of nrgy.f90, flattened: lambda, total, q (bonded), el/vdw of qx, qq, qp, qw, restraint.
ENERGY_TERMS = ("lambda", "total", "q_bond", "q_angle", "q_torsion", "q_improper", "qx_el", "qx_vdw",
                "qq_el", "qq_vdw", "qp_el", "qp_vdw", "qw_el", "qw_vdw", "restraint")
_STATE_BYTES = 4 + 8 * len(ENERGY_TERMS)
_OFFDIAG_FIELDS = [("i", "<i4"), ("j", "<i4"), ("Hij", "<f8"), ("rkl", "<f8")]
_OFFDIAG_BYTES = 24


class Energies:
    """The energy records of one run, frame by frame.

    q_energies: structured array, one row per (frame, state): ``frame`` (from 1), ``state``
        (the state number qdyn wrote), then the q_energies terms of nrgy.f90 (ENERGY_TERMS).
    offdiag:    structured array, one row per (frame, off-diagonal element): ``frame``, then
        offdiag_save of nrgy.f90: ``i``, ``j``, ``Hij``, ``rkl``. Empty without EVB couplings.
    n_states, n_offdiag: records per frame.

    to_bytes() gives the .en file qdyn wrote, byte for byte.
    """

    def __init__(self, q_energies, offdiag, n_states, n_offdiag):
        self.q_energies, self.offdiag = q_energies, offdiag
        self.n_states, self.n_offdiag = int(n_states), int(n_offdiag)
        if self.n_states and len(q_energies) % self.n_states:
            raise InputError(f"{len(q_energies)} energy rows are not whole frames of {self.n_states} states.")
        if len(offdiag) != self.n_frames * self.n_offdiag:
            raise InputError(f"Expected {self.n_frames * self.n_offdiag} off-diagonal rows, got {len(offdiag)}.")

    @property
    def n_frames(self):
        return len(self.q_energies) // self.n_states if self.n_states else 0

    def to_bytes(self):
        import numpy as np

        n_frames = self.n_frames
        frames = np.zeros(n_frames, _frame_dtype(self.n_states, self.n_offdiag))
        rows = self.q_energies.reshape(n_frames, self.n_states)
        for k in range(self.n_states):
            record = frames[f"s{k}"]
            record["size"] = record["size_end"] = _STATE_BYTES
            for term in ("state", *ENERGY_TERMS):
                record[term] = rows[term][:, k]
        frames["offdiag"]["size"] = frames["offdiag"]["size_end"] = _OFFDIAG_BYTES * self.n_offdiag
        if self.n_offdiag:
            cells = self.offdiag.reshape(n_frames, self.n_offdiag)
            for name, _ in _OFFDIAG_FIELDS:
                frames["offdiag"]["data"][name] = cells[name]
        return frames.tobytes()

    def __eq__(self, other):
        import numpy as np

        return (isinstance(other, Energies) and (self.n_states, self.n_offdiag) == (other.n_states, other.n_offdiag)
                and np.array_equal(self.q_energies, other.q_energies) and np.array_equal(self.offdiag, other.offdiag))

    def __repr__(self):
        return f"Energies({self.n_frames} frames, {self.n_states} states, {self.n_offdiag} off-diagonal)"


def q_energies_dtype():
    """Row type of Energies.q_energies."""
    return [("frame", "<i4"), ("state", "<i4"), *((term, "<f8") for term in ENERGY_TERMS)]


def offdiag_dtype():
    """Row type of Energies.offdiag."""
    return [("frame", "<i4"), *_OFFDIAG_FIELDS]


def _frame_dtype(n_states, n_offdiag):
    """One frame of a .en file, record markers included (packed, little-endian)."""
    state = [("size", "<i4"), ("state", "<i4"), *((term, "<f8") for term in ENERGY_TERMS), ("size_end", "<i4")]
    offdiag = [("size", "<i4"), ("data", _OFFDIAG_FIELDS, (n_offdiag,)), ("size_end", "<i4")]
    return [*((f"s{k}", state) for k in range(n_states)), ("offdiag", offdiag)]


def parse_energies(data, name="energy file"):
    """The content of a Q .en file as :class:`Energies`.

    Raises InputError unless the content is whole frames that to_bytes() writes back unchanged.
    """
    import numpy as np

    data = bytes(data)
    if not data:  # a run shorter than one energy interval
        return Energies(np.zeros(0, q_energies_dtype()), np.zeros(0, offdiag_dtype()), 0, 0)
    # the layout of the first frame: state records up to the off-diagonal record
    position, n_states = 0, 0
    while True:
        if position + 4 > len(data):
            raise InputError(f"Truncated {name}: no off-diagonal record after the first states.")
        (size,) = struct.unpack_from("<i", data, position)
        if size != _STATE_BYTES:
            break
        n_states += 1
        position += size + 8
    if n_states == 0 or size < 0 or size % _OFFDIAG_BYTES:
        raise InputError(f"Not a Q energy file: {name} (record of {size} bytes).")
    n_offdiag = size // _OFFDIAG_BYTES
    dtype = np.dtype(_frame_dtype(n_states, n_offdiag))
    if len(data) % dtype.itemsize:
        raise InputError(f"Truncated or irregular {name}: {len(data)} bytes are not whole frames "
                         f"of {dtype.itemsize} bytes.")
    frames = np.frombuffer(data, dtype)
    n_frames = len(frames)

    q_energies = np.zeros(n_frames * n_states, q_energies_dtype())
    rows = q_energies.reshape(n_frames, n_states)
    rows["frame"] = np.arange(1, n_frames + 1)[:, None]
    for k in range(n_states):
        record = frames[f"s{k}"]
        for term in ("state", *ENERGY_TERMS):
            rows[term][:, k] = record[term]
    offdiag = np.zeros(n_frames * n_offdiag, offdiag_dtype())
    if n_offdiag:
        cells = offdiag.reshape(n_frames, n_offdiag)
        cells["frame"] = np.arange(1, n_frames + 1)[:, None]
        for field, _ in _OFFDIAG_FIELDS:
            cells[field] = frames["offdiag"]["data"][field]
    energies = Energies(q_energies, offdiag, n_states, n_offdiag)
    if energies.to_bytes() != data:  # record markers, or a layout that changes between frames
        raise InputError(f"Irregular {name}: it cannot be read back losslessly.")
    return energies


def read_energy_file(path):
    """A Q .en file as :class:`Energies`."""
    return parse_energies(Path(path).read_bytes(), str(path))


def write_energy_file(path, energies):
    """Write :class:`Energies` as the .en file qdyn writes."""
    Path(path).write_bytes(energies.to_bytes())

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


# Sections that are not key-value settings, or whose keys are all Config fields.
_LIST_SECTIONS = frozenset({"lambdas", "trajectory_atoms", *_RESTRAINT_SECTIONS.values()})


def read_inp(path):
    """The Config values in a Q .inp file. Settings that the file leaves out are not returned."""
    path = Path(path).expanduser().resolve()
    return _values(_sections(path), path.parent)


def read_sections(sections, folder="."):
    """As read_inp, for .inp content given as {section: {key: value}} (see Config.from_sections)."""
    rows = {}
    for section, content in sections.items():
        name = _SECTION_ALIASES.get(str(section).lower(), str(section).lower())
        if name == "lambdas":
            content = [list(content)]
        elif name == "trajectory_atoms":
            content = [str(content).split()]
        elif isinstance(content, dict):
            content = [[str(key), *_value(value, None).split()] for key, value in content.items()]
        rows.setdefault(name, []).extend([str(x) for x in row] for row in content)
    return _values(rows, Path(folder).expanduser().resolve())


def _values(sections, folder):
    keys = {(section, key): (name, kind) for name, (section, key, kind) in _INP_KEYS.items()}
    values, extra = {}, {}
    for section, rows in sections.items():
        if section in _LIST_SECTIONS:
            continue
        for row in rows:
            key = _KEY_ALIASES.get(row[0].lower(), row[0].lower())
            if (section, key) in keys and len(row) > 1:
                name, kind = keys[section, key]
                values[name] = _BOOL[row[1].lower()] if kind is bool else kind(row[1])
            elif section == "md" and key == "constraint_algorithm" and len(row) > 2:
                values["solute_constraint_algorithm"] = row[1].lower()
                values["solvent_constraint_algorithm"] = row[2].lower()
            elif section == "files" and key in _FILE_FIELDS and len(row) > 1:
                values[_FILE_FIELDS[key]] = str((folder / row[1]).resolve())
            elif len(row) > 1:
                # no Config field: kept for Fortran Q (other files in [files] by their full path)
                value = str((folder / row[1]).resolve()) if section == "files" else " ".join(row[1:])
                extra.setdefault(section, {})[key] = value
    if extra:
        values["extra"] = extra
    if sections.get("lambdas"):
        values["lambdas"] = tuple(float(x) for row in sections["lambdas"] for x in row)
    if sections.get("trajectory_atoms"):
        values["trajectory_atoms"] = " ".join(sections["trajectory_atoms"][0])
    for name in _RESTRAINTS:
        values[name] = tuple(_read_restraint(name, row) for row in sections.get(_RESTRAINT_SECTIONS[name], []))
    return values


def normalise_extra(extra):
    """Config.extra with lower-case names; keywords that have a Config field are refused."""
    named = {(section, key) for section, key, _ in _INP_KEYS.values()}
    named |= {("md", "constraint_algorithm"), *(("files", key) for key in _FILE_FIELDS)}
    result = {}
    for section, items in extra.items():
        section = _SECTION_ALIASES.get(str(section).lower(), str(section).lower())
        if section in _LIST_SECTIONS or not isinstance(items, dict):
            raise InputError(f"extra: [{section}] is not a section of key-value settings.")
        for key, value in items.items():
            key = str(key).lower()
            if (section, _KEY_ALIASES.get(key, key)) in named:
                raise InputError(f"extra: [{section}] {key} has its own Config field; set that instead.")
            result.setdefault(section, {})[key] = value
    return result


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
