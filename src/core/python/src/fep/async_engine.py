"""An engine in its own process, with a run() that does not wait.

    engines = [AsyncEngine(config) for _ in range(16)]     # 16 processes, started at the same time
    for engine, state in zip(engines, states):
        engine.set_state(state)
        engine.update_config(steps=500, lambdas=...)
    handles = [engine.run() for engine in engines]         # returns at once
    results = [handle.result() for handle in handles]      # RunResult of each
    for engine in engines:
        engine.close()

An AsyncEngine has the methods of qgpu.Engine. All but run() wait for their answer, and they
first wait for a run of this engine that is not over yet. run() returns a
concurrent.futures.Future. At most MAX_RUNNING runs (of all engines together) are under way at
a time; the others wait in a queue.

Each engine is a process of its own, so a crash of one does not take the main process or the
other engines with it: the call that was waiting raises an error. On one GPU the processes
overlap only with NVIDIA MPS running; without it they take turns.
"""
import multiprocessing
import os
import traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np

MAX_RUNNING = 10
_runner = None      # the threads that wait for the runs under way; made at the first run()


def _queue():
    global _runner
    if _runner is None:
        _runner = ThreadPoolExecutor(max_workers=MAX_RUNNING, thread_name_prefix="engine-run")
    return _runner


def _serve(pipe, config, device, precision, log):
    """The process of one engine: build it, then answer calls until "close" or the pipe ends."""
    # The engine prints its energies to standard output: send that to the log file (or nowhere).
    target = os.open(log if log else os.devnull, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.dup2(target, 1)
    try:
        from qgpu import Engine
        engine = Engine(config, device=device, precision=precision)
        pipe.send(("ok", engine.n_atoms))
    except BaseException as error:
        pipe.send(("error", error, traceback.format_exc()))
        return
    while True:
        try:
            name, args, kwargs = pipe.recv()
        except EOFError:        # the main process is gone
            return
        if name == "close":
            return
        try:
            if name == "_coordinates":     # only the atoms asked for: the whole state is large
                pipe.send(("ok", np.asarray(engine.get_state()["coords"])[args[0]].astype(np.float32)))
                continue
            value = getattr(engine, name)
            pipe.send(("ok", value(*args, **kwargs) if callable(value) else value))
        except BaseException as error:
            pipe.send(("error", error, traceback.format_exc()))


class EngineDied(RuntimeError):
    """The process of an engine ended without an answer (for example a crash in the engine)."""


class AsyncEngine:
    def __init__(self, config, device=0, precision="spfp", log=None):
        """Start the process. This returns at once; the first call waits until the engine is built.

        log: file for what the engine prints; None throws it away.
        """
        spawn = multiprocessing.get_context("spawn")    # a GPU connection does not survive a fork
        self._pipe, theirs = spawn.Pipe()
        self._process = spawn.Process(target=_serve, args=(theirs, config, device, precision, log),
                                      daemon=True)      # ends with the main process
        self._process.start()
        theirs.close()
        self._built = False
        self._pending = None    # the run that may still be under way
        self._n_atoms = None

    def _answer(self):
        try:
            status, *value = self._pipe.recv()
        except (EOFError, ConnectionError):
            raise EngineDied(f"engine process {self._process.pid} died") from None
        if status == "error":
            error, text = value
            raise error from RuntimeError("in the engine process:\n" + text)
        return value[0]

    def _send(self, name, *args, **kwargs):
        """Call a method of the engine in its process and wait for the answer."""
        if not self._built:
            self._n_atoms = self._answer()
            self._built = True
        try:
            self._pipe.send((name, args, kwargs))
        except (BrokenPipeError, OSError):
            raise EngineDied(f"engine process {self._process.pid} died") from None
        return self._answer()

    def _call(self, name, *args, **kwargs):
        if self._pending is not None:
            pending, self._pending = self._pending, None
            pending.exception()     # wait for it; its error belongs to whoever holds the handle
        return self._send(name, *args, **kwargs)

    def run(self):
        """Start a run and return a handle; handle.result() waits and gives the RunResult."""
        previous, self._pending = self._pending, None

        def job():
            if previous is not None:
                previous.exception()
            return self._send("run")

        self._pending = _queue().submit(job)
        return self._pending

    def update_config(self, **changes):
        return self._call("update_config", **changes)

    def set_state(self, state):
        return self._call("set_state", state)

    def get_state(self):
        return self._call("get_state")

    def load_restart(self, path):
        return self._call("load_restart", str(path))

    def coordinates(self, atoms):
        """Coordinates of the atoms (0-based indices) in the current state, (n, 3) float32."""
        return self._call("_coordinates", np.asarray(atoms, dtype=np.int64))

    @property
    def config(self):
        return self._call("config")

    @property
    def n_atoms(self):
        if self._n_atoms is None:
            self._call("n_atoms")
        return self._n_atoms

    def close(self):
        """End the process. A run under way is finished first."""
        if self._process.is_alive():
            try:
                if self._pending is not None:
                    self._pending.exception()
                self._pipe.send(("close", (), {}))
            except (EngineDied, OSError):
                pass
            self._process.join(timeout=10)
            if self._process.is_alive():
                self._process.terminate()
        self._pipe.close()

    def __enter__(self):
        return self

    def __exit__(self, *error):
        self.close()
