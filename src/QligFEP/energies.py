"""All lambda windows of a replicate in one energy file, instead of one ``.en`` file per window.

The plan runner (``qligfep_run``) gets the energies of every window from the engine in
memory and writes them to ``energies.parquet`` in the replicate folder: one row per
(window, frame, state), with the window's name and lambdas, the frame (from 1), the state,
and then the energy terms qdyn writes (``qgpu.fortran.ENERGY_TERMS``)::

    import pyarrow.parquet as pq
    table = pq.read_table("energies.parquet").to_pandas()
    table[table.window == "md_0500_0500"]

(``pandas.read_parquet`` works as well; on Habrok's qligfep_new environment it aborts
when Python exits, after reading.)

The file is lossless: ``qligfep_energies export energies.parquet --to DIR`` writes the
``md_XXXX_YYYY.en`` files back, byte for byte as qdyn wrote them, so ``qfep`` can be run on
them at any time (``qligfep_energies qfep`` does that in a temporary folder).

Off-diagonal elements (EVB couplings; QligFEP's FEP files have none) are kept in the
column ``offdiag``, on the row of a frame's first state.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

FORMAT_VERSION = 1
METADATA_KEY = b"qligfep.energies"
#: Name of the energy file in a replicate folder.
FILE_NAME = "energies.parquet"


class EnergyStoreError(ValueError):
    """An energy file cannot be written or read as asked."""


@dataclass(frozen=True)
class Window:
    """The energies of one lambda window."""

    name: str  # md_0500_0500: the .en file name without suffix
    lambdas: tuple[float, ...]
    energies: object  # qgpu.fortran.Energies


def _fortran():
    from .run_plan import _import_qgpu

    _import_qgpu()
    from qgpu import fortran

    return fortran


def _pyarrow():
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        raise EnergyStoreError("The energy file needs pyarrow: pip install pyarrow") from None
    return pa, pq


def _table(windows: list[Window]):
    import numpy as np

    pa, _ = _pyarrow()
    fortran = _fortran()
    n_lambdas = max((len(w.lambdas) for w in windows), default=0)
    with_offdiag = any(w.energies.n_offdiag for w in windows)
    columns: dict[str, list] = {"window": [], **{f"lambda{k + 1}": [] for k in range(n_lambdas)}}
    columns.update({name: [] for name, _ in fortran.q_energies_dtype()})
    offdiag: list = []
    for window in windows:
        rows = window.energies.q_energies
        columns["window"].append(np.full(len(rows), window.name, dtype=object))
        for k in range(n_lambdas):
            value = window.lambdas[k] if k < len(window.lambdas) else np.nan
            columns[f"lambda{k + 1}"].append(np.full(len(rows), value))
        for name, _ in fortran.q_energies_dtype():
            columns[name].append(rows[name])
        if with_offdiag:
            energies = window.energies
            cells = energies.offdiag.reshape(energies.n_frames, energies.n_offdiag)
            for position, row in enumerate(rows):
                if position % energies.n_states == 0:  # a frame's first state
                    frame = cells[row["frame"] - 1]
                    offdiag.append([{n: cell[n].item() for n in ("i", "j", "Hij", "rkl")} for cell in frame])
                else:
                    offdiag.append(None)
    arrays = {}
    for name, parts in columns.items():
        values = np.concatenate(parts) if parts else np.zeros(0)
        if name == "window":
            arrays[name] = pa.array(values.tolist() if len(values) else [], pa.string()).dictionary_encode()
        else:
            arrays[name] = pa.array(values)
    if with_offdiag:
        kind = pa.list_(pa.struct([("i", pa.int32()), ("j", pa.int32()), ("Hij", pa.float64()), ("rkl", pa.float64())]))
        arrays["offdiag"] = pa.array(offdiag, kind)
    metadata = {
        "version": FORMAT_VERSION,
        "terms": list(fortran.ENERGY_TERMS),
        "windows": [
            {"name": w.name, "lambdas": list(w.lambdas), "frames": w.energies.n_frames,
             "states": w.energies.n_states, "offdiag": w.energies.n_offdiag}
            for w in windows
        ],  # fmt: skip
    }
    return pa.table(arrays).replace_schema_metadata({METADATA_KEY: json.dumps(metadata).encode()})


def write_energies(path, windows: Iterable[Window]) -> Path:
    """Write the windows, in the given order, to one parquet file.

    The file is written next to ``path`` under a temporary name and then renamed, so
    ``path`` is always either absent, the previous complete file, or the new complete one.
    """
    _, pq = _pyarrow()
    windows = list(windows)
    names = [w.name for w in windows]
    if len(set(names)) != len(names):
        raise EnergyStoreError("Two windows have the same name.")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial")
    pq.write_table(_table(windows), partial, compression="zstd")
    os.replace(partial, path)
    return path


def _metadata(table) -> dict:
    raw = (table.schema.metadata or {}).get(METADATA_KEY)
    if raw is None:
        raise EnergyStoreError("Not a QligFEP energy file (no window index in its metadata).")
    metadata = json.loads(raw)
    if metadata.get("version") != FORMAT_VERSION:
        raise EnergyStoreError(f"Energy file format {metadata.get('version')} is not supported here.")
    return metadata


def read_energies(path, names: Iterable[str] | None = None) -> list[Window]:
    """The windows stored in an energy file, in the order they were run (or of ``names``)."""
    import numpy as np

    _, pq = _pyarrow()
    fortran = _fortran()
    table = pq.read_table(path)
    metadata = _metadata(table)
    index = {w["name"]: w for w in metadata["windows"]}
    names = list(index) if names is None else list(names)
    unknown = [name for name in names if name not in index]
    if unknown:
        raise EnergyStoreError(f"{path} has no window {', '.join(unknown)}.")

    window_column = np.asarray(table.column("window").to_pylist(), dtype=object)
    data = {name: table.column(name).to_numpy() for name, _ in fortran.q_energies_dtype()}
    offdiag = table.column("offdiag").to_pylist() if "offdiag" in table.column_names else None
    windows = []
    for name in names:
        info = index[name]
        selected = np.flatnonzero(window_column == name)
        rows = np.zeros(len(selected), fortran.q_energies_dtype())
        for column, values in data.items():
            rows[column] = values[selected]
        cells = np.zeros(info["frames"] * info["offdiag"], fortran.offdiag_dtype())
        if info["offdiag"]:
            flat = [cell for i in selected if offdiag[i] is not None for cell in offdiag[i]]
            for k, cell in enumerate(flat):
                cells[k] = (k // info["offdiag"] + 1, cell["i"], cell["j"], cell["Hij"], cell["rkl"])
        energies = fortran.Energies(rows, cells, info["states"], info["offdiag"])
        if energies.n_frames != info["frames"]:
            raise EnergyStoreError(f"{path}: window {name} has {energies.n_frames} frames, its index says {info['frames']}.")
        windows.append(Window(name, tuple(info["lambdas"]), energies))
    return windows


def export(path, folder, names: Iterable[str] | None = None) -> list[Path]:
    """Write ``<window>.en`` files into ``folder``, byte for byte as qdyn wrote them."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for window in read_energies(path, names):
        target = folder / f"{window.name}.en"
        target.write_bytes(window.energies.to_bytes())
        written.append(target)
    return written


def run_qfep(path, qfep, qfep_input, outdir, timeout: float | None = 180) -> None:
    """Run qfep on the windows of an energy file and put ``qfep.out`` and ``energies.csv`` in ``outdir``.

    The ``.en`` files exist only in a temporary folder while qfep runs. ``qfep_input`` names
    them as qligfep wrote it (inputfiles/qfep.inp). As in the run scripts, a qfep that does
    not finish within ``timeout`` seconds is stopped and what it wrote is kept.
    """
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    qfep_input = Path(qfep_input).resolve()
    with tempfile.TemporaryDirectory(prefix="qligfep_qfep_") as scratch:
        export(path, scratch)
        with open(qfep_input) as stdin, open(Path(scratch) / "qfep.out", "w") as stdout:
            try:
                subprocess.run([str(qfep)], stdin=stdin, stdout=stdout, cwd=scratch, timeout=timeout, check=True)
            except subprocess.TimeoutExpired:
                print(f"WARNING: qfep did not finish within {timeout:g} s; keeping its output so far.", file=sys.stderr)
        for name in ("qfep.out", "energies.csv"):
            if (Path(scratch) / name).exists():
                shutil.move(Path(scratch) / name, outdir / name)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="qligfep_energies",
        description="Work with the energy file of a replicate (energies.parquet): all lambda windows in one file.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="List the windows: name, lambdas, frames.")
    listing.add_argument("file")
    exporting = commands.add_parser("export", help="Write the md_XXXX_YYYY.en files, byte for byte as qdyn wrote them.")
    exporting.add_argument("file")
    exporting.add_argument("--to", required=True, help="Folder for the .en files (created if needed).")
    exporting.add_argument("--window", nargs="+", help="Only these windows (names without .en).")
    qfep = commands.add_parser("qfep", help="Run qfep on the windows; writes qfep.out and energies.csv.")
    qfep.add_argument("file")
    qfep.add_argument("--qfep", default="qfep", help="The qfep executable. Defaults to qfep on PATH.")
    qfep.add_argument("--input", required=True, help="The qfep input, e.g. ../../../inputfiles/qfep.inp.")
    qfep.add_argument("--out", default=".", help="Folder for qfep.out and energies.csv. Defaults to here.")
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            for window in read_energies(args.file):
                lambdas = " ".join(f"{x:g}" for x in window.lambdas)
                print(f"{window.name}  lambdas {lambdas}  {window.energies.n_frames} frames")
        elif args.command == "export":
            written = export(args.file, args.to, args.window)
            print(f"Wrote {len(written)} energy files to {args.to}")
        else:
            run_qfep(args.file, args.qfep, args.input, args.out)
    except (EnergyStoreError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


def main_exe() -> None:
    sys.exit(main())


if __name__ == "__main__":
    main_exe()
