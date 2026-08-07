"""Read native Q energy files and consolidate them into one CSV file.

Q writes ``.en`` files as Fortran sequential-unformatted records.  QGPU's
``NativeOutput`` intentionally uses the same layout.  Each sampled frame is:

* one state record per lambda state: ``int32 state`` followed by the fifteen
  ``float64`` fields of ``q_energies``;
* one record containing zero or more ``offdiag_save`` values, each represented
  by two ``int32`` and two ``float64`` values.

The converter uses a long-form CSV so no binary data is discarded.  State and
off-diagonal records share identifying columns and leave unrelated columns
blank.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import struct
from collections.abc import Iterator, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO

_STATE_FORMAT = "i15d"
_OFFDIAG_FORMAT = "ii2d"
_STATE_SIZE = struct.calcsize("<" + _STATE_FORMAT)
_OFFDIAG_SIZE = struct.calcsize("<" + _OFFDIAG_FORMAT)

CSV_FIELDS = [
    "source_file",
    "frame",
    "record_type",
    "state",
    "lambda",
    "total",
    "q_bond",
    "q_angle",
    "q_torsion",
    "q_improper",
    "qx_el",
    "qx_vdw",
    "qq_el",
    "qq_vdw",
    "qp_el",
    "qp_vdw",
    "qw_el",
    "qw_vdw",
    "restraint",
    "offdiag_index",
    "offdiag_i",
    "offdiag_j",
    "offdiag_hij",
    "offdiag_rkl",
]

_STATE_VALUE_FIELDS = CSV_FIELDS[3:19]


class EnergyFormatError(ValueError):
    """Raised when a native energy file is malformed or unsupported."""


@dataclass(frozen=True)
class EnergyFileSummary:
    """Validation information collected while converting one energy file."""

    source_file: str
    frames: int
    states: int
    lambdas: tuple[float, ...]
    state_rows: int
    offdiag_rows: int
    byte_order: str
    record_marker_bytes: int


@dataclass(frozen=True)
class ConversionSummary:
    """Summary of an atomic multi-file CSV conversion."""

    output_file: str
    files: tuple[EnergyFileSummary, ...]

    @property
    def frames(self) -> int:
        return sum(item.frames for item in self.files)

    @property
    def rows(self) -> int:
        return sum(item.state_rows + item.offdiag_rows for item in self.files)

    def to_dict(self) -> dict:
        result = asdict(self)
        result["frames"] = self.frames
        result["rows"] = self.rows
        return result


@dataclass(frozen=True)
class _RecordFormat:
    endian: str
    marker_bytes: int

    @property
    def marker_format(self) -> str:
        return self.endian + ("i" if self.marker_bytes == 4 else "q")

    @property
    def byte_order_name(self) -> str:
        return "little" if self.endian == "<" else "big"


def _detect_record_format(stream: BinaryIO, source: Path) -> _RecordFormat:
    """Detect record-marker width and byte order from the first state record."""
    prefix = stream.read(8)
    stream.seek(0)
    if len(prefix) < 4:
        raise EnergyFormatError(f"{source}: empty or truncated energy file")

    file_size = source.stat().st_size
    candidates = (
        _RecordFormat("<", 4),
        _RecordFormat(">", 4),
        _RecordFormat("<", 8),
        _RecordFormat(">", 8),
    )
    for candidate in candidates:
        if len(prefix) < candidate.marker_bytes:
            continue
        length = struct.unpack(candidate.marker_format, prefix[: candidate.marker_bytes])[0]
        if length != _STATE_SIZE or candidate.marker_bytes * 2 + length > file_size:
            continue

        stream.seek(candidate.marker_bytes)
        payload = stream.read(length)
        trailer = stream.read(candidate.marker_bytes)
        stream.seek(0)
        if len(payload) != length or len(trailer) != candidate.marker_bytes:
            continue
        trailing_length = struct.unpack(candidate.marker_format, trailer)[0]
        if trailing_length != length:
            continue
        state = struct.unpack_from(candidate.endian + "i", payload)[0]
        if state == 1:
            return candidate

    raise EnergyFormatError(
        f"{source}: unsupported energy record layout; expected a Fortran record "
        f"containing int32 + 15 float64 values"
    )


def _iter_fortran_records(
    stream: BinaryIO,
    source: Path,
    record_format: _RecordFormat,
) -> Iterator[bytes]:
    """Yield validated sequential-unformatted record payloads."""
    marker_size = record_format.marker_bytes
    marker_format = record_format.marker_format
    record_number = 0

    while True:
        marker = stream.read(marker_size)
        if not marker:
            return
        record_number += 1
        if len(marker) != marker_size:
            raise EnergyFormatError(f"{source}: truncated leading marker for record {record_number}")

        length = struct.unpack(marker_format, marker)[0]
        if length < 0:
            raise EnergyFormatError(
                f"{source}: split Fortran subrecords are not supported (record {record_number})"
            )
        payload = stream.read(length)
        trailer = stream.read(marker_size)
        if len(payload) != length or len(trailer) != marker_size:
            raise EnergyFormatError(f"{source}: truncated record {record_number}")
        trailing_length = struct.unpack(marker_format, trailer)[0]
        if trailing_length != length:
            raise EnergyFormatError(
                f"{source}: mismatched record markers for record {record_number}: "
                f"{length} != {trailing_length}"
            )
        yield payload


def _record_rows(
    source: Path,
) -> tuple[Iterator[dict], EnergyFileSummary]:
    """Parse one file and return rows plus a summary.

    Rows are materialized per source file before being returned.  This ensures a
    malformed trailing frame cannot result in a successful-looking partial
    source conversion.
    """
    rows: list[dict] = []
    with source.open("rb") as stream:
        record_format = _detect_record_format(stream, source)
        records = iter(_iter_fortran_records(stream, source, record_format))

        frames = 0
        expected_states: int | None = None
        expected_lambdas: tuple[float, ...] | None = None
        state_rows = 0
        offdiag_rows = 0

        while True:
            try:
                payload = next(records)
            except StopIteration:
                break

            frames += 1
            frame_states = 0
            frame_lambdas: list[float] = []
            while len(payload) == _STATE_SIZE:
                values = struct.unpack(record_format.endian + _STATE_FORMAT, payload)
                state = values[0]
                if state != frame_states + 1:
                    raise EnergyFormatError(
                        f"{source}: frame {frames} expected state {frame_states + 1}, found {state}"
                    )
                row = {field: "" for field in CSV_FIELDS}
                row.update(
                    {
                        "source_file": source.name,
                        "frame": frames,
                        "record_type": "state",
                    }
                )
                row.update(zip(_STATE_VALUE_FIELDS, values, strict=True))
                rows.append(row)
                frame_lambdas.append(values[1])
                frame_states += 1
                state_rows += 1
                try:
                    payload = next(records)
                except StopIteration as exc:
                    raise EnergyFormatError(
                        f"{source}: frame {frames} has no trailing off-diagonal record"
                    ) from exc

            if frame_states == 0:
                raise EnergyFormatError(f"{source}: frame {frames} contains no state records")
            if expected_states is None:
                expected_states = frame_states
                expected_lambdas = tuple(frame_lambdas)
            elif frame_states != expected_states:
                raise EnergyFormatError(
                    f"{source}: frame {frames} has {frame_states} states; expected {expected_states}"
                )
            elif tuple(frame_lambdas) != expected_lambdas:
                raise EnergyFormatError(
                    f"{source}: frame {frames} lambdas differ from the first frame"
                )

            if len(payload) % _OFFDIAG_SIZE != 0:
                raise EnergyFormatError(
                    f"{source}: frame {frames} off-diagonal record has invalid size {len(payload)}"
                )
            for offset in range(0, len(payload), _OFFDIAG_SIZE):
                index = offset // _OFFDIAG_SIZE + 1
                i, j, hij, rkl = struct.unpack_from(
                    record_format.endian + _OFFDIAG_FORMAT,
                    payload,
                    offset,
                )
                row = {field: "" for field in CSV_FIELDS}
                row.update(
                    {
                        "source_file": source.name,
                        "frame": frames,
                        "record_type": "offdiag",
                        "offdiag_index": index,
                        "offdiag_i": i,
                        "offdiag_j": j,
                        "offdiag_hij": hij,
                        "offdiag_rkl": rkl,
                    }
                )
                rows.append(row)
                offdiag_rows += 1

    if frames == 0 or expected_states is None or expected_lambdas is None:
        raise EnergyFormatError(f"{source}: energy file contains no complete frames")

    summary = EnergyFileSummary(
        source_file=source.name,
        frames=frames,
        states=expected_states,
        lambdas=expected_lambdas,
        state_rows=state_rows,
        offdiag_rows=offdiag_rows,
        byte_order=record_format.byte_order_name,
        record_marker_bytes=record_format.marker_bytes,
    )
    return iter(rows), summary


def convert_energy_files(
    sources: Sequence[str | Path],
    output: str | Path,
) -> ConversionSummary:
    """Atomically consolidate native Q ``.en`` files into one CSV.

    The source files are never modified.  The destination is first written to a
    sibling ``.partial`` file and only published after every source has parsed
    and validated successfully.
    """
    source_paths = tuple(Path(source) for source in sources)
    if not source_paths:
        raise ValueError("at least one energy file is required")
    for source in source_paths:
        if not source.is_file():
            raise FileNotFoundError(source)

    output_path = Path(output)
    resolved_output = output_path.resolve()
    if any(source.resolve() == resolved_output for source in source_paths):
        raise ValueError("output CSV must not replace an input energy file")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial_path = output_path.with_name(output_path.name + ".partial")
    summaries: list[EnergyFileSummary] = []
    open_output = gzip.open if output_path.suffix == ".gz" else Path.open

    try:
        with open_output(partial_path, "wt" if output_path.suffix == ".gz" else "w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for source in source_paths:
                rows, summary = _record_rows(source)
                if summaries and summary.states != summaries[0].states:
                    raise EnergyFormatError(
                        f"{source}: has {summary.states} states; expected {summaries[0].states}"
                    )
                writer.writerows(rows)
                summaries.append(summary)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(partial_path, output_path)
    except Exception:
        # Leave no publishable partial conversion.  Source energy files remain
        # untouched and can still be consumed by qfep or retried.
        if partial_path.exists():
            partial_path.unlink()
        raise

    return ConversionSummary(str(output_path), tuple(summaries))


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="qen2csv",
        description="Consolidate native Q/QGPU binary energy files into one validated CSV.",
    )
    parser.add_argument("energy_files", nargs="+", help="Input .en files in lambda-window order")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output .csv path; use .csv.gz for transparent gzip compression",
    )
    parser.add_argument(
        "--summary-json",
        help="Optional path for machine-readable conversion/validation metadata",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> ConversionSummary:
    args = parse_arguments(argv)
    summary = convert_energy_files(args.energy_files, args.output)
    if args.summary_json:
        summary_path = Path(args.summary_json)
        summary_path.write_text(json.dumps(summary.to_dict(), indent=2) + "\n")
    print(json.dumps(summary.to_dict(), indent=2))
    return summary


def main_exe() -> None:
    main()


if __name__ == "__main__":
    main_exe()
