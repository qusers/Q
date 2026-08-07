"""Tests for native Q/QGPU energy-file consolidation."""

import csv
import gzip
import struct
from pathlib import Path

import pytest

from QligFEP.energy_converter import EnergyFormatError, convert_energy_files

STATE_FORMAT = "i15d"
OFFDIAG_FORMAT = "ii2d"


def _write_record(stream, payload: bytes, endian: str = "<", marker_bytes: int = 4) -> None:
    marker_format = endian + ("i" if marker_bytes == 4 else "q")
    marker = struct.pack(marker_format, len(payload))
    stream.write(marker)
    stream.write(payload)
    stream.write(marker)


def _write_energy_file(
    path: Path,
    frames: list[list[tuple]],
    offdiagonals: list[list[tuple]] | None = None,
    endian: str = "<",
    marker_bytes: int = 4,
) -> None:
    offdiagonals = offdiagonals or [[] for _ in frames]
    with path.open("wb") as stream:
        for states, frame_offdiagonals in zip(frames, offdiagonals, strict=True):
            for values in states:
                _write_record(
                    stream,
                    struct.pack(endian + STATE_FORMAT, *values),
                    endian,
                    marker_bytes,
                )
            payload = b"".join(
                struct.pack(endian + OFFDIAG_FORMAT, *values) for values in frame_offdiagonals
            )
            _write_record(stream, payload, endian, marker_bytes)


def _state(state: int, lambda_value: float, base: float) -> tuple:
    # state index plus the 15 float64 values in nrgy.q_energies order
    return (state, lambda_value, *(base + i for i in range(14)))


def test_consolidates_windows_and_preserves_source_filename(tmp_path):
    first = tmp_path / "md_1000_0000.en"
    second = tmp_path / "md_0000_1000.en"
    _write_energy_file(
        first,
        frames=[
            [_state(1, 1.0, 10.0), _state(2, 0.0, 20.0)],
            [_state(1, 1.0, 30.0), _state(2, 0.0, 40.0)],
        ],
    )
    _write_energy_file(
        second,
        frames=[[_state(1, 0.0, 50.0), _state(2, 1.0, 60.0)]],
    )

    output = tmp_path / "energies.csv"
    summary = convert_energy_files([first, second], output)

    with output.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert summary.frames == 3
    assert summary.rows == 6
    assert [item.frames for item in summary.files] == [2, 1]
    assert summary.files[0].lambdas == (1.0, 0.0)
    assert summary.files[1].lambdas == (0.0, 1.0)
    assert {row["source_file"] for row in rows} == {first.name, second.name}
    assert [row["frame"] for row in rows if row["source_file"] == first.name] == ["1", "1", "2", "2"]
    assert rows[0]["record_type"] == "state"
    assert rows[0]["state"] == "1"
    assert float(rows[0]["lambda"]) == 1.0
    assert float(rows[0]["total"]) == 10.0
    assert float(rows[0]["restraint"]) == 23.0


def test_writes_transparently_readable_gzip_csv(tmp_path):
    source = tmp_path / "md_0500_0500.en"
    _write_energy_file(source, frames=[[_state(1, 0.5, 1.0), _state(2, 0.5, 2.0)]])
    output = tmp_path / "energies.csv.gz"

    summary = convert_energy_files([source], output)

    with gzip.open(output, "rt", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert summary.rows == 2
    assert [row["state"] for row in rows] == ["1", "2"]


def test_preserves_offdiagonal_records_in_same_csv(tmp_path):
    source = tmp_path / "evb.en"
    _write_energy_file(
        source,
        frames=[[_state(1, 0.5, 1.0), _state(2, 0.5, 2.0)]],
        offdiagonals=[[(1, 2, 12.5, 3.25), (2, 3, -4.0, 8.5)]],
    )

    output = tmp_path / "energies.csv"
    summary = convert_energy_files([source], output)
    with output.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert summary.files[0].offdiag_rows == 2
    offdiag = [row for row in rows if row["record_type"] == "offdiag"]
    assert len(offdiag) == 2
    assert offdiag[0]["source_file"] == source.name
    assert offdiag[0]["frame"] == "1"
    assert offdiag[0]["state"] == ""
    assert offdiag[0]["offdiag_i"] == "1"
    assert offdiag[0]["offdiag_j"] == "2"
    assert float(offdiag[0]["offdiag_hij"]) == 12.5
    assert float(offdiag[0]["offdiag_rkl"]) == 3.25


@pytest.mark.parametrize(("endian", "marker_bytes"), [("<", 4), (">", 4), ("<", 8)])
def test_detects_supported_fortran_record_variants(tmp_path, endian, marker_bytes):
    source = tmp_path / f"variant-{marker_bytes}.en"
    _write_energy_file(
        source,
        frames=[[_state(1, 1.0, 1.0), _state(2, 0.0, 2.0)]],
        endian=endian,
        marker_bytes=marker_bytes,
    )

    summary = convert_energy_files([source], tmp_path / "energies.csv")

    assert summary.files[0].byte_order == ("little" if endian == "<" else "big")
    assert summary.files[0].record_marker_bytes == marker_bytes


def test_rejects_inconsistent_state_count_without_publishing_csv(tmp_path):
    source = tmp_path / "broken.en"
    _write_energy_file(
        source,
        frames=[
            [_state(1, 1.0, 1.0), _state(2, 0.0, 2.0)],
            [_state(1, 1.0, 3.0)],
        ],
    )
    output = tmp_path / "energies.csv"

    with pytest.raises(EnergyFormatError, match="has 1 states; expected 2"):
        convert_energy_files([source], output)

    assert source.exists()
    assert not output.exists()
    assert not (tmp_path / "energies.csv.partial").exists()


def test_rejects_truncated_file_and_preserves_sources(tmp_path):
    source = tmp_path / "truncated.en"
    _write_energy_file(source, frames=[[_state(1, 1.0, 1.0)]])
    source.write_bytes(source.read_bytes()[:-2])

    with pytest.raises(EnergyFormatError, match="truncated"):
        convert_energy_files([source], tmp_path / "energies.csv")

    assert source.exists()
