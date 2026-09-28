"""Utilities for analysing Q energy files."""

import struct
from pathlib import Path

import numpy as np
import pandas as pd


ENERGY_COLUMNS = (
    "frame",
    "state",
    "lambda",
    "total",
    "bond",
    "angle",
    "torsion",
    "improper",
    "qx_coul",
    "qx_vdw",
    "qq_coul",
    "qq_vdw",
    "qp_coul",
    "qp_vdw",
    "qw_coul",
    "qw_vdw",
    "restraint",
)

_STATE_VALUE_COUNT = 15
_STATE_RECORD_SIZE = struct.calcsize("<i15d")


def _energy_dataframe(rows: list[tuple[int | float, ...]]) -> pd.DataFrame:
    dataframe = pd.DataFrame.from_records(rows, columns=ENERGY_COLUMNS)
    column_types = {"frame": "int64", "state": "int64"}
    column_types.update({column: "float64" for column in ENERGY_COLUMNS[2:]})
    return dataframe.astype(column_types)


def read_en(file_path: str | Path) -> pd.DataFrame:
    """Read a binary QGPU ``.en`` file into a tidy DataFrame.

    Each row represents one lambda state in one energy frame. ``frame`` is a
    zero-based index, while ``state`` preserves the one-based state number
    stored in the file.

    Q energy files use Fortran-style records: a 32-bit byte count, the record
    payload, and the same byte count again. QGPU writes one 124-byte state
    record per lambda state and an empty record to terminate each frame.

    Args:
        file_path: Path to the QGPU energy file.

    Returns:
        A DataFrame with the columns listed in :data:`ENERGY_COLUMNS`.

    Raises:
        ValueError: If the file is truncated, corrupt, or does not use the
            QGPU energy-record layout.
    """
    path = Path(file_path)
    rows: list[tuple[int | float, ...]] = []

    with path.open("rb") as energy_file:
        marker = energy_file.read(4)
        if not marker:
            return _energy_dataframe(rows)
        if len(marker) != 4:
            raise ValueError(
                f"Invalid Q energy file {path}: truncated first record marker")

        if struct.unpack("<i", marker)[0] == _STATE_RECORD_SIZE:
            byte_order = "<"
        elif struct.unpack(">i", marker)[0] == _STATE_RECORD_SIZE:
            byte_order = ">"
        else:
            little_size = struct.unpack("<i", marker)[0]
            raise ValueError(
                f"Invalid Q energy file {path} at byte 0: expected a "
                f"{_STATE_RECORD_SIZE}-byte state record, found {little_size}"
            )

        marker_format = f"{byte_order}i"
        state_format = f"{byte_order}i{_STATE_VALUE_COUNT}d"
        frame = 0
        states_in_frame = 0
        expected_states: int | None = None

        while True:
            record_offset = energy_file.tell() - 4
            record_size = struct.unpack(marker_format, marker)[0]
            if record_size not in (0, _STATE_RECORD_SIZE):
                raise ValueError(
                    f"Invalid Q energy file {path} at byte {record_offset}: "
                    f"unsupported record size {record_size}"
                )

            payload = energy_file.read(record_size)
            if len(payload) != record_size:
                raise ValueError(
                    f"Invalid Q energy file {path} at byte {record_offset}: "
                    "truncated record payload"
                )

            trailing_marker = energy_file.read(4)
            if len(trailing_marker) != 4:
                raise ValueError(
                    f"Invalid Q energy file {path} at byte {record_offset}: "
                    "missing trailing record marker"
                )
            trailing_size = struct.unpack(marker_format, trailing_marker)[0]
            if trailing_size != record_size:
                raise ValueError(
                    f"Invalid Q energy file {path} at byte {record_offset}: "
                    f"record markers disagree ({record_size} != {trailing_size})"
                )

            if record_size == 0:
                if states_in_frame == 0:
                    raise ValueError(
                        f"Invalid Q energy file {path} at byte {record_offset}: empty frame"
                    )
                if expected_states is None:
                    expected_states = states_in_frame
                elif states_in_frame != expected_states:
                    raise ValueError(
                        f"Invalid Q energy file {path} at byte {record_offset}: "
                        f"frame {frame} has {states_in_frame} states; expected {expected_states}"
                    )
                frame += 1
                states_in_frame = 0
            else:
                state_record = struct.unpack(state_format, payload)
                state = state_record[0]
                expected_state = states_in_frame + 1
                if state != expected_state:
                    raise ValueError(
                        f"Invalid Q energy file {path} at byte {record_offset}: "
                        f"frame {frame} contains state {state}; expected {expected_state}"
                    )
                rows.append((frame, *state_record))
                states_in_frame += 1

            marker = energy_file.read(4)
            if not marker:
                if states_in_frame:
                    raise ValueError(
                        f"Invalid Q energy file {path}: frame {frame} is missing "
                        "its terminating empty record"
                    )
                break
            if len(marker) != 4:
                raise ValueError(
                    f"Invalid Q energy file {path}: truncated record marker")

    return _energy_dataframe(rows)


def read_en_order(fep_dir: str | Path) -> list[Path]:
    """Return energy files in the order specified by qfep.inp."""
    fep_dir = Path(fep_dir)
    qfep_path = fep_dir / "qfep.inp"
    lines = qfep_path.read_text().splitlines()

    number_of_files = int(lines[0].split()[0])

    try:
        start = next(
            index
            for index, line in enumerate(lines)
            if line.strip().upper() == "!ENERGY_FILES"
        ) + 1
    except StopIteration as error:
        raise ValueError(f"No !ENERGY_FILES section in {qfep_path}") from error

    energy_files = []
    for line in lines[start:]:
        line = line.strip()
        if not line or line.startswith(("!", "#")):
            continue

        energy_files.append(fep_dir / line.split()[0])
        if len(energy_files) == number_of_files:
            break

    if len(energy_files) != number_of_files:
        raise ValueError(
            f"{qfep_path} declares {number_of_files} energy files, "
            f"but only {len(energy_files)} were listed"
        )

    missing = [path for path in energy_files if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing energy file: {missing[0]}")

    return energy_files


def read_fep_directory(fep_dir: str | Path) -> pd.DataFrame:
    tables = []

    for window_order, energy_path in enumerate(read_en_order(fep_dir)):
        table = read_en(energy_path)
        table.insert(0, "window_order", window_order)
        table.insert(1, "file", energy_path.name)
        tables.append(table)

    return pd.concat(tables, ignore_index=True)


def _log_mean_exp(values: np.ndarray) -> float:
    """Compute log(mean(exp(values))) without overflowing."""
    maximum = float(np.max(values))
    return maximum + float(np.log(np.mean(np.exp(values - maximum))))


def _fermi(values: np.ndarray) -> np.ndarray:
    """Compute 1 / (1 + exp(values)) without overflowing."""
    return np.exp(-np.logaddexp(0.0, values))


def calculate_bar_dg(
    energies: pd.DataFrame,
    kT: float,
    skip: int = 0,
    tolerance: float = 0.001,
    max_iterations: int = 10_000,
    alpha: dict[int, float] | None = None,
) -> pd.DataFrame:
    """Calculate adjacent-window BAR free energies using QFEP's equations.

    ``bar_dg`` is the free-energy difference from the preceding window and
    ``bar_dg_sum`` is the cumulative free energy from the first window.
    """
    required_columns = {
        "window_order",
        "file",
        "frame",
        "state",
        "lambda",
        "total",
    }
    missing_columns = required_columns.difference(energies.columns)
    if missing_columns:
        raise ValueError(f"Missing energy columns: {sorted(missing_columns)}")
    if kT <= 0:
        raise ValueError("kT must be positive")
    if skip < 0:
        raise ValueError("skip must be non-negative")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    if max_iterations <= 0:
        raise ValueError("max_iterations must be positive")

    alpha = alpha or {}
    windows = []

    for window_order in sorted(energies["window_order"].unique()):
        complete_window = energies[energies["window_order"] == window_order]
        sampled_window = complete_window[complete_window["frame"] >= skip]
        if sampled_window.empty:
            raise ValueError(
                f"Window {window_order} has no frames after skipping {skip}")
        if sampled_window.duplicated(["frame", "state"]).any():
            raise ValueError(
                f"Window {window_order} contains duplicate frame/state rows")

        state_ids = sorted(complete_window["state"].unique())
        lambda_counts = complete_window.groupby("state")["lambda"].nunique()
        if not lambda_counts.eq(1).all():
            raise ValueError(
                f"Window {window_order} has inconsistent lambda values")

        lambdas = (
            complete_window.groupby("state")["lambda"]
            .first()
            .reindex(state_ids)
            .to_numpy(dtype=float)
        )
        total_by_frame = sampled_window.pivot(
            index="frame", columns="state", values="total")
        total_by_frame = total_by_frame.reindex(columns=state_ids).sort_index()
        if total_by_frame.isna().any().any():
            raise ValueError(
                f"Window {window_order} has incomplete energy frames")

        offsets = np.array([alpha.get(int(state), 0.0) for state in state_ids])
        state_energies = total_by_frame.to_numpy(dtype=float) + offsets
        if not np.isfinite(state_energies).all() or not np.isfinite(lambdas).all():
            raise ValueError(
                f"Window {window_order} contains non-finite values")

        files = complete_window["file"].unique()
        if len(files) != 1:
            raise ValueError(
                f"Window {window_order} refers to multiple energy files")

        windows.append(
            {
                "window_order": int(window_order),
                "file": files[0],
                "states": state_ids,
                "lambdas": lambdas,
                "energies": state_energies,
            }
        )

    if not windows:
        raise ValueError("No energy windows were provided")

    def result_row(window, bar_dg, bar_dg_sum, iterations):
        row = {
            "window_order": window["window_order"],
            "file": window["file"],
        }
        row.update(
            {
                f"lambda_{state}": value
                for state, value in zip(window["states"], window["lambdas"])
            }
        )
        row.update(
            {
                "n_samples": len(window["energies"]),
                "bar_dg": bar_dg,
                "bar_dg_sum": bar_dg_sum,
                "bar_iterations": iterations,
            }
        )
        return row

    results = [result_row(windows[0], 0.0, 0.0, 0)]
    cumulative_dg = 0.0

    for previous, current in zip(windows, windows[1:]):
        if previous["states"] != current["states"]:
            raise ValueError(
                f"Windows {previous['window_order']} and {current['window_order']} "
                "contain different states"
            )

        delta_lambda = current["lambdas"] - previous["lambdas"]
        # pre cur
        
        """
        pre -> cur
        dU = U(cur) - U(pre)
           = delta_lambda * previous
           
        cur -> pre
        dU = U(pre) - U(cur)
           = -delta_lambda * current
        """
        dv_forward = previous["energies"] @ delta_lambda
        dv_reverse = current["energies"] @ delta_lambda

        # QFEP uses overlap sampling as the initial Bennett constant.
        log_sum_forward = _log_mean_exp(-dv_forward / (2.0 * kT))
        log_sum_reverse = _log_mean_exp(dv_reverse / (2.0 * kT))
        constant = -kT * (log_sum_forward - log_sum_reverse)

        n_forward = len(dv_forward)
        n_reverse = len(dv_reverse)
        forward_reverse_ratio = n_forward / n_reverse
        reverse_forward_ratio = n_reverse / n_forward

        for iteration in range(1, max_iterations + 1):
            sum_forward = np.mean(
                _fermi(
                    np.log(forward_reverse_ratio)
                    + (dv_forward - constant) / kT
                )
            )
            sum_reverse = np.mean(
                _fermi(
                    np.log(reverse_forward_ratio)
                    + (-dv_reverse + constant) / kT
                )
            )
            new_constant = -kT * (
                np.log(sum_forward)
                - np.log(sum_reverse)
                - constant / kT
                + np.log(forward_reverse_ratio)
            )
            difference = abs(constant - new_constant)
            constant = float(new_constant)
            if difference <= tolerance:
                break
        else:
            raise RuntimeError(
                f"BAR did not converge for windows {previous['window_order']} and "
                f"{current['window_order']} after {max_iterations} iterations"
            )

        cumulative_dg += constant
        results.append(result_row(current, constant, cumulative_dg, iteration))

    return pd.DataFrame(results)


def main():
    fep_1h1s_31_dir = "/home/mcpi-02/code/qligfepv2-BenchmarkExperiments-shen/perturbations/cdk2/2.protein/FEP_1h1s_31/FEP1/298/1"

    # df = read_en(fep_1h1s_31_dir + "/md_0500_0500.en")
    # print(df)

    # en_oder = read_en_order(fep_1h1s_31_dir)
    # # print(en_oder)
    energies = read_fep_directory(fep_1h1s_31_dir)
    bar = calculate_bar_dg(energies, kT=0.592, skip=100)
    print(bar)
    print("Final BAR dG:", bar.iloc[-1]["bar_dg_sum"])


if __name__ == "__main__":
    main()
