"""BAR analysis for Q's exact foreign-Hamiltonian work (.fhe) sidecars."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

KB_KCAL = 0.00198720425864083
LAMBDA_ATOL = 1e-12
SUPPORTED_FORMATS = {"Q_FOREIGN_HAMILTONIANS_V1", "Q_FOREIGN_WORKS_V2"}


@dataclass
class ForeignHamiltonianSamples:
    path: Path
    format: str
    source_lambda1: float
    source_lambda2: float
    targets: np.ndarray
    values: np.ndarray

    def _target_column(self, target_lambda2: float, atol: float = LAMBDA_ATOL) -> np.ndarray:
        matches = np.flatnonzero(np.isclose(self.targets, target_lambda2, atol=atol, rtol=0.0))
        if matches.size != 1:
            raise ValueError(
                f"{self.path}: expected one target column for lambda2={target_lambda2}, "
                f"found {matches.size}"
            )
        return self.values[:, matches[0]]

    def energy_at(self, target_lambda2: float, atol: float = LAMBDA_ATOL) -> np.ndarray:
        """Return a V1 total-Hamiltonian column."""
        if self.format != "Q_FOREIGN_HAMILTONIANS_V1":
            raise ValueError(f"{self.path}: {self.format} stores work, not total energy")
        return self._target_column(target_lambda2, atol=atol)

    def work_to(self, target_lambda2: float, atol: float = LAMBDA_ATOL) -> np.ndarray:
        """Return exact H(target)-H(source) for every retained source frame."""
        target = self._target_column(target_lambda2, atol=atol)
        if self.format == "Q_FOREIGN_WORKS_V2":
            return target
        source = self._target_column(self.source_lambda2, atol=atol)
        return target - source



def read_fhe(path: str | Path) -> ForeignHamiltonianSamples:
    """Read and validate a supported Q foreign-work sidecar."""
    path = Path(path)
    targets = None
    rows = []
    format_name = None
    units = None
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("# Q_FOREIGN_"):
                candidate = stripped[2:]
                if candidate in SUPPORTED_FORMATS:
                    format_name = candidate
            elif stripped.startswith("# target_state2_weights"):
                targets = np.asarray([float(value) for value in stripped.split()[2:]], dtype=float)
            elif stripped.startswith("# units "):
                units = stripped.removeprefix("# units ")
            elif not stripped.startswith("#"):
                rows.append([float(value) for value in stripped.split()])

    if format_name is None:
        raise ValueError(f"{path}: unsupported or missing .fhe format version")
    if units != "kcal/mol":
        raise ValueError(f"{path}: expected units kcal/mol, found {units!r}")
    if targets is None or targets.size == 0:
        raise ValueError(f"{path}: no target lambda columns")
    if not np.isfinite(targets).all() or np.any((targets < 0) | (targets > 1)):
        raise ValueError(f"{path}: target lambdas must be finite and in [0, 1]")
    if np.any(np.diff(np.sort(targets)) <= LAMBDA_ATOL):
        raise ValueError(f"{path}: duplicate target lambdas")
    if not rows:
        raise ValueError(f"{path}: no sampled frames")

    data = np.asarray(rows, dtype=float)
    if data.shape[1] != 3 + targets.size:
        raise ValueError(f"{path}: expected {3 + targets.size} columns, found {data.shape[1]}")
    if not np.isfinite(data).all():
        raise ValueError(f"{path}: non-finite foreign work value")
    if np.any((data[:, 1:3] < 0) | (data[:, 1:3] > 1)) or not np.allclose(
        data[:, 1:3].sum(axis=1), 1.0, atol=LAMBDA_ATOL, rtol=0.0
    ):
        raise ValueError(f"{path}: source lambdas must be in [0, 1] and sum to one")
    if not np.allclose(data[:, 1:3], data[0, 1:3], atol=LAMBDA_ATOL, rtol=0.0):
        raise ValueError(f"{path}: source lambda changed within one equilibrium sidecar")
    if np.any(data[:, 0] < 1) or np.any(data[:, 0] != np.floor(data[:, 0])) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f"{path}: frame numbers must be positive, strictly increasing integers")
    source_column = np.flatnonzero(np.isclose(targets, data[0, 2], atol=LAMBDA_ATOL, rtol=0.0))
    if source_column.size != 1:
        raise ValueError(f"{path}: missing source lambda target column")
    if format_name == "Q_FOREIGN_WORKS_V2" and not np.allclose(
        data[:, 3 + source_column[0]], 0.0, atol=1e-7, rtol=0.0
    ):
        raise ValueError(f"{path}: source-to-source work must be zero")

    return ForeignHamiltonianSamples(
        path=path,
        format=format_name,
        source_lambda1=float(data[0, 1]),
        source_lambda2=float(data[0, 2]),
        targets=targets,
        values=data[:, 3:],
    )



def _fermi_weights(beta_work: np.ndarray) -> np.ndarray:
    """Fermi weights, evaluated without exponential overflow."""
    result = np.empty_like(beta_work)
    positive = beta_work >= 0
    result[positive] = np.exp(-beta_work[positive]) / (1.0 + np.exp(-beta_work[positive]))
    result[~positive] = 1.0 / (1.0 + np.exp(beta_work[~positive]))
    return result



def _weight_ess(weights: np.ndarray) -> float:
    denominator = float(np.dot(weights, weights))
    return 0.0 if denominator == 0.0 else float(weights.sum() ** 2 / denominator)


def _bar_delta_f(forward: np.ndarray, reverse: np.ndarray, beta: float) -> float:
    """Count-corrected BAR in log space, including very large endpoint works.

    Uses the BAR equation documented by pymbar.other_estimators.bar_zero:
    https://pymbar.readthedocs.io/en/stable/other_estimators.html
    Keep this independent of the legacy NEQ estimator and its count convention.
    """
    wf, wr = beta * forward, beta * reverse
    if not np.isfinite(wf).all() or not np.isfinite(wr).all():
        raise ValueError("reduced foreign work is non-finite")
    count_offset = np.log(forward.size / reverse.size)

    def objective(df):
        return float(
            logsumexp(-np.logaddexp(0.0, wf - df + count_offset))
            - logsumexp(-np.logaddexp(0.0, wr + df - count_offset))
        )

    lo = min(float(wf.min()), float(-wr.max()))
    hi = max(float(wf.max()), float(-wr.min()))
    if lo == hi:
        return lo / beta
    return float(brentq(objective, lo, hi)) / beta



def analyze_fhe(
    paths: list[str | Path],
    temperature: float = 298.0,
    skip: int = 0,
    stride: int = 1,
) -> dict:
    """Analyze adjacent windows and return BAR/overlap diagnostics.

    ``skip`` discards initial records from every window as burn-in. ``stride``
    then selects every Nth retained record; it controls analysis density but
    does not make temporally correlated configurations independent.
    """
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be positive")
    if skip < 0:
        raise ValueError("skip must be non-negative")
    if stride < 1:
        raise ValueError("stride must be positive")
    beta = 1.0 / (KB_KCAL * temperature)

    grouped: dict[float, list[ForeignHamiltonianSamples]] = {}
    formats = set()
    seen = set()
    for path in paths:
        resolved = Path(path).resolve()
        if resolved in seen:
            raise ValueError(f"duplicate sidecar: {path}")
        seen.add(resolved)
        sample = read_fhe(path)
        sample.values = sample.values[skip::stride]
        if sample.values.shape[0] == 0:
            raise ValueError(f"{sample.path}: no frames remain after skip={skip}, stride={stride}")
        grouped.setdefault(sample.source_lambda2, []).append(sample)
        formats.add(sample.format)
    if len(grouped) < 2:
        raise ValueError("BAR needs at least two source lambda windows")

    windows = sorted(grouped)
    for samples in grouped.values():
        for sample in samples:
            for target in sample.targets:
                if not np.any(np.isclose(windows, target, atol=LAMBDA_ATOL, rtol=0.0)):
                    raise ValueError(f"{sample.path}: missing source window for target lambda2={target}")
    results = []
    cumulative = 0.0
    for lambda_i, lambda_j in zip(windows, windows[1:]):
        samples_i = grouped[lambda_i]
        samples_j = grouped[lambda_j]
        work_forward = np.concatenate([sample.work_to(lambda_j) for sample in samples_i])
        work_reverse = np.concatenate([sample.work_to(lambda_i) for sample in samples_j])
        delta_g = _bar_delta_f(work_forward, work_reverse, beta)
        cumulative += delta_g
        gap_kt = float(beta * (work_forward.mean() + work_reverse.mean()))
        count_offset = np.log(work_forward.size / work_reverse.size)
        weights_forward = _fermi_weights(beta * (work_forward - delta_g) + count_offset)
        weights_reverse = _fermi_weights(beta * (work_reverse + delta_g) - count_offset)
        acceptance_forward = float(weights_forward.mean())
        acceptance_reverse = float(weights_reverse.mean())
        ess_forward = _weight_ess(weights_forward)
        ess_reverse = _weight_ess(weights_reverse)
        bar_residual = float((weights_forward.sum() - weights_reverse.sum()) / (work_forward.size + work_reverse.size))
        poor = abs(gap_kt) > 2.0 or min(acceptance_forward, acceptance_reverse) < 0.05
        results.append(
            {
                "lambda2_i": lambda_i,
                "lambda2_j": lambda_j,
                "frames": int(work_forward.size),
                "frames_forward": int(work_forward.size),
                "frames_reverse": int(work_reverse.size),
                "delta_g_kcal": delta_g,
                "cumulative_delta_g_kcal": cumulative,
                "mean_forward_work_kcal": float(work_forward.mean()),
                "mean_reverse_work_kcal": float(work_reverse.mean()),
                "overlap_gap_kT": gap_kt,
                "fermi_acceptance_forward": acceptance_forward,
                "fermi_acceptance_reverse": acceptance_reverse,
                "fermi_ess_forward": ess_forward,
                "fermi_ess_reverse": ess_reverse,
                "bar_residual": bar_residual,
                "poor_overlap": poor,
            }
        )

    return {
        "format": next(iter(formats)) if len(formats) == 1 else "mixed",
        "temperature_K": temperature,
        "beta_per_kcal": beta,
        "frame_selection": {"skip": skip, "stride": stride},
        "windows": windows,
        "total_delta_g_kcal": cumulative,
        "poor_windows": sum(row["poor_overlap"] for row in results),
        "adjacent": results,
    }



def write_results(result: dict, output_prefix: str | Path) -> None:
    prefix = Path(output_prefix)
    prefix.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    rows = result["adjacent"]
    with prefix.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)



def main_exe() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help=".fhe files or directories containing .fhe files")
    parser.add_argument("-T", "--temperature", type=float, default=298.0)
    parser.add_argument("--skip", type=int, default=0, help="initial records to discard from every window")
    parser.add_argument("--stride", type=int, default=1, help="retain every Nth record after burn-in")
    parser.add_argument("-o", "--output-prefix", default="foreign_bar")
    args = parser.parse_args()

    paths = []
    for value in args.paths:
        path = Path(value)
        if path.is_dir():
            paths.extend(sorted(path.glob("*.fhe")))
        else:
            paths.append(path)
    result = analyze_fhe(paths, temperature=args.temperature, skip=args.skip, stride=args.stride)
    write_results(result, args.output_prefix)
    print(f"Delta G = {result['total_delta_g_kcal']:.6f} kcal/mol")
    print(f"Poor-overlap windows: {result['poor_windows']} / {len(result['adjacent'])}")


if __name__ == "__main__":
    main_exe()
