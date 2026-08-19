"""Collect and combine the results of amino-acid mutation FEPs (QresFEP).

The ligand workflow needs a mapping file to know which perturbations exist,
because the edges of a ligand network are an arbitrary choice made upstream. A
mutation set is not arbitrary: every perturbation is fully described by its own
name, so the results are discovered by walking the directory tree that
``qresfep`` created.

A mutation is only meaningful once both legs have run::

    ddG_fold = dG_protein - dG_tripeptide

so the two legs are matched by mutation name and combined here. The layout this
walks is the one the run script produces::

    <leg>/FEP_<WT><POS><MUT>/FEP<stage>/<temperature>/<replicate>/qfep.out
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .logger import logger, setup_logger

#: Free energy estimators read from each qfep.out, in reporting order. BAR is the
#: primary result; the one-directional Zwanzig sums are kept because the gap
#: between them is the most direct signal that a stage is underconverged.
ESTIMATORS = ("dG", "dG_forward", "dG_reverse", "dG_bar")

#: Directory name pattern written by `qresfep`.
_FEP_DIR = re.compile(r"^FEP_([A-Z]{3})(\d+)([A-Z]{3})$")

#: Section labels in qfep output, keyed by the fourth token of each heading.
_QFEP_BLOCKS = {
    "Free": "zwanzig",
    "Termodynamic": "ti",
    "Overlap": "overlap",
    "BAR": "bar",
    "Reaction": None,
}


@dataclass
class LegResult:
    """One leg (protein or reference) of one mutation.

    Attributes:
        leg: ``protein`` or ``tripeptide``.
        replicates: Per-replicate total free energies, keyed by estimator; each
            total is the sum over the two FEP stages.
        crashed: Replicates that did not produce a readable result, with a reason.
    """

    leg: str
    replicates: dict[str, list[float]] = field(default_factory=dict)
    crashed: list[str] = field(default_factory=list)

    def mean(self, estimator: str = "dG_bar") -> float:
        values = [v for v in self.replicates.get(estimator, []) if not math.isnan(v)]
        return float(np.mean(values)) if values else math.nan

    def sem(self, estimator: str = "dG_bar") -> float:
        """Standard error over replicates; NaN with fewer than two."""
        values = [v for v in self.replicates.get(estimator, []) if not math.isnan(v)]
        if len(values) < 2:
            return math.nan
        return float(np.std(values, ddof=1) / math.sqrt(len(values)))

    @property
    def n_replicates(self) -> int:
        values = self.replicates.get("dG_bar", [])
        return sum(1 for v in values if not math.isnan(v))


@dataclass
class MutationResult:
    """Both legs of one mutation, and the folding free energy shift they imply."""

    mutation: str
    wild_type: str
    position: int
    mutant: str
    legs: dict[str, LegResult] = field(default_factory=dict)

    def ddG(self, estimator: str = "dG_bar") -> float:
        """Return ddG_fold = dG_protein - dG_tripeptide, or NaN if a leg is missing."""
        if "protein" not in self.legs or "tripeptide" not in self.legs:
            return math.nan
        return self.legs["protein"].mean(estimator) - self.legs["tripeptide"].mean(estimator)

    def ddG_sem(self, estimator: str = "dG_bar") -> float:
        """Propagate the two legs' standard errors, which are independent."""
        if "protein" not in self.legs or "tripeptide" not in self.legs:
            return math.nan
        protein = self.legs["protein"].sem(estimator)
        reference = self.legs["tripeptide"].sem(estimator)
        if math.isnan(protein) or math.isnan(reference):
            return math.nan
        return math.sqrt(protein**2 + reference**2)


def read_qfep_stage(qfep_out: Path, reverse_direction: bool) -> dict[str, float]:
    """Read one stage's free energies from a ``qfep.out``.

    Args:
        qfep_out: The qfep output file.
        reverse_direction: True when the stage was run from lambda 0 to 1 rather
            than 1 to 0, in which case the reported values carry the opposite
            sign. Mutations out of glycine are run this way, because the
            transformation has to start from the larger side chain.

    Returns:
        One value per entry of :data:`ESTIMATORS`; NaN where qfep reported
        ``-Infinity`` or the block was absent.

    Raises:
        OSError: If qfep recorded an error instead of results.
    """
    start, end = ("0.000000", "1.000000") if reverse_direction else ("1.000000", "0.000000")
    sign = -1.0 if reverse_direction else 1.0
    results = dict.fromkeys(ESTIMATORS, math.nan)

    def value(token: str) -> float:
        return math.nan if token == "-Infinity" else sign * float(token)

    block = None
    for raw in qfep_out.read_text(errors="ignore").splitlines():
        if "ERROR:" in raw:
            raise OSError(f"qfep reported an error in {qfep_out}: {raw.strip()}")
        tokens = raw.split()
        if len(tokens) > 3:
            block = _QFEP_BLOCKS.get(tokens[3], block)
        if len(tokens) < 2:
            continue

        if block == "zwanzig":
            if tokens[0] == start and len(tokens) > 4:
                results["dG_reverse"] = value(tokens[4])
            elif tokens[0] == end and len(tokens) > 5:
                results["dG_forward"] = value(tokens[2])
                results["dG"] = value(tokens[5])
        elif block == "bar" and tokens[0] == end:
            results["dG_bar"] = value(tokens[2])

    return results


def _stage_direction(wild_type: str, start: str, stage: int) -> bool:
    """Return whether a stage ran from lambda 0 to 1.

    Mutations out of glycine are set up in the reverse direction throughout. A
    run started at mid-lambda splits: its first stage goes one way, its second
    the other.
    """
    if wild_type == "GLY":
        return True
    if start == "0.5":
        return stage == 1
    return False


def read_mutation_leg(fep_dir: Path, leg: str, temperature: str) -> LegResult:
    """Sum both FEP stages for every replicate of one leg.

    A replicate only counts when both stages produced a result: the stages are
    consecutive halves of one transformation, so half of one is not a free energy.

    Args:
        fep_dir: The ``FEP_<mutation>`` directory.
        leg: Name to record for this leg.
        temperature: Temperature subdirectory to read.
    """
    result = LegResult(leg=leg)
    match = _FEP_DIR.match(fep_dir.name)
    wild_type = match.group(1) if match else ""

    start = "1"
    config = fep_dir / "inputfiles" / "resfep_config.json"
    if config.exists():
        start = str(json.loads(config.read_text()).get("start", "1"))

    stage_dirs = [fep_dir / "FEP1", fep_dir / "FEP2"]
    existing_stages = [stage for stage in stage_dirs if stage.is_dir()]
    if not existing_stages:
        logger.debug(f"{fep_dir.name} ({leg}): no FEP stage directories yet")
        return result

    replicates = sorted(
        {
            path.name
            for stage in existing_stages
            for path in (stage / temperature).glob("*")
            if path.is_dir()
        },
        key=lambda name: int(name) if name.isdigit() else 0,
    )

    for replicate in replicates:
        totals = dict.fromkeys(ESTIMATORS, 0.0)
        failed = None
        for stage_dir in stage_dirs:
            stage = int(stage_dir.name[-1])
            qfep_out = stage_dir / temperature / replicate / "qfep.out"
            if not qfep_out.exists():
                failed = f"replicate {replicate}: {stage_dir.name} has no qfep.out"
                break
            try:
                values = read_qfep_stage(
                    qfep_out, _stage_direction(wild_type, start, stage)
                )
            except OSError as exc:
                failed = f"replicate {replicate}: {exc}"
                break
            for estimator in ESTIMATORS:
                totals[estimator] += values[estimator]

        if failed:
            result.crashed.append(failed)
            logger.debug(failed)
            continue
        for estimator in ESTIMATORS:
            result.replicates.setdefault(estimator, []).append(totals[estimator])

    return result


def collect(
    protein_dir: Path,
    tripeptide_dir: Path,
    temperature: str = "298",
) -> list[MutationResult]:
    """Find every mutation that has been set up and read whatever has finished.

    Mutations present in only one leg are still reported, with a NaN ddG, so that
    an incomplete set is visible rather than silently dropped.

    Args:
        protein_dir: Directory holding the protein leg's ``FEP_*`` directories.
        tripeptide_dir: Same for the reference peptide leg.
        temperature: Temperature subdirectory to read.
    """
    found: dict[str, MutationResult] = {}
    for leg, directory in (("protein", Path(protein_dir)), ("tripeptide", Path(tripeptide_dir))):
        if not directory.is_dir():
            logger.warning(f"{directory} does not exist; skipping the {leg} leg")
            continue
        for fep_dir in sorted(directory.glob("FEP_*")):
            match = _FEP_DIR.match(fep_dir.name)
            if not match:
                logger.debug(f"Skipping {fep_dir.name}: not a mutation directory")
                continue
            wild_type, position, mutant = match.group(1), int(match.group(2)), match.group(3)
            name = f"{wild_type}{position}{mutant}"
            record = found.setdefault(
                name,
                MutationResult(
                    mutation=name, wild_type=wild_type, position=position, mutant=mutant
                ),
            )
            record.legs[leg] = read_mutation_leg(fep_dir, leg, temperature)

    return [found[name] for name in sorted(found, key=lambda n: found[n].position)]


def to_frame(results: list[MutationResult], estimator: str = "dG_bar") -> pd.DataFrame:
    """Lay the results out as one row per mutation.

    Args:
        results: Output of :func:`collect`.
        estimator: Which estimator the reported values come from.

    Returns:
        A DataFrame with the per-leg free energies, the folding shift and the
        replicate counts.
    """
    rows = []
    for record in results:
        protein = record.legs.get("protein")
        reference = record.legs.get("tripeptide")
        rows.append(
            {
                "mutation": record.mutation,
                "wild_type": record.wild_type,
                "position": record.position,
                "mutant": record.mutant,
                "dG_protein": protein.mean(estimator) if protein else math.nan,
                "dG_protein_sem": protein.sem(estimator) if protein else math.nan,
                "n_protein": protein.n_replicates if protein else 0,
                "dG_tripeptide": reference.mean(estimator) if reference else math.nan,
                "dG_tripeptide_sem": reference.sem(estimator) if reference else math.nan,
                "n_tripeptide": reference.n_replicates if reference else 0,
                "ddG_fold": record.ddG(estimator),
                "ddG_fold_sem": record.ddG_sem(estimator),
                "crashes": (len(protein.crashed) if protein else 0)
                + (len(reference.crashed) if reference else 0),
            }
        )
    return pd.DataFrame(rows)


def compare_to_experiment(frame: pd.DataFrame, experimental_csv: Path) -> pd.DataFrame:
    """Join experimental values onto the results and report the deviations.

    Args:
        frame: Output of :func:`to_frame`.
        experimental_csv: CSV with a ``mutation`` column and a ``ddG_exp`` column.

    Returns:
        The frame with ``ddG_exp`` and ``error`` columns added, restricted to
        mutations that carry an experimental value.

    Raises:
        ValueError: If the CSV lacks either required column.
    """
    experimental = pd.read_csv(experimental_csv)
    missing = {"mutation", "ddG_exp"} - set(experimental.columns)
    if missing:
        raise ValueError(f"{experimental_csv} is missing column(s): {', '.join(sorted(missing))}")

    merged = frame.merge(experimental[["mutation", "ddG_exp"]], on="mutation", how="left")
    merged["error"] = merged["ddG_fold"] - merged["ddG_exp"]
    return merged


def summarize(frame: pd.DataFrame) -> str:
    """Return a short report of agreement with experiment, for logging."""
    scored = frame.dropna(subset=["ddG_fold", "ddG_exp"]) if "ddG_exp" in frame else pd.DataFrame()
    if scored.empty:
        return "No mutations have both a computed and an experimental value."

    error = scored["error"]
    rmse = float(np.sqrt((error**2).mean()))
    mue = float(error.abs().mean())
    correlation = float(scored["ddG_fold"].corr(scored["ddG_exp"]))
    return (
        f"{len(scored)} mutation(s) with experimental data: "
        f"RMSE {rmse:.2f}, MUE {mue:.2f}, R {correlation:.2f} kcal/mol"
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="qresfep_analyze",
        description=(
            "Collect QresFEP results and report the folding free energy shift of each "
            "mutation. Mutations are discovered from the FEP_<WT><POS><MUT> directory "
            "names, so no mapping file is needed."
        ),
    )
    parser.add_argument(
        "-p",
        "--protein",
        dest="protein_dir",
        default="protein",
        help="Directory holding the protein leg's FEP_* directories. Defaults to `protein`.",
    )
    parser.add_argument(
        "-t",
        "--tripeptide",
        dest="tripeptide_dir",
        default="tripeptide",
        help="Directory holding the reference leg's FEP_* directories. Defaults to `tripeptide`.",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        default="298",
        help="Temperature subdirectory to read. Defaults to 298.",
    )
    parser.add_argument(
        "-e",
        "--estimator",
        dest="estimator",
        default="dG_bar",
        choices=list(ESTIMATORS),
        help="Free energy estimator to report. Defaults to dG_bar (Bennett acceptance ratio).",
    )
    parser.add_argument(
        "-exp",
        "--experimental",
        dest="experimental",
        default=None,
        help="Optional CSV with `mutation` and `ddG_exp` columns, to score the results against.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="resfep_results.csv",
        help="Where to write the results table. Defaults to resfep_results.csv.",
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log",
        default="info",
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Log level. Defaults to info.",
    )
    return parser.parse_args()


def main(args: argparse.Namespace | None = None) -> pd.DataFrame:
    """Collect the results, write the table and report what was found."""
    args = args or parse_arguments()
    setup_logger(level=args.log.upper())

    results = collect(Path(args.protein_dir), Path(args.tripeptide_dir), args.temperature)
    if not results:
        logger.warning(
            f"No FEP_* directories found under {args.protein_dir} or {args.tripeptide_dir}"
        )
        return pd.DataFrame()

    frame = to_frame(results, args.estimator)
    if args.experimental:
        frame = compare_to_experiment(frame, Path(args.experimental))
        logger.info(summarize(frame))

    frame.to_csv(args.output, index=False)
    logger.info(f"{len(frame)} mutation(s) written to {args.output}")

    incomplete = frame[frame["ddG_fold"].isna()]
    if not incomplete.empty:
        logger.warning(
            f"{len(incomplete)} mutation(s) have no ddG yet (a leg is missing or unfinished): "
            + ", ".join(incomplete["mutation"])
        )
    crashed = frame[frame["crashes"] > 0]
    if not crashed.empty:
        logger.warning(
            f"{int(crashed['crashes'].sum())} replicate(s) failed across "
            f"{len(crashed)} mutation(s): " + ", ".join(crashed["mutation"])
        )
    return frame


def main_exe():
    main()


if __name__ == "__main__":
    main_exe()
