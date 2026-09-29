"""Quality-control utilities for QligFEP analysis."""

from collections.abc import Sequence

import numpy as np


def replicate_statistics(energies: Sequence[float]) -> dict:
    """Calculate descriptive statistics for replicate FEP energies.

    NaN values are treated as failed or unavailable replicate results and are
    excluded from the descriptive statistics.

    Args:
        energies: Energy values from individual FEP replicates.

    Returns:
        Dictionary containing replicate counts and descriptive statistics.
    """
    values = np.asarray(energies, dtype=float)
    valid_values = values[~np.isnan(values)]

    n_replicates = len(values)
    n_valid = len(valid_values)
    n_failed = n_replicates - n_valid

    if n_valid == 0:
        return {
            "n_replicates": n_replicates,
            "n_valid": 0,
            "n_failed": n_failed,
            "mean": None,
            "std": None,
            "sem": None,
            "range": None,
        }

    std = float(np.std(valid_values))

    return {
        "n_replicates": n_replicates,
        "n_valid": n_valid,
        "n_failed": n_failed,
        "mean": float(np.mean(valid_values)),
        "std": std,
        "sem": float(std / np.sqrt(n_valid)),
        "range": float(np.ptp(valid_values)),
    }


def analyze_fep_edge(
    data: dict,
    system: str,
    fep: str,
    method: str = "dGbar",
) -> dict:
    """Calculate replicate QC statistics for a single FEP edge.

    Args:
        data: FepReader data dictionary.
        system: System containing the FEP edge, e.g. ``2.protein``.
        fep: Name of the FEP edge, e.g. ``FEP_lig1_lig2``.
        method: Energy estimator to analyze. Defaults to ``dGbar``.

    Returns:
        Dictionary containing edge metadata and replicate statistics.
    """
    energies = data[system][fep]["FEP_result"][method]["energies"]
    statistics = replicate_statistics(energies)

    return {
        "system": system,
        "fep": fep,
        "method": method,
        **statistics,
    }


def analyze_fep_system(
    data: dict,
    system: str,
    method: str = "dGbar",
) -> list[dict]:
    """Calculate replicate QC statistics for all FEP edges in a system.

    Args:
        data: FepReader data dictionary.
        system: System containing the FEP edges, e.g. ``2.protein``.
        method: Energy estimator to analyze. Defaults to ``dGbar``.

    Returns:
        List of QC dictionaries, one for each FEP edge.
    """
    results = []

    for fep in sorted(data[system]):
        if "FEP_result" not in data[system][fep]:
            continue

        results.append(
            analyze_fep_edge(
                data=data,
                system=system,
                fep=fep,
                method=method,
            )
        )

    return results
