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


def check_fep_pair_consistency(
    data: dict,
    fep: str,
    water_sys: str = "1.water",
    protein_sys: str = "2.protein",
) -> dict:
    """Check whether matching water and protein FEPs use consistent settings."""
    water = data[water_sys][fep]
    protein = data[protein_sys][fep]

    checks = {
        "fep_stage_match": water["fep_stage"] == protein["fep_stage"],
        "temperature_match": water["temperature"] == protein["temperature"],
        "lambda_sum_match": water["lambda_sum"] == protein["lambda_sum"],
    }

    return {
        "fep": fep,
        **checks,
        "consistent": all(checks.values()),
    }


def analyze_ddg_edge(
    data: dict,
    fep: str,
    method: str = "ddGbar",
) -> dict:
    """Extract QC-relevant statistics for a calculated ddG edge."""
    result = data["result"][method][fep]

    return {
        "fep": fep,
        "method": method,
        "from": result["from"],
        "to": result["to"],
        "avg": result[f"{method}_avg"],
        "sem": result[f"{method}_sem"],
        "std": result[f"{method}_std"],
    }


def summarize_fep_edge_qc(
    data: dict,
    fep: str,
    method: str = "dGbar",
    water_sys: str = "1.water",
    protein_sys: str = "2.protein",
) -> dict:
    """Summarize QC information for one complete FEP edge."""
    protein_qc = analyze_fep_edge(
        data,
        system=protein_sys,
        fep=fep,
        method=method,
    )

    water_qc = analyze_fep_edge(
        data,
        system=water_sys,
        fep=fep,
        method=method,
    )

    consistency = check_fep_pair_consistency(
        data,
        fep=fep,
        water_sys=water_sys,
        protein_sys=protein_sys,
    )

    ddg_method = f"d{method}"
    ddg_qc = analyze_ddg_edge(
        data,
        fep=fep,
        method=ddg_method,
    )

    return {
        "fep": fep,
        "from": ddg_qc["from"],
        "to": ddg_qc["to"],
        "ddg": ddg_qc["avg"],
        "ddg_sem": ddg_qc["sem"],
        "ddg_std": ddg_qc["std"],
        "protein_n_valid": protein_qc["n_valid"],
        "protein_n_failed": protein_qc["n_failed"],
        "protein_std": protein_qc["std"],
        "protein_range": protein_qc["range"],
        "water_n_valid": water_qc["n_valid"],
        "water_n_failed": water_qc["n_failed"],
        "water_std": water_qc["std"],
        "water_range": water_qc["range"],
        "consistent": consistency["consistent"],
    }


def summarize_fep_system_qc(
    data: dict,
    method: str = "dGbar",
    water_sys: str = "1.water",
    protein_sys: str = "2.protein",
) -> list[dict]:
    """Summarize QC information for all FEP edges in a calculation."""
    protein_feps = sorted(data[protein_sys])
    water_feps = sorted(data[water_sys])

    if protein_feps != water_feps:
        raise ValueError("FEPs do not match between protein and water systems.")

    return [
        summarize_fep_edge_qc(
            data=data,
            fep=fep,
            method=method,
            water_sys=water_sys,
            protein_sys=protein_sys,
        )
        for fep in protein_feps
    ]


def cycle_closure_error(
    data: dict,
    cycle: list[str],
    method: str = "ddGbar",
) -> float:
    """Calculate the thermodynamic closure error for a ligand cycle.

    Args:
        data: FepReader data dictionary containing calculated ddG results.
        cycle: Ordered ligand names forming a closed cycle, e.g.
            ["lig1", "lig2", "lig3", "lig1"].
        method: Calculated ddG method. Defaults to ``ddGbar``.

    Returns:
        Sum of the directed ddG values around the cycle.

    Raises:
        ValueError: If the cycle is not closed or an edge cannot be found.
    """
    if len(cycle) < 4 or cycle[0] != cycle[-1]:
        raise ValueError("Cycle must contain at least three ligands and be closed.")

    edges = data["result"][method]

    edge_lookup = {}

    for edge in edges.values():
        source = edge["from"]
        target = edge["to"]
        value = edge[f"{method}_avg"]

        edge_lookup[(source, target)] = value
        edge_lookup[(target, source)] = -value

    closure_error = 0.0

    for source, target in zip(cycle[:-1], cycle[1:]):
        try:
            closure_error += edge_lookup[(source, target)]
        except KeyError as exc:
            raise ValueError(f"Missing FEP edge: {source} -> {target}") from exc

    return float(closure_error)
