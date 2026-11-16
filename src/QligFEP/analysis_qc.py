"""Quality-control utilities for QligFEP analysis."""

from collections.abc import Sequence

import numpy as np


def replicate_statistics(energies: Sequence[float]) -> dict:
    """Calculate descriptive statistics for replicate FEP energies.

    Non-finite values are treated as failed or unavailable replicate results and are
    excluded from the descriptive statistics.

    Args:
        energies: Energy values from individual FEP replicates.

    Returns:
        Dictionary containing replicate counts and descriptive statistics.
    """
    values = np.asarray(energies, dtype=float)
    valid_values = values[np.isfinite(values)]

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
    """Compare stage, temperature, and independently observed window counts.

    Missing input files leave the window-count check unverified. The legacy
    ``lambda_sum`` field may reflect a user override or a cached count in older
    result files, so it is not evidence of matching inputs.
    """
    water = data[water_sys][fep]
    protein = data[protein_sys][fep]

    water_count = water.get("input_n_lambdas")
    protein_count = protein.get("input_n_lambdas")
    checks = {
        "fep_stage_match": water["fep_stage"] == protein["fep_stage"],
        "temperature_match": water["temperature"] == protein["temperature"],
        "lambda_sum_match": (
            water_count == protein_count if water_count is not None and protein_count is not None else None
        ),
    }

    if False in checks.values():
        consistent = False
        status = "mismatch"
    elif None in checks.values():
        consistent = None
        status = "unverified"
    else:
        consistent = True
        status = "consistent"

    return {
        "fep": fep,
        **checks,
        "protein_input_n_lambdas": protein_count,
        "water_input_n_lambdas": water_count,
        "consistent": consistent,
        "consistency_status": status,
    }


def analyze_ddg_edge(
    data: dict,
    fep: str,
    method: str = "ddGbar",
) -> dict:
    """Extract QC-relevant statistics for a calculated ddG edge."""
    result = data["result"][method][fep]

    statistics = {}
    for name in ("avg", "sem", "std"):
        value = result[f"{method}_{name}"]
        statistics[name] = value if value is not None and np.isfinite(value) else None

    return {
        "fep": fep,
        "method": method,
        "from": result["from"],
        "to": result["to"],
        **statistics,
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
        **{key: value for key, value in consistency.items() if key != "fep"},
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
        if value is None or not np.isfinite(value):
            continue
        edge_lookup[(source, target)] = value
        edge_lookup[(target, source)] = -value

    closure_error = 0.0

    for source, target in zip(cycle[:-1], cycle[1:]):
        try:
            closure_error += edge_lookup[(source, target)]
        except KeyError as exc:
            raise ValueError(f"Missing FEP edge: {source} -> {target}") from exc

    return float(closure_error)


def summarize_cycle_closure_qc(
    data: dict,
    method: str = "ddGbar",
) -> list[dict]:
    """Report basis cycles, marking unavailable closures without aborting QC."""
    cycles = find_cycle_basis(data, method=method)
    results = []

    for cycle in cycles:
        try:
            closure_error = cycle_closure_error(data, cycle=cycle, method=method)
        except ValueError as exc:
            closure_error = None
            status = "unavailable"
            reason = str(exc)
        else:
            status = "ok"
            reason = ""

        results.append(
            {
                "cycle": cycle,
                "n_edges": len(cycle) - 1,
                "closure_error": closure_error,
                "abs_closure_error": abs(closure_error) if closure_error is not None else None,
                "status": status,
                "reason": reason,
            }
        )

    return results


def _canonicalize_cycle(cycle: list[str]) -> list[str]:
    """Return a deterministic representation of a closed cycle."""
    nodes = cycle[:-1]

    rotations = []

    for sequence in (nodes, list(reversed(nodes))):
        for index in range(len(sequence)):
            rotated = sequence[index:] + sequence[:index]
            rotations.append(tuple(rotated))

    canonical = min(rotations)

    return list(canonical) + [canonical[0]]


def find_cycle_basis(
    data: dict,
    method: str = "ddGbar",
) -> list[list[str]]:
    """Find an independent cycle basis for the FEP network."""
    edges = data["result"][method]

    adjacency: dict[str, set[str]] = {}

    for edge in edges.values():
        source = edge["from"]
        target = edge["to"]

        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set()).add(source)

    visited = set()
    parent = {}
    depth = {}
    tree_edges = set()
    back_edges = []

    def edge_key(a: str, b: str) -> tuple[str, str]:
        return tuple(sorted((a, b)))

    def dfs(node: str, node_parent: str | None, node_depth: int) -> None:
        visited.add(node)
        parent[node] = node_parent
        depth[node] = node_depth

        for neighbor in sorted(adjacency[node]):
            if neighbor == node_parent:
                continue

            key = edge_key(node, neighbor)

            if neighbor not in visited:
                tree_edges.add(key)
                dfs(neighbor, node, node_depth + 1)
            elif key not in tree_edges and depth[neighbor] < depth[node]:
                back_edges.append((node, neighbor))

    for start in sorted(adjacency):
        if start not in visited:
            dfs(start, None, 0)

    cycles = []

    for source, target in back_edges:
        path_source = []
        node = source

        while node is not None:
            path_source.append(node)
            node = parent[node]

        path_target = []
        node = target

        while node is not None:
            path_target.append(node)
            node = parent[node]

        source_ancestors = set(path_source)
        common = next(node for node in path_target if node in source_ancestors)

        left = path_source[: path_source.index(common) + 1]
        right = path_target[: path_target.index(common)]

        cycle = left + list(reversed(right)) + [source]
        cycles.append(_canonicalize_cycle(cycle))

    return cycles
