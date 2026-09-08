"""Reproduce four archived CMET reductions without changing the raw reference.

Run from an environment importing this worktree's QligFEP; see RESULTS.md.
The downloaded raw data are deliberately not included in Git.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess

from QligFEP import endpoint_trim


ROOT = Path(__file__).resolve().parent


def digest(path):
    return endpoint_trim._hash(path)


def verify_download(selection):
    """Reconstruct the exact remote sha256sum inventory without a network call."""
    files = []
    for item in selection["selection"]:
        for local, remote in [("inputfiles", "remote_input"), ("energy", "remote_energy")]:
            for path in sorted((ROOT / "raw" / item["label"] / local).iterdir()):
                if not path.is_file():
                    raise ValueError(f"Unexpected non-file in reference: {path}")
                files.append((item[remote]+"/"+path.name, path))
    for path in sorted((ROOT / "raw/engine").iterdir()):
        files.append((selection["campaign"]+"/"+selection["reference_engine_subdir"]+"/"+path.name, path))
    for path in sorted((ROOT / "raw/metadata").iterdir()):
        prefix = "/scripts/" if path.name == "analyze_cmet_results.py" else "/"
        files.append((selection["campaign"]+prefix+path.name, path))
    inventory = "".join(f"{digest(path)}  {remote}\n" for remote, path in sorted(files))
    got = hashlib.sha256(inventory.encode()).hexdigest()
    if got != selection["remote_inventory_sha256"]:
        raise ValueError(f"Downloaded reference differs from remote inventory: {got}")
    validation = json.loads((ROOT / "raw/metadata/qfep_safebar_validation.json").read_text())
    archived_binary_md5 = hashlib.md5((ROOT / "raw/engine/qfep").read_bytes()).hexdigest()
    if archived_binary_md5 != validation["qfep_safebar_md5"]:
        raise ValueError("Archived executable differs from historical validation record")
    return {"sha256": got, "file_count": len(files),
            "bytes": sum(path.stat().st_size for _, path in files),
            "archived_binary_md5_matches_validation": True}


def bar_rows(path):
    rows, active = [], False
    for line in path.read_text().splitlines():
        if line.startswith("# Part 6: BAR"):
            if rows or active:
                raise ValueError("Multiple BAR sections")
            active = True
        elif line.startswith("# Part"):
            active = False
        elif active and line.strip() and not line.lstrip().startswith("#"):
            # Native QFEP format: 3x,f8.6,5f9.3 (three columns used for BAR).
            fields = [line[3:11], line[11:20], line[20:29]]
            rows.append(tuple(float(value) for value in fields))
    if not rows:
        raise ValueError(f"No BAR rows in {path}")
    return rows


def records(path):
    with path.open("rb") as stream:
        while marker := stream.read(4):
            if len(marker) != 4:
                raise ValueError("Truncated Fortran record")
            size = struct.unpack("<i", marker)[0]
            if size < 0 or size > path.stat().st_size:
                raise ValueError("Unsupported Fortran record format")
            data = stream.read(size)
            if len(data) != size or stream.read(4) != marker:
                raise ValueError("Damaged Fortran record")
            yield data


def restart_offsets(path):
    blocks = list(records(path))
    if len(blocks) != 3 or len(blocks[2]) != 16:
        raise ValueError("Expected spherical restart with three single-precision shell offsets")
    n, *offsets = struct.unpack("<i3f", blocks[2])
    if n != 3 or not all(math.isfinite(value) for value in offsets):
        raise ValueError("Invalid shell offsets")
    return offsets


def check_energy_mappings(manifest):
    counts = []
    for window in manifest["retained"]:
        expected = [float(value) for value in window["lambdas"]]
        frames, block = 0, []
        for payload in records(Path(window["source_energy"])):
            block.append(payload)
            if len(block) == 3:
                if block[2] != b"":
                    raise ValueError("Unexpected off-diagonal energy record")
                for state, data in enumerate(block[:2], 1):
                    if len(data) != 124:
                        raise ValueError("Unexpected Q energy record size")
                    values = struct.unpack("<i15d", data)
                    if values[0] != state or abs(values[1]-expected[state-1]) > 1e-12:
                        raise ValueError("Stored mapping disagrees with original MD input")
                    if not all(math.isfinite(value) for value in values[1:]):
                        raise ValueError("Nonfinite retained energy component")
                frames += 1
                block = []
        if block or frames == 0:
            raise ValueError("Incomplete energy file")
        counts.append(frames)
    return {"frames_per_retained_window_min": min(counts),
            "frames_per_retained_window_max": max(counts),
            "all_stored_lambda_pairs_match_md_inputs": True}


def boundary_parameters(path):
    """Independently implement the archived rounded-log Born convention."""
    result = {}
    for line in path.read_text().splitlines():
        if "Eff. solvent radius" in line:
            result["radius"] = float(line.split("=")[-1])
        elif "q_from" not in result and line.lstrip().startswith("SUM"):
            result["q_from"], result["q_to"] = map(float, line.split()[1:3])
        elif "Total charge of non-Q atoms" in line:
            result["q_env"] = float(line.split("=")[-1])
    if set(result) != {"radius", "q_from", "q_to", "q_env"}:
        raise ValueError("Incomplete archived boundary log")
    return result


def run(output, binary):
    selection = json.loads((ROOT / "selection.json").read_text())
    verified = verify_download(selection)
    if output.exists():
        raise ValueError("Use a new output directory; previous results are never overwritten")
    binary = binary.resolve(strict=True)
    if binary == (ROOT / "raw/engine/qfep").resolve():
        raise ValueError("Rebuild the archived source locally; do not run the archived Linux binary")
    output.mkdir(parents=True)
    reports = []
    for item in selection["selection"]:
        original = ROOT / "raw" / item["label"]
        analysis = output / item["label"]
        source_input = original / "inputfiles/qfep.inp"
        header = [endpoint_trim._tokens(line) for line in source_input.read_text().splitlines()
                  if endpoint_trim._tokens(line)]
        control = header[2]
        config = json.loads((original / "inputfiles/fep_config.json").read_text())
        if config["perstate_born"] or not config["perstate_polarization"]:
            raise ValueError("Selected campaign is not the expected post-hoc-Born baseline")
        manifest = endpoint_trim.prepare(source_input, original / "inputfiles", analysis,
                                        selection["lambda1_end"], selection["lambda1_start"],
                                        energy_dir=original / "energy")
        mapping_check = check_energy_mappings(manifest)
        with (analysis / "qfep.inp").open() as source, (analysis / "qfep.out").open("w") as dest:
            process = subprocess.run([str(binary)], stdin=source, stdout=dest, stderr=subprocess.PIPE,
                                     text=True, cwd=analysis, timeout=30)
        (analysis / "qfep.stderr").write_text(process.stderr)
        if process.returncode:
            raise ValueError(f"QFEP failed for {item['label']}: {process.stderr}")
        summary = endpoint_trim.summarize(analysis)
        archived = bar_rows(original / "energy/qfep.out")
        if len(archived) != int(header[0][0]):
            raise ValueError("Archived full-ladder BAR row count is incomplete")
        interior = [row for row in archived if float(selection["lambda1_end"]) <= row[0] <=
                    float(selection["lambda1_start"])]
        current = bar_rows(analysis / "qfep.out")
        if len(interior) != len(current):
            raise ValueError("Archived and reproduced intervals differ")
        errors = []
        for old, new in zip(interior, current):
            if old[0] != new[0] or not all(math.isfinite(v) for v in old+new):
                raise ValueError("Invalid or mismatched reference BAR row")
            errors.append(abs(new[2]-(old[2]-interior[0][2])))
        if max(errors) > selection["reference_tolerance_kcal_mol"]:
            raise ValueError(f"Archived BAR reproduction failed for {item['label']}: {max(errors)}")
        diag = boundary_parameters(original / "energy/md_0000_1000.log")
        convention = selection["historical_posthoc_born"]
        coefficient = convention["ke"]*(1-1/convention["epsilon"])/(2*diag["radius"])
        full_born = -coefficient*((diag["q_env"]+diag["q_to"])**2 -
                                 (diag["q_env"]+diag["q_from"])**2)
        multiplier = float(summary["state_constant_gap_multiplier"])
        offsets = restart_offsets(original / "energy/eq5.re")
        report = {"label": item["label"], "direction": item["direction"], "leg": item["leg"],
                  "bar_kcal_mol": summary["bar_kcal_mol"],
                  "archived_cumulative_difference_kcal_mol": interior[-1][2]-interior[0][2],
                  "maximum_cumulative_row_error_kcal_mol": max(errors),
                  "compared_rows": len(current), "retained_pairs": summary["retained_pairs"],
                  "kT_kcal_mol": float(control[0]), "discarded_frames_per_window": int(control[1]),
                  "qfep_alpha2_kcal_mol": float(header[5][0]),
                  "qfep_alpha_interval_contribution_kcal_mol": multiplier*float(header[5][0]),
                  "historical_born_parameters": diag, "born_full_state_gap_kcal_mol": full_born,
                  "born_interval_shift_kcal_mol": multiplier*full_born,
                  "corrected_interval_kcal_mol": summary["bar_kcal_mol"]+multiplier*full_born,
                  "eq5_shell_offsets_radians": offsets,
                  "configured_softcore_method": config["softcore_method"],
                  "endpoint_manifest_sha256": digest(analysis / "endpoint-trim.json"),
                  "archived_qfep_output_sha256": digest(original / "energy/qfep.out"),
                  "reproduced_qfep_output_sha256": summary["qfep_output_sha256"], **mapping_check}
        reports.append(report)
        print(f"{item['label']}: BAR {summary['bar_kcal_mol']:.3f}; max archived-row error {max(errors):.6f}", flush=True)
    by_label = {record["label"]: record for record in reports}
    directional = {}
    for direction in ("fwd", "rev"):
        water, protein = (by_label[f"{direction}-{leg}"] for leg in ("water", "protein"))
        if water["qfep_alpha_interval_contribution_kcal_mol"] != protein["qfep_alpha_interval_contribution_kcal_mol"]:
            raise ValueError("QFEP analysis constants do not cancel across legs")
        directional[direction] = {
            "protein_minus_water_bar_kcal_mol": protein["bar_kcal_mol"]-water["bar_kcal_mol"],
            "protein_minus_water_born_kcal_mol": protein["born_interval_shift_kcal_mol"]-water["born_interval_shift_kcal_mol"],
            "corrected_protein_minus_water_kcal_mol": protein["corrected_interval_kcal_mol"]-water["corrected_interval_kcal_mol"]}
    offsets_match = all(by_label[f"fwd-{leg}"]["eq5_shell_offsets_radians"] ==
                        by_label[f"rev-{leg}"]["eq5_shell_offsets_radians"] for leg in ("water", "protein"))
    result = {"schema_version": 1, "gate": "archived_interior_BAR_reproduction_passed",
              "full_endpoint_free_energy": False, "independent_replicates": 1,
              "lambda1_start": selection["lambda1_start"], "lambda1_end": selection["lambda1_end"],
              "download_verification": verified, "selection_sha256": digest(ROOT / "selection.json"),
              "reproducer_sha256": digest(Path(__file__)),
              "endpoint_utility_sha256": digest(Path(endpoint_trim.__file__)),
              "local_archived_solver_binary_sha256": digest(binary),
              "archived_solver_source_sha256": {p.name: digest(p) for p in sorted((ROOT / "raw/engine").glob("*.f90"))},
              "reductions": reports, "directional_results": directional,
              "forward_reverse_offsets_identical": offsets_match,
              "interpretation": "Numerical reproduction only. No equilibrium/overlap qualification or full-endpoint estimate. Different forward/reverse offsets prevent interpreting their sum as same-Hamiltonian closure."}
    (output / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--qfep", type=Path, default=ROOT / "engine-build/qfep")
    args = parser.parse_args()
    run(args.output_dir.resolve(), args.qfep)
