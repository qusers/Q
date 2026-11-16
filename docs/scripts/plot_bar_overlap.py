"""Figure for docs/source/concepts/softcore.md: BAR overlap per lambda step.

The figure is drawn from ``docs/scripts/data/bar_overlap_acceptance.tsv``, a
summary of the per-step data of the charge-changing campaign (26 edges, two
legs, 5 replicas). The summary is rebuilt from the dashboard data, which is
not part of the repository, with:

    python docs/scripts/plot_bar_overlap.py --windows-json /path/to/bar-overlap-dashboard/windows.json

Without the option the committed summary is used:

    python docs/scripts/plot_bar_overlap.py
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed

SUMMARY = Path(__file__).resolve().parent / "data" / "bar_overlap_acceptance.tsv"
LEGS = (("1.water", "water"), ("2.protein", "protein"))
THRESHOLD = 0.05  # acceptance threshold used in the campaign
FIELDS = ["lambda_start", "lambda_end"] + [f"{name}_{stat}" for _, name in LEGS for stat in ("p10", "median", "p90")]


def summarise(windows_json: Path) -> None:
    """Median and 10th/90th percentile over all edges and replicas, per step and leg."""
    data = json.loads(windows_json.read_text())
    windows = data["windows"]
    columns = {}
    for leg, name in LEGS:
        values = np.array([replica["acc"] for edge in data["edges"] for replica in edge["legs"][leg]])
        print(f"{name}: {values.shape[0]} edge-replica series, {values.shape[1]} steps")
        low, median, high = np.percentile(values, [10, 50, 90], axis=0)
        columns.update({f"{name}_p10": low, f"{name}_median": median, f"{name}_p90": high})
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        for i in range(len(windows) - 1):
            writer.writerow([f"{windows[i]:.3f}", f"{windows[i + 1]:.3f}"] + [f"{columns[k][i]:.4f}" for k in FIELDS[2:]])
    print(f"wrote {SUMMARY}")


def load() -> dict[str, np.ndarray]:
    with SUMMARY.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {field: np.array([float(row[field]) for row in rows]) for field in FIELDS}


def steps(start, end, values):
    """x and y of a line that is flat across each step."""
    return np.column_stack([start, end]).ravel(), np.repeat(values, 2)


def draw(theme):
    table = load()
    start, end = table["lambda_start"], table["lambda_end"]
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 5.2), sharex=True, sharey=True)
    for ax, (_, name), colour in zip(axes, LEGS, theme.series):
        x, low = steps(start, end, table[f"{name}_p10"])
        _, high = steps(start, end, table[f"{name}_p90"])
        _, median = steps(start, end, table[f"{name}_median"])
        ax.fill_between(x, low, high, color=colour, alpha=0.25, lw=0, label="10th–90th percentile")
        ax.plot(x, median, color=colour, lw=1.6, label="median")
        ax.axhline(THRESHOLD, color=theme.muted, lw=1, ls=":")
        ax.text(0.995, THRESHOLD + 0.012, "threshold used in the campaign, 0.05", ha="right", va="bottom",
                fontsize=8.5, color=theme.muted)
        ax.set_ylim(0, 0.52)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_xlim(0, 1)
        ax.grid(axis="y")
        ax.set_ylabel("Fermi acceptance")
        ax.text(0.0, 1.02, f"{name} leg", transform=ax.transAxes, ha="left", va="bottom", fontsize=10)
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.12), ncols=2)
    axes[1].set_xlabel("λ, weight of state 2; each segment spans one step between adjacent windows")
    fig.suptitle("BAR acceptance per λ step, 26 charge-changing edges", x=0.012, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 1.0))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--windows-json", type=Path, help="rebuild the summary from the dashboard data first")
    args = parser.parse_args()
    if args.windows_json:
        summarise(args.windows_json)
    save_themed(draw, "bar-overlap")
    table = load()
    for _, name in LEGS:
        median, low = table[f"{name}_median"], table[f"{name}_p10"]
        i = int(median.argmin())
        print(
            f"{name}: lowest median {median[i]:.3f} at step {table['lambda_start'][i]:.3f}-{table['lambda_end'][i]:.3f}; "
            f"lowest 10th percentile {low.min():.3f}; steps with median >= 0.45: {(median >= 0.45).sum()} of {median.size}; "
            f"median < 0.40 between {table['lambda_start'][median < 0.40].min():.3f} and {table['lambda_end'][median < 0.40].max():.3f}"
        )


if __name__ == "__main__":
    main()
