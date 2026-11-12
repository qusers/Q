"""Figure for docs/source/concepts/softcore.md: a lambda-weighted Lennard-Jones wall.

Run from the repository root:
    python docs/scripts/plot_endpoint_wall.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EPSILON = 0.10  # kcal/mol, typical heavy-atom pair
R_MIN = 3.6  # A
KT = 0.596  # kcal/mol at 300 K
WEIGHTS = [1.0, 0.1, 0.01, 0.001]
SHADES = ["#08306b", "#2171b5", "#6baed6", "#9ecae1"]  # one hue, dark = full weight
OUT = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures" / "endpoint-wall.svg"


def lennard_jones(r):
    x = (R_MIN / r) ** 6
    return EPSILON * (x * x - 2.0 * x)


def wall_radius(weight):
    """Distance at which weight * V_LJ equals kT (bisection on the repulsive branch)."""
    lo, hi = 0.5, R_MIN
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if weight * lennard_jones(mid) > KT:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    r = np.linspace(1.2, 5.0, 800)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.axhline(KT, color="#888888", lw=1, ls="--")
    ax.text(4.95, KT + 0.08, "kT at 300 K", ha="right", va="bottom", fontsize=9, color="#555555")
    for weight, shade in zip(WEIGHTS, SHADES):
        ax.plot(r, weight * lennard_jones(r), color=shade, lw=2, label=f"λ = {weight:g}")
        ax.plot(wall_radius(weight), KT, "o", color=shade, ms=6)
    ax.set_ylim(-0.3, 4.0)
    ax.set_xlim(1.2, 5.0)
    ax.set_xlabel("distance r (Å)")
    ax.set_ylabel("λ · V_LJ(r)  (kcal/mol)")
    ax.set_title("A small state weight moves the repulsive wall, it does not remove it", fontsize=10.5, loc="left")
    ax.legend(frameon=False, title="state weight", fontsize=9, title_fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e6e6e6", lw=0.8)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    for weight in WEIGHTS:
        print(f"lambda={weight:g}  wall at {wall_radius(weight):.2f} A")


if __name__ == "__main__":
    main()
