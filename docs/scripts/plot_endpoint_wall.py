"""Figure for docs/source/concepts/softcore.md: a Lennard-Jones pair scaled by a state weight.

    python docs/scripts/plot_endpoint_wall.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed

EPSILON = 0.10  # kcal/mol, typical heavy-atom pair
R_MIN = 3.6  # A
KT = 0.596  # kcal/mol at 300 K
WEIGHTS = (1.0, 0.1, 0.01, 0.001)


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


def draw(theme):
    r = np.linspace(1.2, 5.0, 800)
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.axhline(KT, color=theme.muted, lw=1, ls="--")
    ax.text(4.95, KT + 0.08, "kT at 300 K", ha="right", va="bottom", fontsize=9, color=theme.muted)
    for weight, shade in zip(WEIGHTS, theme.ramp):
        ax.plot(r, weight * lennard_jones(r), color=shade, label=f"λ = {weight:g}")
        ax.plot(wall_radius(weight), KT, "o", color=shade, ms=6, mec=theme.surface, mew=0.8)
    ax.set_ylim(-0.3, 4.0)
    ax.set_xlim(1.2, 5.0)
    ax.set_xlabel("distance r (Å)")
    ax.set_ylabel("λ · V$_{LJ}$(r)  (kcal/mol)")
    ax.set_title("Lennard-Jones pair potential at state weights 1, 0.1, 0.01 and 0.001")
    ax.legend(title="state weight")
    ax.grid(axis="y")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    save_themed(draw, "endpoint-wall")
    for weight in WEIGHTS:
        print(f"lambda={weight:g}  wall at {wall_radius(weight):.2f} A")
