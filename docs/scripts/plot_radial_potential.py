"""Figure for docs/source/concepts/spherical-boundary.md: the radial solvent restraint.

The potential and the functions for its depth and width are those of
``water_sphere`` in ``src/q6/md.f90`` (equations on
docs/source/manual/reference/boundaries.md). The effective solvent radius is
the one qprep reports for the Tyk2 sphere of the benchmark set.

    python docs/scripts/plot_radial_potential.py
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed

R_W = 23.29  # A, effective solvent radius of the Tyk2 sphere (25 A nominal radius)
FORCE_CONSTANT = 60.0  # kcal/mol/A^2, radial_force
TEMPERATURE = 298.0  # K
BOLTZMANN = 0.0019872  # kcal/mol/K

R_0 = R_W - math.sqrt(BOLTZMANN * TEMPERATURE / FORCE_CONSTANT)
DEPTH = 0.26 * math.exp(-0.19 * (R_W - 15.0)) + 0.74  # kcal/mol
WIDTH = 0.2 / (1.0 + math.exp(0.4 * (R_W - 25.0))) + 0.3  # 1/A


def wall(r):
    """Half-harmonic wall alone."""
    return np.where(r > R_0, 0.5 * FORCE_CONSTANT * (r - R_0) ** 2, 0.0)


def solvent_potential(r):
    decay = np.exp(WIDTH * np.minimum(r - R_0, 0.0))
    return np.where(r > R_0, 0.5 * FORCE_CONSTANT * (r - R_0) ** 2 - DEPTH, DEPTH * (decay**2 - 2.0 * decay))


def draw(theme):
    fig, (wide, zoom) = plt.subplots(1, 2, figsize=(7.4, 3.9), sharey=True, gridspec_kw={"width_ratios": [1.5, 1]})
    for ax, (low, high) in ((wide, (12.0, 24.0)), (zoom, (22.8, 23.6))):
        r = np.linspace(low, high, 2401)
        ax.axhline(0.0, color=theme.muted, lw=0.8)
        ax.plot(r, wall(r), color=theme.muted, lw=1.4, ls=(0, (4, 2)), label="half-harmonic wall alone")
        ax.plot(r, solvent_potential(r), color=theme.series[0], label="$V_{solvent}(r)$")
        ax.set_xlim(low, high)
        ax.set_ylim(-1.0, 2.0)
        ax.grid(axis="y")
    wide.annotate(
        "", xy=(R_0 - 0.5, -DEPTH), xytext=(R_0 - 0.5, 0.0),
        arrowprops={"arrowstyle": "<|-|>", "color": theme.text, "lw": 1.0, "shrinkA": 0, "shrinkB": 0},
    )
    wide.text(R_0 - 0.5, 0.1, f"$D_e$ = {DEPTH:.2f} kcal/mol", ha="right", va="bottom", fontsize=9)
    wide.set_ylabel("potential (kcal/mol)")
    wide.legend(loc="upper left")
    wide.text(0.5, 1.0, "12 to 24 Å", transform=wide.transAxes, ha="center", va="bottom", fontsize=10)
    zoom.text(0.5, 1.0, "22.8 to 23.6 Å", transform=zoom.transAxes, ha="center", va="bottom", fontsize=10)
    for x, y, label, align in ((R_0, 1.9, f"$r_0$ = {R_0:.2f} Å", "right"), (R_W, -0.75, f"$r_w$ = {R_W:.2f} Å", "left")):
        zoom.axvline(x, color=theme.muted, lw=1, ls=":")
        zoom.text(x + (-0.015 if align == "right" else 0.015), y, label, ha=align, va="center", fontsize=9)
    fig.supxlabel("distance of the water oxygen from the sphere centre, r (Å)", fontsize=10)
    fig.suptitle("Radial restraint on water for an effective solvent radius of 23.29 Å", x=0.012, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 1.0))
    return fig


if __name__ == "__main__":
    save_themed(draw, "radial-potential")
    print(f"r_0 = {R_0:.3f} A, D_e = {DEPTH:.3f} kcal/mol, alpha = {WIDTH:.3f} 1/A")
    for target in (0.1, 0.01):
        # distance inside r_0 at which the attraction has decayed to this fraction of D_e
        x = np.linspace(-15.0, 0.0, 150001)
        value = -(np.exp(2 * WIDTH * x) - 2 * np.exp(WIDTH * x))
        print(f"attraction is {target:.0%} of D_e at {-x[np.argmin(np.abs(value - target))]:.1f} A inside r_0")
    print(f"V = +kT ({BOLTZMANN * TEMPERATURE:.3f}) at r = {R_0 + math.sqrt(2 * (BOLTZMANN * TEMPERATURE + DEPTH) / FORCE_CONSTANT):.2f} A")
