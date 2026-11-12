"""Figure for docs/source/concepts/electrostatics.md: the ranges of the electrostatics.

A schematic drawn to scale for the QligFEP defaults: a sphere of 25 A radius
and cut-offs of 10 A. No data.

    python docs/scripts/plot_electrostatics_ranges.py
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
from figstyle import save_themed
from matplotlib.patches import Circle, Rectangle

SPHERE = 25.0  # A
CUTOFF = 10.0  # A
GROUP = (-7.0, -3.0)  # position of the charge group in the left panel
LIMIT = 33.0
ANGLE = 40.0  # direction of the drawn sphere radius, degrees
COS, SIN = math.cos(math.radians(ANGLE)), math.sin(math.radians(ANGLE))


def draw(theme):
    explicit, lrf = theme.series[0], theme.series[2]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.1))
    for ax in axes:
        ax.set_xlim(-LIMIT, LIMIT)
        ax.set_ylim(-LIMIT, LIMIT)
        ax.set_aspect("equal")
        ax.axis("off")
        # outside the sphere
        ax.add_patch(Rectangle((-LIMIT, -LIMIT), 2 * LIMIT, 2 * LIMIT, fc=theme.muted, alpha=0.13, lw=0))
        ax.add_patch(Rectangle((-LIMIT, -LIMIT), 2 * LIMIT, 2 * LIMIT, fc="none", ec=theme.muted, alpha=0.5,
                               hatch="....", lw=0))
        ax.text(-LIMIT + 1.5, LIMIT - 1.5, "outside the sphere:\nnot in the energy", ha="left", va="top",
                fontsize=8.5, linespacing=1.3)

    left, right = axes
    for ax in axes:
        ax.plot([0, SPHERE * COS], [0, SPHERE * SIN], color=theme.text, lw=0.8, zorder=5)
        ax.text(SPHERE * 0.62 * COS - 1.0, SPHERE * 0.62 * SIN + 1.0, "25 Å", ha="center", va="center", fontsize=8.5,
                rotation=ANGLE, zorder=5)
    left.add_patch(Circle((0, 0), SPHERE, fc=theme.surface, ec="none"))
    left.add_patch(Circle((0, 0), SPHERE, fc=lrf, alpha=0.30, ec="none"))
    left.add_patch(Circle(GROUP, CUTOFF, fc=theme.surface, ec="none"))
    left.add_patch(Circle(GROUP, CUTOFF, fc=explicit, alpha=0.40, ec="none"))
    left.add_patch(Circle(GROUP, CUTOFF, fc="none", ec=theme.text, lw=1.2, ls=(0, (4, 3))))
    left.add_patch(Circle((0, 0), SPHERE, fc="none", ec=theme.text, lw=1.4))
    left.plot(*GROUP, "o", color=theme.text, ms=5)
    left.text(GROUP[0], GROUP[1] + 4.6, "explicit\nCoulomb", ha="center", va="center", fontsize=9, linespacing=1.3)
    left.plot([GROUP[0], GROUP[0]], [GROUP[1], GROUP[1] - CUTOFF], color=theme.text, lw=0.8)
    left.text(GROUP[0] + 0.8, GROUP[1] - CUTOFF / 2, "10 Å", ha="left", va="center", fontsize=8.5, linespacing=1.3)
    left.text(1.0, -16.8, "local reaction field\n(multipole expansion)", ha="center", va="center", fontsize=9,
              linespacing=1.3)
    left.set_title("Charge group that is not a Q-atom", loc="center", fontsize=10)

    right.add_patch(Circle((0, 0), SPHERE, fc=theme.surface, ec="none"))
    right.add_patch(Circle((0, 0), SPHERE, fc=explicit, alpha=0.40, ec="none"))
    right.add_patch(Circle((0, 0), SPHERE, fc="none", ec=theme.text, lw=1.4))
    right.plot(0, 0, "o", color=theme.text, ms=5)
    right.text(-1.2, 0, "Q-atom", ha="right", va="center", fontsize=8.5)
    right.text(0, -9.0, "explicit Coulomb\nwith every atom in the sphere", ha="center", va="center", fontsize=9,
               linespacing=1.3)
    right.set_title("Q-atom", loc="center", fontsize=10)

    fig.suptitle("Treatment of electrostatic interactions by distance, QligFEP defaults", x=0.012, ha="left",
                 fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 1.0))
    return fig


if __name__ == "__main__":
    save_themed(draw, "electrostatics-ranges")
