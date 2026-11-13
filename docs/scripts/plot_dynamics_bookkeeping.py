"""Figure for docs/source/concepts/fep-in-q.md: dynamics and bookkeeping.

A schematic with no data: what qdyn does with the state energies of one
configuration, and what is stored and re-weighted by qfep.

    python docs/scripts/plot_dynamics_bookkeeping.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from figstyle import save_themed
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def box(ax, theme, x, y, width, height, text, edge=None):
    ax.add_patch(
        FancyBboxPatch(
            (x - width / 2, y - height / 2), width, height, boxstyle="round,pad=0,rounding_size=0.012",
            fc="none", ec=edge or theme.muted, lw=1.2,
        )
    )
    ax.text(x, y, text, ha="center", va="center", fontsize=9, linespacing=1.5)


def arrow(ax, theme, start, end, style="arc3,rad=0"):
    ax.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=11, lw=1.2, color=theme.text, connectionstyle=style,
            shrinkA=0, shrinkB=0,
        )
    )


def draw(theme):
    fig = plt.figure(figsize=(7.4, 4.1))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    dynamics, bookkeeping = theme.series[0], theme.series[1]
    top, bottom, height = 0.64, 0.25, 0.17

    ax.text(0.02, 0.955, "One configuration $x$ of window $m$: dynamics and bookkeeping", fontsize=10.5, va="center")
    ax.text(0.02, 0.865, "Dynamics: qdyn, window $m$", fontsize=9.5, fontweight="bold", va="center")
    ax.text(0.02, 0.085, "Bookkeeping: energy file and qfep", fontsize=9.5, fontweight="bold", va="center")

    box(ax, theme, 0.10, top, 0.16, height, "configuration\n$x$", dynamics)
    box(ax, theme, 0.34, top, 0.19, height, "state energies\n" r"$\varepsilon_1(x),\ \varepsilon_2(x)$", dynamics)
    box(
        ax, theme, 0.62, top, 0.27, height,
        "mapping potential\n" r"$V_m = \lambda_1^{(m)}\varepsilon_1 + \lambda_2^{(m)}\varepsilon_2$", dynamics,
    )
    box(ax, theme, 0.89, top, 0.14, height, "forces", dynamics)
    arrow(ax, theme, (0.18, top), (0.245, top))
    arrow(ax, theme, (0.435, top), (0.485, top))
    arrow(ax, theme, (0.755, top), (0.82, top))
    loop = top + height / 2 + 0.065
    ax.plot([0.89, 0.89, 0.10], [top + height / 2, loop, loop], color=theme.text, lw=1.2, solid_joinstyle="miter")
    arrow(ax, theme, (0.10, loop), (0.10, top + height / 2))
    ax.text(0.495, loop + 0.012, "next step", ha="center", va="bottom", fontsize=8.5, color=theme.muted)

    box(
        ax, theme, 0.19, bottom, 0.30, height,
        ".en file\n" r"$\varepsilon_1(x),\ \varepsilon_2(x)$, weights of window $m$", bookkeeping,
    )
    box(
        ax, theme, 0.595, bottom, 0.39, height,
        "qfep: re-weighting to window $m+1$\n"
        r"$V_{m+1}(x) - V_m(x) = \sum_s (\lambda_s^{(m+1)} - \lambda_s^{(m)})\,\varepsilon_s(x)$",
        bookkeeping,
    )
    box(ax, theme, 0.905, bottom, 0.15, height, "Zwanzig, TI,\noverlap\nsampling, BAR", bookkeeping)
    arrow(ax, theme, (0.30, top - height / 2), (0.30, bottom + height / 2))
    ax.text(0.315, 0.445, "written at every\nenergy-output step", ha="left", va="center", fontsize=8.5,
            color=theme.muted, linespacing=1.4)
    arrow(ax, theme, (0.34, bottom), (0.40, bottom))
    arrow(ax, theme, (0.79, bottom), (0.83, bottom))
    ax.text(
        0.98, 0.085,
        r"re-weighting uses the stored $\varepsilon_s(x)$: $\varepsilon_s$ must not depend on the weights",
        ha="right", va="center", fontsize=8.5, color=theme.muted,
    )
    return fig


if __name__ == "__main__":
    save_themed(draw, "dynamics-bookkeeping")
