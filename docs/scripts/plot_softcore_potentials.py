"""Figures for docs/source/concepts/softcore.md: weight-dependent soft-core pair potentials.

The formulas are those of ``src/q6/qsoftcore.f90`` on the branch
``feature/softcore-charge-changing`` (Gapsys et al. 2012; Beutler et al. 1994),
in Q's units (angstrom, kcal/mol). Each curve is the state potential
V_s(r; lambda_s), before it is multiplied by the state weight lambda_s.

    python docs/scripts/plot_softcore_potentials.py
    python docs/scripts/plot_softcore_potentials.py --table   # values for checking against Q
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed

WEIGHTS = (1.0, 0.5, 0.1, 0.0)
COULOMB_CONSTANT = 332.0  # kcal A / (mol e^2), Q default

# Carbonyl oxygen (type O) and TIP3P oxygen (type OW), AMBER14sb.prm, arithmetic rule
R_STAR = 1.6612 + 1.7683  # A
EPSILON = (0.21 * 0.152) ** 0.5  # kcal/mol
SIGMA = R_STAR / 2.0 ** (1.0 / 6.0)  # A
C12 = EPSILON * R_STAR**12
C6 = 2.0 * EPSILON * R_STAR**6

# Carbonyl oxygen and TIP3P hydrogen
Q_I, Q_J = -0.55, 0.417  # e
QQ = COULOMB_CONSTANT * Q_I * Q_J  # kcal A / mol

# Defaults of the branch (qatom.f90)
GAPSYS_ALPHA_LJ, GAPSYS_ALPHA_Q, GAPSYS_SIGMA_Q = 0.85, 3.0, 1.0  # -, A, e^-2
BEUTLER_ALPHA_LJ, BEUTLER_ALPHA_C = 0.5, 100.0  # -, A^2


def gapsys_lj_radius(weight):
    return GAPSYS_ALPHA_LJ * ((26.0 / 7.0) * SIGMA**6 * (1.0 - weight)) ** (1.0 / 6.0)


def gapsys_coulomb_radius(weight):
    return GAPSYS_ALPHA_Q * (1.0 + GAPSYS_SIGMA_Q * abs(Q_I * Q_J)) * (1.0 - weight) ** (1.0 / 6.0)


def gapsys_lj(r, weight):
    """Lennard-Jones above the switching radius, quadratic (linear force) below it."""
    rsc = gapsys_lj_radius(weight)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        hard = C12 / r**12 - C6 / r**6
    if rsc <= 0.0:
        return hard
    a = 78.0 * C12 / rsc**14 - 21.0 * C6 / rsc**8
    b = 168.0 * C12 / rsc**13 - 48.0 * C6 / rsc**7
    c = 91.0 * C12 / rsc**12 - 28.0 * C6 / rsc**6
    return np.where(r < rsc, a * r**2 - b * r + c, hard)


def gapsys_coulomb(r, weight):
    rsc = gapsys_coulomb_radius(weight)
    with np.errstate(divide="ignore"):
        hard = QQ / r
    if rsc <= 0.0:
        return hard
    soft = QQ / rsc**3 * r**2 - 3.0 * QQ / rsc**2 * r + 3.0 * QQ / rsc
    return np.where(r < rsc, soft, hard)


def beutler_lj(r, weight):
    """Beutler Eq. 7 with n = 1: r^6 is replaced by r^6 + alpha_LJ sigma^6 (1 - weight)^2."""
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        denominator = r**6 + BEUTLER_ALPHA_LJ * SIGMA**6 * (1.0 - weight) ** 2
        return C12 / denominator**2 - C6 / denominator


def beutler_coulomb(r, weight):
    """Beutler Eq. 8: r^2 is replaced by r^2 + alpha_C (1 - weight)^2."""
    with np.errstate(divide="ignore"):
        return QQ / np.sqrt(r**2 + BEUTLER_ALPHA_C * (1.0 - weight) ** 2)


def make_draw(kind):
    if kind == "lj":
        panels = (("Gapsys", gapsys_lj, gapsys_lj_radius), ("Beutler", beutler_lj, None))
        ylim, ylabel = (-4.0, 110.0), "V$_{LJ}$(r; λ)  (kcal/mol)"
        title = "Soft-core Lennard-Jones potential of an O–OW pair at state weights 1, 0.5, 0.1 and 0"
        legend_loc = "upper right"
    else:
        panels = (("Gapsys", gapsys_coulomb, gapsys_coulomb_radius), ("Beutler", beutler_coulomb, None))
        ylim, ylabel = (-85.0, 3.0), "V$_{el}$(r; λ)  (kcal/mol)"
        title = "Soft-core Coulomb potential of a −0.55 e / +0.417 e pair at state weights 1, 0.5, 0.1 and 0"
        legend_loc = "lower right"

    def draw(theme):
        r = np.linspace(0.0, 6.0, 1201)
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9), sharey=True)
        for ax, (name, potential, radius) in zip(axes, panels):
            ax.axhline(0.0, color=theme.muted, lw=0.8)
            for weight, shade in zip(WEIGHTS, theme.ramp):
                ax.plot(r, potential(r, weight), color=shade, label=f"λ = {weight:g}")
                if radius is not None and radius(weight) > 0.0:
                    rsc = radius(weight)
                    ax.plot(rsc, potential(np.array([rsc]), weight)[0], "o", color=shade, ms=6,
                            mec=theme.surface, mew=0.8)
            ax.set_xlim(0.0, 6.0)
            ax.set_ylim(*ylim)
            ax.set_xlabel("distance r (Å)")
            ax.grid(axis="y")
            ax.text(0.5, 1.0, name, transform=ax.transAxes, ha="center", va="bottom", fontsize=10)
        axes[0].set_ylabel(ylabel)
        axes[1].legend(title="state weight", loc=legend_loc)
        fig.suptitle(title, x=0.012, ha="left", fontsize=10.5)
        fig.tight_layout(rect=(0, 0, 1, 1.0))
        return fig

    return draw


def table():
    """Values on a short grid, for comparison with the Fortran module."""
    r = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 6.0])
    for weight in WEIGHTS:
        columns = (gapsys_lj(r, weight), gapsys_coulomb(r, weight), beutler_lj(r, weight), beutler_coulomb(r, weight))
        for i, distance in enumerate(r):
            print(f"{weight:4.1f} {distance:4.1f} " + " ".join(f"{column[i]:16.8e}" for column in columns))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", action="store_true", help="print values instead of writing the figures")
    args = parser.parse_args()
    if args.table:
        table()
        return
    save_themed(make_draw("lj"), "softcore-lj")
    save_themed(make_draw("coulomb"), "softcore-coulomb")
    zero = np.array([0.0])
    print(f"pair: R* {R_STAR:.4f} A, epsilon {EPSILON:.4f} kcal/mol, sigma {SIGMA:.3f} A, qq {QQ:.2f} kcal A/mol")
    for weight in WEIGHTS[1:]:
        print(
            f"lambda={weight:g}: Gapsys r_LJ {gapsys_lj_radius(weight):.2f} A, V_LJ(0) {gapsys_lj(zero, weight)[0]:.2f}; "
            f"r_Q {gapsys_coulomb_radius(weight):.2f} A, V_el(0) {gapsys_coulomb(zero, weight)[0]:.1f}; "
            f"Beutler V_LJ(0) {beutler_lj(zero, weight)[0]:.2f}, V_el(0) {beutler_coulomb(zero, weight)[0]:.1f}"
        )


if __name__ == "__main__":
    main()
