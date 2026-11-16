"""Figure for docs/source/concepts/electrostatics.md: the local reaction field expansion.

A constructed example, not a simulation: seven point charges beyond the 10 A
cut-off of a charge group at the origin. Their potential near the group is
compared with its Taylor expansion around the group centre, the expansion
that the local reaction field uses.

    python docs/scripts/plot_lrf_expansion.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Circle, Rectangle

COULOMB_CONSTANT = 332.0  # kcal A / (mol e^2)
CUTOFF = 10.0  # A
# x, y (A) and charge (e) of the distant charges, all in the plane of the figure
CHARGES = np.array(
    [
        (12.5, 3.0, -0.60),
        (14.0, -9.0, 0.45),
        (-11.0, 7.5, 0.55),
        (-16.0, -4.0, -0.50),
        (-2.5, 16.5, 0.40),
        (-5.0, -13.0, -0.45),
        (10.5, 13.0, -0.35),
    ]
)
ATOMS = np.array([(-0.61, 0.0), (0.61, 0.0)])  # two atoms of the group, 1.22 A apart
HALF = 3.0  # A, half-width of the enlarged region
DIVERGING = ("#184f95", "#3987e5", "#9ec5f4", "#e2e2df", "#f3aaa9", "#e66767", "#a92525")


def exact(x, y):
    """Potential of the distant charges, in kcal/mol per unit charge."""
    total = np.zeros(np.broadcast(x, y).shape)
    for cx, cy, q in CHARGES:
        total += COULOMB_CONSTANT * q / np.hypot(x - cx, y - cy)
    return total


def coefficients():
    """Potential, gradient and second derivatives at the group centre (the origin)."""
    phi0, gradient, hessian = 0.0, np.zeros(2), np.zeros((2, 2))
    for cx, cy, q in CHARGES:
        s = -np.array([cx, cy])  # from the charge to the centre
        d = np.linalg.norm(s)
        phi0 += COULOMB_CONSTANT * q / d
        gradient += -COULOMB_CONSTANT * q * s / d**3
        hessian += COULOMB_CONSTANT * q * (3.0 * np.outer(s, s) - d**2 * np.eye(2)) / d**5
    return phi0, gradient, hessian


def expansion(x, y, order):
    phi0, gradient, hessian = coefficients()
    value = np.full(np.broadcast(x, y).shape, phi0)
    if order >= 1:
        value = value + gradient[0] * x + gradient[1] * y
    if order >= 2:
        value = value + 0.5 * (hessian[0, 0] * x * x + 2.0 * hessian[0, 1] * x * y + hessian[1, 1] * y * y)
    return value


def draw(theme):
    fig = plt.figure(figsize=(7.4, 6.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.38, wspace=0.3)
    geometry, field, profile = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])

    geometry.add_patch(Circle((0, 0), CUTOFF, fc="none", ec=theme.text, lw=1.2, ls=(0, (4, 3))))
    geometry.add_patch(Rectangle((-HALF, -HALF), 2 * HALF, 2 * HALF, fc="none", ec=theme.muted, lw=1.0))
    geometry.plot(0, 0, "o", color=theme.text, ms=4)
    for cx, cy, q in CHARGES:
        colour = DIVERGING[-2] if q > 0 else DIVERGING[1]
        geometry.plot(cx, cy, "o", ms=13, color=colour, mec=theme.surface, mew=0.8)
        geometry.text(cx, cy, "+" if q > 0 else "−", ha="center", va="center", fontsize=10, color="#ffffff",
                      fontweight="bold")
        geometry.text(cx, cy - 2.3, f"{q:+.2f}".replace("-", "−"), ha="center", va="top", fontsize=7.5)
    geometry.text(0, HALF + 0.9, "cut-off, 10 Å", ha="center", va="bottom", fontsize=8.5)
    geometry.set_xlim(-21, 21)
    geometry.set_ylim(-19, 21)
    geometry.set_aspect("equal")
    geometry.set_xlabel("x (Å)")
    geometry.set_ylabel("y (Å)")
    geometry.set_title("Charges (e) beyond the cut-off of a group", fontsize=9.5)

    axis = np.linspace(-HALF, HALF, 241)
    gx, gy = np.meshgrid(axis, axis)
    phi0 = coefficients()[0]
    difference = exact(gx, gy) - phi0
    limit = float(np.abs(difference).max())
    cmap = LinearSegmentedColormap.from_list("blue_red", DIVERGING)
    mesh = field.pcolormesh(gx, gy, difference, cmap=cmap, norm=TwoSlopeNorm(0.0, -limit, limit), shading="gouraud",
                            rasterized=True)
    field.contour(gx, gy, difference, levels=np.arange(-12, 13, 2.0), colors="#3a3a37", linewidths=0.5)
    field.plot(ATOMS[:, 0], ATOMS[:, 1], "-", color="#0b0b0b", lw=1.6)
    field.plot(ATOMS[:, 0], ATOMS[:, 1], "o", color="#0b0b0b", ms=6, mec="#ffffff", mew=0.8)
    field.plot(0, 0, "+", color="#0b0b0b", ms=8, mew=1.4)
    field.set_aspect("equal")
    field.set_xlabel("x (Å)")
    field.set_ylabel("y (Å)")
    field.set_title("Potential relative to the group centre", fontsize=9.5)
    bar = fig.colorbar(mesh, ax=field, fraction=0.046, pad=0.04)
    bar.set_label("kcal/mol per e", fontsize=8.5)
    bar.outline.set_edgecolor(theme.muted)

    line = np.linspace(-HALF, HALF, 601)
    profile.axvspan(ATOMS[0, 0], ATOMS[1, 0], color=theme.muted, alpha=0.15, lw=0)
    profile.plot(line, exact(line, 0.0), color=theme.text, lw=3.5, alpha=0.35, label="exact")
    styles = ((0, ":", theme.series[3]), (1, "--", theme.series[1]), (2, "-", theme.series[0]))
    names = ("value at the centre", "first order", "second order")
    for (order, style, colour), name in zip(styles, names):
        profile.plot(line, expansion(line, 0.0, order), color=colour, ls=style, lw=1.8, label=name)
    profile.set_xlim(-HALF, HALF)
    profile.set_xlabel("x (Å), along the line through the two atoms; the shaded band is the extent of the group")
    profile.set_ylabel("potential (kcal/mol per e)")
    profile.set_title("Potential along y = 0: exact and expanded around the centre", fontsize=9.5)
    profile.grid(axis="y")
    profile.legend(ncols=4, loc="upper center", bbox_to_anchor=(0.5, -0.24), handlelength=2.4)
    fig.suptitle("Taylor expansion of the potential of distant charges around a charge-group centre", x=0.012, ha="left",
                 fontsize=10.5)
    fig.subplots_adjust(left=0.09, right=0.92, top=0.9, bottom=0.14)
    return fig


if __name__ == "__main__":
    save_themed(draw, "lrf-expansion")
    phi0, gradient, _ = coefficients()
    print(f"potential at the centre {phi0:.2f} kcal/mol/e, gradient {np.hypot(*gradient):.2f} kcal/mol/e/A")
    for reach in (0.61, 1.5, 3.0):
        line = np.linspace(-reach, reach, 2001)
        errors = [np.abs(expansion(line, 0.0, order) - exact(line, 0.0)).max() for order in (0, 1, 2)]
        print(f"|x| <= {reach}: max error order 0 {errors[0]:.3f}, order 1 {errors[1]:.4f}, order 2 {errors[2]:.5f}")
