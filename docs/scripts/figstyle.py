"""Shared style for the documentation figures.

Every plot is written twice, as ``<name>-light.svg`` and ``<name>-dark.svg``,
with a transparent background. The pages include the light file; the
``themed_figures`` Sphinx extension adds the dark one and the Furo classes
that show the file matching the active theme.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGURES = Path(__file__).resolve().parents[1] / "source" / "_static" / "figures"


@dataclass(frozen=True)
class Theme:
    name: str
    text: str
    muted: str
    grid: str
    surface: str  # page background of the Furo theme; used for halos, not painted
    series: tuple[str, ...]  # categorical, fixed order
    ramp: tuple[str, ...]  # one hue; first entry has the most contrast with the page


LIGHT = Theme(
    name="light",
    text="#0b0b0b",
    muted="#52514e",
    grid="#e2e2df",
    surface="#ffffff",
    series=("#2a78d6", "#eb6834", "#008300", "#4a3aa7"),
    ramp=("#0d366b", "#256abf", "#5598e7", "#9ec5f4"),
)
DARK = Theme(
    name="dark",
    text="#ffffff",
    muted="#c3c2b7",
    grid="#3a3a37",
    surface="#131416",
    series=("#3987e5", "#d95926", "#199e70", "#9085e9"),
    ramp=("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab"),
)
THEMES = (LIGHT, DARK)


def rc(theme: Theme) -> dict:
    return {
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.transparent": True,
        "text.color": theme.text,
        "axes.labelcolor": theme.text,
        "axes.edgecolor": theme.muted,
        "axes.titlecolor": theme.text,
        "xtick.color": theme.muted,
        "ytick.color": theme.muted,
        "xtick.labelcolor": theme.text,
        "ytick.labelcolor": theme.text,
        "grid.color": theme.grid,
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 10.5,
        "axes.titlelocation": "left",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.title_fontsize": 9,
        "legend.frameon": False,
        "lines.linewidth": 2,
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "svg.hashsalt": "q-docs",  # stable ids, so regenerated files diff cleanly
    }


def save_themed(draw, name: str) -> list[Path]:
    """Call ``draw(theme)`` for each theme and save the returned figure as SVG."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    written = []
    for theme in THEMES:
        with plt.rc_context(rc(theme)):
            fig = draw(theme)
            out = FIGURES / f"{name}-{theme.name}.svg"
            fig.savefig(out, metadata={"Date": None})
            plt.close(fig)
        written.append(out)
        print(f"wrote {out.relative_to(FIGURES.parents[3])}")
    return written
