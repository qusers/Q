"""Figure for docs/source/concepts/fep-in-q.md: the two ligands in the two states.

Each ligand is drawn once per state: in colour where it has its real atom
types and charges, and as a grey shade where its atoms are ``DUM`` with charge 0.

The two ligand images are read in place from a folder that is not part of the
repository and are embedded, downscaled, in the SVG:

    python docs/scripts/plot_dual_topology_states.py \
        --ligand-a /path/to/memprep/4l.png --ligand-b /path/to/memprep/4r.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed
from matplotlib.colors import to_rgb
from plot_lambda_schedule import load_ligand

SHADE_ALPHA = 0.55


def shade(image: np.ndarray, colour: str) -> np.ndarray:
    """Replace the atom colours by one grey and make the image translucent."""
    out = image.copy()
    out[..., :3] = np.round(np.array(to_rgb(colour)) * 255)
    out[..., 3] = np.round(image[..., 3] * SHADE_ALPHA)
    return out


def make_draw(ligand_a: np.ndarray, ligand_b: np.ndarray):
    def draw(theme):
        fig = plt.figure(figsize=(6.6, 4.5))
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.text(0.02, 0.955, "Dual topology: ligand A and ligand B in state 1 and state 2", fontsize=10.5, va="center")

        columns = (0.36, 0.78)
        rows = (0.62, 0.22)
        for x, state in zip(columns, ("State 1", "State 2")):
            ax.text(x, 0.86, state, ha="center", va="center", fontsize=10.5, fontweight="bold")
        for y, name in zip(rows, ("Ligand A", "Ligand B")):
            ax.text(0.02, y + 0.02, name, ha="left", va="center", fontsize=10.5, fontweight="bold")
        ax.plot([0.02, 0.98], [0.82, 0.82], color=theme.muted, lw=0.8)
        ax.plot([0.02, 0.98], [0.42, 0.42], color=theme.grid, lw=0.8)
        ax.plot([0.57, 0.57], [0.03, 0.89], color=theme.grid, lw=0.8)

        real = "real atom types and charges"
        dummy = "atom type DUM, charge 0"
        cells = (
            (0, 0, ligand_a, real),
            (0, 1, shade(ligand_a, theme.muted), dummy),
            (1, 0, shade(ligand_b, theme.muted), dummy),
            (1, 1, ligand_b, real),
        )
        for row, column, image, label in cells:
            x, y = columns[column], rows[row]
            inset = ax.inset_axes((x - 0.17, y - 0.115, 0.34, 0.30), transform=ax.transData)
            inset.imshow(image, interpolation="none")
            inset.axis("off")
            ax.text(x, y - 0.15, label, ha="center", va="center", fontsize=9, color=theme.text)
        return fig

    return draw


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ligand-a", type=Path, required=True, help="PNG of the ligand that is real in state 1")
    parser.add_argument("--ligand-b", type=Path, required=True, help="PNG of the ligand that is real in state 2")
    args = parser.parse_args()
    save_themed(make_draw(load_ligand(args.ligand_a), load_ligand(args.ligand_b)), "dual-topology-states")


if __name__ == "__main__":
    main()
