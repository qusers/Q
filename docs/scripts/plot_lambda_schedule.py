"""Figure for docs/source/qligfep/equilibrium.md: the default lambda schedule.

The schedule is the one written by ``QligFEP.get_lambdas`` in
``src/QligFEP/qligfep.py`` for ``--windows 100 --sampling sigmoidal``. The
formula is repeated here so the script runs without installing QligFEP.

The two ligand images are optional. They are read in place from a folder that
is not part of the repository and are embedded, downscaled, in the SVG:

    python docs/scripts/plot_lambda_schedule.py \
        --ligand-a /path/to/memprep/4l.png --ligand-b /path/to/memprep/4r.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed

WINDOWS = 100
K_SIGMOIDAL = -1.1


def rational_sigmoid(t: float, k: float) -> float:
    return k * t / (k - t + 1.0)


def sigmoidal_schedule(windows: int = WINDOWS) -> np.ndarray:
    """Weight of state 1 in each window, from 1.000 to 0.000 (windows + 1 values)."""
    half = windows // 2
    upper = [round(0.5 * (rational_sigmoid(i / half, K_SIGMOIDAL) + 1.0), 3) for i in range(half + 1)]
    lower = [round(1.0 - value, 3) for value in upper[1:]]
    return np.array(upper[::-1] + lower)


def load_ligand(path: Path, size: int = 440) -> np.ndarray:
    """Crop to the molecule and downscale, resampling with premultiplied alpha."""
    from PIL import Image

    image = Image.open(path).convert("RGBA")
    image = image.crop(image.getbbox()).convert("RGBa")
    image.thumbnail((size, size), Image.LANCZOS)
    return np.asarray(image.convert("RGBA"))


def make_draw(ligand_a: np.ndarray | None, ligand_b: np.ndarray | None):
    lambdas = sigmoidal_schedule()
    index = np.arange(lambdas.size)
    start = int(np.argmin(np.abs(lambdas - 0.5)))

    def draw(theme):
        fig, (ax, rug) = plt.subplots(
            1, 2, figsize=(6.6, 4.2), sharey=True, gridspec_kw={"width_ratios": [1, 0.045], "wspace": 0.04}
        )
        colour = theme.series[0]
        ax.plot(index, lambdas, color=colour, lw=1.2, zorder=2)
        ax.plot(index, lambdas, "o", color=colour, ms=3, zorder=3, label="window")
        ax.plot(
            index[start], lambdas[start], "o", ms=9, mfc=theme.surface, mec=theme.text, mew=1.6, zorder=4,
            label="first window run (λ = 0.5)",
        )
        for target in (start - 9, start + 9):
            ax.annotate(
                "", xy=(target, 0.5), xytext=(index[start] + np.sign(target - start) * 2.2, 0.5),
                arrowprops={"arrowstyle": "-|>", "color": theme.text, "lw": 1.2}, zorder=4,
            )
        ax.text(start - 10, 0.5, "towards λ = 1", ha="right", va="center", fontsize=9, color=theme.text)
        ax.text(start + 10, 0.5, "towards λ = 0", ha="left", va="center", fontsize=9, color=theme.text)

        ax.set_xlim(-2, WINDOWS + 2)
        ax.set_ylim(-0.03, 1.03)
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_xlabel("window index")
        ax.set_ylabel("λ, weight of state 1 (ligand A)")
        ax.grid(axis="y")
        ax.set_title("Default QligFEP λ schedule: 101 windows with sigmoidal spacing")
        ax.legend(loc="upper right", handlelength=1.2)

        for image, box, label in (
            (ligand_a, (1, 0.57, 29, 0.35), "ligand A"),
            (ligand_b, (71, 0.09, 29, 0.35), "ligand B"),
        ):
            if image is None:
                continue
            inset = ax.inset_axes(box, transform=ax.transData)
            inset.imshow(image, interpolation="none")
            inset.axis("off")
            below = label == "ligand A"
            ax.text(
                box[0] + box[2] / 2, box[1] - 0.01 if below else box[1] + box[3] + 0.01, label,
                ha="center", va="top" if below else "bottom", fontsize=9,
            )

        rug.hlines(lambdas, 0, 1, color=theme.text, lw=0.5, alpha=0.75)
        rug.set_xlim(0, 1)
        rug.axis("off")
        rug.set_title("λ values", fontsize=8, loc="center", color=theme.muted)
        fig.subplots_adjust(left=0.105, right=0.975, top=0.92, bottom=0.125)
        return fig

    return draw


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ligand-a", type=Path, help="PNG of the ligand that is real in state 1")
    parser.add_argument("--ligand-b", type=Path, help="PNG of the ligand that is real in state 2")
    args = parser.parse_args()
    ligand_a = load_ligand(args.ligand_a) if args.ligand_a else None
    ligand_b = load_ligand(args.ligand_b) if args.ligand_b else None
    save_themed(make_draw(ligand_a, ligand_b), "lambda-schedule")
    lambdas = sigmoidal_schedule()
    steps = -np.diff(lambdas)
    print(f"{lambdas.size} windows; first values {lambdas[:4]}, around the start {lambdas[49:52]}")
    print(f"smallest step {steps.min():.3f}, largest step {steps.max():.3f}")
    print(f"windows with lambda >= 0.9: {(lambdas >= 0.9).sum()}, <= 0.1: {(lambdas <= 0.1).sum()}")


if __name__ == "__main__":
    main()
