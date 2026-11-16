"""Figure for docs/source/qligfep/equilibrium.md: the four ``--sampling`` options.

The schedules are those written by ``QligFEP.get_lambdas`` in
``src/QligFEP/qligfep.py`` for ``--windows 100``. The formulas are repeated
here so the script runs without installing QligFEP.

    python docs/scripts/plot_lambda_spacings.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from figstyle import save_themed
from plot_lambda_schedule import WINDOWS, rational_sigmoid, sigmoidal_schedule

# k of the rational sigmoid for the options that are not "sigmoidal"
K = {"linear": 1000.0, "exponential": -1.1, "reverse_exponential": 1.1}
ORDER = ("linear", "sigmoidal", "exponential", "reverse_exponential")
LINESTYLES = ("-", "-", "--", ":")


def schedule(sampling: str, windows: int = WINDOWS) -> np.ndarray:
    """Weight of state 1 in each window, from 1.000 to 0.000 (windows + 1 values)."""
    if sampling == "sigmoidal":
        return sigmoidal_schedule(windows)
    values = [round(rational_sigmoid(i / windows, K[sampling]), 3) for i in range(windows + 1)]
    return np.array(values[::-1])


def draw(theme):
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    index = np.arange(WINDOWS + 1)
    for name, colour, style in zip(ORDER, theme.series, LINESTYLES):
        ax.plot(index, schedule(name), color=colour, ls=style, label=name)
    ax.set_xlim(-2, WINDOWS + 2)
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_xlabel("window index")
    ax.set_ylabel("λ, weight of state 1 (ligand A)")
    ax.grid(axis="y")
    ax.set_title("λ schedules of the four --sampling options, 101 windows")
    ax.legend(title="--sampling", loc="lower left", handlelength=2.6)
    fig.tight_layout()
    return fig


def main() -> None:
    save_themed(draw, "lambda-spacings")
    for name in ORDER:
        lambdas = schedule(name)
        steps = -np.diff(lambdas)
        print(
            f"{name:20s} step min {steps.min():.3f} max {steps.max():.3f}; "
            f"windows at lambda >= 0.9: {(lambdas >= 0.9).sum():2d}, <= 0.1: {(lambdas <= 0.1).sum():2d}; "
            f"has 0.500: {bool((lambdas == 0.5).any())}"
        )


if __name__ == "__main__":
    main()
