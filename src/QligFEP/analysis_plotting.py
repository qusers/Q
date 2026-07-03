"""Shared plotting and result-shaping helpers for FEP analysis.

These are used by both the equilibrium (analyze_FEP) and non-equilibrium (analyze_neq)
analyzers to turn a per-edge results table into the experimental-vs-calculated ddG plot.
"""

from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import kendalltau

from .logger import logger


def prepare_df(json_dict, experimental_data: bool = True):
    pref = "dg" if "dg_error" in json_dict["edges"][0] else "ddg"
    df = pd.DataFrame(json_dict["edges"])
    if experimental_data:
        df = (
            df.assign(
                ddg_value=lambda x: x[pref + "_value"],
                residual=lambda x: x[pref + "_value"] - x["Q_ddG_avg"],
                residual_abs=lambda x: x["residual"].abs(),
            )
            .sort_values("residual_abs", ascending=False)
            .drop(columns="residual_abs")
        )
    df = df.assign(
        fep_name=lambda x: "FEP_" + x["from"] + "_" + x["to"],
    )
    return df


def bootstrap_statistic(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    statistic: str,
    ci: float = 0.95,
    nbootstrap: int = 1000,
    rng=None,
) -> dict:
    """Point estimate and bootstrap confidence interval for a comparison statistic.

    Resamples the paired ``(y_true, y_pred)`` values with replacement ``nbootstrap``
    times and reports the statistic computed on the full data (``mle``) together with
    the ``ci`` confidence-interval bounds (``low``/``high``). Supported statistics are
    "RMSE", "MUE" and "KTAU".
    """
    if rng is None:
        rng = np.random.default_rng(12345)
    compute = {
        "RMSE": lambda a, b: float(np.sqrt(np.mean((a - b) ** 2))),
        "MUE": lambda a, b: float(np.mean(np.abs(a - b))),
        "KTAU": lambda a, b: float(kendalltau(a, b)[0]),
    }[statistic]

    n = len(y_true)
    s_n = np.empty(nbootstrap)
    for replicate in range(nbootstrap):
        idx = rng.choice(n, size=n, replace=True)
        s_n[replicate] = compute(y_true[idx], y_pred[idx])
    s_n.sort()

    low_frac = (1.0 - ci) / 2.0
    low_idx = int(np.floor(nbootstrap * low_frac))
    high_idx = min(int(np.ceil(nbootstrap * (1.0 - low_frac))), nbootstrap - 1)
    return {
        "mle": compute(y_true, y_pred),
        "low": float(s_n[low_idx]),
        "high": float(s_n[high_idx]),
    }


def create_ddG_plot(
    results_df: pd.DataFrame,
    margin: float = 1.0,
    xylims: tuple | None = None,
    output_path: str | None = None,
    target_name: str | None = None,
    savefig: bool = False,
    font: str | None = None,
):
    """Creates the ddG plot for the FEP that has already been analyzed. The plot will
    show the experimental (X axis) vs mean predicted values (Y axis), with error bars
    representing the standard error of the mean (SEM). Points are colored based on
    their deviation from experimental values (blue = 0 deviation, red = 3+ kcal/mol).

    Args:
        reuslts_df: pd.DataFrame with the results from the FEP, output from `prepare_df`.
        margin: margin value to be added/subtracted to the max/min values obtained. Defaults to 1.0.
        xylims: if values are passed, x&y min will be xylims[0] and max will be [1]. Defaults to None.
        output_path: path to save the plot. If None, the plot will not be saved. Defaults to None.
        target_name: name of the target protein to be added in the plot. Defaults to None.
        savefig: if True, will save the plot to the output_path. Defaults to False.

    Returns:
        the matplotlib figure and axis objects (fig, ax).
    """
    fep_names = results_df["fep_name"].values
    avg_values = results_df["Q_ddG_avg"].values
    sem_values = results_df["Q_ddG_sem"].values
    exp_values = results_df["ddg_value"].values
    nan_val_idxs = np.where(np.isnan(avg_values))[0]
    if len(nan_val_idxs) > 0:
        logger.warning(f"Dropping FEPs with nan values: {fep_names[nan_val_idxs]}")
        mask = ~np.isin(np.arange(len(avg_values)), nan_val_idxs)
        avg_values = avg_values[mask]
        sem_values = sem_values[mask]
        exp_values = exp_values[mask]

    # Calculate absolute deviations for coloring
    deviations = np.abs(avg_values - exp_values)

    # Create colormap from blue to red, with max at 3 kcal/mol
    norm = mcolors.Normalize(vmin=0, vmax=4)
    cmap = cm.RdBu_r

    ## CALCULATE STATISTICS
    def result_to_latex(res, latexify_each=False):  # TODO: move this out of this method?
        """Round the statistic's output to one decimal case and return a LaTeX string."""
        mle = round(res["mle"], 2)
        low = round(res["low"], 2)
        high = round(res["high"], 2)

        if latexify_each:
            return f"${mle:.2f}_{{{low}}}^{{{high}}}$"
        else:
            return f"{mle:.2f}_{{{low}}}^{{{high}}}"

    statistics = ["RMSE", "MUE", "KTAU"]
    stats_dict = {}
    for stat in statistics:
        boot = bootstrap_statistic(avg_values, exp_values, statistic=stat)
        stats_dict[stat] = result_to_latex(boot)

    if xylims is not None:
        assert len(xylims) == 2, "xylims must be a tuple with 2 elements."
        assert xylims[0] < xylims[1], "xylims[0] must be smaller than xylims[1]."
        min_val = xylims[0]
        max_val = xylims[1]
    else:
        all_values = np.concatenate([avg_values, exp_values])
        min_val = min(all_values) - margin
        max_val = max(all_values) + margin

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Plot colored points with error bars
    scatter = plt.errorbar(
        exp_values,
        avg_values,
        fmt="o",
        yerr=sem_values,
        ecolor="gray",
        elinewidth=1.5,
        capsize=2,
        zorder=4,
        linestyle="None",
        markersize=8,
    )

    # Remove the default markers and add colored ones
    scatter[0].remove()

    # Add the colored scatter points
    scatter_points = plt.scatter(
        exp_values,
        avg_values,
        c=deviations,
        cmap=cmap,
        norm=norm,
        s=45,
        zorder=5,
        edgecolors="black",
        linewidths=0.5,
        alpha=0.8,
    )

    plt.plot([min_val, max_val], [min_val, max_val], "k-", linewidth=1.5, zorder=3)  # Black identity line

    # Highlight predictions within 1 and 2 kcal/mol of the experimental affinity
    ax.fill_between(
        [min_val, max_val],
        [min_val - 1, max_val - 1],
        [min_val + 1, max_val + 1],
        color="darkgray",
        alpha=0.3,
        zorder=2,
    )
    ax.fill_between(
        [min_val, max_val],
        [min_val - 2, max_val - 2],
        [min_val + 2, max_val + 2],
        color="lightgray",
        alpha=0.3,
        zorder=1,
    )

    # Add colorbar
    cbar = plt.colorbar(scatter_points, ax=ax, shrink=0.40, aspect=10, anchor=(0.0, 0.85), pad=0.05)

    cbar.set_label("|Deviation| (kcal/mol)", rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_ticks([0, 1, 2, 3, 4])
    cbar.set_ticklabels(["0", "1", "2", "3", "≥4"])
    cbar.ax.tick_params(labelsize=10)

    # set labels, make it square and add legend
    plt.title(
        f"{(target_name + ' ' if target_name is not None else '')}"
        r"$\Delta\Delta \text{G}_{\text{BAR}}$ ($\mathrm{N}="
        f"{len(exp_values)}$)"
    )
    plt.xlabel(r"$\Delta\Delta G_{exp} (kcal/mol)$")
    plt.ylabel(r"$\Delta\Delta G_{calc} (kcal/mol)$")
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    ax.set_aspect("equal", adjustable="box")

    # add statistics to the plot
    unit = r"\frac{kcal}{mol}"
    text_body = (
        f"$\\tau = {stats_dict['KTAU']}$",
        f"RMSE = ${stats_dict['RMSE']}  {unit}$",
        f"MUE = ${stats_dict['MUE']}  {unit}$",
    )
    logger.info(f"Stats: {' '.join(text_body)}")
    hori_height = 0.35
    spacing = 0.085
    txt_positions = (
        (1.04, hori_height),
        (1.04, hori_height - spacing),
        (1.04, hori_height - spacing * 2),
    )
    for txt_position, body in zip(txt_positions, text_body):
        plt.text(
            *txt_position,
            body,
            fontsize=12,
            verticalalignment="bottom",
            horizontalalignment="left",
            transform=ax.transAxes,
            fontproperties=font,
        )

    legend_elements = [
        Line2D([0], [0], color="k", linestyle="-", label="Identity line"),
        Patch(facecolor="darkgray", alpha=0.3, label="Within 1 kcal/mol"),
        Patch(facecolor="lightgray", alpha=0.3, label="Within 2 kcal/mol"),
    ]

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.04, 0),
        loc="lower left",
        borderaxespad=0,
        frameon=False,
    )
    # Remove top and right spines using matplotlib
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)

    if savefig:
        if output_path is None:
            output_path = Path().cwd()
            logger.info("Using default name to save the plot at the current working directory...")
            fig.savefig(f"{target_name}_ddG_plot.png", dpi=300, bbox_inches="tight")
            return fig, ax
        if isinstance(output_path, str):
            output_path = Path(output_path)
        assert isinstance(output_path, Path), "output_path must be a string or a Path object."
        if output_path.is_dir():
            output_path = output_path / f"{target_name}_ddG_plot.png"
            logger.info(f"Using default name to save the plot at {output_path}")
        elif output_path.exists():
            logger.warning(f"File {output_path} already exists. Overwriting...")
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig, ax
