import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import glob
import sys
from scipy.stats import gaussian_kde

#this is a script to plot cumulative work vs time for forward and reverse processes with sideways histogram + KDE of final work values
# Usage: python analysis_histo_plot.py <log_folder>
# The log_folder should contain the log files from the nonequilibrium simulations
# use this tool to visualize the work distributions and their evolution over time, and check if the work values are reasonable
# and if the choice of parameters for the simulations was appropriate.

def plot_forward_reverse_work_with_hist_kde(Wf, Wr, outfile="cumulative_work_with_hist_kde.png", time_step=50):
    """
    Plot cumulative work vs time for forward and reverse processes with
    sideways histogram + KDE of final work values.
    """

    # === Transform reverse work ===
    Wr_shifted = Wr - Wr[:, -1][:, None]  # make final value zero
    Wr_rev = Wr_shifted[:, ::-1]          # reverse time

    # === Compute mean ± SEM ===
    wf_mean = Wf.mean(axis=0)
    wf_sem = Wf.std(axis=0, ddof=1) / np.sqrt(Wf.shape[0])
    wr_mean = Wr_rev.mean(axis=0)
    wr_sem = Wr_rev.std(axis=0, ddof=1) / np.sqrt(Wr_rev.shape[0])

    # === Time axis ===
    n_steps = Wf.shape[1]
    time = np.arange(n_steps) * time_step

    # === Figure with gridspec for side histogram/KDE ===
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4,1], wspace=0.05)
    ax_main = fig.add_subplot(gs[0])
    ax_side = fig.add_subplot(gs[1], sharey=ax_main)

    # === Plot ensemble trajectories (light for context) ===
    for w in Wf:
        ax_main.plot(time, w, color='C0', alpha=0.2, lw=1)
    for w in Wr_rev:
        ax_main.plot(time, w, color='C1', alpha=0.2, lw=1)

    # === Plot mean ± SEM ===
    ax_main.plot(time, wf_mean, color='C0', lw=0.5, label='Forward')
    ax_main.fill_between(time, wf_mean - wf_sem, wf_mean + wf_sem, color='C0', alpha=0.3)
    ax_main.plot(time, wr_mean, color='C1', lw=0.5, label='Reverse (reversed, sign-flipped)')
    ax_main.fill_between(time, wr_mean - wr_sem, wr_mean + wr_sem, color='C1', alpha=0.3)

    ax_main.set_xlabel("Time (simulation steps)")
    ax_main.set_ylabel("Cumulative work")
    ax_main.set_title("Work Accumulation Over Time")
    ax_main.legend(frameon=False)
    ax_main.grid(True, alpha=0.2)

    # === Side histogram + KDE ===
    wf_final = Wf[:, -1]
    wr_final = Wr_rev[:, -1]

    # Histogram (horizontal)
    bins = 30
    ax_side.hist(wf_final, bins=bins, orientation='horizontal', color='C0', alpha=0.4, label='Forward', density = True)
    ax_side.hist(wr_final, bins=bins, orientation='horizontal', color='C1', alpha=0.4, label='Reverse', density = True)
    
    

    # KDE
    y_vals = np.linspace(min(wf_final.min(), wr_final.min()),
                         max(wf_final.max(), wr_final.max()), 200)
    kde_wf = gaussian_kde(wf_final)
    kde_wr = gaussian_kde(wr_final)
    ax_side.plot(kde_wf(y_vals), y_vals, color='C0', lw=2)
    ax_side.plot(kde_wr(y_vals), y_vals, color='C1', lw=2)
    ax_side.fill_betweenx(y_vals, 0, kde_wf(y_vals), color='C0', alpha=0.2)
    ax_side.fill_betweenx(y_vals, 0, kde_wr(y_vals), color='C1', alpha=0.2)

    ax_side.set_xlabel("Count / Density")
    ax_side.tick_params(axis="y", labelleft=False)
    ax_side.legend(frameon=False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved plot with histogram + KDE to {outfile}")



def get_works_from_logs(folder):
    Wf = []
    Wr = []
    wf_singrep = []
    wr_singrep = []

    forward_logs = glob.glob(os.path.join(folder, "neq_1_*.log"))
    reverse_logs = glob.glob(os.path.join(folder, "neq_0_*.log"))
    file_default_array = []
    # Process forward logs
    for fname in forward_logs:
        wf_singrep = []
      
        with open(fname, "r") as f:
            lines = f.readlines()
            if not lines:
                continue
            if lines[-1].strip() == "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" or lines[-1].strip() != "================================================================================":
                continue
            for line in lines:
                if "work" in line:
                    try:
                        wf_singrep.append(float(line.split()[6]))
                    except Exception:
                        print(f"Failed to parse work in {fname}")
        wf_singrep=np.array(wf_singrep)                
        wf_singrep=list(wf_singrep)
        Wf.append(wf_singrep)


    # Process reverse logs
    for fname in reverse_logs:

        wr_singrep = []
        with open(fname, "r") as f:
            lines = f.readlines()
            if not lines:
                continue
            if lines[-1].strip() == "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" or lines[-1].strip() != "================================================================================":
                continue
            for line in lines:
                if "work" in line:
                    try:
                        wr_singrep.append(float(line.split()[6]))
                    except Exception:
                        print(f"Failed to parse work in {fname}")
        wr_singrep=np.array(wr_singrep)
        wr_singrep=list(wr_singrep)
        Wr.append(wr_singrep)

    print(f'work values in {folder}: Forward={len(forward_logs)}, Reverse={len(reverse_logs)}')
    print(len(Wf),len(Wr))
    Wr = np.array(Wr)
    Wf = np.array(Wf)
    print(Wf.T[-1])
    return Wf, Wr

Wf, Wr = get_works_from_logs(sys.argv[1])
print(Wf.T[-1])
plot_forward_reverse_work_with_hist_kde(Wf, Wr, "cumulative_work_with_hist_costmatch_protein.png", time_step=50)

