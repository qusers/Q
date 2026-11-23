from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class RunResult:
    """Energies of one Engine.run(), in kcal/mol.

    q_energy[f, s] is the Q energy of state s at frame f, unweighted (as in the .en file).
    The potential at any lambda vector L is  E_rest + L @ q_energy[f]; E_rest is the same
    for every L, so differences between windows (exchange, BAR, MBAR) need only q_energy.
    """
    lambdas: np.ndarray     # (n_states,)
    # (n_frames,) one frame every config.energy steps; empty if energy = 0
    steps: np.ndarray
    q_energy: np.ndarray    # (n_frames, n_states)
    upot: np.ndarray        # (n_frames,)
    ukin: np.ndarray
    utot: np.ndarray
    # (n_states,) at the coordinates after the last step
    final_q: np.ndarray
    final_upot: float
    final_ukin: float
    final_utot: float
