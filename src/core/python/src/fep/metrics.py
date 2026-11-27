"""Numbers computed from the work between windows."""
import numpy as np

KB = 0.0019872041          # Boltzmann constant, kcal/(mol K)


def bar(work_forward, work_backward, temperature):
    """Free energy F_b - F_a in kcal/mol, by the Bennett acceptance ratio.

    work_forward: U_b - U_a at the frames of window a; work_backward: U_a - U_b at the frames
    of window b; both in kcal/mol.
    """
    kT = KB * temperature
    forward = np.asarray(work_forward, dtype=float) / kT
    backward = np.asarray(work_backward, dtype=float) / kT
    ratio = np.log(len(forward) / len(backward))

    def excess(df):
        """Zero at the estimate, and rising with df."""
        from_a = np.logaddexp.reduce(-np.logaddexp(0.0, ratio + forward - df))
        from_b = np.logaddexp.reduce(-np.logaddexp(0.0, -ratio + backward + df))
        return from_a - from_b

    low = min(forward.min(), -backward.max()) - 1.0
    high = max(forward.max(), -backward.min()) + 1.0
    for _ in range(100):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if excess(middle) < 0 else (low, middle)
    return 0.5 * (low + high) * kT
