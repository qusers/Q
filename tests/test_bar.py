"""Tests for the BAR estimator used by the non-equilibrium FEP analysis.

The forward and reverse switching work distributions of a process satisfy the Crooks
fluctuation theorem, P_f(W) / P_r(-W) = exp(beta * (W - dF)). For Gaussian works of equal
variance sigma^2 this fixes the means to mu_f = dF + sigma^2/2 (forward) and
mu_r = -dF + sigma^2/2 (reverse), so BAR should recover the known dF. Working in kT units
(beta = 1) keeps the synthetic data unit-agnostic.
"""

import numpy as np
import pytest

from QligFEP.analyze_neq import bar_delta_f, bar_with_uncertainty, work_overlap


def crooks_consistent_works(true_dF, sigma, n, seed):
    rng = np.random.default_rng(seed)
    work_forward = rng.normal(true_dF + sigma**2 / 2, sigma, n)
    work_reverse = rng.normal(-true_dF + sigma**2 / 2, sigma, n)
    return work_forward, work_reverse


@pytest.mark.parametrize("true_dF", [-5.0, 0.0, 2.5, 8.0])
def test_bar_recovers_known_free_energy(true_dF):
    work_forward, work_reverse = crooks_consistent_works(true_dF, sigma=2.0, n=5000, seed=1)
    estimate = bar_delta_f(work_forward, work_reverse, beta=1.0)
    assert estimate == pytest.approx(true_dF, abs=0.2)


def test_bar_is_antisymmetric_in_swapping_directions():
    work_forward, work_reverse = crooks_consistent_works(3.0, sigma=2.0, n=5000, seed=2)
    forward = bar_delta_f(work_forward, work_reverse, beta=1.0)
    # Swapping forward<->reverse must flip the sign of the estimated free energy.
    reverse = bar_delta_f(work_reverse, work_forward, beta=1.0)
    assert forward == pytest.approx(-reverse, abs=1e-6)


def test_bootstrap_ci_brackets_the_estimate():
    work_forward, work_reverse = crooks_consistent_works(2.5, sigma=2.0, n=4000, seed=3)
    rng = np.random.default_rng(0)
    dF, dF_err, overlap = bar_with_uncertainty(work_forward, work_reverse, beta=1.0, n_bootstrap=500, rng=rng)
    assert dF == pytest.approx(2.5, abs=0.25)
    assert dF_err > 0
    assert abs(dF - 2.5) < 5 * dF_err
    assert 0.0 <= overlap <= 1.0


def test_overlap_high_for_identical_distributions():
    rng = np.random.default_rng(4)
    # Wf and -Wr drawn from the same distribution -> near-complete overlap.
    work_forward = rng.normal(0.0, 1.0, 4000)
    work_reverse = -rng.normal(0.0, 1.0, 4000)
    assert work_overlap(work_forward, work_reverse) > 0.8


def test_bar_raises_without_both_directions():
    with pytest.raises(ValueError):
        bar_delta_f([], [1.0, 2.0], beta=1.0)
