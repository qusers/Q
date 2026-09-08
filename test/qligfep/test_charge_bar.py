"""Known Gaussian/harmonic free energies and conservative uncertainty handling."""
import math
import subprocess
import sys

import numpy as np
import pytest
from scipy.special import ndtri

from QligFEP import charge_bar as cb
from QligFEP.charge_analysis import native_boltzmann, aligned_observables


def test_cli_requires_explicit_discard_choice(tmp_path):
    result = subprocess.run(
        [sys.executable, '-m', 'QligFEP.charge_analysis', str(tmp_path/'plan.json')],
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 2
    assert '--discard-frames' in result.stderr
    assert not result.stdout


@pytest.mark.parametrize('counts', [(2000, 2000), (2000, 7000), (7000, 2000)])
def test_bar_gaussian_work_with_known_free_energy_and_unequal_counts(counts):
    delta, sigma = 2.3, 1.2
    forward = delta+.5*sigma**2+sigma*ndtri((np.arange(counts[0])+.5)/counts[0])
    reverse = -delta+.5*sigma**2+sigma*ndtri((np.arange(counts[1])+.5)/counts[1])
    result = cb.bar(forward, reverse)
    assert result['delta_f'] == pytest.approx(delta, abs=1e-4)
    assert 0 < result['normalized_overlap'] < 1
    assert result['sample_counts'] == list(counts)
    assert abs(result['log_equation_residual']) < 1e-10


def test_bar_direction_and_constant_shift_identities():
    rng = np.random.default_rng(73)
    forward, reverse = rng.normal(2, 1, 300), rng.normal(-1, 1, 500)
    original = cb.bar(forward, reverse)
    reversed_result = cb.bar(reverse, forward)
    shifted = cb.bar(forward+1e6, reverse-1e6)
    assert reversed_result['delta_f'] == pytest.approx(-original['delta_f'], abs=1e-11)
    assert shifted['delta_f']-1e6 == pytest.approx(original['delta_f'], abs=2e-9)
    assert shifted['normalized_overlap'] == pytest.approx(original['normalized_overlap'], abs=1e-10)
    assert reversed_result['normalized_overlap'] == pytest.approx(original['normalized_overlap'], abs=1e-12)


def test_constant_gap_has_exact_bar_and_full_overlap():
    result = cb.bar(np.full(100, 123.), np.full(731, -123.))
    assert result['delta_f'] == pytest.approx(123., abs=1e-12)
    assert result['normalized_overlap'] == pytest.approx(1., abs=1e-12)
    assert cb.correlation(np.full(300, .1))['constant_trace'] is True


def test_exact_discrete_partition_functions_and_overlap():
    # Exact p_A=(1/4,3/4), p_B=(3/5,2/5), with unequal exact population counts.
    # The mixture here is (1/2,1/2); normalized overlap sums to 0.9.
    delta = 1.37
    gap = delta+np.log(np.array([.25, .75])/np.array([.6, .4]))
    result = cb.bar(np.repeat(gap, [100, 300]), -np.repeat(gap, [600, 400]))
    assert result['delta_f'] == pytest.approx(delta, abs=1e-12)
    assert result['normalized_overlap'] == pytest.approx(.9, abs=1e-12)


def test_disjoint_work_does_not_become_zero_uncertainty_evidence():
    result = cb.bar(np.full(300, 10000.), np.full(300, 10000.))
    assert result['normalized_overlap'] == 0
    assert result['overlap_resolved'] is False
    ladder = cb.ladder([np.full(300, 10000.), np.full(300, -10000.)], [0, 1], beta=1, bootstrap=50)
    assert ladder['gap_statistical_gates_passed'] is False
    assert ladder['conditional_interval_95'] is None
    assert ladder['delta_g_0_to_sign'] is None
    assert any('overlap' in failure for failure in ladder['failures'])


@pytest.mark.parametrize('bad', [[], [math.nan], [math.inf], [[1., 2.]]])
def test_invalid_energy_series_rejected(bad):
    with pytest.raises(ValueError, match='one-dimensional'):
        cb.bar(bad, [0.])


def harmonic_ladder(n=1000):
    # U_w(x)=x^2/2+w*(d*x+c), x~N(-w*d,1), beta=1.
    # Exact F(1)-F(0)=c-d^2/2, no force-field or charged-result fitting.
    rng = np.random.default_rng(751)
    weights, d, c = [0., .5, 1.], .5, 1.3
    return [d*rng.normal(-w*d, 1., n)+c for w in weights], weights, c-.5*d*d


def test_harmonic_ladder_bootstrap_preserves_shared_edge_covariance():
    gaps, weights, truth = harmonic_ladder()
    result = cb.ladder(gaps, weights, beta=1., block_length=10, bootstrap=150, seed=8)
    assert result['gap_statistical_gates_passed'], result['failures']
    assert result['delta_g_0_to_sign'] == pytest.approx(truth, abs=.035)
    assert result['conditional_interval_95'][0] < truth < result['conditional_interval_95'][1]
    covariance = np.asarray(result['edge_covariance'])
    assert covariance.shape == (2, 2) and covariance[0, 1] > 0
    assert covariance.sum() == pytest.approx(result['conditional_standard_error']**2, rel=1e-12)
    assert result['production_ready'] is False


def test_reversed_ladder_reports_same_canonical_transformation():
    gaps, weights, _ = harmonic_ladder(500)
    forward = cb.ladder(gaps, weights, beta=1., block_length=10, bootstrap=50)
    reverse = cb.ladder(gaps[::-1], weights[::-1], beta=1., block_length=10, bootstrap=50)
    assert reverse['delta_g_0_to_sign'] == pytest.approx(forward['delta_g_0_to_sign'], abs=1e-12)
    assert reverse['path_direction'] == -1


def test_ar1_correlation_scale_and_short_block_rejection():
    rng, phi = np.random.default_rng(901), .9
    data = np.zeros(50000)
    innovations = rng.normal(0, math.sqrt(1-phi*phi), len(data))
    for i in range(1, len(data)):
        data[i] = phi*data[i-1]+innovations[i]
    diagnostic = cb.correlation(data)
    assert diagnostic['g'] == pytest.approx((1+phi)/(1-phi), rel=.2)
    result = cb.ladder([data, data-.1], [0, 1], beta=.1, block_length=1, bootstrap=50)
    assert result['conditional_interval_95'] is None
    assert any('block shorter' in failure for failure in result['failures'])


def test_short_native_like_trace_does_not_get_confidence_interval():
    gaps, weights, _ = harmonic_ladder(9)
    result = cb.ladder(gaps, weights, beta=1., bootstrap=50)
    assert result['conditional_interval_95'] is None
    assert not result['gap_statistical_gates_passed']


def test_slow_nonenergy_observable_controls_blocks_without_changing_bar():
    gaps, weights, _ = harmonic_ladder(1000)
    rng = np.random.default_rng(922)
    slow = np.zeros(1000)
    for index in range(1, len(slow)):
        slow[index] = .98*slow[index-1]+rng.normal()
    panel = [{'polarization': slow, 'empty_density_bin': np.zeros(1000)} for _ in gaps]
    gap_only = cb.ladder(gaps, weights, beta=1., bootstrap=50)
    result = cb.ladder(gaps, weights, beta=1., bootstrap=50, observables=panel)
    assert gap_only['statistical_gates_passed']
    assert result['delta_g_0_to_sign'] == gap_only['delta_g_0_to_sign']
    assert not result['statistical_gates_passed']
    assert result['conditional_interval_95'] is None
    for metric in result['windows']:
        assert metric['slowest_observable'] == 'polarization'
        assert metric['g'] > metric['gap_g']
        assert metric['block_length'] >= 5*metric['g']
        assert metric['observable_correlations']['empty_density_bin']['constant_trace']


def test_explicit_blocks_must_cover_slow_observable_too():
    gaps, weights, _ = harmonic_ladder(1000)
    signal = np.sin(np.arange(1000)/30.)
    result = cb.ladder(gaps, weights, beta=1., block_length=10, bootstrap=50,
                       observables=[{'slow': signal} for _ in gaps])
    assert result['conditional_interval_95'] is None
    assert any('block shorter' in failure for failure in result['failures'])


@pytest.mark.parametrize('panels', [[{}]*3, [{'energy_gap': [0]*100}]*3,
                                  [{'temp': [0]*99}]*3, [{'temp': [0]*100}]*2])
def test_missing_misaligned_or_reserved_observables_rejected(panels):
    gaps, weights, _ = harmonic_ladder(100)
    with pytest.raises(ValueError, match='observable|Observable'):
        cb.ladder(gaps, weights, beta=1., bootstrap=50, observables=panels)


def test_native_observables_match_energy_steps_not_trace_array_indices(monkeypatch, tmp_path):
    from QligFEP import charge_diagnostics
    records = [{'step': step, 'values': [float(step)]*13, 'density': [step]*21,
                'shells': [[step, .1*step, .2*step]]} for step in range(0, 41, 10)]
    monkeypatch.setattr(charge_diagnostics, 'trace', lambda *a, **k: records)
    definition = {'signature': {'md': {'steps': '40'}, 'intervals': {'output': '10', 'energy': '10'}}}
    native = {'meta': [0, 0, 100, 0, 0, 1], 'parameters': [10.]}
    panel, matching = aligned_observables(tmp_path/'native.log', definition, native, 1, 2)
    assert panel['temperature_total_kelvin'] == [20., 30.]
    assert panel['q_radius_angstrom'] == [20., 30.]
    assert panel['shell_1_radial_sum'] == [2., 3.]
    assert matching['first_retained_energy_step'] == 20
    assert matching['last_retained_energy_step'] == 30
    assert matching['matched_frames'] == 2
    records.pop(2)  # An absent middle sample must not be interpolated or skipped.
    with pytest.raises(ValueError, match='every retained energy step'):
        aligned_observables(tmp_path/'native.log', definition, native, 1, 2)


@pytest.mark.parametrize('weights', [[.001, .5, .999], [0, 1, .5], [0, 0, 1]])
def test_partial_or_nonmonotonic_endpoint_ladder_rejected(weights):
    gaps, _, _ = harmonic_ladder(100)
    with pytest.raises(ValueError, match='full monotonic'):
        cb.ladder(gaps, weights, beta=1, bootstrap=50)


@pytest.mark.parametrize('declaration,valid', [
    ('real, parameter :: boltz = 0.001986', True),
    ('real(8), parameter :: boltz = 0.001986', False),
    ('real, parameter :: boltz = -0.001986', False),
    ('real, parameter :: boltz = 0.001986\nreal, parameter :: boltz = 0.1', False),
])
def test_beta_comes_from_supported_native_precision(tmp_path, declaration, valid):
    source = tmp_path/'src/q6'
    source.mkdir(parents=True)
    (source/'md.f90').write_text(declaration+'\n')
    if valid:
        result = native_boltzmann(tmp_path/'build.json')
        assert result == float(np.float32(.001986))
        assert result != .001986
    else:
        with pytest.raises(ValueError, match='Boltzmann'):
            native_boltzmann(tmp_path/'build.json')
