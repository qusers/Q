"""BAR and conditional block uncertainty for a fixed charge-only lambda ladder.

This is an analysis estimator, not a sampler. Neither overlap nor a bootstrap
interval establishes equilibrium sampling by the underlying MD/thermostat.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, log_expit, logsumexp


def _array(values):
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or not len(result) or not np.isfinite(result).all():
        raise ValueError('Require a nonempty finite one-dimensional energy series')
    return result


def bar(forward_work, reverse_work):
    """Dimensionless two-state Bennett acceptance ratio with actual sample counts.

Forward work is u_B-u_A sampled at A; reverse is u_A-u_B sampled at B.
The normalized two-state overlap is (1/N_A+1/N_B)*sum p*(1-p), where
p is the pooled posterior B probability at the BAR optimum. Range: 0..1.
"""
    forward, reverse = _array(forward_work), _array(reverse_work)
    nforward, nreverse = len(forward), len(reverse)
    ratio = math.log(nforward/nreverse)
    with np.errstate(over='raise', invalid='raise'):
        pooled = np.concatenate((forward, -reverse))
        center = float(np.median(pooled))
        left, right = forward-center, reverse+center
        def score(value):
            return float(logsumexp(log_expit(value-ratio-left))-
                         logsumexp(log_expit(ratio-value-right)))
        low = float(min(np.min(left), np.min(-right))-1.)
        high = float(max(np.max(left), np.max(-right))+1.)
        root = brentq(score, low, high, xtol=1e-12, rtol=1e-14, maxiter=200)
        z = root-ratio-(pooled-center)
        overlap = float((1/nforward+1/nreverse)*np.sum(expit(z)*expit(-z)))
        delta = root+center
    if not math.isfinite(delta) or not -1e-12 <= overlap <= 1+1e-12:
        raise ValueError('Invalid BAR result')
    return {'delta_f': delta, 'normalized_overlap': min(1., max(0., overlap)),
            'sample_counts': [nforward, nreverse], 'log_equation_residual': score(root),
            'overlap_resolved': overlap > 1e-12}


def correlation(values):
    """Conservative paired-positive, monotone autocorrelation-sum diagnostic.

Use the biased finite-series autocorrelation, truncate nonpositive successive
lag-pair sums and enforce nonincreasing positive sums. Clamp g to at least one.
This cannot detect an unsampled slow mode or prove stationarity.
    """
    data = _array(values)
    if np.all(data == data[0]):
        return {'g': 1., 'effective_samples': float(len(data)), 'constant_trace': True}
    centered = data-data.mean()
    if not np.isfinite(centered).all():
        raise ValueError('Energy centering overflow')
    scale = float(np.max(np.abs(centered)))
    if scale == 0:
        return {'g': 1., 'effective_samples': float(len(data)), 'constant_trace': True}
    centered = centered/scale
    nfft = 1 << (2*len(data)-1).bit_length()
    transform = np.fft.rfft(centered, n=nfft)
    acf = np.fft.irfft(transform*transform.conjugate(), n=nfft)[:len(data)]
    acf /= acf[0]
    previous, total = math.inf, 0.
    for lag in range(1, len(data)-1, 2):
        pair = float(acf[lag]+acf[lag+1])
        if pair <= 0:
            break
        previous = min(previous, pair)
        total += previous
    g = max(1., 1+2*total)
    return {'g': g, 'effective_samples': len(data)/g, 'constant_trace': False}


def block_indices(length, block_length, rng):
    """Circular moving blocks, retaining the original number of observations."""
    starts = rng.integers(0, length, size=math.ceil(length/block_length))
    return ((starts[:, None]+np.arange(block_length)) % length).ravel()[:length]


def ladder(gaps, weights, *, beta, block_length=None, bootstrap=1000, seed=112):
    """Analyze full endpoint ladders; canonical result is always 0 to signed charge.

Resample each window once per bootstrap draw, reusing it for both neighboring
BAR terms. This preserves their shared-frame covariance. Residual dependence
across chained windows and between replicas is NOT modeled by this interval.
"""
    data = [_array(gap) for gap in gaps]
    weights = _array(weights)
    if len(data) != len(weights) or len(data) < 2 or not math.isfinite(beta) or beta <= 0:
        raise ValueError('Require matching windows/weights and positive beta')
    direction = 1 if (weights[0], weights[-1]) == (0., 1.) else -1
    if ((weights[0], weights[-1]) not in ((0., 1.), (1., 0.)) or
            np.any(direction*np.diff(weights) <= 0) or np.any((weights < 0) | (weights > 1))):
        raise ValueError('Require a full monotonic ladder with exact endpoints')
    if block_length is not None and (type(block_length) is not int or block_length < 1):
        raise ValueError('Block length must be a positive integer')
    if type(bootstrap) is not int or bootstrap < 50 or type(seed) is not int or seed < 0:
        raise ValueError('Require at least 50 bootstrap draws and a nonnegative integer seed')
    def estimate(arrays):
        return [bar(beta*(b-a)*left, -beta*(b-a)*right)
                for a, b, left, right in zip(weights, weights[1:], arrays, arrays[1:])]
    pairs = estimate(data)
    metrics, failures, lengths = [], [], []
    for index, values in enumerate(data):
        metric = correlation(values)
        length = block_length if block_length is not None else max(1, math.ceil(5*metric['g']))
        lengths.append(length)
        metrics.append({**metric, 'frames': len(values), 'block_length': length,
                        'full_blocks': len(values)//length})
        if metric['effective_samples'] < 100 or len(values)//length < 20:
            failures.append(f'window {index}: insufficient effective samples/blocks')
        if length < 5*metric['g']:
            failures.append(f'window {index}: block shorter than five estimated correlation factors')
        if metric['constant_trace']:
            failures.append(f'window {index}: constant gap cannot diagnose mixing')
    for index, pair in enumerate(pairs):
        if pair['normalized_overlap'] < .03:
            failures.append(f'pair {index}: normalized overlap below 0.03')
    point = (direction*sum(pair['delta_f'] for pair in pairs)/beta
             if all(pair['overlap_resolved'] for pair in pairs) else None)
    draws, covariance, interval, uncertainty = None, None, None, None
    evaluated = 0
    if not failures:
        rng = np.random.default_rng(seed)
        estimates = []
        for _ in range(bootstrap):
            # One draw per window, not independent draws for its two BAR edges.
            resampled = [values[block_indices(len(values), length, rng)] for values, length in zip(data, lengths)]
            replicate = estimate(resampled)
            evaluated += 1
            if any(pair['normalized_overlap'] < .03 for pair in replicate):
                failures.append('bootstrap draw has inadequate overlap; no draws discarded or selectively retained')
                break
            estimates.append([direction*pair['delta_f']/beta for pair in replicate])
        if not failures:
            draws = np.asarray(estimates)
            covariance = np.atleast_2d(np.cov(draws, rowvar=False, ddof=1)).tolist()
            totals = draws.sum(axis=1)
            interval = np.quantile(totals, [.025, .975]).tolist()
            uncertainty = float(totals.std(ddof=1))
    return {'delta_g_0_to_sign': point, 'path_direction': direction,
            'pairs': pairs, 'windows': metrics, 'conditional_interval_95': interval,
            'conditional_standard_error': uncertainty, 'edge_covariance': covariance,
            'bootstrap_draws_evaluated': evaluated,
            'bootstrap_draws_used_for_interval': 0 if draws is None else len(draws),
            'bootstrap_draws_requested': bootstrap, 'bootstrap_seed': seed,
            'gap_statistical_gates_passed': not failures, 'failures': failures,
            'production_ready': False,
            'limitations': ['gap-based diagnostics only, not equilibration or stationarity proof',
                            'conditional within-window block uncertainty, not between-replica uncertainty',
                            'residual dependence across chained windows is not included',
                            'existing MD/thermostat equilibrium assumptions remain unverified']}
