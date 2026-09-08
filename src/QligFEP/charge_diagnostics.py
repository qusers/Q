"""Read-only native charge trajectory diagnostics, not equilibrium qualification."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from . import charge_protocol as cp


def trace(path, *, steps, interval, waters, shells):
    """Read Q's version-1 geometry/temperature records at force-evaluation times."""
    expected = list(range(0, steps, interval))+[steps]
    records, active = [], None
    with path.open() as stream:
        for line in stream:
            if 'WARNING: hot atom' in line:
                raise ValueError('Native hot-atom velocity reset occurred; temperature summaries cannot qualify this run')
            tokens = line.split()
            if not tokens or not (tokens[0].startswith('Q_CHARGE_TRACE') or tokens[0].startswith('QCT_')):
                continue
            label, fields = tokens[0], tokens[1:]
            if label == 'Q_CHARGE_TRACE_V1':
                if len(fields) != 16:
                    raise ValueError('Invalid charge trace field count')
                step, evaluations, invalid = map(int, fields[:3])
                values = list(map(float, fields[3:]))
                if (evaluations != step+1 or invalid != 0 or not all(map(math.isfinite, values)) or
                        any(value < 0 for value in values) or min(values[3:]) <= 0):
                    raise ValueError('Invalid/nonfinite geometry or temperature in native charge trace')
                active = {'step': step, 'evaluations': evaluations, 'values': values, 'shells': []}
                records.append(active)
            elif label == 'QCT_DENSITY':
                if active is None or len(fields) != 22 or 'density' in active:
                    raise ValueError('Missing/duplicate/invalid density trace')
                step, *counts = map(int, fields)
                if step != active['step'] or min(counts) < 0 or sum(counts) != waters:
                    raise ValueError('Density trace does not account for every water')
                active['density'] = counts
            elif label == 'QCT_SHELL':
                if active is None or len(fields) != 5:
                    raise ValueError('Invalid shell trace')
                step, index, count = map(int, fields[:3])
                moment, square = map(float, fields[3:])
                if (step != active['step'] or index != len(active['shells'])+1 or not 0 <= count <= waters or
                        not math.isfinite(moment) or not math.isfinite(square) or
                        abs(moment) > count+1e-9 or not 0 <= square <= count+1e-9 or
                        moment*moment > count*square+1e-8):
                    raise ValueError('Inconsistent shell population/orientation moments')
                active['shells'].append([count, moment, square])
            else:
                raise ValueError('Unsupported native charge trace version/record')
    if [row['step'] for row in records] != expected:
        raise ValueError('Missing/duplicate/out-of-order native charge trace snapshots')
    for index, row in enumerate(records):
        if 'density' not in row or len(row['shells']) != shells or sum(s[0] for s in row['shells']) > waters:
            raise ValueError('Incomplete native density/shell trace')
        values = row['values']
        if (values[1] > values[2] or values[2] > values[0] or values[3] > values[4] or values[5] > values[6] or
                not values[9] <= values[7] <= values[10] or not values[11] <= values[8] <= values[12]):
            raise ValueError('Invalid cumulative native trace bounds')
        if index:
            previous = records[index-1]['values']
            if (any(values[j] < previous[j] for j in (0, 2, 4, 6, 10, 12)) or
                    any(values[j] > previous[j] for j in (3, 5, 9, 11))):
                raise ValueError('Native trace cumulative extrema are not monotonic')
    return records


def _summary(values):
    data = np.asarray(values, dtype=float)
    half = max(1, len(data)//2)
    return {'minimum': float(data.min()), 'maximum': float(data.max()), 'mean': float(data.mean()),
            'early_mean': float(data[:half].mean()), 'late_mean': float(data[half:].mean())}


def assess(window, native, log_path):
    """Validate evaluation coverage; report slow observables without fitting gates."""
    settings = window['signature']
    steps, interval = int(settings['md']['steps']), int(settings['intervals']['output'])
    waters, shells = native['meta'][2], native['meta'][5]
    if waters != int(waters) or shells != int(shells) or min(waters, shells) < 1:
        raise ValueError('Require positive integer native water/shell counts')
    waters, shells = int(waters), int(shells)
    rows = trace(log_path, steps=steps, interval=interval, waters=waters, shells=shells)
    first, last = rows[0]['values'], rows[-1]['values']
    cutoff = min(native['cutoffs'][:4])
    if 2*last[0] >= cutoff:
        raise ValueError('All-evaluation radius bound does not establish untruncated pair coverage')
    # A predeclared geometry-drift alarm, not a water-model calibration. Absolute
    # bond lengths and topology identity remain separately inspectable.
    drift = max(first[3]-last[3], last[4]-first[4], first[5]-last[5], last[6]-first[6])
    if drift > .005:
        raise ValueError('Water distance extrema drift by more than 0.005 angstrom from initialized geometry')
    radius = native['parameters'][0]
    edges = np.linspace(0, radius, 21)
    volumes = 4*math.pi/3*np.diff(edges**3)
    densities = np.asarray([row['density'][:20] for row in rows])/volumes
    return {'schema_version': 1, 'gate': 'native_trace_consistency_passed', 'production_ready': False,
            'force_geometries_observed': rows[-1]['evaluations'], 'expected_force_geometries': steps+1,
            'snapshot_count': len(rows), 'snapshot_steps': [r['step'] for r in rows],
            'maximum_all_atom_radius_angstrom': last[0], 'maximum_q_atom_radius_angstrom': last[2],
            'all_evaluation_pair_bound_angstrom': 2*last[0], 'smallest_active_cutoff_angstrom': cutoff,
            'cutoff_margin_angstrom': cutoff-2*last[0],
            'oh_distance_range_angstrom': last[3:5], 'hh_distance_range_angstrom': last[5:7],
            'water_distance_extrema_drift_angstrom': drift,
            'total_temperature_all_evaluation_range_kelvin': last[9:11],
            'free_temperature_all_evaluation_range_kelvin': last[11:13],
            'total_temperature_snapshot_kelvin': _summary([r['values'][7] for r in rows]),
            'free_temperature_snapshot_kelvin': _summary([r['values'][8] for r in rows]),
            'q_radius_snapshot_angstrom': _summary([r['values'][1] for r in rows]),
            'radial_bin_edges_angstrom': edges.tolist(),
            'oxygen_density_per_angstrom_cubed': [_summary(densities[:, i]) for i in range(20)],
            'outside_effective_radius_water_count': _summary([r['density'][20] for r in rows]),
            'shells': [{'index': i+1, 'population': _summary([r['shells'][i][0] for r in rows]),
                        'radial_orientation_sum': _summary([r['shells'][i][1] for r in rows]),
                        'radial_orientation_square_sum': _summary([r['shells'][i][2] for r in rows])}
                       for i in range(shells)],
            'limitations': ['temperature, density and orientation stationarity not established',
                            'positive temperature is not thermostat ensemble validation',
                            'geometry drift check is relative to initialization, not a force-field identity check',
                            'shell moments and density are snapshots; no bound on unobserved slow modes',
                            'no energy clipping or boundary-parameter fitting']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--window', type=int, required=True)
    args = parser.parse_args()
    from . import charge_chain
    try:
        plan = charge_chain.inspect_plan(args.plan)
        completed, _ = charge_chain.progress(plan)
        if not 0 <= args.window < completed:
            raise ValueError('Require a completed, verified window index')
        window = plan['series']['windows'][args.window]
        receipt = json.loads((Path(window['input']).parent/'charge-completed.json').read_text())
        native = receipt['result']['native_initialization']['native']
        report = assess(window, native, Path(window['input']).parent/'charge-native.log')
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Charge diagnostics failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
