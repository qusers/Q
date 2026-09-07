"""Read-only charge-only chain analysis with explicit Born and uncertainty views."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import struct

import numpy as np
import scipy

from . import charge_bar, charge_chain, charge_completion, charge_protocol as cp


def native_boltzmann(build_report):
    """Supported build dialect declares boltz as default real(4), not real(8)."""
    source = (build_report.parent/'src/q6/md.f90').read_text()
    pattern = r'^\s*real\s*,\s*parameter\s*::\s*boltz\s*=\s*([0-9.eEdD+\-]+)\s*(?:!.*)?$'
    values = re.findall(pattern, source, flags=re.MULTILINE | re.IGNORECASE)
    if len(values) != 1:
        raise ValueError('Unsupported/ambiguous native Boltzmann declaration; do not guess beta')
    value = struct.unpack('<f', struct.pack('<f', float(cp._number(values[0]))))[0]
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Invalid native Boltzmann constant')
    return value


def analyze_chain(path, *, discard_frames, block_length=None, bootstrap=1000, seed=112):
    if type(discard_frames) is not int or discard_frames < 0:
        raise ValueError('Declare nonnegative integer discarded frames per window')
    analysis_sources = {name: cp.fingerprint(Path(__file__).with_name(name))
                        for name in ('charge_analysis.py', 'charge_bar.py')}
    plan = charge_chain.inspect_plan(path)
    completed, _ = charge_chain.progress(plan)
    if completed != len(plan['series']['windows']):
        raise ValueError('Analyze only a completely verified chain')
    series, gaps, weights, provenance, constants = plan['series'], [], [], [], None
    for index, definition in enumerate(series['windows']):
        directory = Path(definition['input']).parent
        receipt_path = directory/'charge-completed.json'
        receipt = json.loads(receipt_path.read_text())
        native = receipt['result']['native_initialization']
        candidate_constants = native['born_state_constants_kcal_mol']
        if constants is None:
            constants = candidate_constants
        elif constants != candidate_constants:
            raise ValueError('Born constants differ across windows')
        window_weights = list(map(float, definition['lambdas']))
        energy_path = Path(definition['paths']['energy'])
        values = np.array([frame[1][1]-frame[0][1]
                           for frame in charge_completion.frames(energy_path, window_weights)])
        if discard_frames >= len(values):
            raise ValueError(f'Discard removes all frames in window {index}')
        values = values[discard_frames:]
        # Work in the raw finite-sphere view first. A force-free state constant
        # does not require new trajectories. Restore it exactly once below.
        born_gap = constants[1]-constants[0]
        if series['born_mode'] == 'integrated':
            values = values-born_gap
        gaps.append(values)
        weights.append(window_weights[1])
        provenance.append({'input': definition['input'], 'input_sha256': definition['input_sha256'],
                           'energy_sha256': cp.fingerprint(energy_path),
                           'completion_sha256': cp.fingerprint(receipt_path),
                           'retained_frames': len(values), 'discarded_frames': discard_frames})
    temperature = float(series['windows'][0]['signature']['md']['temperature'])
    boltzmann = native_boltzmann(Path(plan['build_report']))
    beta = 1/(temperature*boltzmann)
    result = charge_bar.ladder(gaps, weights, beta=beta, block_length=block_length, bootstrap=bootstrap, seed=seed)
    raw = result['delta_g_0_to_sign']
    born_gap = constants[1]-constants[0]
    corrected = None if raw is None else raw+born_gap
    raw_interval = result['conditional_interval_95']
    corrected_interval = None if raw_interval is None else [value+born_gap for value in raw_interval]
    declared_corrected = series['born_mode'] != 'control'
    # Recheck file and driver identities after the statistical calculation too.
    if charge_chain.inspect_plan(path) != plan:
        raise ValueError('Plan/build/input changed during analysis')
    for item, definition in zip(provenance, series['windows']):
        directory = Path(definition['input']).parent
        if (cp.fingerprint(Path(definition['paths']['energy'])) != item['energy_sha256'] or
                cp.fingerprint(directory/'charge-completed.json') != item['completion_sha256']):
            raise ValueError('Saved chain data changed during analysis')
    if any(cp.fingerprint(Path(__file__).with_name(name)) != checksum for name, checksum in analysis_sources.items()):
        raise ValueError('Analysis source changed during calculation')
    return {'schema_version': 1, 'gate': 'audited_chain_analysis', 'production_ready': False,
            'series': {key: value for key, value in series.items() if key != 'windows'},
            'plan_sha256': plan['plan_sha256'], 'build_report_sha256': plan['build_report_sha256'],
            'temperature_kelvin': temperature, 'native_boltzmann_kcal_mol_kelvin': boltzmann, 'beta_mol_per_kcal': beta,
            'windows': provenance, 'raw_delta_g_kcal_mol': raw, 'with_born_delta_g_kcal_mol': corrected,
            'born_delta_0_to_sign_kcal_mol': born_gap,
            'declared_result_kcal_mol': corrected if declared_corrected else raw,
            'control_corrected_view_is_comparison_only': not declared_corrected,
            'raw_conditional_interval_95_kcal_mol': raw_interval,
            'with_born_conditional_interval_95_kcal_mol': corrected_interval,
            'analysis': result, 'runtime_versions': {'numpy': np.__version__, 'scipy': scipy.__version__},
            'estimate_status': 'conditional_gap_statistics_only' if result['gap_statistical_gates_passed'] else 'insufficient_sampling',
            'analysis_sources_sha256': analysis_sources,
            'limitations': ['one ladder only; no between-replica, direction or radius conclusion',
                            'explicit discard is a choice, not an equilibration proof',
                            'raw and Born views share data and uncertainty; they are not independent experiments',
                            'no inference about ghost endpoints or experimental binding accuracy']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--discard-frames', type=int, required=True)
    parser.add_argument('--block-length', type=int)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=112)
    args = parser.parse_args()
    path = args.plan
    del args.plan
    try:
        report = analyze_chain(path, **vars(args))
    except (OSError, ValueError, KeyError, TypeError, OverflowError, FloatingPointError, RuntimeError) as error:
        parser.exit(2, f'Charge analysis failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
