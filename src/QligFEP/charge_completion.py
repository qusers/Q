"""Read-only completion checks for one staged charge-only Q MD window.

Normal exit and consistent outputs are prerequisites for restart chaining, not
proof of equilibration, full-trajectory interaction coverage or physical validity.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct

from . import boundary_native as bn, charge_protocol as cp


def _record(stream):
    marker = stream.read(4)
    if not marker:
        return None
    if len(marker) != 4:
        raise ValueError('Truncated energy record marker')
    size = struct.unpack('<i', marker)[0]
    if size not in (0, 124):
        raise ValueError('Unsupported two-state energy record size')
    payload = stream.read(size)
    if len(payload) != size or stream.read(4) != marker:
        raise ValueError('Damaged energy record')
    return payload


def frames(path, weights):
    """Stream Q's 2*(integer + 15 doubles) plus empty off-diagonal records.

The double fields are lambda, total, four bonded terms, qx(el,vdW),
qq(el,vdW), qp(el,vdW), qw(el,vdW), restraint. qx is qq+qp+qw,
not an additional interaction. Integrated Born is in total, not restraint.
"""
    if (len(weights) != 2 or not all(math.isfinite(w) and 0 <= w <= 1 for w in weights)
            or not math.isclose(sum(weights), 1., abs_tol=1e-12, rel_tol=0)):
        raise ValueError('Invalid expected two-state lambda mapping')
    with path.open('rb') as stream:
        while (first := _record(stream)) is not None:
            second = _record(stream)
            states = []
            for state, payload in enumerate((first, second), 1):
                if payload is None or len(payload) != 124 or struct.unpack('<i', payload[:4])[0] != state:
                    raise ValueError('Missing or incorrect saved state mapping')
                values = struct.unpack('<15d', payload[4:])
                if not all(map(math.isfinite, values)):
                    raise ValueError('Nonfinite saved energy value')
                if not math.isclose(values[0], weights[state-1], abs_tol=1e-12, rel_tol=0):
                    raise ValueError('Incorrect saved lambda')
                states.append(values)
            if _record(stream) != b'':
                raise ValueError('Missing or nonempty off-diagonal energy record')
            yield tuple(states)


def validate_window(window, born_mode, log_path):
    native = bn.validate_window(window, born_mode, log_path)
    log = log_path.read_text()
    if log.count('terminated normally.') != 1 or 'terminated abnormally' in log:
        raise ValueError('Require one normally completed native MD run')
    if log.count('FINAL  Energy summary') != 1:
        raise ValueError('Missing or duplicate final energy summary')
    settings = window['signature']
    expected = (int(settings['md']['steps'])-1)//int(settings['intervals']['energy'])
    energy_path, final_path = (Path(window['paths'][key]) for key in ('energy', 'final'))
    constants = [row[-1] for row in native['native']['state']]
    count, max_residual = 0, 0.
    minimum, maximum = math.inf, -math.inf
    for states in frames(energy_path, list(map(float, window['lambdas']))):
        count += 1
        for values, constant in zip(states, constants):
            for got, wanted, name in (
                (values[6], math.fsum(values[i] for i in (8, 10, 12)), 'electrostatic subtotal'),
                (values[7], math.fsum(values[i] for i in (9, 11, 13)), 'LJ subtotal'),
                (values[1], math.fsum((*values[2:8], values[14], constant)), 'pure-state total including applied Born'),
            ):
                if not math.isclose(got, wanted, abs_tol=1e-8, rel_tol=1e-12):
                    raise ValueError(f'Saved {name} is inconsistent at frame {count}')
                max_residual = max(max_residual, abs(got-wanted))
        gap = states[1][1]-states[0][1]
        if not math.isfinite(gap):
            raise ValueError('Nonfinite pure-state energy difference')
        minimum, maximum = min(minimum, gap), max(maximum, gap)
    if not count or count != expected:
        raise ValueError(f'Incomplete/extra energy frames: got {count}, expected {expected}')
    final = cp.restart_offsets(final_path)
    initial = window['restart']
    if final['atoms'] != initial['atoms'] or final['offset_record_sha256'] != initial['offset_record_sha256']:
        raise ValueError('Final restart dimensions or frozen offsets differ from input')
    return {'schema_version': 1, 'gate': 'completed_window_consistency_passed', 'production_ready': False,
            'native_initialization': native, 'saved_frames': count,
            'energy_sha256': cp.fingerprint(energy_path), 'final_sha256': cp.fingerprint(final_path),
            'final_restart': final, 'maximum_accounting_residual_kcal_mol': max_residual,
            'energy_gap_range_kcal_mol': [minimum, maximum],
            'limitations': ['output consistency, not source-build or launch attribution',
                            'normal exit and finite outputs do not establish trajectory stability or equilibration',
                            'no whole-trajectory cutoff, solvent-geometry or temperature audit',
                            'uncertainty, overlap, convergence and physical correction remain unqualified']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('preflight', type=Path)
    parser.add_argument('--series', required=True)
    parser.add_argument('--window', type=int, required=True)
    parser.add_argument('--log', type=Path, required=True)
    args = parser.parse_args()
    try:
        window, mode = bn.load_window(args.preflight, args.series, args.window)
        report = validate_window(window, mode, args.log)
        report['preflight_sha256'] = cp.fingerprint(args.preflight)
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Charge completion audit failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
