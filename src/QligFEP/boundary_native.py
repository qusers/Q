"""Check native Q initialization records against one staged charge-only window."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct

from . import charge_protocol as cp


def parse(text):
    begin, end = 'Q_BOUNDARY_AUDIT_V3 BEGIN', 'Q_BOUNDARY_AUDIT_V3 END'
    lines = text.splitlines()
    if lines.count(begin) != 1 or lines.count(end) != 1 or lines.index(begin) >= lines.index(end):
        raise ValueError('Require one complete native boundary audit V3 block')
    records = {}
    for line in lines[lines.index(begin)+1:lines.index(end)]:
        fields = line.split()
        if not fields or not fields[0].startswith('QBA_'):
            raise ValueError('Unexpected content inside native boundary audit')
        values = [float(v) for v in fields[1:]]
        if not all(map(math.isfinite, values)):
            raise ValueError('Nonfinite native boundary audit value')
        records.setdefault(fields[0][4:].lower(), []).append(values)
    singles = {'convention': 1, 'meta': 8, 'flags': 6, 'parameters': 12, 'center': 3,
               'solute_boundary': 3, 'cutoffs': 5, 'water': 2, 'water_compatibility': 2,
               'restraint_counts': 6}
    if set(records)-{'position'} != set(singles) | {'state', 'qatom', 'shell', 'water_atom'}:
        raise ValueError('Missing or unsupported native audit records')
    result = {}
    for key, size in singles.items():
        if len(records[key]) != 1 or len(records[key][0]) != size:
            raise ValueError(f'Invalid native {key} record')
        result[key] = records[key][0]
    if any(v != int(v) for v in result['meta']) or any(v not in (0, 1) for v in result['flags']):
        raise ValueError('Invalid native integer/flag values')
    natom, nsolute, nwater, nqat, nstates, nshell, _, _ = result['meta']
    if min(natom, nsolute, nwater, nqat, nshell) <= 0 or nstates != 2 or nsolute+3*nwater != natom:
        raise ValueError('Unsupported native system dimensions')
    counts = result['restraint_counts']
    if any(v < 0 or v != int(v) for v in counts) or counts[-1] not in (0, 1):
        raise ValueError('Invalid native restraint counts')
    positions = records.get('position', [])
    if len(positions) != counts[1] or any(len(row) != 9 for row in positions):
        raise ValueError('Invalid native position dimensions')
    if [row[0] for row in positions] != list(range(1, int(counts[1])+1)):
        raise ValueError('Duplicate/unordered native position indices')
    if any(row[1] != int(row[1]) or not 1 <= row[1] <= natom or
           row[8] != int(row[8]) or not 0 <= row[8] <= nstates for row in positions):
        raise ValueError('Invalid native position atom/state indices')
    result['position'] = positions
    for key, count, size in [('state', nstates, 6), ('qatom', nqat, 11+nstates),
                             ('shell', nshell, 4+nstates), ('water_atom', 3, 10)]:
        rows = records[key]
        if len(rows) != count or any(len(row) != size for row in rows):
            raise ValueError(f'Invalid native {key} dimensions')
        if [row[0] for row in rows] != list(range(1, int(count)+1)):
            raise ValueError(f'Duplicate/unordered native {key} indices')
        result[key] = rows
    return result


def _close(got, expected, label, *, tolerance=1e-12):
    if not math.isclose(got, expected, abs_tol=tolerance, rel_tol=1e-12):
        raise ValueError(f'Native {label} mismatch: {got} != {expected}')


def _single(value):
    return struct.unpack('=f', struct.pack('=f', float(value)))[0]


def validate_window(window, born_mode, log_path):
    """Initialization consistency only; neither a job-completion nor physics gate."""
    if born_mode not in {'integrated', 'posthoc', 'control'}:
        raise ValueError('Unknown Born accounting mode')
    if cp.fingerprint(Path(window['input'])) != window['input_sha256']:
        raise ValueError('MD input changed after preflight')
    for key in ('topology', 'fep', 'restart'):
        if cp.fingerprint(Path(window['paths'][key])) != window['assets_sha256'][key]:
            raise ValueError(f'{key} changed after preflight')
    if cp.inspect_window(Path(window['input']), born_mode) != window:
        raise ValueError('Staged window description differs from the actual inputs')
    audit = parse(log_path.read_text())
    if audit['convention'] != [1]:
        raise ValueError('Unsupported native angular charge convention')
    natom, nsolute, nwater, nqat, _, nshell, rule, solvent_type = audit['meta']
    if natom != window['restart']['atoms'] or nqat != len(window['states']):
        raise ValueError('Native atom dimensions disagree with staged assets')
    if audit['flags'] != [1, int(born_mode == 'integrated'), 0, 0, 0, 0]:
        raise ValueError('Native flags disagree with frozen direct charge-only protocol')
    uniform, zero_h_lj = audit['water_compatibility']
    if uniform != 1 or zero_h_lj not in (0, 1) or solvent_type not in (0, 1):
        raise ValueError('Unsupported or nonuniform native water model')
    if solvent_type == 0 and not zero_h_lj:
        raise ValueError('Optimized SPC-like water code requires zero hydrogen LJ coefficients')
    for site in audit['water_atom']:
        if site[1] <= 0 or site[1] != int(site[1]) or site[2] <= 0 or min(site[4:]) < 0:
            raise ValueError('Invalid native water type, mass or LJ coefficients')
    actual_zero_h = all(v == 0 for site in audit['water_atom'][1:] for v in site[4:])
    if actual_zero_h != bool(zero_h_lj):
        raise ValueError('Native hydrogen LJ compatibility flag disagrees with coefficients')
    for site in audit['water_atom'][1:]:
        _close(site[3], _single(-0.5*audit['water_atom'][0][3]), 'symmetric neutral water charges')
    radius, requested_radius, ke, eps, override, env, excluded_env, kpol, krad, depth, width, max_radius = audit['parameters']
    config = window['signature']
    expected_positions = config['atom_restraints']
    if audit['restraint_counts'] != [0, len(expected_positions), 0, 0, 0, 0]:
        raise ValueError('Native extra restraint counts disagree with supported shared positions')
    for row, expected in zip(audit['position'], expected_positions):
        for got, value in zip(row[1:], expected):
            _close(got, float(value), 'shared position restraint')
    _close(radius, float(config['solvent']['radius']), 'effective radius')
    _close(requested_radius, radius, 'requested/effective radius')
    _close(eps, 80, 'dielectric')
    _close(kpol, float(config['solvent']['polarization_force']), 'angular force constant')
    _close(krad, float(config['solvent']['radial_force']), 'radial force constant')
    inner, exclusion_radius, solute_force = audit['solute_boundary']
    declared_inner = float(config['sphere']['shell_radius'])
    _close(inner, declared_inner if declared_inner > 1 else declared_inner*exclusion_radius,
           'solute shell radius')
    _close(solute_force, float(config['sphere']['shell_force']), 'solute shell force constant')
    if not 0 < inner <= exclusion_radius:
        raise ValueError('Invalid native solute boundary radii')
    if min(radius, ke, kpol, krad, depth, width, *audit['water']) <= 0 or override > 0:
        raise ValueError('Invalid native constants or unsupported Born override')
    if max_radius < 0:
        raise ValueError('Invalid native coordinate extent')
    cutoff_keys = ('solute_solute', 'solute_solvent', 'solvent_solvent', 'q_atom')
    for key, value in zip(cutoff_keys, audit['cutoffs'][:4]):
        _close(value, float(config['cut-offs'][key]), key+' cutoff')
    # The native report uses zero as an explicit disabled sentinel. The engine
    # does not read this input key when LRF is off; its variable may be unset.
    _close(audit['cutoffs'][4], 0., 'disabled LRF cutoff')
    # Raw FEP charges are stored as real(4) by Q. Compare against that explicit
    # representation, not an arbitrary loose charge tolerance.
    for row, expected in zip(audit['qatom'], window['states']):
        _, atom, excluded, atom_type, mass, *parameters_and_charges = row
        if atom != expected[0] or atom != int(atom) or atom > nsolute or excluded != 0:
            raise ValueError('Native Q mapping differs or contains excluded/non-solute Q atoms')
        if atom_type <= 0 or atom_type != int(atom_type) or mass <= 0:
            raise ValueError('Invalid native Q atom type/mass')
        avdw, bvdw = parameters_and_charges[:3], parameters_and_charges[3:6]
        if avdw[0] <= 0 or bvdw[0] <= 0:
            raise ValueError('Charge-only validation requires real nonzero Q-atom LJ parameters')
        for got, charge in zip(parameters_and_charges[6:], expected[1:]):
            _close(got, _single(charge), 'Q-atom charge')
    coefficient = ke*(1-1/eps)/(2*radius)
    state_charges, state_constants = [], []
    for s, row in enumerate(audit['state']):
        _, weight, q_all, q_in, q_out, applied = row
        expected_q = sum(_single(atom[s+1]) for atom in window['states'])
        _close(weight, float(window['lambdas'][s]), 'lambda')
        _close(q_all, expected_q, 'Q-region charge')
        _close(q_in, q_all, 'included Q-region charge')
        _close(q_out, 0., 'excluded Q-region charge')
        constant = -coefficient*(env+q_all)**2
        _close(applied, constant if born_mode == 'integrated' else 0., 'applied Born constant')
        state_charges.append(env+q_all)
        state_constants.append(constant)
    offsets = window['restart']['offsets_radians']
    if len(offsets) != nshell:
        raise ValueError('Native shell count differs from staged restart')
    outer = radius
    for row, offset in zip(audit['shell'], offsets):
        _, rout, dr, actual_offset, *strengths = row
        # wshell radii, unlike the global effective radius, are stored real(4).
        _close(rout, _single(outer), 'shell outer radius')
        if not 0 < dr < rout:
            raise ValueError('Invalid native shell width')
        if actual_offset != offset:
            raise ValueError('Native frozen offset differs from staged restart')
        outer -= dr
    return {'gate': 'native_initialization_consistency_passed', 'production_ready': False,
            'input_sha256': window['input_sha256'], 'log_sha256': cp.fingerprint(log_path),
            'angular_charge_convention': 'legacy_q_region', 'native': audit,
            'effective_radius_angstrom': radius, 'included_non_q_charge': env,
            'excluded_non_q_topology_charge': excluded_env, 'sphere_state_charges': state_charges,
            'born_coefficient_kcal_mol_per_e2': coefficient, 'born_state_constants_kcal_mol': state_constants,
            'born_mode': born_mode,
            'initial_conservative_cutoff_bound_passed': min(audit['cutoffs'][:4]) > 2*max_radius,
            'limitations': ['initialization only, before initial bond-constraint projection',
                            'input hashes do not prove which executable produced a supplied log',
                            'initial coordinate bound does not prove whole-trajectory interaction coverage',
                            'unchanged real topology LJ verified only for the restricted charge-only FEP dialect',
                            'job completion, saved energies, final offsets, sampling and physical validity remain separate']}


def load_window(preflight, series_id, index):
    """Select an immutable staged window; file fingerprints are not build proof."""
    report = json.loads(preflight.read_text())
    if report['gate'] != 'staged_input_consistency_passed':
        raise ValueError('Expected a staged-input preflight report')
    if cp.fingerprint(Path(report['manifest_path'])) != report['manifest_sha256']:
        raise ValueError('Manifest changed after preflight')
    if cp.fingerprint(Path(report['engine']['binary'])) != report['engine']['sha256']:
        raise ValueError('Engine changed after preflight')
    manifest_path = Path(report['manifest_path'])
    manifest = json.loads(manifest_path.read_text())
    expected_engine = {**manifest['engine'],
                       'binary': str((manifest_path.parent/manifest['engine']['binary']).resolve())}
    if report['engine'] != expected_engine:
        raise ValueError('Preflight engine declaration differs from retained manifest')
    selected = [s for s in report['series'] if s['id'] == series_id]
    declared = [s for s in manifest['series'] if s['id'] == series_id]
    if len(selected) != 1 or len(declared) != 1 or not 0 <= index < len(selected[0]['windows']):
        raise ValueError('No unique declared series/window')
    series, spec = selected[0], declared[0]
    if ({k: v for k, v in series.items() if k != 'windows'} !=
            {k: v for k, v in spec.items() if k != 'windows'} or
            len(series['windows']) != len(spec['windows'])):
        raise ValueError('Preflight series declaration differs from retained manifest')
    window, item = series['windows'][index], spec['windows'][index]
    if (window['input'] != str((manifest_path.parent/item['input']).resolve()) or
            window['input_sha256'] != item['sha256'] or window['assets_sha256'] != item['assets_sha256']):
        raise ValueError('Preflight window declaration differs from retained manifest')
    return window, series['born_mode']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('preflight', type=Path, help='Saved staged-input preflight JSON report')
    parser.add_argument('--series', required=True)
    parser.add_argument('--window', required=True, type=int, help='Zero-based index in the declared series')
    parser.add_argument('--log', required=True, type=Path)
    args = parser.parse_args()
    try:
        window, mode = load_window(args.preflight, args.series, args.window)
        result = validate_window(window, mode, args.log)
        result['preflight_sha256'] = cp.fingerprint(args.preflight)
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Native boundary audit failed: {error}\n')
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
