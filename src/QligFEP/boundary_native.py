"""Check native Q initialization records against one staged charge-only window."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct

from . import charge_protocol as cp


def parse(text):
    begin, end = 'Q_BOUNDARY_AUDIT_V1 BEGIN', 'Q_BOUNDARY_AUDIT_V1 END'
    lines = text.splitlines()
    if lines.count(begin) != 1 or lines.count(end) != 1 or lines.index(begin) >= lines.index(end):
        raise ValueError('Require one complete native boundary audit block')
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
               'solute_boundary': 3, 'cutoffs': 5, 'water': 2}
    if set(records) != set(singles) | {'state', 'qatom', 'shell'}:
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
    for key, count, size in [('state', nstates, 6), ('qatom', nqat, 11+nstates), ('shell', nshell, 4+nstates)]:
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
    radius, requested_radius, ke, eps, override, env, excluded_env, kpol, krad, depth, width, max_radius = audit['parameters']
    config = window['signature']
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('preflight', type=Path, help='Saved staged-input preflight JSON report')
    parser.add_argument('--series', required=True)
    parser.add_argument('--window', required=True, type=int, help='Zero-based index in the declared series')
    parser.add_argument('--log', required=True, type=Path)
    args = parser.parse_args()
    try:
        report = json.loads(args.preflight.read_text())
        if report['gate'] != 'staged_input_consistency_passed':
            raise ValueError('Expected a staged-input preflight report')
        if cp.fingerprint(Path(report['manifest_path'])) != report['manifest_sha256']:
            raise ValueError('Manifest changed after preflight')
        if cp.fingerprint(Path(report['engine']['binary'])) != report['engine']['sha256']:
            raise ValueError('Engine changed after preflight')
        selected = [s for s in report['series'] if s['id'] == args.series]
        if len(selected) != 1 or not 0 <= args.window < len(selected[0]['windows']):
            raise ValueError('No unique declared series/window')
        result = validate_window(selected[0]['windows'][args.window], selected[0]['born_mode'], args.log)
        result['preflight_sha256'] = cp.fingerprint(args.preflight)
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Native boundary audit failed: {error}\n')
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
