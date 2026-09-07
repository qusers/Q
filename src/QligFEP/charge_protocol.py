"""Read-only preflight for staged, charge-only Q validation inputs.

This deliberately accepts a small explicit input dialect, not every Q feature.
Passing is an input-consistency gate, never physical or production qualification.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import struct

from .endpoint_trim import _number, _tokens


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sections(path):
    result, active = {}, None
    for line in path.read_text().splitlines():
        tokens = _tokens(line)
        if not tokens:
            continue
        if tokens[0].startswith('['):
            if len(tokens) != 1 or not tokens[0].endswith(']'):
                raise ValueError(f'Malformed section in {path}')
            active = tokens[0][1:-1].lower()
            if active in result:
                raise ValueError(f'Duplicate section {active} in {path}')
            result[active] = []
        elif active is None:
            raise ValueError(f'Content before first section in {path}')
        else:
            result[active].append(tokens)
    return result


def keyed(rows):
    result = {}
    for row in rows:
        key = row[0].lower()
        if key in result or len(row) != 2:
            raise ValueError(f'Duplicate or non-scalar key: {key}')
        result[key] = row[1]
    return result


def restart_offsets(path):
    """Supported dialect: little-endian, 4-byte markers/integers, 8-byte x/v."""
    blocks = []
    with path.open('rb') as stream:
        while marker := stream.read(4):
            if len(marker) != 4:
                raise ValueError('Truncated restart record marker')
            size = struct.unpack('<i', marker)[0]
            if size < 0 or size > path.stat().st_size:
                raise ValueError('Unsupported restart record format')
            block = stream.read(size)
            if len(block) != size or stream.read(4) != marker:
                raise ValueError('Damaged restart record')
            blocks.append(block)
    if len(blocks) != 3 or any(len(b) < 4 for b in blocks):
        raise ValueError('Expected spherical coordinate, velocity and offset records')
    nat3 = struct.unpack('<i', blocks[0][:4])[0]
    for block in blocks[:2]:
        if nat3 <= 0 or nat3 % 3 or len(block) != 4+8*nat3 or struct.unpack('<i', block[:4])[0] != nat3:
            raise ValueError('Invalid restart coordinate/velocity dimensions')
        if not all(math.isfinite(v[0]) for v in struct.iter_unpack('<d', block[4:])):
            raise ValueError('Nonfinite restart coordinates/velocities')
    count = struct.unpack('<i', blocks[2][:4])[0]
    if count <= 0 or len(blocks[2]) != 4+4*count:
        raise ValueError('Invalid restart offset dimensions')
    offsets = list(struct.unpack(f'<{count}f', blocks[2][4:]))
    if not all(map(math.isfinite, offsets)):
        raise ValueError('Nonfinite restart offsets')
    return {'atoms': nat3//3, 'offsets_radians': offsets,
            'offset_record_sha256': hashlib.sha256(blocks[2]).hexdigest()}


def charge_states(path, atom_count):
    fep = sections(path)
    if set(fep) != {'fep', 'atoms', 'change_charges'} or keyed(fep['fep']) != {'states': '2'}:
        raise ValueError('Charge-only FEP requires exactly states, atoms and change_charges; no softcore/type/bond changes')
    mapping, charges = {}, {}
    for row in fep['atoms']:
        if len(row) != 2:
            raise ValueError('Invalid Q-atom mapping')
        q, atom = map(int, row)
        if q in mapping or atom in mapping.values() or not 1 <= atom <= atom_count:
            raise ValueError('Duplicate or out-of-range Q-atom mapping')
        mapping[q] = atom
    if not mapping or sorted(mapping) != list(range(1, len(mapping)+1)):
        raise ValueError('Q-atom indices must be contiguous from 1')
    for row in fep['change_charges']:
        if len(row) != 3:
            raise ValueError('Expected exactly two charges per Q atom')
        q = int(row[0])
        if q in charges or q not in mapping:
            raise ValueError('Duplicate or unmapped charge entry')
        charges[q] = [_number(v) for v in row[1:]]
    if set(charges) != set(mapping):
        raise ValueError('Every Q atom needs two explicit charges')
    canonical = [(mapping[q], *charges[q]) for q in sorted(mapping)]
    totals = [sum((v[s] for v in charges.values()), Decimal(0)) for s in (0, 1)]
    return canonical, totals


def _required(values, expected):
    for key, value in expected.items():
        if values.get(key, '').lower() != value:
            raise ValueError(f'Require explicit {key} {value}')


def _positive(values, keys):
    for key in keys:
        if key not in values or _number(values[key]) <= 0:
            raise ValueError(f'Require positive explicit {key}')


def inspect_window(path, born_mode):
    raw = sections(path)
    required = {'md', 'cut-offs', 'sphere', 'solvent', 'intervals', 'files', 'lambdas'}
    if set(raw) != required:
        raise ValueError(f'Unsupported/missing MD sections: {set(raw)^required}')
    config = {k: keyed(v) for k, v in raw.items() if k != 'lambdas'}
    md, solvent, files = config['md'], config['solvent'], config['files']
    allowed = {
        'md': {'steps', 'stepsize', 'temperature', 'bath_coupling', 'random_seed', 'initial_temperature',
               'shake_solvent', 'shake_hydrogens', 'shake_solute', 'lrf'},
        'cut-offs': {'solute_solvent', 'solute_solute', 'solvent_solvent', 'q_atom', 'lrf'},
        'sphere': {'shell_force', 'shell_radius'},
        'solvent': {'radius', 'radial_force', 'polarization', 'polarization_force', 'charge_correction',
                    'perstate_polarization', 'polarization_adaptation', 'perstate_born_correction', 'born_dielectric'},
        'intervals': {'output', 'non_bond', 'energy'},
        'files': {'topology', 'fep', 'restart', 'final', 'energy'},
    }
    for section, entries in config.items():
        if set(entries)-allowed[section]:
            raise ValueError(f'Unsupported {section} keys: {set(entries)-allowed[section]}')
    _required(md, {'lrf': 'off', 'shake_solvent': 'on', 'shake_hydrogens': 'on', 'shake_solute': 'off'})
    _positive(md, ['steps', 'stepsize', 'temperature', 'bath_coupling', 'random_seed', 'initial_temperature'])
    _positive(config['sphere'], ['shell_force', 'shell_radius'])
    _positive(solvent, ['radius', 'radial_force', 'polarization_force'])
    _required(solvent, {'polarization': 'on', 'charge_correction': 'on', 'perstate_polarization': 'on',
                        'polarization_adaptation': 'off',
                        'perstate_born_correction': 'on' if born_mode == 'integrated' else 'off'})
    # The angular dielectric factor is hardcoded to epsilon 80 in this branch.
    if _number(solvent.get('born_dielectric', 'nan')) != 80:
        raise ValueError('Use explicit dielectric 80 to match the current angular target')
    _positive(config['cut-offs'], allowed['cut-offs'])
    _positive(config['intervals'], allowed['intervals'])
    for key in ['steps', 'random_seed']:
        if _number(md[key]) != int(_number(md[key])):
            raise ValueError(f'{key} must be an integer')
    for key, value in config['intervals'].items():
        if _number(value) != int(_number(value)):
            raise ValueError(f'{key} interval must be an integer')
    if int(config['intervals']['energy']) >= int(md['steps']):
        raise ValueError('Energy interval must retain production records')
    weights = [_number(v) for row in raw['lambdas'] for v in row]
    if len(weights) != 2 or sum(weights) != 1 or any(v < 0 or v > 1 for v in weights):
        raise ValueError('Require normalized two-state lambdas')
    if set(files) != allowed['files']:
        raise ValueError('Declare topology, fep, restart, final and energy explicitly')
    if any(len(v.encode()) > 80 for v in files.values()):
        raise ValueError('Q filename exceeds its 80-byte field')
    paths = {key: (path.parent/value).resolve() for key, value in files.items()}
    assets = {key: fingerprint(paths[key]) for key in ('topology', 'fep', 'restart')}
    restart = restart_offsets(paths['restart'])
    states, totals = charge_states(paths['fep'], restart['atoms'])
    # Paths, mixing weights and random seeds are not potential parameters.
    signature = {k: dict(v) for k, v in config.items() if k != 'files'}
    signature['md'].pop('random_seed')
    signature['solvent'].pop('perstate_born_correction')
    return {'input': str(path), 'input_sha256': fingerprint(path), 'assets_sha256': assets,
            'paths': {k: str(v) for k, v in paths.items()}, 'restart': restart,
            'lambdas': list(map(str, weights)), 'q_region_charges': list(map(str, totals)),
            'states': [[atom, str(q0), str(q1)] for atom, q0, q1 in states],
            'signature': signature, 'random_seed': int(md['random_seed'])}


def validate(manifest_path):
    spec = json.loads(manifest_path.read_text())
    if set(spec) != {'schema_version', 'engine', 'series'} or spec['schema_version'] != 1:
        raise ValueError('Unsupported manifest schema')
    engine = spec['engine']
    if set(engine) != {'binary', 'sha256', 'source_commit'}:
        raise ValueError('Declare engine binary, sha256 and source_commit')
    binary = (manifest_path.parent/engine['binary']).resolve()
    if fingerprint(binary) != engine['sha256']:
        raise ValueError('Engine hash mismatch')
    if len(engine['source_commit']) != 40 or any(c not in '0123456789abcdef' for c in engine['source_commit']):
        raise ValueError('Declare a full source commit hash')
    if not spec['series']:
        raise ValueError('No series declared')
    systems, ids, outputs, protected, results = {}, set(), set(), {str(binary), str(manifest_path.resolve())}, []
    for series in spec['series']:
        if set(series) != {'id', 'system', 'sign', 'direction', 'replica', 'born_mode', 'apply_born_posthoc', 'windows'}:
            raise ValueError('Unsupported series declaration')
        if not series['id'] or series['id'] in ids or not series['system']:
            raise ValueError('Series IDs must be nonempty and unique; system is required')
        ids.add(series['id'])
        mode = series['born_mode']
        if mode not in {'integrated', 'posthoc', 'control'}:
            raise ValueError('Unknown Born accounting mode')
        if type(series['apply_born_posthoc']) is not bool or series['apply_born_posthoc'] != (mode == 'posthoc'):
            raise ValueError('Born accounting would omit or double-count the declared correction')
        if type(series['sign']) is not int or series['sign'] not in (-1, 1):
            raise ValueError('Charge sign must be -1 or +1')
        if type(series['replica']) is not int or series['replica'] < 1:
            raise ValueError('Replica must be a positive integer')
        if series['direction'] not in ('forward', 'reverse') or len(series['windows']) < 2:
            raise ValueError('Declare forward/reverse direction and at least two windows')
        windows = []
        for item in series['windows']:
            if set(item) != {'input', 'sha256', 'assets_sha256'}:
                raise ValueError('Each window requires input path, sha256 and assets_sha256')
            path = (manifest_path.parent/item['input']).resolve()
            if fingerprint(path) != item['sha256']:
                raise ValueError('MD input hash mismatch')
            window = inspect_window(path, mode)
            if window['assets_sha256'] != item['assets_sha256']:
                raise ValueError('Topology/FEP/restart asset hash mismatch')
            windows.append(window)
            protected.update([str(path), *(window['paths'][k] for k in ('topology', 'fep', 'restart'))])
            for key in ('energy', 'final'):
                target = window['paths'][key]
                if target in outputs or Path(target).exists():
                    raise ValueError('Output collision or existing output; refuse overwrite')
                outputs.add(target)
        first = windows[0]
        # Same system means same boundary and common state-1 atom charges, even across signs.
        identity = [first['assets_sha256']['topology'], first['signature'],
                    first['restart']['offset_record_sha256'],
                    [[atom, q0] for atom, q0, q1 in first['states']]]
        if systems.setdefault(series['system'], identity) != identity:
            raise ValueError('System boundary, topology, Q mapping, neutral state or frozen offsets differ')
        for window in windows:
            if (window['signature'] != first['signature'] or window['states'] != first['states'] or
                    window['assets_sha256']['topology'] != first['assets_sha256']['topology'] or
                    window['restart']['offset_record_sha256'] != first['restart']['offset_record_sha256']):
                raise ValueError('Within-series Hamiltonian or simulation settings differ')
            q0, q1 = map(Decimal, window['q_region_charges'])
            if q0 != 0 or q1 != series['sign']:
                raise ValueError('Validation requires Q-region charge 0 to declared +/-1')
        weights = [Decimal(w['lambdas'][1]) for w in windows]
        expected_ends = (0, 1) if series['direction'] == 'forward' else (1, 0)
        orientation = 1 if series['direction'] == 'forward' else -1
        if (weights[0], weights[-1]) != expected_ends or any(orientation*(b-a) <= 0 for a,b in zip(weights, weights[1:])):
            raise ValueError('Require complete monotonic charge-only ladder including both endpoints')
        results.append({**series, 'windows': windows})
    if outputs & protected:
        raise ValueError('Output would overwrite an input asset')
    return {'schema_version': 1, 'gate': 'staged_input_consistency_passed', 'production_ready': False,
            'manifest_path': str(manifest_path.resolve()), 'manifest_sha256': fingerprint(manifest_path),
            'engine': {**engine, 'binary': str(binary)}, 'series': results,
            'unverified': ['binary-to-source build provenance', 'native effective radius and included/excluded charge',
                           'native nonzero unchanged LJ parameters and interaction coverage',
                           'both-sign/control/replica campaign completeness and independence',
                           'final restart offsets and saved energy mappings',
                           'physical charge/background convention, convergence, uncertainty and overlap',
                           'compute budget and HPC approval']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    args = parser.parse_args()
    try:
        report = validate(args.manifest.resolve())
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(2, f'Charge protocol preflight failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
