"""Assemble the declared two-radius charge pilot; write plans, never launch MD."""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil

from . import charge_chain as chain, charge_diagnostics as diagnostics
from . import charge_probe as probe, charge_protocol as cp

# Unmodified repository force-field assets at the charged-boundary baseline
# 55aa555c, checked before pilot outcomes. Not a parameter fit or free-energy target.
FORCE_FIELD_HASHES = {
    'water.lib': 'bcce8b17a6159e7a2fce6e04e706ece9a9b4d838ad44086bc42f16169220420d',
    'water.prm': '32868461147cd815c274bf104981886f49e324baf703629b1129163139896d0b',
}


def protocol(profile, timestep):
    if timestep not in (.5, 1.):
        raise ValueError('Require a declared 0.5/1 fs timestep')
    if profile == 'feasibility_pilot':
        preparation, window, discard, weights = 100., 30., 10., [i/10 for i in range(11)]
    elif profile == 'software_smoke':
        preparation, window, discard, weights = .1, .1, 0., [0., .5, 1.]
    else:
        raise ValueError('Unsupported campaign profile')
    return {'preparation_ps_per_ladder': preparation, 'window_ps': window,
            'discard_ps_per_window': discard, 'canonical_weights': weights,
            'ladders': 16, 'total_ps': 16*(preparation+len(weights)*window),
            'total_steps': int(Decimal(str(16*(preparation+len(weights)*window)))*1000/Decimal(str(timestep))),
            'physical_qualification_possible_from_profile_alone': False}


def topology_identity(path):
    """Restricted Qprep probe: compare all topology content except date/coordinates.

The one solute coordinate must be zero; all remaining atoms must form waters.
This exception is not appropriate for arbitrary protein reference restraints.
"""
    lines = path.read_text().splitlines()
    markers = [i for i, line in enumerate(lines) if '= Total no. of atoms, no. of solute atoms.' in line]
    codes = [i for i, line in enumerate(lines) if '= No. of integer atom codes.' in line]
    if len(markers) != 1 or len(codes) != 1 or codes[0] <= markers[0]:
        raise ValueError('Unsupported probe topology coordinate layout')
    start, stop = markers[0], codes[0]
    atoms, solute = map(int, lines[start].split('=')[0].split())
    coordinates = [float(v) for line in lines[start+1:stop] for v in line.split()]
    if (solute != 1 or atoms < 4 or (atoms-1) % 3 or len(coordinates) != 3*atoms or
            coordinates[:3] != [0., 0., 0.] or not all(map(math.isfinite, coordinates))):
        raise ValueError('Require one centered synthetic probe and complete finite water coordinates')
    retained = [line for line in lines[:start+1]+lines[stop:] if not line.startswith('DATE ')]
    normalized = '\n'.join(' '.join(line.split()) for line in retained)
    radii = [line for line in lines if '= Exclusion, solvent radii' in line]
    if len(radii) != 1:
        raise ValueError('Require the native topology exclusion/solvent radii')
    exclusion, effective = map(float, radii[0].split('=')[0].split())
    return {'parameters_sha256': hashlib.sha256(normalized.encode()).hexdigest(),
            'coordinates_sha256': hashlib.sha256(json.dumps(coordinates).encode()).hexdigest(),
            'atoms': atoms, 'exclusion_radius_angstrom': exclusion, 'effective_radius_angstrom': effective}


def _endpoint(path, checksum):
    if cp.fingerprint(path) != checksum:
        raise ValueError('Endpoint plan hash mismatch')
    # Reject an origin pointing at another ladder, before any recursive inspection.
    raw = json.loads(path.read_text())
    if raw.get('schema_version') != 2 or raw.get('purpose') != 'endpoint_preparation':
        raise ValueError('Require a fixed-endpoint preparation plan, not another ladder')
    report = chain.inspect_plan(path)
    if not report['native_trace_required']:
        raise ValueError('Campaign requires a build with all-evaluation native diagnostics')
    completed, final_hash = chain.progress(report)
    return report, completed, final_hash


def validate_transfer(path, spec, windows):
    origin = spec['endpoint_origin']
    if set(origin) != {'plan', 'sha256', 'final_sha256'}:
        raise ValueError('Unsupported endpoint-to-ladder origin declaration')
    endpoint_path = (path.parent/origin['plan']).resolve()
    report, completed, final_hash = _endpoint(endpoint_path, origin['sha256'])
    source = report['series']
    if completed != len(source['windows']):
        raise ValueError('Endpoint preparation must be completely verified before ladder transfer')
    last, first = source['windows'][-1], windows[0]
    if (final_hash != origin['final_sha256'] or final_hash != spec['initial_restart_sha256'] or
            (path.parent/spec['initial_restart']).resolve() != Path(last['paths']['final'])):
        raise ValueError('Ladder does not consume the exact verified endpoint final restart')
    engine = {**spec['engine'], 'binary': str((path.parent/spec['engine']['binary']).resolve())}
    if engine != report['engine'] or spec['build_report_sha256'] != report['build_report_sha256']:
        raise ValueError('Endpoint and ladder builds differ')
    for key in ('system', 'sign', 'direction', 'replica', 'born_mode', 'apply_born_posthoc'):
        if spec['series'][key] != source[key]:
            raise ValueError('Endpoint and ladder series identities differ')
    signature = json.loads(json.dumps(first['signature']))
    signature['md']['steps'] = last['signature']['md']['steps']
    if (signature != last['signature'] or first['states'] != last['states'] or
            first['assets_sha256'] != last['assets_sha256'] or first['lambdas'] != last['lambdas']):
        raise ValueError('Endpoint-to-ladder Hamiltonian, initial weight or velocity mode differs')
    return {'plan': str(endpoint_path), 'sha256': origin['sha256'], 'final_sha256': final_hash,
            'equilibrated': False, 'independent_replica_established': False,
            'preparation_total_ps': report['preparation_origin']['total_ps']}


def inspect(path):
    path = path.resolve()
    spec = json.loads(path.read_text())
    if set(spec) != {'schema_version', 'profile', 'timestep_fs', 'budget', 'endpoints'} or spec['schema_version'] != 1:
        raise ValueError('Unsupported campaign schema')
    budget = protocol(spec['profile'], spec['timestep_fs'])
    if spec['budget'] != budget or len(spec['endpoints']) != 16:
        raise ValueError('Campaign matrix or budget differs from its declared profile')
    expected = set(itertools.product((10., 14.), (-1, 1), ('forward', 'reverse'), (1, 2)))
    seen, seeds, coordinates, restarts, paths, identities, files, rows = set(), set(), set(), set(), set(), {}, None, []
    velocity_seeds = set()
    for entry in spec['endpoints']:
        if set(entry) != {'plan', 'sha256', 'ladder_directory'}:
            raise ValueError('Unsupported campaign endpoint entry')
        endpoint_path = (path.parent/entry['plan']).resolve()
        endpoint, completed, final_hash = _endpoint(endpoint_path, entry['sha256'])
        series, origin = endpoint['series'], endpoint['preparation_origin']
        key = (origin['grid_radius_angstrom'], series['sign'], series['direction'], series['replica'])
        if key not in expected or key in seen:
            raise ValueError('Missing/duplicate/unexpected radius-sign-direction-replica cell')
        seen.add(key)
        if (origin['total_ps'] != budget['preparation_ps_per_ladder'] or
                float(series['windows'][0]['signature']['md']['stepsize']) != spec['timestep_fs']):
            raise ValueError('Endpoint preparation duration/timestep differs from the campaign budget')
        label = f'r{int(key[0])}-q{key[1]}-{key[2]}-rep{key[3]}'
        if entry['ladder_directory'] != 'ladders/'+label:
            raise ValueError('Use the unique deterministic ladder destination for each campaign cell')
        prepared = Path(origin['prepared_report'])
        if any(cp.fingerprint(prepared.parent/name) != checksum for name, checksum in FORCE_FIELD_HASHES.items()):
            raise ValueError('Campaign requires the unchanged predeclared repository force-field assets')
        topology = topology_identity(prepared.parent/'system.top')
        orientation_lines = [line.split() for line in (prepared.parent/'prepare.inp').read_text().splitlines()
                             if line.startswith('set random_seed_solvent ')]
        if (orientation_lines != [['set', 'random_seed_solvent', str(origin['orientation_seed'])]] or
                topology['exclusion_radius_angstrom'] != key[0] or
                topology['effective_radius_angstrom'] != origin['effective_radius_angstrom']):
            raise ValueError('Campaign radius/orientation labels differ from the actual prepared assets')
        if (origin['orientation_seed'] in seeds or origin['velocity_seed'] in velocity_seeds or
                topology['coordinates_sha256'] in coordinates or endpoint['initial_restart_sha256'] in restarts or
                prepared in paths):
            raise ValueError('Reused seed, prepared coordinates or initial restart cannot represent another replica')
        seeds.add(origin['orientation_seed']); velocity_seeds.add(origin['velocity_seed'])
        coordinates.add(topology['coordinates_sha256']); restarts.add(endpoint['initial_restart_sha256']); paths.add(prepared)
        identity = (topology['parameters_sha256'], origin['effective_radius_angstrom'],
                    endpoint['initial_offsets']['offset_record_sha256'])
        if identities.setdefault(key[0], identity) != identity:
            raise ValueError('Same-radius probe Hamiltonian/topology or offsets differ across replicas')
        shared = (endpoint['build_report_sha256'], endpoint['engine'],
                  {name: cp.fingerprint(prepared.parent/name) for name in ('water.lib', 'water.prm', 'probe.lib', 'probe.pdb')})
        if files is not None and shared != files:
            raise ValueError('Campaign engine/build or force-field/probe assets differ')
        files = shared
        seed_path = Path(origin['seed_report'])
        seed = json.loads(seed_path.read_text())
        seed_window = {'signature': {'md': {'steps': '20'}, 'intervals': {'output': '10'}}}
        diagnostics.assess(seed_window, seed['native'], seed_path.parent/'native.log')
        ladder_path = path.parent/entry['ladder_directory']/'plan.json'
        ladder_completed = 0
        if ladder_path.exists():
            ladder = chain.inspect_plan(ladder_path)
            if ladder['endpoint_origin'] is None or ladder['endpoint_origin']['plan'] != str(endpoint_path):
                raise ValueError('Campaign ladder belongs to a different endpoint origin')
            weights = budget['canonical_weights'][::1 if series['direction'] == 'forward' else -1]
            if (list(map(float, (w['lambdas'][1] for w in ladder['series']['windows']))) != weights or
                    any(int(w['signature']['md']['steps']) != int(budget['window_ps']*1000/spec['timestep_fs'])
                        for w in ladder['series']['windows'])):
                raise ValueError('Generated ladder differs from the campaign profile')
            ladder_completed, _ = chain.progress(ladder)
        elif ladder_path.parent.exists():
            raise ValueError('Incomplete existing ladder directory; preserve it, do not overwrite or silently resume generation')
        rows.append({'cell': list(key), 'endpoint_plan': str(endpoint_path), 'endpoint_sha256': entry['sha256'],
                     'completed_preparation_segments': completed, 'preparation_segments': len(series['windows']),
                     'ladder_path': str(ladder_path), 'completed_ladder_windows': ladder_completed,
                     'orientation_seed': origin['orientation_seed'], 'velocity_seed': origin['velocity_seed'],
                     'topology_identity': topology, 'effective_radius_angstrom': origin['effective_radius_angstrom']})
    if seen != expected:
        raise ValueError('Incomplete campaign matrix')
    return {'schema_version': 1, 'gate': 'campaign_plan_consistency_passed', 'production_ready': False,
            'profile': spec['profile'], 'timestep_fs': spec['timestep_fs'], 'budget': budget,
            'campaign_path': str(path), 'campaign_sha256': cp.fingerprint(path), 'cells': rows,
            'equilibrated_endpoints_established': False, 'independent_sampling_established': False,
            'limitations': ['distinct seeds/initial files are necessary controls, not equilibrium independence proof',
                            'prospective MD step budget, not an HPC allocation or permission to launch',
                            'pilot cannot qualify the general physical correction',
                            'between-replica/direction/radius statistical analysis remains required']}


def assemble(endpoint_plans, directory, *, profile='feasibility_pilot', timestep_fs=1.):
    budget = protocol(profile, timestep_fs)
    if len(endpoint_plans) != 16:
        raise ValueError('Require the full sixteen-cell endpoint matrix')
    entries = []
    for endpoint_path in endpoint_plans:
        endpoint_path = endpoint_path.resolve(strict=True)
        report, _, _ = _endpoint(endpoint_path, cp.fingerprint(endpoint_path))
        series, radius = report['series'], report['preparation_origin']['grid_radius_angstrom']
        label = f'r{int(radius)}-q{series["sign"]}-{series["direction"]}-rep{series["replica"]}'
        entries.append({'plan': str(endpoint_path), 'sha256': cp.fingerprint(endpoint_path), 'ladder_directory': 'ladders/'+label})
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    path = directory/'campaign.json'
    probe._json(path, {'schema_version': 1, 'profile': profile, 'timestep_fs': timestep_fs,
                       'budget': budget, 'endpoints': entries})
    return inspect(path)


def stage_ladder(campaign_path, index):
    campaign = inspect(campaign_path)
    if type(index) is not int or not 0 <= index < len(campaign['cells']):
        raise ValueError('Invalid campaign cell index')
    cell = campaign['cells'][index]
    endpoint, completed, final_hash = _endpoint(Path(cell['endpoint_plan']), cell['endpoint_sha256'])
    if completed != len(endpoint['series']['windows']):
        raise ValueError('Complete and verify endpoint preparation before staging its ladder')
    last = endpoint['series']['windows'][-1]
    directory = Path(cell['ladder_path']).parent
    directory.mkdir(parents=True, exist_ok=False)
    weights = campaign['budget']['canonical_weights'][::1 if endpoint['series']['direction'] == 'forward' else -1]
    steps = int(campaign['budget']['window_ps']*1000/campaign['timestep_fs'])
    windows = []
    for i, weight in enumerate(weights):
        run = directory/f'w{i:03d}'
        run.mkdir()
        for key, name in [('topology', 'system.top'), ('fep', 'charge.fep')]:
            shutil.copyfile(last['paths'][key], run/name)
        restart = os.path.relpath(last['paths']['final'], run) if i == 0 else f'../w{i-1:03d}/final.re'
        text = probe.md_input(cell['effective_radius_angstrom'], steps, campaign['timestep_fs'], 0, weight)
        (run/'run.inp').write_text(text.replace('[files]\n', f'[files]\nrestart {restart}\n'))
        windows.append({'input': f'w{i:03d}/run.inp', 'sha256': cp.fingerprint(run/'run.inp'),
                        'assets_sha256': {'topology': cp.fingerprint(run/'system.top'), 'fep': cp.fingerprint(run/'charge.fep')}})
    series = {**endpoint['series'], 'id': endpoint['series']['id'].replace('endpoint-', 'ladder-', 1), 'windows': windows}
    spec = {'schema_version': 3, 'purpose': 'charge_ladder', 'engine': endpoint['engine'],
            'build_report': endpoint['build_report'], 'build_report_sha256': endpoint['build_report_sha256'],
            'initial_restart': last['paths']['final'], 'initial_restart_sha256': final_hash, 'series': series,
            'endpoint_origin': {'plan': cell['endpoint_plan'], 'sha256': cell['endpoint_sha256'], 'final_sha256': final_hash}}
    path = directory/'plan.json'
    probe._json(path, spec)
    return chain.inspect_plan(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    assembly = commands.add_parser('assemble')
    assembly.add_argument('directory', type=Path)
    assembly.add_argument('endpoint_plans', type=Path, nargs='+')
    assembly.add_argument('--profile', default='feasibility_pilot', choices=('feasibility_pilot', 'software_smoke'))
    assembly.add_argument('--timestep-fs', type=float, default=1.)
    for action in ('inspect', 'stage-ladder'):
        command = commands.add_parser(action)
        command.add_argument('campaign_path', type=Path)
        if action == 'stage-ladder':
            command.add_argument('--index', type=int, required=True)
    args = vars(parser.parse_args())
    action = args.pop('action')
    try:
        report = assemble(**args) if action == 'assemble' else stage_ladder(**args) if action == 'stage-ladder' else inspect(args['campaign_path'])
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Campaign {action} failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
