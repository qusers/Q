"""Plan bounded fixed-endpoint preparation using existing Q restart dynamics.

This command writes inputs only. It neither launches MD nor qualifies an
endpoint as equilibrated. A separately retained 20-step grid-start seed is
included in the preparation budget; it is never treated as production data.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import json
import os
from pathlib import Path
import shutil

from . import charge_build, charge_probe as probe, charge_protocol as cp


def _assets(report_path):
    report = json.loads(report_path.read_text())
    for name, checksum in report['assets_sha256'].items():
        if Path(name).name != name or cp.fingerprint(report_path.parent/name) != checksum:
            raise ValueError('Preparation/seed asset changed or has an invalid path')
    return report


def validate_origin(path, spec, build_path, build, windows):
    """Bind the grid seed, Qprep, Qdyn, endpoint identity and full step budget."""
    declaration = spec['preparation']
    if set(declaration) != {'prepared_report', 'prepared_sha256', 'seed_report', 'seed_sha256',
                            'total_steps', 'segment_steps'}:
        raise ValueError('Unsupported endpoint origin declaration')
    reports = []
    for label in ('prepared', 'seed'):
        report_path = (path.parent/declaration[label+'_report']).resolve()
        if cp.fingerprint(report_path) != declaration[label+'_sha256']:
            raise ValueError('Endpoint origin report hash mismatch')
        reports.append((report_path, _assets(report_path)))
    (prepared_path, prepared), (seed_path, seed) = reports
    qprep = (build_path.parent/build['binaries']['qprep']['path']).resolve()
    qdyn = (build_path.parent/build['binaries']['qdyn']['path']).resolve()
    if (prepared['gate'] != 'probe_preparation_only' or prepared['production_ready'] is not False or
            prepared['qprep'] != {'binary': str(qprep), 'sha256': cp.fingerprint(qprep)} or
            seed['gate'] != 'bounded_md_timing_only' or seed['production_ready'] is not False or
            seed['engine'] != {'binary': str(qdyn), 'sha256': cp.fingerprint(qdyn)} or
            seed['prepared_sha256'] != cp.fingerprint(prepared_path)):
        raise ValueError('Endpoint seed/preparation does not match the retained isolated build')
    series = spec['series']
    weight = 0. if series['direction'] == 'forward' else 1.
    timestep = seed['timestep_fs']
    if (seed['steps'] != 20 or timestep not in (.5, 1.) or seed['sign'] != series['sign'] or
            seed['state2_weight'] != weight or series['born_mode'] != 'integrated'):
        raise ValueError('Require a 20-step integrated-Born seed at the intended starting endpoint')
    inp = seed_path.parent/'run.inp'
    velocity_seed = int(cp.keyed(cp.sections(inp)['md'])['random_seed'])
    if not 1 <= velocity_seed < 100_000_000:
        raise ValueError('Grid seed requires a positive explicit velocity seed')
    if inp.read_text() != probe.md_input(prepared['effective_radius_angstrom'], 20, timestep, velocity_seed, weight):
        raise ValueError('Grid seed input differs from the declared probe protocol')
    expected_fep = prepared_path.parent/('positive.fep' if series['sign'] == 1 else 'negative.fep')
    for seed_file, original in [('system.top', prepared_path.parent/'system.top'), ('charge.fep', expected_fep)]:
        if cp.fingerprint(seed_path.parent/seed_file) != cp.fingerprint(original):
            raise ValueError('Grid seed changed its prepared topology or charge definition')
    initial = (path.parent/spec['initial_restart']).resolve()
    if initial != seed_path.parent/'final.re' or cp.fingerprint(initial) != seed['assets_sha256']['final.re']:
        raise ValueError('Endpoint initial restart is not the recorded grid-seed output')
    offsets = cp.restart_offsets(initial)
    if offsets['atoms'] != prepared['atoms'] or any(offsets['offsets_radians']):
        raise ValueError('Require the prepared atom count and zero frozen seed offsets')
    log = (seed_path.parent/'native.log').read_text()
    if (log.count('terminated normally.') != 1 or 'terminated abnormally' in log or
            probe.bn.parse(log) != seed['native']):
        raise ValueError('Grid seed native log disagrees with its retained report')
    probe._check_energies(seed_path.parent/'states.en', [1-weight, weight], 1)
    total, segment = declaration['total_steps'], declaration['segment_steps']
    if (type(total) is not int or type(segment) is not int or segment < 20 or segment % 10 or
            total <= 20 or total % 10 or total*timestep > 100000):
        raise ValueError('Invalid endpoint step budget; maximum total is 100 ps including grid seed')
    remaining, sizes = total-20, []
    while remaining:
        size = min(segment, remaining)
        if size < 20:
            raise ValueError('Preparation remainder needs at least 20 steps for saved energy checks')
        sizes.append(size)
        remaining -= size
    if [int(w['signature']['md']['steps']) for w in windows] != sizes:
        raise ValueError('Endpoint segments do not match the declared preparation budget')
    for window in windows:
        if (window['assets_sha256'] != {'topology': cp.fingerprint(prepared_path.parent/'system.top'),
                                        'fep': cp.fingerprint(expected_fep)} or
                float(window['signature']['md']['stepsize']) != timestep or
                float(window['signature']['solvent']['radius']) != prepared['effective_radius_angstrom']):
            raise ValueError('Endpoint continuation differs from its grid seed system/timestep')
        # Reconstruct the prospective protocol, not just a self-consistent set of
        # modified windows. No heating, restraint or solvent-parameter retuning.
        expected = probe.md_input(prepared['effective_radius_angstrom'],
                                  int(window['signature']['md']['steps']), timestep, 0, weight)
        actual_path = Path(window['input'])
        restart_name = cp.keyed(cp.sections(actual_path)['files'])['restart']
        expected = expected.replace('[files]\n', f'[files]\nrestart {restart_name}\n')
        if actual_path.read_text() != expected:
            raise ValueError('Endpoint input differs from the frozen continuation protocol')
    return {'prepared_report': str(prepared_path), 'seed_report': str(seed_path),
            'orientation_seed': prepared['orientation_seed'], 'velocity_seed': velocity_seed,
            'grid_seed_steps': 20, 'continuation_steps': total-20, 'total_steps': total,
            'total_ps': total*timestep/1000, 'segment_steps': sizes,
            'grid_radius_angstrom': prepared['grid_radius_angstrom'],
            'effective_radius_angstrom': prepared['effective_radius_angstrom'],
            'equilibrated': False, 'independent_replica_established': False}


def generate(prepared_report, seed_report, build_report, directory, *, sign, direction, replica,
             total_ps=100., segment_ps=20.):
    """Write a new endpoint schedule; all calculation remains a separate action."""
    if type(sign) is not int or sign not in (-1, 1) or direction not in ('forward', 'reverse'):
        raise ValueError('Require a sign and forward/reverse starting-endpoint identity')
    if type(replica) is not int or replica < 1:
        raise ValueError('Require a positive replica identifier')
    prepared_report, seed_report, build_report = (p.resolve(strict=True) for p in
                                                 (prepared_report, seed_report, build_report))
    prepared, seed = _assets(prepared_report), _assets(seed_report)
    build = charge_build.validate(build_report)
    timestep = seed['timestep_fs']
    if timestep not in (.5, 1.):
        raise ValueError('Require a supported 0.5/1 fs grid seed')
    sizes = [Decimal(str(value))*1000/Decimal(str(timestep)) for value in (total_ps, segment_ps)]
    if any(not v.is_finite() or v != int(v) or v < 20 or int(v) % 10 for v in sizes):
        raise ValueError('Durations must give integer multiples of ten steps, at least twenty each')
    total, segment = map(int, sizes)
    if not 20 < total <= 100000/timestep or (total-20) % segment == 10:
        raise ValueError('Invalid total preparation cap or ten-step remainder')
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    weight, windows = (0. if direction == 'forward' else 1.), []
    remaining, index = total-20, 0
    while remaining:
        steps = min(remaining, segment)
        run = directory/f'p{index:03d}'
        run.mkdir()
        shutil.copyfile(prepared_report.parent/'system.top', run/'system.top')
        shutil.copyfile(prepared_report.parent/('positive.fep' if sign == 1 else 'negative.fep'), run/'charge.fep')
        # Native filenames have an 80-byte limit; a local symlink is deliberately
        # not used to conceal an overlong external path. The input checker fails.
        restart = os.path.relpath(seed_report.parent/'final.re', run) if index == 0 else f'../p{index-1:03d}/final.re'
        content = probe.md_input(prepared['effective_radius_angstrom'], steps, timestep, 0, weight)
        (run/'run.inp').write_text(content.replace('[files]\n', f'[files]\nrestart {restart}\n'))
        windows.append({'input': f'p{index:03d}/run.inp', 'sha256': cp.fingerprint(run/'run.inp'),
                        'assets_sha256': {'topology': cp.fingerprint(run/'system.top'), 'fep': cp.fingerprint(run/'charge.fep')}})
        remaining -= steps
        index += 1
    qdyn = (build_report.parent/build['binaries']['qdyn']['path']).resolve()
    initial = seed_report.parent/'final.re'
    spec = {'schema_version': 2, 'purpose': 'endpoint_preparation',
            'engine': {'binary': str(qdyn), 'sha256': cp.fingerprint(qdyn), 'source_commit': build['source_commit']},
            'build_report': str(build_report), 'build_report_sha256': cp.fingerprint(build_report),
            'initial_restart': str(initial), 'initial_restart_sha256': cp.fingerprint(initial),
            'preparation': {'prepared_report': str(prepared_report), 'prepared_sha256': cp.fingerprint(prepared_report),
                            'seed_report': str(seed_report), 'seed_sha256': cp.fingerprint(seed_report),
                            'total_steps': total, 'segment_steps': segment},
            'series': {'id': f'endpoint-{sign}-{direction}-{replica}',
                       'system': f'probe-r{prepared["grid_radius_angstrom"]}', 'sign': sign,
                       'direction': direction, 'replica': replica, 'born_mode': 'integrated',
                       'apply_born_posthoc': False, 'windows': windows}}
    path = directory/'plan.json'
    probe._json(path, spec)
    from . import charge_chain
    return charge_chain.inspect_plan(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('prepared_report', 'seed_report', 'build_report', 'directory'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--sign', type=int, required=True)
    parser.add_argument('--direction', choices=('forward', 'reverse'), required=True)
    parser.add_argument('--replica', type=int, required=True)
    parser.add_argument('--total-ps', type=float, default=100.)
    parser.add_argument('--segment-ps', type=float, default=20.)
    args = parser.parse_args()
    try:
        report = generate(**vars(args))
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as error:
        parser.exit(2, f'Endpoint plan failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
