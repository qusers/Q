"""Bounded real-target compatibility checks, not a free-energy campaign.

Preserve reference dual-topology files and restraints; test the current spherical
Hamiltonian at lambda 0.5 with integrated/post-hoc Born accounting. This is
separate from the deliberately restricted charge-only campaign validator.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import time

from . import boundary_native as bn, charge_build, charge_diagnostics as diag, charge_protocol as cp
from .charge_completion import frames
from .charge_probe import _json

EDGES = {'cmet': ('CHEMBL3402742_23', 'CHEMBL3402744_300'),
         'eg5': ('CHEMBL1085666', 'CHEMBL1089056')}
STEPS = 20


def zero_offsets(source, destination):
    """Copy coordinates/velocities exactly; declare a new zero-offset Hamiltonian."""
    info = cp.restart_offsets(source)
    data = source.read_bytes()
    offset = 0
    for _ in range(2):
        offset += struct.unpack('<i', data[offset:offset+4])[0]+8
    count = len(info['offsets_radians'])
    payload = struct.pack('<i', count)+bytes(4*count)
    marker = struct.pack('<i', len(payload))
    with destination.open('xb') as stream:
        stream.write(data[:offset]+marker+payload+marker)
    assert destination.read_bytes()[:offset] == data[:offset]
    assert cp.restart_offsets(destination)['offsets_radians'] == [0.]*count
    return info


def render(source, radius, integrated, *, steps=STEPS, weight=.5):
    if type(steps) is not int or not 20 <= steps <= 2000 or steps % 10:
        raise ValueError('Target checks require 20..2000 steps in multiples of ten')
    if not math.isfinite(weight) or not 0 < weight < 1:
        raise ValueError('Target checks exclude exact endpoints')
    raw = cp.sections(source)
    required = {'md', 'cut-offs', 'sphere', 'solvent', 'intervals', 'files', 'lambdas'}
    extra = {'trajectory_atoms', 'correction', 'sequence_restraints', 'distance_restraints', 'wall_restraints'}
    if not required <= raw.keys() or set(raw)-required-extra:
        raise ValueError('Unsupported target MD sections')
    if [v for row in raw['lambdas'] for v in row] != ['0.500', '0.500']:
        raise ValueError('Select the existing midpoint template explicitly')
    original = cp.keyed(raw['md'])
    if original.get('shake_solvent') != 'on' or original.get('shake_hydrogens') != 'on':
        raise ValueError('Target smoke expects existing rigid water/hydrogen constraints')
    if original.get('shake_solute') != 'off':
        raise ValueError('Unexpected solute constraint selection')
    md = {'steps': str(steps), 'stepsize': '1.0', 'temperature': '298',
          'bath_coupling': original['bath_coupling'], 'random_seed': '0',
          'initial_temperature': '298', 'shake_solvent': 'on', 'shake_hydrogens': 'on',
          'shake_solute': 'off', 'constraint_algorithm': 'shake shake', 'lrf': 'off'}
    if 'separate_scaling' in original:
        md['separate_scaling'] = original['separate_scaling']
    # Refuse extra dynamics options rather than silently dropping them.
    if set(original)-set(md):
        raise ValueError('Unsupported original dynamics options')
    solvent = cp.keyed(raw['solvent'])
    if set(solvent)-{'radial_force', 'polarisation', 'polarisation_force'}:
        raise ValueError('Unsupported original solvent options')
    if solvent.get('polarisation') != 'on':
        raise ValueError('Expected original solvent polarization')
    config = {
        'MD': md, 'cut-offs': {key: '99' for key in ('solute_solute', 'solute_solvent', 'solvent_solvent', 'q_atom', 'lrf')},
        'sphere': cp.keyed(raw['sphere']),
        'solvent': {'radius': str(radius), 'radial_force': solvent['radial_force'],
                    'polarization': 'on', 'polarization_force': solvent['polarisation_force'],
                    'charge_correction': 'on', 'perstate_polarization': 'on',
                    'polarization_adaptation': 'off', 'born_dielectric': '80',
                    'perstate_born_correction': 'on' if integrated else 'off'},
        'intervals': {'output': '10', 'energy': '1', 'non_bond': '1'},
        'files': {'topology': 'dualtop.top', 'fep': 'FEP1.fep', 'restart': 'start.re',
                  'final': 'final.re', 'energy': 'states.en'},
    }
    result = ''.join('['+section+']\n'+''.join(k+' '+v+'\n' for k, v in values.items()) for section, values in config.items())
    result += f'[lambdas]\n{1-weight:.12g} {weight:.12g}\n'
    for section in ('sequence_restraints', 'distance_restraints', 'wall_restraints'):
        if raw.get(section):
            result += '['+section+']\n'+''.join(' '.join(row)+'\n' for row in raw[section])
    return result


def without_softcore(source):
    """Remove only the archived softcore settings, preserving all physical tables."""
    raw = cp.sections(source)
    if set(raw) != {'fep', 'atoms', 'change_charges', 'atom_types', 'softcore', 'change_atoms'}:
        raise ValueError('Unsupported target FEP dialect')
    if cp.keyed(raw['fep']) != {'states': '2', 'softcore_use_max_potential': 'on', 'softcore_method': 'gapsys'}:
        raise ValueError('Expected historical Gapsys source for explicit no-softcore conversion')
    result = '[FEP]\nstates 2\nsoftcore_use_max_potential off\nsoftcore_method standard\n'
    for section, rows in raw.items():
        if section not in {'fep', 'softcore'}:
            result += '['+section+']\n'+''.join(' '.join(row)+'\n' for row in rows)
    return result


def stage_case(input_directory, restart, destination, identity, *, no_softcore=False, steps=STEPS, weight=.5):
    destination.mkdir(parents=True, exist_ok=False)
    reference = destination/'reference'
    reference.mkdir()
    originals = {name: input_directory/name for name in ('dualtop.top', 'FEP1.fep', 'md_0500_0500.inp')}
    originals['eq5.re'] = restart
    provenance = {}
    for name, source in originals.items():
        checksum = cp.fingerprint(source)
        shutil.copyfile(source, reference/name)
        if cp.fingerprint(source) != checksum or cp.fingerprint(reference/name) != checksum:
            raise ValueError('Reference changed while copying')
        provenance[name] = {'path': str(source.resolve()), 'sha256': checksum}
    radii = [line for line in (reference/'dualtop.top').read_text().splitlines() if '= Exclusion, solvent radii' in line]
    if len(radii) != 1:
        raise ValueError('Missing topology radius')
    radius = float(radii[0].split('=')[0].split()[1])
    fep = cp.sections(reference/'FEP1.fep')
    if set(fep) != {'fep', 'atoms', 'change_charges', 'atom_types', 'softcore', 'change_atoms'}:
        raise ValueError('Unsupported target FEP dialect')
    if cp.keyed(fep['fep']) != {'states': '2', 'softcore_use_max_potential': 'on', 'softcore_method': 'gapsys'}:
        raise ValueError('Expected unchanged historical Gapsys softcore')
    old_offsets = None
    for mode in ('integrated', 'posthoc'):
        run = destination/mode
        run.mkdir()
        for name in ('dualtop.top', 'FEP1.fep'):
            shutil.copyfile(reference/name, run/name)
        if no_softcore:
            (run/'FEP1.fep').write_text(without_softcore(reference/'FEP1.fep'))
        old_offsets = zero_offsets(reference/'eq5.re', run/'start.re')
        (run/'run.inp').write_text(render(reference/'md_0500_0500.inp', radius, mode == 'integrated',
                                        steps=steps, weight=weight))
    plan = {'schema_version': 1, 'identity': identity, 'originals': provenance,
            'source_restart': old_offsets, 'radius': radius, 'steps_per_run': steps,
            'timestep_fs': 1., 'lambda': [1-weight, weight], 'production_ready': False,
            'softcore': 'none' if no_softcore else 'historical_gapsys',
            'files': {str(f.relative_to(destination)): cp.fingerprint(f)
                      for f in sorted(destination.rglob('*')) if f.is_file()}}
    _json(destination/'plan.json', plan)
    return plan


def stage(root, destination, *, no_softcore=False, weights=(.5,), steps=STEPS):
    destination.mkdir(parents=True, exist_ok=False)
    cases = []
    for target, (left, right) in EDGES.items():
        for direction in ('fwd', 'rev'):
            edge = f'FEP_{left}_{right}' if direction == 'fwd' else f'FEP_{right}_{left}'
            for leg in ('1.water', '2.protein'):
                origin = root/f'{target}_pilot__none__r20__{direction}__rep1'/leg/edge
                identity = {'target': target, 'direction': direction, 'leg': leg, 'edge': edge, 'replica': 1}
                for weight in weights:
                    name = f'{target}-{direction}-{leg}'
                    if tuple(weights) != (.5,):
                        name += f'-w{weight:.4f}'
                    stage_case(origin/'inputfiles', origin/'FEP1/298/1/eq5.re', destination/name, identity,
                               no_softcore=no_softcore, steps=steps, weight=weight)
                    cases.append(name)
    _json(destination/'matrix.json', {'cases': cases, 'total_steps': len(cases)*2*steps,
                                    'aggregate_ps': len(cases)*2*steps/1000,
                                    'weights': list(weights), 'softcore': 'none' if no_softcore else 'historical_gapsys',
                                    'production_ready': False})


def check_run(directory, plan, integrated):
    text = (directory/'native.log').read_text()
    audit = bn.parse(text)
    # Preserve initialization evidence even when a later diagnostic rejects a run.
    _json(directory/'initialization.json', audit)
    if text.count('terminated normally.') != 1 or 'terminated abnormally' in text:
        raise ValueError('Target did not complete normally')
    if plan.get('softcore') == 'none':
        active = cp.sections(directory/'FEP1.fep')
        if 'softcore' in active or cp.keyed(active['fep']) != {
                'states': '2', 'softcore_use_max_potential': 'off', 'softcore_method': 'standard'}:
            raise ValueError('No-softcore input contract violated')
        if 'No softcore section found. Using normal LJ potentials.' not in text:
            raise ValueError('Native executable did not confirm normal LJ potentials')
    if audit['flags'] != [1, int(integrated), 0, 1, 0, 0] or audit.get('constraint_algorithms') != ['shake', 'shake']:
        raise ValueError('Target flags/solver mismatch')
    if audit['water_compatibility'] != [1, 1]:
        raise ValueError('Reference water is incompatible with optimized water interactions')
    if any(row[2] != 0 for row in audit['qatom']):
        raise ValueError('Excluded Q atoms require a separate physical definition')
    radius, _, ke, eps, _, env, *_ = audit['parameters']
    bn._close(radius, plan['radius'], 'effective radius')
    bn._close(eps, 80., 'dielectric')
    constants = [-ke*(1-1/eps)*(env+row[2])**2/(2*radius) for row in audit['state']]
    weights = plan['lambda']
    steps = plan['steps_per_run']
    for row, constant, weight in zip(audit['state'], constants, weights):
        bn._close(row[1], weight, 'state weight')
        bn._close(row[-1], constant if integrated else 0., 'Born self energy')
    observations = diag.assess({'signature': {'md': {'steps': str(steps)}, 'intervals': {'output': '10'}}},
                               audit, directory/'native.log')
    final = cp.restart_offsets(directory/'final.re')
    initial = cp.restart_offsets(directory/'start.re')
    if final['atoms'] != initial['atoms'] or final['offset_record_sha256'] != initial['offset_record_sha256']:
        raise ValueError('Target dimensions/frozen offsets changed')
    saved = list(frames(directory/'states.en', weights))
    if len(saved) != steps-1:
        raise ValueError('Wrong target energy-frame count')
    for frame in saved:
        for values, row in zip(frame, audit['state']):
            for got, wanted in ((values[6], sum(values[i] for i in (8, 10, 12))),
                                (values[7], sum(values[i] for i in (9, 11, 13))),
                                (values[1], sum(values[2:8])+values[14]+row[-1])):
                if not math.isclose(got, wanted, abs_tol=1e-8, rel_tol=1e-12):
                    raise ValueError('Target pure-state bookkeeping mismatch')
    gaps = [frame[1][1]-frame[0][1] for frame in saved]
    return {'native': audit, 'born_constants': constants, 'diagnostics': observations,
            'energy_gap_range_kcal_mol': [min(gaps), max(gaps)],
            'pure_state_total_ranges_kcal_mol': [[min(f[s][1] for f in saved), max(f[s][1] for f in saved)]
                                                for s in range(2)],
            'saved_frames': len(saved), 'production_ready': False}, saved


def run_case(directory, binary, build_report):
    plan = json.loads((directory/'plan.json').read_text())
    binary = binary.resolve(strict=True)
    build = charge_build.validate(build_report)
    if (build_report.parent/build['binaries']['qdyn']['path']).resolve() != binary:
        raise ValueError('Target executable is not the recorded isolated build')
    pinned = {name: cp.fingerprint(path) for name, path in
              [('binary', binary), ('build_report', build_report), ('driver', Path(__file__)),
               ('plan', directory/'plan.json')]}
    for name, expected in plan['files'].items():
        if cp.fingerprint(directory/name) != expected:
            raise ValueError('Staged input changed')
    drivers = {name: cp.fingerprint(Path(module.__file__)) for name, module in
               [('boundary', bn), ('build', charge_build), ('diagnostics', diag), ('protocol', cp)]}
    _json(directory/'started.json', {'hashes': pinned, 'drivers': drivers,
                                     'binary': str(binary), 'build_report': str(build_report)})
    try:
        outcomes, saved = {}, {}
        for mode in ('integrated', 'posthoc'):
            run = directory/mode
            started = time.perf_counter()
            with (run/'native.log').open('x') as log:
                result = subprocess.run([str(binary), 'run.inp'], cwd=run, stdout=log,
                                        stderr=subprocess.STDOUT,
                                        timeout=900 if plan['steps_per_run'] > STEPS else 120)
            if result.returncode:
                raise ValueError(f'Native {mode} failed with exit {result.returncode}')
            outcomes[mode], saved[mode] = check_run(run, plan, mode == 'integrated')
            outcomes[mode]['elapsed_seconds'] = time.perf_counter()-started
        if (directory/'integrated/final.re').read_bytes() != (directory/'posthoc/final.re').read_bytes():
            raise ValueError('Force-free Born changed the trajectory')
        largest = 0.
        for a, b in zip(saved['integrated'], saved['posthoc']):
            for state, (on, off) in enumerate(zip(a, b)):
                constant = outcomes['integrated']['born_constants'][state]
                for i, (x, y) in enumerate(zip(on, off)):
                    residual = abs((x-y)-(constant if i == 1 else 0.))
                    largest = max(largest, residual)
                    if residual > 1e-8:
                        raise ValueError('Integrated/post-hoc state energies disagree')
        if cp.fingerprint(binary) != pinned['binary'] or cp.fingerprint(Path(__file__)) != pinned['driver']:
            raise ValueError('Executable or driver changed')
        charge_build.validate(build_report)
        if cp.fingerprint(build_report) != pinned['build_report'] or cp.fingerprint(directory/'plan.json') != pinned['plan']:
            raise ValueError('Build report or plan changed')
        for name, module in [('boundary', bn), ('build', charge_build), ('diagnostics', diag), ('protocol', cp)]:
            if cp.fingerprint(Path(module.__file__)) != drivers[name]:
                raise ValueError('Validation implementation changed')
        for name, expected in plan['files'].items():
            if cp.fingerprint(directory/name) != expected:
                raise ValueError('Staged input changed during run')
        report = {'gate': 'real_target_compatibility_passed', 'identity': plan['identity'],
                  'lambda': plan['lambda'], 'softcore': plan.get('softcore', 'historical_gapsys'),
                  'runs': outcomes, 'maximum_born_difference_residual': largest,
                  'outputs': {str(f.relative_to(directory)): cp.fingerprint(f) for f in sorted(directory.rglob('*')) if f.is_file()},
                  'production_ready': False,
                  'limitations': ['short fixed-lambda check; no free-energy or endpoint qualification',
                                  'no-softcore dual topology' if plan.get('softcore') == 'none' else 'archived dual-topology Gapsys softcore',
                                  'new zero-offset/direct-electrostatics Hamiltonian is not equilibrated by old restart',
                                  'legacy Q-region-only angular target; charged-background adequacy untested']}
        _json(directory/'completed.json', report)
        return report
    except Exception as error:
        _json(directory/'failed.json', {'error': str(error), 'production_ready': False})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    stage_parser = sub.add_parser('stage')
    stage_parser.add_argument('root', type=Path)
    stage_parser.add_argument('destination', type=Path)
    run_parser = sub.add_parser('run')
    run_parser.add_argument('directory', type=Path)
    run_parser.add_argument('binary', type=Path)
    run_parser.add_argument('build_report', type=Path)
    args = vars(parser.parse_args())
    command = args.pop('command')
    if command == 'stage':
        stage(**args)
    else:
        print(json.dumps(run_case(**args), allow_nan=False))


if __name__ == '__main__':
    main()
