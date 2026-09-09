"""Stage archived 101-window inputs with QligFEP's standard Snellius array renderer.

No topology/force-field regeneration, new sampler, or automatic resubmission.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from types import SimpleNamespace

from QligFEP import charge_protocol as cp, boundary_native as bn, charge_diagnostics as diag
from QligFEP.charge_completion import frames
from QligFEP.charge_target_smoke import EDGES, without_softcore
from QligFEP.charge_probe import _json


def serialize(raw):
    return ''.join('['+key+']\n'+''.join(' '.join(row)+'\n' for row in rows) for key, rows in raw.items())


def adapt_input(source, radius):
    raw = cp.sections(source)
    md = cp.keyed(raw['md'])
    # Preserve physical duration and the original heating/restraint schedule.
    timestep = float(md['stepsize'])
    if timestep == 2.:
        md['steps'] = str(2*int(md['steps']))
        md['stepsize'] = '1.0'
    elif timestep != .2:
        raise ValueError('Unexpected archived timestep')
    md['constraint_algorithm'] = 'shake shake'
    md['lrf'] = 'off'
    raw['md'] = [line.split() for line in (k+' '+v for k, v in md.items())]
    raw['cut-offs'] = [[key, '120'] for key in ('solute_solute', 'solute_solvent', 'solvent_solvent', 'q_atom', 'lrf')]
    solvent = cp.keyed(raw['solvent'])
    assert solvent['polarisation'] == 'on'
    solvent.update(radius=str(radius), charge_correction='on', perstate_polarization='on',
                   polarization_adaptation='off', perstate_born_correction='on', born_dielectric='80')
    raw['solvent'] = [[key, value] for key, value in solvent.items()]
    intervals = cp.keyed(raw['intervals'])
    intervals.update(output='100', non_bond='1')
    intervals.pop('trajectory', None)
    raw['intervals'] = [[key, value] for key, value in intervals.items()]
    files = cp.keyed(raw['files'])
    files.pop('trajectory', None)
    raw['files'] = [[key, value] for key, value in files.items()]
    # The old correction logger belongs to a different accounting convention.
    raw.pop('correction', None)
    raw.pop('trajectory_atoms', None)
    return serialize(raw)


def trim_qfep(source):
    lines = source.read_text().splitlines()
    names = [line.strip() for line in lines[9:] if line.strip()]
    assert int(lines[0]) == len(names) == 101
    assert names[0] == 'md_1000_0000.en' and names[-1] == 'md_0000_1000.en'
    assert len(set(names)) == 101
    lines[0] = '99'
    return '\n'.join(lines[:9]+names[1:-1])+'\n'


def stage(reference_root, destination, remote_root, cluster):
    if cluster != 'SNELLIUS':
        raise ValueError('This protocol is scoped to SNELLIUS')
    from QligFEP.qligfep import QligFEP
    from QligFEP.settings.settings import CLUSTER_DICT
    destination.mkdir(parents=True, exist_ok=False)
    original_profile = CLUSTER_DICT['SNELLIUS']
    CLUSTER_DICT['SNELLIUS'] = {**original_profile, 'ACCOUNT': 'ugsei19097', 'TIME': '1-00:00:00',
        'MODULES': 'module purge\nmodule load 2023 OpenMPI/4.1.5-GCC-12.3.0\n',
        'QDYN': f'qdyn={remote_root}/build/src/q6/qdynp',
        'QFEP': f'{remote_root}/build/src/q6/qfep'}
    all_cases = []
    try:
        for target, (left, right) in EDGES.items():
            for direction in ('fwd', 'rev'):
                lig1, lig2 = (left, right) if direction == 'fwd' else (right, left)
                for leg in ('1.water', '2.protein'):
                    relative = Path(f'{target}-{direction}')/leg/f'FEP_{lig1}_{lig2}'
                    source = reference_root/relative/'inputfiles'
                    edge = destination/relative
                    original = edge/'reference'
                    shutil.copytree(source, original)
                    inputs = edge/'inputfiles'
                    inputs.mkdir()
                    for name in ('dualtop.top',):
                        shutil.copyfile(original/name, inputs/name)
                    (inputs/'FEP1.fep').write_text(without_softcore(original/'FEP1.fep'))
                    lines = [line for line in (inputs/'dualtop.top').read_text().splitlines()
                             if '= Exclusion, solvent radii' in line]
                    assert len(lines) == 1
                    radius = float(lines[0].split('=')[0].split()[1])
                    names = sorted(original.glob('eq*.inp'))+sorted(original.glob('md*.inp'))
                    assert len(names) == 106
                    for source_input in names:
                        (inputs/source_input.name).write_text(adapt_input(source_input, radius))
                    mds = sorted(inputs.glob('md*.inp'))
                    weights = {path.name: list(map(float, cp.sections(path)['lambdas'][0])) for path in mds}
                    assert len(weights) == len({tuple(w) for w in weights.values()}) == 101
                    restraints = [json.dumps({k: cp.sections(path).get(k, []) for k in
                                  ('sequence_restraints', 'distance_restraints', 'wall_restraints')}, sort_keys=True)
                                  for path in mds]
                    assert len(set(restraints)) == 1
                    upper = sorted([name for name, w in weights.items() if w[0] > .5], key=lambda n: weights[n][0])
                    lower = sorted([name for name, w in weights.items() if w[0] < .5], key=lambda n: -weights[n][0])
                    assert len(upper) == len(lower) == 50
                    order = ['eq1.inp', 'eq2.inp', 'eq3.inp', 'eq4.inp', 'eq5.inp', 'md_0500_0500.inp']
                    order += [name for pair in zip(upper, lower) for name in pair]
                    available = set()
                    for name in order:
                        files = cp.keyed(cp.sections(inputs/name)['files'])
                        if 'restart' in files:
                            assert files['restart'] in available, (name, files['restart'])
                        available.add(files['final'])
                    (inputs/'qfep.inp').write_text(trim_qfep(original/'qfep.inp'))
                    shutil.copyfile(original/'qfep.inp', original/'qfep-full-endpoints.inp')
                    # Invoke the same renderer used by qligfep --cluster SNELLIUS.
                    obj = SimpleNamespace(cluster=cluster, start='0.5', replicates=3, temperature='298',
                                          seeds=[2924, 25360, 21448], system=leg.split('.')[1],
                                          lig1=lig1, lig2=lig2, to_clean=None)
                    QligFEP.write_runfile(obj, str(inputs), [[], upper, lower])
                    QligFEP.write_submitfile(obj, str(edge))
                    script = inputs/'runSNELLIUS.sh'
                    text = script.read_text()
                    text = text.replace('#SBATCH --nodes=1', '#SBATCH --nodes=1\n#SBATCH --no-requeue')
                    # Keep standard per-edge working-directory and replica-array semantics.
                    text = text.replace('workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"',
                                        'workdir="${SLURM_SUBMIT_DIR:?Submit from the FEP edge directory}"')
                    text = text.replace('mkdir -p $rundir', 'mkdir "$rundir"')
                    driver = f'{remote_root}/standard_101.py'
                    python = '/gpfs/home3/dvidal/bash/envs/qligfep_new/bin/python'
                    source_python = '/projects/prjs2157/astra-charge-change-perturbation/releases/no-softcore-787e76eb-20260908/source/src'
                    exports = (f'export PYTHONPATH={source_python}\nexport PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1\n'
                               'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n')
                    text = text.replace('set -eo pipefail', 'set -eo pipefail\n'+exports+
                                        f'{python} {driver} preflight "$workdir" {remote_root}\n')
                    text = re.sub(r'^(time mpirun .* (\S+)\.inp > (\S+)\.log)$',
                                  lambda m: m[0]+f'\n{python} {driver} check {m[2]}.inp {m[3]}.log', text, flags=re.M)
                    assert text.count(f'{driver} check ') == 106
                    # A timeout is not a successful analysis. Keep raw files regardless.
                    text = text.replace(' || [ $? -eq 124 ]', '')
                    script.write_text(text)
                    _json(edge/'analysis-scope.json', {'sampled_windows': 101, 'analyzed_windows': 99,
                          'excluded_energy_files': ['md_1000_0000.en', 'md_0000_1000.en'],
                          'state_2_interval': [.001, .999], 'born_mode': 'integrated',
                          'apply_born_posthoc': False, 'free_energy_qualification': 'pending cross-lambda and sampling checks'})
                    production_ps = sum(int(cp.keyed(cp.sections(p)['md'])['steps'])*
                                        float(cp.keyed(cp.sections(p)['md'])['stepsize'])/1000 for p in mds)
                    all_ps = sum(int(cp.keyed(cp.sections(inputs/n)['md'])['steps'])*
                                 float(cp.keyed(cp.sections(inputs/n)['md'])['stepsize'])/1000 for n in order)
                    assert production_ps == 1010.
                    _json(edge/'protocol.json', {'cluster': cluster, 'target': target, 'direction': direction,
                          'leg': leg, 'replicas': 3, 'seeds': obj.seeds, 'radius': radius,
                          'production_ps_per_replica': production_ps, 'total_ps_per_replica': all_ps,
                          'order': order, 'source_reference': str(source.resolve()),
                          'input_hashes': {p.name: cp.fingerprint(p) for p in inputs.iterdir() if p.is_file()},
                          'reference_hashes': {p.name: cp.fingerprint(p) for p in original.iterdir() if p.is_file()},
                          'submission_sha256': cp.fingerprint(edge/'FEP_submit.sh')})
                    subprocess.run(['bash', '-n', str(script)], check=True)
                    subprocess.run(['bash', '-n', str(edge/'FEP_submit.sh')], check=True)
                    all_cases.append(str(relative))
    finally:
        CLUSTER_DICT['SNELLIUS'] = original_profile
    _json(destination/'matrix.json', {'cluster': cluster, 'cases': all_cases, 'arrays': 8,
                                     'replica_jobs': 24, 'windows_per_replica': 101,
                                     'source_native_commit': '787e76ebf40c37b3a68dd63d1132df4ce78bba9f'})


def preflight(edge, root):
    gate = json.loads((root/'mpi-gate-passed.json').read_text())
    assert cp.fingerprint(root/'build/build.json') == gate['build_sha256']
    build = json.loads((root/'build/build.json').read_text())
    for binary in build['binaries'].values():
        assert cp.fingerprint(Path(binary['path'])) == binary['sha256']
    manifest = json.loads((root/'driver-hashes.json').read_text())
    for name, checksum in manifest.items():
        assert cp.fingerprint(root/name) == checksum
    plan = json.loads((edge/'protocol.json').read_text())
    for name, checksum in plan['input_hashes'].items():
        assert cp.fingerprint(edge/'inputfiles'/name) == checksum
    assert cp.fingerprint(edge/'FEP_submit.sh') == plan['submission_sha256']
    print('Pinned engine, driver, and standard inputs verified.', flush=True)


def check(source, log):
    raw = cp.sections(source)
    md, files = cp.keyed(raw['md']), cp.keyed(raw['files'])
    text = log.read_text()
    if text.count('terminated normally.') != 1 or 'terminated abnormally' in text or 'WARNING: hot atom' in text:
        raise ValueError('Native completion or hot-atom gate failed')
    assert 'No softcore section found. Using normal LJ potentials.' in text
    audit = bn.parse(text)
    assert audit['flags'] == [1, 1, 0, 1, 0, 0]
    assert audit['constraint_algorithms'] == ['shake', 'shake']
    assert audit['water_compatibility'] == [1, 1]
    observations = diag.assess({'signature': {'md': md, 'intervals': cp.keyed(raw['intervals'])}}, audit, log)
    final = cp.restart_offsets(source.parent/files['final'])
    assert all(value == 0 for value in final['offsets_radians'])
    if 'restart' in files:
        initial = cp.restart_offsets(source.parent/files['restart'])
        assert final['offset_record_sha256'] == initial['offset_record_sha256']
    weights = list(map(float, raw['lambdas'][0]))
    count = 0
    exact_endpoint = min(weights) == 0.
    if 'energy' in files and not exact_endpoint:
        for frame in frames(source.parent/files['energy'], weights):
            count += 1
            for values, state in zip(frame, audit['state']):
                assert math.isclose(values[1], sum(values[2:8])+values[14]+state[-1], abs_tol=1e-8, rel_tol=1e-12)
        assert count == (int(md['steps'])-1)//int(cp.keyed(raw['intervals'])['energy'])
    _json(source.with_suffix('.checked.json'), {'normal_completion': True, 'diagnostics': observations,
          'native': audit, 'checked_energy_frames': count, 'exact_endpoint_excluded_from_energy_analysis': exact_endpoint,
          'input_sha256': cp.fingerprint(source), 'log_sha256': cp.fingerprint(log),
          'final_restart_sha256': cp.fingerprint(source.parent/files['final'])})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    setup = commands.add_parser('stage')
    setup.add_argument('reference_root', type=Path)
    setup.add_argument('destination', type=Path)
    setup.add_argument('--remote-root', required=True)
    setup.add_argument('-c', '--cluster', required=True, choices=['SNELLIUS'])
    pre = commands.add_parser('preflight')
    pre.add_argument('edge', type=Path)
    pre.add_argument('root', type=Path)
    chk = commands.add_parser('check')
    chk.add_argument('source', type=Path)
    chk.add_argument('log', type=Path)
    args = vars(parser.parse_args())
    command = args.pop('command')
    globals()[command](**args)
