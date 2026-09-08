"""One capped Snellius build/probe/real-target smoke job; never production."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile

from QligFEP import charge_build, charge_probe as probe, charge_target_smoke as targets
from QligFEP import charge_protocol as cp, charge_diagnostics as diag


def main(release):
    release = release.resolve(strict=True)
    project = Path('/projects/prjs2157/astra-charge-change-perturbation').resolve(strict=True)
    if not release.is_relative_to(project) or release == project:
        raise ValueError('Release must be a dedicated directory in the authorized project')
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Native cluster checks must run in a scheduler allocation')
    manifest = json.loads((release/'release.json').read_text())
    if cp.fingerprint(release/'source.tar') != manifest['source_archive_sha256']:
        raise ValueError('Source package hash mismatch')
    sources = {}
    with tarfile.open(release/'source.tar') as archive:
        if archive.pax_headers.get('comment') != manifest['commit']:
            raise ValueError('Source package commit mismatch')
        import hashlib
        for member in archive.getmembers():
            if member.isdir():
                continue
            name = Path(member.name)
            if not member.isfile() or name.is_absolute() or '..' in name.parts:
                raise ValueError('Unsafe source package member')
            expected = hashlib.sha256(archive.extractfile(member).read()).hexdigest()
            if cp.fingerprint(release/'source'/name) != expected:
                raise ValueError('Unpacked source differs from package')
            sources[str(name)] = expected
    probe._json(release/'started.json', {'job': os.environ['SLURM_JOB_ID'],
                                       'host': platform.node(), 'python': sys.version,
                                       'manifest_sha256': cp.fingerprint(release/'release.json'),
                                       'source_files': sources})
    build = charge_build.build_archive(
        release/'native-source.tar', release/'build', '/usr/bin/gfortran',
        expected_commit=manifest['commit'], expected_sha256=manifest['native_archive_sha256'])
    charge_build.validate(release/'build/build.json')
    qdyn = release/'build'/build['binaries']['qdyn']['path']
    qprep = release/'build'/build['binaries']['qprep']['path']
    with (release/'solver-tests.log').open('x') as log:
        subprocess.run(['make', '-C', str(release/'build/src/q6'), 'test-lincs', 'test-settle',
                        'FC=/usr/bin/gfortran'], stdout=log, stderr=subprocess.STDOUT, check=True, timeout=60)
    # Four short, fresh charge-only checks before using the real target adapter.
    for radius in (10., 14.):
        prepared = release/f'probe-r{int(radius)}'
        probe.prepare(prepared, qprep, radius, 758971)
        for sign in (-1, 1):
            run = release/f'probe-r{int(radius)}-q{sign}'
            report = probe.timing(prepared, run, qdyn, sign=sign, weight=1., steps=20)
            observed = diag.assess({'signature': {'md': {'steps': '20'}, 'intervals': {'output': '10'}}},
                                   report['native'], run/'native.log')
            probe._json(run/'geometry.json', observed)
    targets.stage(Path('/projects/prjs2157/charge-change/runs'), release/'targets')
    matrix = json.loads((release/'targets/matrix.json').read_text())
    results = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(targets.run_case, release/'targets'/name, qdyn, release/'build/build.json'): name
                   for name in matrix['cases']}
        for future in as_completed(futures):
            name = futures[future]
            try:
                report = future.result()
                native = report['runs']['integrated']['native']
                results[name] = {'passed': True, 'included_non_q_charge': native['parameters'][5],
                                 'excluded_non_q_charge': native['parameters'][6],
                                 'q_state_charges': [row[2] for row in native['state']],
                                 'effective_radius': native['parameters'][0],
                                 'maximum_born_difference_residual': report['maximum_born_difference_residual']}
            except Exception as error:
                results[name] = {'passed': False, 'error': str(error)}
            print(name, json.dumps(results[name]), flush=True)
    if any(cp.fingerprint(release/'source'/name) != expected for name, expected in sources.items()):
        raise ValueError('Source package changed during job')
    charge_build.validate(release/'build/build.json')
    report = {'job': os.environ['SLURM_JOB_ID'], 'commit': manifest['commit'], 'cases': results,
              'all_target_checks_passed': all(r['passed'] for r in results.values()),
              'production_ready': False, 'free_energy_estimated': False,
              'maximum_native_steps': 400, 'maximum_aggregate_ps': .4,
              'requested_cpus': 4, 'requested_walltime_minutes': 20}
    probe._json(release/'summary.json', report)
    print(json.dumps(report, indent=2), flush=True)
    return 0 if report['all_target_checks_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main(Path(sys.argv[1])))
