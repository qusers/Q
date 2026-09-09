"""Build the pinned engine and qualify standard 16-rank Snellius execution."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

import numpy as np

from QligFEP import charge_target_smoke as target, charge_protocol as cp
from QligFEP.charge_probe import _json


def main(root):
    root = root.resolve(strict=True)
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('An allocation is required')
    old = Path('/projects/prjs2157/astra-charge-change-perturbation/releases/no-softcore-787e76eb-20260908')
    archive = old/'native-source.tar'
    expected = '3e1dcd376a906e893bb15dffb44b3072771e43fed86ad5e3245cabbbdac58c52'
    assert cp.fingerprint(archive) == expected
    build = root/'build'
    build.mkdir()
    shutil.copyfile(archive, build/'native-source.tar')
    with tarfile.open(archive) as tar:
        assert tar.pax_headers['comment'] == '787e76ebf40c37b3a68dd63d1132df4ce78bba9f'
        sources = {}
        for member in tar.getmembers():
            path = Path(member.name)
            assert not path.is_absolute() and '..' not in path.parts
            assert member.isdir() or member.isfile()
            if member.isfile():
                sources[member.name] = hashlib.sha256(tar.extractfile(member).read()).hexdigest()
        tar.extractall(build)
    compiler = shutil.which('gfortran')
    mpifc = shutil.which('mpif90')
    command = ['make', '-j1', 'qdyn', 'qprep', 'qfep', 'qdynp', 'FC='+compiler, 'MPIFC='+mpifc]
    with (build/'build.log').open('x') as log:
        subprocess.run(command, cwd=build/'src/q6', stdout=log, stderr=subprocess.STDOUT,
                       timeout=180, check=True)
    binaries = {name: {'path': str(build/'src/q6'/name), 'sha256': cp.fingerprint(build/'src/q6'/name)}
                for name in ('qdyn', 'qdynp', 'qprep', 'qfep')}
    _json(build/'build.json', {'source_commit': '787e76ebf40c37b3a68dd63d1132df4ce78bba9f',
                             'archive_sha256': expected, 'source_files': sources,
                             'binaries': binaries, 'command': command,
                             'compiler_version': subprocess.check_output([compiler, '--version'], text=True),
                             'mpi_version': subprocess.check_output(['mpirun', '--version'], text=True)})
    reports = {}
    for name in ('cmet-fwd-2.protein-w0.0001', 'eg5-fwd-2.protein-w0.0001'):
        reference = old/'targets'/name/'reference'
        case = root/'mpi-gate'/name
        plan = target.stage_case(reference, reference/'eq5.re', case,
                                 {'target': name.split('-')[0]}, no_softcore=True, weight=.0001)
        arrays = {}
        elapsed = {}
        for mode in ('integrated', 'posthoc', 'serial'):
            directory = case/mode
            if mode == 'serial':
                shutil.copytree(case/'posthoc', directory,
                                ignore=shutil.ignore_patterns('native.log', 'initialization.json', 'gate-check.json',
                                                              'states.en', 'final.re'))
            cmd = ([binaries['qdyn']['path']] if mode == 'serial' else
                   ['mpirun', '-n', '16', '--map-by', 'core', '--bind-to', 'core', binaries['qdynp']['path']])
            start = time.perf_counter()
            with (directory/'native.log').open('x') as log:
                subprocess.run(cmd+['run.inp'], cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=120)
            report, saved = target.check_run(directory, plan, mode == 'integrated')
            arrays[mode] = np.asarray(saved)
            elapsed[mode] = time.perf_counter()-start
            _json(directory/'gate-check.json', report)
        assert (case/'integrated/final.re').read_bytes() == (case/'posthoc/final.re').read_bytes()
        delta = arrays['integrated']-arrays['posthoc']
        delta[:, :, 1] -= report['born_constants']
        assert np.max(abs(delta)) < 1e-8
        serial_delta = float(np.max(abs(arrays['serial']-arrays['posthoc'])))
        assert serial_delta < 1e-6, serial_delta
        reports[name] = {'born_residual': float(np.max(abs(delta))), 'serial_mpi_energy_residual': serial_delta,
                         'elapsed_seconds': elapsed}
    assert all(cp.fingerprint(build/name) == sha for name, sha in sources.items())
    assert all(cp.fingerprint(Path(v['path'])) == v['sha256'] for v in binaries.values())
    _json(root/'mpi-gate-passed.json', {'job': os.environ['SLURM_JOB_ID'], 'cases': reports,
                                      'build_sha256': cp.fingerprint(build/'build.json')})
    print(json.dumps(reports, indent=2), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
