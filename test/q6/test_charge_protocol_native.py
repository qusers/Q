"""Exercise staged preflight through actual short charge-only Q MD runs."""
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from QligFEP import charge_protocol as cp
from test_born_serialization import QDYN, PROJECT_ROOT, _energy_state_records, _read_fortran_records
from test_frozen_offset_restart import _run


def test_staged_both_signs_and_directions_run_in_native_q(tmp_path):
    assert QDYN.is_file(), 'Build serial Qdyn first'
    assert QDYN.stat().st_mtime >= max(p.stat().st_mtime for p in (PROJECT_ROOT/'src/q6').glob('*.f90')), 'Rebuild stale Qdyn'
    seed = tmp_path/'seed'
    result = _run(seed)
    assert result.returncode == 0, result.stdout+result.stderr
    offset_record = _read_fortran_records(seed/'final.re')[2]
    series = []
    for sign in (-1, 1):
        for direction in ('forward', 'reverse'):
            windows = []
            for index, weight in enumerate((1., .5, 0.) if direction == 'forward' else (0., .5, 1.)):
                run = tmp_path/f'{sign}-{direction}-{index}'
                run.mkdir()
                shutil.copyfile(seed/'system.top', run/'system.top')
                shutil.copyfile(seed/'final.re', run/'start.re')
                (run/'charge.fep').write_text(f'[FEP]\nstates 2\n[atoms]\n1 1\n[change_charges]\n1 0 {sign}\n')
                inp = (seed/'run.inp').read_text().replace('steps 2\n', 'steps 4\n')
                inp = inp.replace('shake_solvent off', 'shake_solvent on').replace('shake_hydrogens off', 'shake_hydrogens on')
                inp = inp.replace('[solvent]\n', '[solvent]\nradius 20.1\nperstate_born_correction on\nborn_dielectric 80\n')
                inp = inp.replace('[files]\n', '[files]\nrestart start.re\n')
                inp = inp.replace('0.500000 0.500000', f'{weight:.6f} {1-weight:.6f}')
                path = run/'run.inp'
                path.write_text(inp)
                window = cp.inspect_window(path, 'integrated')
                windows.append({'input': str(path), 'sha256': cp.fingerprint(path),
                                'assets_sha256': window['assets_sha256']})
            series.append({'id': f'{sign}-{direction}', 'system': 'native-fixture', 'sign': sign,
                           'direction': direction, 'replica': 1, 'born_mode': 'integrated',
                           'apply_born_posthoc': False, 'windows': windows})
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=PROJECT_ROOT,
                            check=True, capture_output=True, text=True, timeout=10).stdout.strip()
    manifest = tmp_path/'manifest.json'
    manifest.write_text(json.dumps({'schema_version': 1, 'engine': {'binary': str(QDYN),
                                   'sha256': cp.fingerprint(QDYN), 'source_commit': commit}, 'series': series}))
    report = cp.validate(manifest)
    assert report['production_ready'] is False
    for series in report['series']:
        for window in series['windows']:
            run = Path(window['input']).parent
            result = subprocess.run([str(QDYN), 'run.inp'], cwd=run,
                                    capture_output=True, text=True, timeout=30)
            assert result.returncode == 0, result.stdout+result.stderr
            assert _read_fortran_records(run/'final.re')[2] == offset_record
            frames = _energy_state_records(run/'audit.en')
            assert len(frames) == 3
            for frame in frames:
                np.testing.assert_array_equal([frame[1][0], frame[2][0]], list(map(float, window['lambdas'])))
                assert np.isfinite([frame[1], frame[2]]).all()
