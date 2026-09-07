"""Exercise staged preflight through actual short charge-only Q MD runs."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from QligFEP import charge_protocol as cp
from QligFEP import boundary_native as bn
from test_born_serialization import DATA, QDYN, PROJECT_ROOT, _energy_state_records, _read_fortran_records
from test_frozen_offset_restart import _run


@pytest.fixture(scope='module')
def completed_windows(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp('native-protocol')
    assert QDYN.is_file(), 'Build serial Qdyn first'
    assert QDYN.stat().st_mtime >= max(p.stat().st_mtime for p in (PROJECT_ROOT/'src/q6').glob('*.f90')), 'Rebuild stale Qdyn'
    seed = tmp_path/'seed'
    # The archived topology incorrectly labels nonzero hydrogen LJ as SPC-like.
    # Use the existing general three-site kernel in this generated smoke copy;
    # retain every coefficient/charge and never alter the archived input.
    topology = tmp_path/'general-water.top'
    text = (DATA/'topology/Na-benzene-water.top').read_text()
    marker = '       0 = solvent type (0=SPC,1=3-atom,2=general)'
    assert text.count(marker) == 1
    topology.write_text(text.replace(marker, marker.replace('       0', '       1', 1)))
    result = _run(seed, topology=topology)
    assert result.returncode == 0, (result.stdout+result.stderr)[-6000:]
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
    (tmp_path/'preflight.json').write_text(json.dumps(report))
    completed = []
    for series in report['series']:
        for window in series['windows']:
            run = Path(window['input']).parent
            result = subprocess.run([str(QDYN), 'run.inp'], cwd=run,
                                    capture_output=True, text=True, timeout=30)
            assert result.returncode == 0, result.stdout+result.stderr
            log = run/'run.log'
            log.write_text(result.stdout)
            completed.append((window, series['born_mode'], log))
            assert _read_fortran_records(run/'final.re')[2] == offset_record
            frames = _energy_state_records(run/'audit.en')
            assert len(frames) == 3
            for frame in frames:
                np.testing.assert_array_equal([frame[1][0], frame[2][0]], list(map(float, window['lambdas'])))
                assert np.isfinite([frame[1], frame[2]]).all()
    return completed


def test_staged_both_signs_and_directions_run_in_native_q(completed_windows):
    assert len(completed_windows) == 12
    for window, mode, log in completed_windows:
        native = bn.validate_window(window, mode, log)
        assert native['production_ready'] is False
        assert native['effective_radius_angstrom'] == 20.1
        assert native['included_non_q_charge'] == pytest.approx(1.115, abs=1e-8)
        assert native['initial_conservative_cutoff_bound_passed']
        assert native['angular_charge_convention'] == 'legacy_q_region'


@pytest.mark.parametrize('index', [0, -1])
def test_native_audit_command_line(completed_windows, index):
    _, _, log = completed_windows[0]
    preflight = log.parent.parent/'preflight.json'
    result = subprocess.run([sys.executable, '-m', 'QligFEP.boundary_native', str(preflight),
                             '--series=-1-forward', '--window', str(index), '--log', str(log)],
                            capture_output=True, text=True, timeout=15)
    if index == 0:
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['gate'] == 'native_initialization_consistency_passed'
    else:
        assert result.returncode == 2
        assert 'No unique declared series/window' in result.stderr


@pytest.mark.parametrize('mode', ['posthoc', 'control'])
def test_native_nonintegrated_born_accounting(completed_windows, tmp_path, mode):
    original, _, _ = completed_windows[0]
    for name in ('system.top', 'charge.fep', 'start.re'):
        shutil.copyfile(Path(original['input']).parent/name, tmp_path/name)
    path = tmp_path/'run.inp'
    path.write_text(Path(original['input']).read_text().replace('perstate_born_correction on',
                                                               'perstate_born_correction off'))
    window = cp.inspect_window(path, mode)
    result = subprocess.run([str(QDYN), 'run.inp'], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr
    log = tmp_path/'run.log'
    log.write_text(result.stdout)
    report = bn.validate_window(window, mode, log)
    assert all(row[-1] == 0 for row in report['native']['state'])
    assert all(value < 0 for value in report['born_state_constants_kcal_mol'])
    assert report['born_mode'] == mode


def test_archived_incompatible_water_type_is_rejected(completed_windows, tmp_path):
    original, mode, _ = completed_windows[0]
    for name in ('charge.fep', 'start.re', 'run.inp'):
        shutil.copyfile(Path(original['input']).parent/name, tmp_path/name)
    shutil.copyfile(DATA/'topology/Na-benzene-water.top', tmp_path/'system.top')
    window = cp.inspect_window(tmp_path/'run.inp', mode)
    result = subprocess.run([str(QDYN), 'run.inp'], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr
    log = tmp_path/'run.log'
    log.write_text(result.stdout)
    with pytest.raises(ValueError, match='requires zero hydrogen LJ'):
        bn.validate_window(window, mode, log)


@pytest.mark.parametrize('mutation,message', [
    ('missing_end', 'complete native'), ('duplicate_block', 'complete native'),
    ('duplicate_meta', 'meta record'), ('nan_parameter', 'Nonfinite'),
    ('wrong_radius', 'effective radius'), ('wrong_convention', 'charge convention'),
    ('excluded_q', 'excluded/non-solute'), ('zero_lj', 'nonzero Q-atom LJ'),
    ('wrong_charge', 'Q-atom charge'), ('wrong_born', 'applied Born constant'),
    ('wrong_offset', 'frozen offset'), ('wrong_lambda', 'lambda mismatch'),
    ('old_version', 'complete native'), ('water_nonuniform', 'nonuniform'),
    ('water_optimized', 'zero hydrogen LJ'), ('water_charge', 'neutral water charges'),
    ('water_lj_flag', 'compatibility flag'),
])
def test_native_report_rejects_mismatches(completed_windows, tmp_path, mutation, message):
    window, mode, source = completed_windows[0]
    lines = source.read_text().splitlines()
    if mutation == 'missing_end':
        lines.remove('Q_BOUNDARY_AUDIT_V2 END')
    elif mutation == 'old_version':
        lines = [line.replace('Q_BOUNDARY_AUDIT_V2', 'Q_BOUNDARY_AUDIT_V1') for line in lines]
    elif mutation == 'duplicate_block':
        lines *= 2
    elif mutation == 'duplicate_meta':
        index = next(i for i, line in enumerate(lines) if line.startswith('QBA_META '))
        lines.insert(index, lines[index])
    else:
        record, field, value = {
            'nan_parameter': ('PARAMETERS', 1, 'NaN'),
            'wrong_radius': ('PARAMETERS', 1, '21.1'),
            'wrong_convention': ('CONVENTION', 1, '2'),
            'excluded_q': ('QATOM', 3, '1'),
            'zero_lj': ('QATOM', 6, '0'),
            'wrong_charge': ('QATOM', 12, '0.5'),
            'wrong_born': ('STATE', 6, '0'),
            'wrong_offset': ('SHELL', 4, '0.5'),
            'wrong_lambda': ('STATE', 2, '0.2'),
            'water_nonuniform': ('WATER_COMPATIBILITY', 1, '0'),
            'water_optimized': ('META', 8, '0'),
            'water_charge': ('WATER_ATOM', 4, '-0.75'),
            'water_lj_flag': ('WATER_COMPATIBILITY', 2, '1'),
        }[mutation]
        index = next(i for i, line in enumerate(lines) if line.startswith('QBA_'+record+' '))
        fields = lines[index].split()
        fields[field] = value
        lines[index] = ' '.join(fields)
    log = tmp_path/'mutated.log'
    log.write_text('\n'.join(lines))
    with pytest.raises(ValueError, match=message):
        bn.validate_window(window, mode, log)
