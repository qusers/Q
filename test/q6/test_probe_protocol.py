"""Fresh charge-only probes pass staged and native shared-restraint checks."""
import json
from pathlib import Path
import shutil
import subprocess
import struct
import sys

import pytest

from QligFEP import boundary_native as bn, charge_protocol as cp, charge_probe as probe
from QligFEP import charge_completion as completion
from test_charge_probe import prepared, QDYN, PROJECT_ROOT


@pytest.fixture(scope='module')
def completed_probe_windows(prepared, tmp_path_factory):
    root = tmp_path_factory.mktemp('staged-probe')
    seed = root/'seed'
    # 20 native steps create a real restart with the shared positional restraint
    # and zero offsets. This is not equilibration or independent replica sampling.
    probe.timing(prepared, seed, QDYN, sign=1, weight=0., steps=20)
    radius = json.loads((prepared/'prepared.json').read_text())['effective_radius_angstrom']
    series = []
    for sign in (-1, 1):
        for direction in ('forward', 'reverse'):
            windows = []
            weights = (0., .5, 1.) if direction == 'forward' else (1., .5, 0.)
            for index, weight in enumerate(weights):
                run = root/f'{sign}-{direction}-{index}'
                run.mkdir()
                shutil.copyfile(prepared/'system.top', run/'system.top')
                shutil.copyfile(prepared/('positive.fep' if sign == 1 else 'negative.fep'), run/'charge.fep')
                shutil.copyfile(seed/'final.re', run/'start.re')
                text = probe.md_input(radius, 100, 1., 112+index, weight)
                text = text.replace('[files]\n', '[files]\nrestart start.re\n')
                path = run/'run.inp'
                path.write_text(text)
                window = cp.inspect_window(path, 'integrated')
                windows.append({'input': str(path), 'sha256': cp.fingerprint(path),
                                'assets_sha256': window['assets_sha256']})
            series.append({'id': f'{sign}-{direction}', 'system': 'fresh-probe-r10', 'sign': sign,
                           'direction': direction, 'replica': 1, 'born_mode': 'integrated',
                           'apply_born_posthoc': False, 'windows': windows})
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=PROJECT_ROOT,
                            capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    manifest = root/'manifest.json'
    manifest.write_text(json.dumps({'schema_version': 1, 'engine': {'binary': str(QDYN),
                                   'sha256': cp.fingerprint(QDYN), 'source_commit': commit}, 'series': series}))
    report = cp.validate(manifest)
    (root/'preflight.json').write_text(json.dumps(report))
    completed = []
    for series in report['series']:
        for window in series['windows']:
            run = Path(window['input']).parent
            probe._execute(QDYN, run, ['run.inp'])
            log = run/'native.log'
            assert 'terminated normally.' in log.read_text()
            probe._check_energies(run/'states.en', list(map(float, window['lambdas'])), 9)
            assert cp.restart_offsets(run/'final.re')['offset_record_sha256'] == window['restart']['offset_record_sha256']
            completed.append((window, log))
    return completed


def test_fresh_probe_both_signs_directions_and_endpoints_are_staged(completed_probe_windows):
    assert len(completed_probe_windows) == 12
    for window, log in completed_probe_windows:
        completed = completion.validate_window(window, 'integrated', log)
        assert completed['production_ready'] is False
        assert completed['saved_frames'] == 9
        assert completed['maximum_accounting_residual_kcal_mol'] < 1e-10
        result = completed['native_initialization']
        assert result['production_ready'] is False
        assert result['included_non_q_charge'] == 0
        assert result['effective_radius_angstrom'] == 10.15
        assert result['native']['restraint_counts'] == [0, 1, 0, 0, 0, 0]
        assert result['native']['position'] == [[1, 1, 0, 0, 0, 10, 10, 10, 0]]
        assert result['initial_conservative_cutoff_bound_passed']


@pytest.mark.parametrize('mutation,message', [
    ('missing', 'position dimensions'), ('duplicate', 'position dimensions'),
    ('wrong_index', 'position indices'), ('fractional_atom', 'atom/state indices'),
    ('fractional_state', 'atom/state indices'), ('wrong_atom', 'shared position'),
    ('wrong_state', 'shared position'), ('wrong_center', 'shared position'),
    ('wrong_force', 'shared position'), ('extra_sequence', 'extra restraint counts'),
    ('implicit_file', 'extra restraint counts'), ('negative_count', 'restraint counts'),
])
def test_native_shared_restraint_mismatch_is_rejected(completed_probe_windows, tmp_path, mutation, message):
    window, original = completed_probe_windows[0]
    lines = original.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith('QBA_POSITION '))
    if mutation == 'missing':
        lines.pop(index)
    elif mutation == 'duplicate':
        lines.insert(index, lines[index])
    else:
        record, field, value = {
            'wrong_index': ('POSITION', 1, '2'),
            'fractional_atom': ('POSITION', 2, '1.5'),
            'fractional_state': ('POSITION', 9, '0.5'),
            'wrong_atom': ('POSITION', 2, '2'),
            'wrong_state': ('POSITION', 9, '1'),
            'wrong_center': ('POSITION', 3, '0.1'),
            'wrong_force': ('POSITION', 6, '11'),
            'extra_sequence': ('RESTRAINT_COUNTS', 1, '1'),
            'implicit_file': ('RESTRAINT_COUNTS', 6, '1'),
            'negative_count': ('RESTRAINT_COUNTS', 1, '-1'),
        }[mutation]
        index = next(i for i, line in enumerate(lines) if line.startswith('QBA_'+record+' '))
        fields = lines[index].split()
        fields[field] = value
        lines[index] = ' '.join(fields)
    log = tmp_path/'mutated.log'
    log.write_text('\n'.join(lines))
    with pytest.raises(ValueError, match=message):
        bn.validate_window(window, 'integrated', log)


def test_completion_command_uses_retained_preflight(completed_probe_windows):
    window, log = completed_probe_windows[0]
    root = Path(window['input']).parent.parent
    result = subprocess.run([sys.executable, '-m', 'QligFEP.charge_completion', str(root/'preflight.json'),
                             '--series=-1-forward', '--window=0', '--log', str(log)],
                            capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)['gate'] == 'completed_window_consistency_passed'


@pytest.mark.parametrize('mutation,message', [('window', 'window declaration'),
                                             ('direction', 'series declaration'),
                                             ('source', 'engine declaration')])
def test_retained_manifest_is_not_bypassed_by_report_edits(completed_probe_windows, tmp_path, mutation, message):
    window, _ = completed_probe_windows[0]
    report = json.loads((Path(window['input']).parent.parent/'preflight.json').read_text())
    if mutation == 'window':
        report['series'][0]['windows'][0] = report['series'][0]['windows'][1]
    elif mutation == 'direction':
        report['series'][0]['direction'] = 'reverse'
    else:
        report['engine']['source_commit'] = 'a'*40
    path = tmp_path/'preflight.json'
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match=message):
        bn.load_window(path, '-1-forward', 0)


@pytest.mark.parametrize('mutation,message', [
    ('truncated_energy', 'Damaged energy record'), ('empty_energy', 'energy frames'),
    ('missing_frame', 'energy frames'), ('extra_frame', 'energy frames'),
    ('nonfinite_energy', 'Nonfinite saved energy'), ('wrong_lambda', 'Incorrect saved lambda'),
    ('missing_born', 'pure-state total'), ('double_born', 'pure-state total'),
    ('subtotal', 'electrostatic subtotal'), ('changed_offset', 'frozen offsets'),
    ('nonfinite_coordinates', 'Nonfinite restart'), ('truncated_restart', 'Damaged restart'),
    ('missing_exit', 'normally completed'), ('duplicate_exit', 'normally completed'),
    ('missing_final_summary', 'final energy summary'),
])
def test_completion_rejects_damaged_or_incomplete_outputs(completed_probe_windows, tmp_path, mutation, message):
    original, source_log = completed_probe_windows[0]
    run = tmp_path/'copy'
    shutil.copytree(Path(original['input']).parent, run)
    window = cp.inspect_window(run/'run.inp', 'integrated')
    log = run/'native.log'
    energy, restart = run/'states.en', run/'final.re'
    if mutation in ('truncated_energy', 'empty_energy', 'missing_frame', 'extra_frame'):
        data = energy.read_bytes()
        data = {'truncated_energy': data[:-1], 'empty_energy': b'',
                'missing_frame': data[:-272], 'extra_frame': data+data[:272]}[mutation]
        energy.write_bytes(data)
    elif mutation in ('nonfinite_energy', 'wrong_lambda', 'missing_born', 'double_born', 'subtotal'):
        data = bytearray(energy.read_bytes())
        # Second pure-state record in the first frame, after marker + state index.
        offset = 132+8
        values = list(struct.unpack_from('<15d', data, offset))
        if mutation == 'nonfinite_energy':
            values[1] = float('inf')
        elif mutation == 'wrong_lambda':
            values[0] = .123
        elif mutation == 'subtotal':
            values[6] += 1
        else:
            born = bn.parse(source_log.read_text())['state'][1][-1]
            values[1] += born if mutation == 'double_born' else -born
        struct.pack_into('<15d', data, offset, *values)
        energy.write_bytes(data)
    elif mutation in ('changed_offset', 'nonfinite_coordinates', 'truncated_restart'):
        data = bytearray(restart.read_bytes())
        if mutation == 'changed_offset':
            struct.pack_into('<f', data, len(data)-8, .1)
        elif mutation == 'nonfinite_coordinates':
            struct.pack_into('<d', data, 8, float('nan'))
        else:
            data = data[:-1]
        restart.write_bytes(data)
    else:
        text = log.read_text()
        if mutation == 'missing_exit':
            text = text.replace('terminated normally.', 'exit marker removed')
        elif mutation == 'duplicate_exit':
            text += '\nterminated normally.\n'
        else:
            text = text.replace('FINAL  Energy summary', 'final summary removed')
        log.write_text(text)
    with pytest.raises(ValueError, match=message):
        completion.validate_window(window, 'integrated', log)


@pytest.mark.parametrize('mode', ['posthoc', 'control'])
def test_completion_checks_nonintegrated_accounting(completed_probe_windows, tmp_path, mode):
    original, _ = completed_probe_windows[0]
    run = tmp_path/'run'
    run.mkdir()
    source = Path(original['input']).parent
    for name in ('system.top', 'charge.fep', 'start.re'):
        shutil.copyfile(source/name, run/name)
    (run/'run.inp').write_text((source/'run.inp').read_text().replace('perstate_born_correction on',
                                                                  'perstate_born_correction off'))
    window = cp.inspect_window(run/'run.inp', mode)
    probe._execute(QDYN, run, ['run.inp'])
    result = completion.validate_window(window, mode, run/'native.log')
    assert result['production_ready'] is False
    assert result['saved_frames'] == 9
    assert result['native_initialization']['born_mode'] == mode
