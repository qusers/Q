"""Actual short native restart chains and non-destructive failure handling."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from QligFEP import charge_chain as chain, charge_protocol as cp, charge_probe as probe
from QligFEP import charge_build
from QligFEP import charge_analysis
from test_charge_probe import PROJECT_ROOT


@pytest.fixture(scope='module')
def isolated_build(tmp_path_factory):
    compiler = shutil.which('gfortran-11') or shutil.which('gfortran')
    assert compiler, 'Native chain tests require an installed Fortran compiler'
    directory = tmp_path_factory.mktemp('isolated-native')/'build'
    charge_build.build(PROJECT_ROOT, directory, compiler)
    return directory/'build.json'


@pytest.fixture
def plan_factory(isolated_build, tmp_path):
    build = charge_build.validate(isolated_build)
    qdyn = isolated_build.parent/build['binaries']['qdyn']['path']
    qprep = isolated_build.parent/build['binaries']['qprep']['path']
    prepared = tmp_path/'prepared'
    probe.prepare(prepared, qprep, 10., 758971)
    def make(sign=1, direction='forward', mode='integrated'):
        root = tmp_path/f'{sign}-{direction}-{mode}'
        root.mkdir()
        probe.timing(prepared, root/'seed', qdyn, sign=sign,
                     weight=0. if direction == 'forward' else 1., steps=20)
        radius = json.loads((prepared/'prepared.json').read_text())['effective_radius_angstrom']
        windows = []
        for index, weight in enumerate((0., .5, 1.) if direction == 'forward' else (1., .5, 0.)):
            run = root/f'w{index}'
            run.mkdir()
            shutil.copyfile(prepared/'system.top', run/'system.top')
            shutil.copyfile(prepared/('positive.fep' if sign == 1 else 'negative.fep'), run/'charge.fep')
            restart = '../seed/final.re' if index == 0 else f'../w{index-1}/final.re'
            text = probe.md_input(radius, 100, 1., 112+index, weight).replace(
                '[files]\n', f'[files]\nrestart {restart}\n')
            if mode != 'integrated':
                text = text.replace('perstate_born_correction on', 'perstate_born_correction off')
            (run/'run.inp').write_text(text)
            windows.append({'input': f'w{index}/run.inp', 'sha256': cp.fingerprint(run/'run.inp'),
                            'assets_sha256': {key: cp.fingerprint(run/name) for key, name in
                                              [('topology', 'system.top'), ('fep', 'charge.fep')]}})
        spec = {'schema_version': 1, 'engine': {'binary': str(qdyn), 'sha256': cp.fingerprint(qdyn),
                                               'source_commit': build['source_commit']},
                'build_report': str(isolated_build), 'build_report_sha256': cp.fingerprint(isolated_build),
                'initial_restart': 'seed/final.re', 'initial_restart_sha256': cp.fingerprint(root/'seed/final.re'),
                'series': {'id': f'{sign}-{direction}', 'system': 'fresh-probe', 'sign': sign,
                           'direction': direction, 'replica': 1, 'born_mode': mode,
                           'apply_born_posthoc': mode == 'posthoc', 'windows': windows}}
        path = root/'plan.json'
        path.write_text(json.dumps(spec))
        return path
    return make


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('direction', ['forward', 'reverse'])
def test_native_chain_realizes_and_pins_each_restart(plan_factory, sign, direction):
    path = plan_factory(sign, direction)
    report = chain.inspect_plan(path)
    assert report['production_ready'] is False
    assert report['total_steps'] == 300
    assert not (path.parent/'w0/final.re').exists()
    with pytest.raises(FileNotFoundError):
        cp.inspect_window(path.parent/'w1/run.inp', 'integrated')
    previous_hash = report['initial_restart_sha256']
    for index in range(3):
        result = chain.run_next(path, max_steps=100, timeout=30)
        assert result['completed_windows'] == index+1
        assert result['production_ready'] is False
        run = path.parent/f'w{index}'
        preflight = json.loads((run/'charge-preflight.json').read_text())
        assert preflight['window']['assets_sha256']['restart'] == previous_hash
        receipt = json.loads((run/'charge-completed.json').read_text())
        previous_hash = receipt['result']['final_sha256']
        assert previous_hash == cp.fingerprint(run/'final.re')
        assert receipt['result']['saved_frames'] == 9
    final_receipt = cp.fingerprint(path.parent/'w2/charge-completed.json')
    assert chain.run_next(path)['gate'] == 'chain_outputs_consistent'
    assert cp.fingerprint(path.parent/'w2/charge-completed.json') == final_receipt
    analysis = charge_analysis.analyze_chain(path, discard_frames=0, bootstrap=50)
    assert analysis['production_ready'] is False
    assert analysis['analysis']['gap_statistical_gates_passed'] is False
    assert analysis['analysis']['correlation_basis'] == 'gap_and_aligned_observables'
    assert all(item['matched_frames'] == 9 and item['first_retained_energy_step'] == 10 and
               item['last_retained_energy_step'] == 90 for item in analysis['observable_coverage'])
    assert all('temperature_free_kelvin' in item['observable_names'] for item in analysis['observable_coverage'])
    assert analysis['raw_conditional_interval_95_kcal_mol'] is None
    assert analysis['with_born_delta_g_kcal_mol'] == pytest.approx(
        analysis['raw_delta_g_kcal_mol']+analysis['born_delta_0_to_sign_kcal_mol'], abs=1e-12)
    with pytest.raises(ValueError, match='Discard removes all'):
        charge_analysis.analyze_chain(path, discard_frames=9, bootstrap=50)
    if sign == 1 and direction == 'forward':
        command = subprocess.run(
            [sys.executable, '-m', 'QligFEP.charge_analysis', str(path),
             '--discard-frames', '0', '--bootstrap', '50'],
            capture_output=True, text=True, check=True, timeout=30)
        assert json.loads(command.stdout) == analysis
        assert cp.fingerprint(path.parent/'w2/charge-completed.json') == final_receipt
        diagnostic_command = subprocess.run(
            [sys.executable, '-m', 'QligFEP.charge_diagnostics', str(path), '--window', '0'],
            capture_output=True, text=True, check=True, timeout=30)
        diagnostic = json.loads(diagnostic_command.stdout)
        assert diagnostic['force_geometries_observed'] == 101
        assert diagnostic['production_ready'] is False


def test_changed_predecessor_blocks_successor(plan_factory):
    path = plan_factory()
    chain.run_next(path)
    final = path.parent/'w0/final.re'
    final.write_bytes(final.read_bytes()[:-1])
    with pytest.raises(ValueError, match='Damaged restart'):
        chain.run_next(path)
    assert not (path.parent/'w1/charge-started.json').exists()


def test_analysis_rejects_incomplete_chain(plan_factory):
    path = plan_factory()
    with pytest.raises(ValueError, match='completely verified chain'):
        charge_analysis.analyze_chain(path, discard_frames=0, bootstrap=50)
    assert not (path.parent/'w0/charge-started.json').exists()


def test_source_required_trace_cannot_be_dropped_from_a_completed_window(plan_factory):
    path = plan_factory()
    assert chain.inspect_plan(path)['native_trace_required'] is True
    chain.run_next(path)
    log = path.parent/'w0/charge-native.log'
    log.write_text('\n'.join(line for line in log.read_text().splitlines()
                             if not line.startswith(('Q_CHARGE_TRACE', 'QCT_')))+'\n')
    with pytest.raises(ValueError, match='trace is missing'):
        chain.run_next(path)
    assert not (path.parent/'w1/charge-started.json').exists()


@pytest.mark.parametrize('artifact', ['states.en', 'final.re', 'charge-started.json', 'charge-native.log'])
def test_existing_attempt_or_output_is_never_retried(plan_factory, artifact):
    path = plan_factory()
    target = path.parent/'w0'/artifact
    target.write_bytes(b'preserve this interrupted-attempt evidence')
    before = target.read_bytes()
    with pytest.raises(ValueError, match='Existing attempt/output'):
        chain.run_next(path)
    assert target.read_bytes() == before


def test_step_budget_prevents_launch(plan_factory):
    path = plan_factory()
    with pytest.raises(ValueError, match='step budget'):
        chain.run_next(path, max_steps=99)
    assert not (path.parent/'w0/charge-started.json').exists()


@pytest.mark.parametrize('timeout', [False, True])
def test_failed_native_attempt_is_preserved_and_blocks_retry(plan_factory, monkeypatch, timeout):
    path = plan_factory()
    def fail(command, **kwargs):
        if timeout:
            raise subprocess.TimeoutExpired(command, 1.)
        return subprocess.CompletedProcess(command, 7)
    monkeypatch.setattr(chain.subprocess, 'run', fail)
    with pytest.raises(subprocess.TimeoutExpired if timeout else ValueError):
        chain.run_next(path)
    assert (path.parent/'w0/charge-failed.json').exists()
    assert not (path.parent/'w0/charge-completed.json').exists()
    with pytest.raises(ValueError, match='Existing attempt/output'):
        chain.run_next(path)
    assert not (path.parent/'w1/charge-started.json').exists()


def test_wrong_restart_dependency_rejected_before_launch(plan_factory):
    path = plan_factory()
    spec = json.loads(path.read_text())
    inp = path.parent/'w2/run.inp'
    inp.write_text(inp.read_text().replace('../w1/final.re', '../seed/final.re'))
    spec['series']['windows'][2]['sha256'] = cp.fingerprint(inp)
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match='immediately preceding'):
        chain.inspect_plan(path)


def test_output_cannot_overwrite_source_input(plan_factory):
    path = plan_factory()
    spec = json.loads(path.read_text())
    inp = path.parent/'w0/run.inp'
    inp.write_text(inp.read_text().replace('energy states.en', 'energy system.top'))
    spec['series']['windows'][0]['sha256'] = cp.fingerprint(inp)
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match='overwrite a protected'):
        chain.inspect_plan(path)


def test_isolated_build_retains_exact_git_sources_and_refuses_overwrite(isolated_build):
    report = charge_build.validate(isolated_build)
    assert report['production_ready'] is False
    for name in ('src/q6/md.f90', 'src/q6/boundary_corrections.f90', 'src/q6/makefile'):
        expected = subprocess.run(['git', 'show', report['source_commit']+':'+name], cwd=PROJECT_ROOT,
                                  capture_output=True, check=True, timeout=10).stdout
        assert (isolated_build.parent/name).read_bytes() == expected
    before = cp.fingerprint(isolated_build)
    with pytest.raises(FileExistsError):
        charge_build.build(PROJECT_ROOT, isolated_build.parent, report['compiler']['binary'])
    assert cp.fingerprint(isolated_build) == before


@pytest.mark.parametrize('asset,message', [('src/q6/md.f90', 'source snapshot'),
                                          ('src/q6/qdyn', 'executable'),
                                          ('build.log', 'build log'),
                                          ('native-source.tar', 'source archive')])
def test_mutated_build_artifact_rejected(isolated_build, tmp_path, asset, message):
    root = tmp_path/'build-copy'
    shutil.copytree(isolated_build.parent, root)
    with (root/asset).open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match=message):
        charge_build.validate(root/'build.json')


def test_chain_refuses_unrelated_source_commit_label(plan_factory):
    path = plan_factory()
    spec = json.loads(path.read_text())
    spec['engine']['source_commit'] = 'a'*40
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match='isolated build record'):
        chain.inspect_plan(path)


@pytest.mark.parametrize('mutation,message', [('commit', 'archive metadata'),
                                             ('inventory', 'source inventory')])
def test_build_metadata_must_match_the_archived_source(isolated_build, tmp_path, mutation, message):
    root = tmp_path/'build-copy'
    shutil.copytree(isolated_build.parent, root)
    path = root/'build.json'
    report = json.loads(path.read_text())
    if mutation == 'commit':
        report['source_commit'] = 'a'*40
    else:
        report['source_files_sha256'].pop('src/q6/md.f90')
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match=message):
        charge_build.validate(path)


def test_partial_future_output_blocks_launch_of_earlier_window(plan_factory):
    path = plan_factory()
    (path.parent/'w1/states.en').write_bytes(b'preserve unrelated partial future output')
    with pytest.raises(ValueError, match='beyond an incomplete predecessor'):
        chain.run_next(path)
    assert not (path.parent/'w0/charge-started.json').exists()


def test_native_chain_analysis_applies_born_once_in_all_modes(plan_factory):
    results = {}
    for mode in ('integrated', 'posthoc', 'control'):
        path = plan_factory(mode=mode)
        for _ in range(3):
            chain.run_next(path)
        results[mode] = charge_analysis.analyze_chain(path, discard_frames=0, bootstrap=50)
    for result in results.values():
        assert result['raw_delta_g_kcal_mol'] == pytest.approx(results['control']['raw_delta_g_kcal_mol'], abs=1e-10)
        assert result['with_born_delta_g_kcal_mol'] == pytest.approx(results['integrated']['declared_result_kcal_mol'], abs=1e-10)
        assert result['native_boltzmann_kcal_mol_kelvin'] == pytest.approx(.001986, abs=1e-10)
    assert results['posthoc']['declared_result_kcal_mol'] == pytest.approx(results['integrated']['declared_result_kcal_mol'], abs=1e-10)
    assert results['control']['declared_result_kcal_mol'] == results['control']['raw_delta_g_kcal_mol']
    assert results['control']['control_corrected_view_is_comparison_only'] is True
