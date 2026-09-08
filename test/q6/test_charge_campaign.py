"""Full pilot matrix at software-only durations; no HPC or equilibrium claim."""
import itertools
import json
from pathlib import Path
import subprocess
import sys

import pytest

from QligFEP import charge_analysis, charge_build, charge_campaign as campaign
from QligFEP import charge_chain as chain, charge_endpoint as endpoint, charge_probe as probe, charge_protocol as cp
from test_charge_chain import isolated_build
from test_charge_endpoint import endpoint_factory


@pytest.fixture(scope='module')
def assembled(isolated_build, tmp_path_factory):
    root = tmp_path_factory.mktemp('campaign')
    build = charge_build.validate(isolated_build)
    qdyn = isolated_build.parent/build['binaries']['qdyn']['path']
    qprep = isolated_build.parent/build['binaries']['qprep']['path']
    plans = []
    for index, (radius, sign, direction, replica) in enumerate(itertools.product((10., 14.), (-1, 1), ('forward', 'reverse'), (1, 2))):
        cell = root/f'c{index:02d}'
        prepared = cell/'prepared'
        probe.prepare(prepared, qprep, radius, 758971+index)
        probe.timing(prepared, cell/'seed', qdyn, sign=sign,
                     weight=0. if direction == 'forward' else 1., steps=20, seed=112+index)
        result = endpoint.generate(prepared/'prepared.json', cell/'seed/timing.json', isolated_build,
                                   cell/'endpoint', sign=sign, direction=direction, replica=replica,
                                   total_ps=.1, segment_ps=.03)
        plans.append(Path(result['plan_path']))
    with pytest.MonkeyPatch.context() as patch:
        def no_launch(*args, **kwargs):
            raise AssertionError('Campaign assembly must not launch any process')
        patch.setattr(subprocess, 'run', no_launch)
        result = campaign.assemble(plans, root/'campaign', profile='software_smoke')
    return root/'campaign/campaign.json', result


def test_assembly_records_complete_matrix_without_launch(assembled):
    path, _ = assembled
    result = campaign.inspect(path)
    assert len(result['cells']) == 16
    assert result['budget']['total_ps'] == 6.4
    assert result['budget']['total_steps'] == 6400
    assert all(cell['completed_preparation_segments'] == 0 for cell in result['cells'])
    assert result['production_ready'] is False
    assert result['independent_sampling_established'] is False
    assert len({cell['topology_identity']['coordinates_sha256'] for cell in result['cells']}) == 16
    assert len({cell['topology_identity']['parameters_sha256'] for cell in result['cells']}) == 2


def test_uncompleted_endpoint_cannot_be_transferred(endpoint_factory, tmp_path):
    path, _ = endpoint_factory()
    spec = {'endpoint_origin': {'plan': str(path), 'sha256': cp.fingerprint(path), 'final_sha256': 'a'*64}}
    with pytest.raises(ValueError, match='completely verified'):
        campaign.validate_transfer(tmp_path/'ladder.json', spec, [])


@pytest.fixture(scope='module')
def completed_campaign(assembled):
    path, _ = assembled
    initial = campaign.inspect(path)
    for index, cell in enumerate(initial['cells']):
        endpoint_path = Path(cell['endpoint_plan'])
        for _ in range(cell['preparation_segments']):
            chain.run_next(endpoint_path, max_steps=30, timeout=30)
        result = campaign.stage_ladder(path, index)
        assert result['endpoint_origin']['equilibrated'] is False
        assert result['endpoint_origin']['independent_replica_established'] is False
        ladder_path = Path(result['plan_path'])
        for _ in range(3):
            chain.run_next(ladder_path, max_steps=100, timeout=30)
    return path, campaign.inspect(path)


def test_all_matrix_cells_transfer_run_and_analyze(completed_campaign):
    _, result = completed_campaign
    assert all(cell['completed_ladder_windows'] == 3 for cell in result['cells'])
    for cell in result['cells']:
        analysis = charge_analysis.analyze_chain(Path(cell['ladder_path']), discard_frames=0, bootstrap=50)
        assert analysis['estimate_status'] == 'insufficient_sampling'
        assert analysis['raw_conditional_interval_95_kcal_mol'] is None
        assert analysis['production_ready'] is False


def test_campaign_inspection_cli_and_staging_overwrite_refusal(completed_campaign):
    path, result = completed_campaign
    command = subprocess.run([sys.executable, '-m', 'QligFEP.charge_campaign', 'inspect', str(path)],
                             capture_output=True, text=True, check=True, timeout=60)
    assert json.loads(command.stdout) == result
    first = Path(result['cells'][0]['ladder_path'])
    before = cp.fingerprint(first)
    with pytest.raises(FileExistsError):
        campaign.stage_ladder(path, 0)
    assert cp.fingerprint(first) == before


@pytest.mark.parametrize('mutation,message', [('duplicate', 'duplicate'), ('budget', 'budget'),
                                             ('missing', 'matrix'), ('destination', 'unique deterministic')])
def test_campaign_matrix_cannot_be_relabelled(assembled, tmp_path, mutation, message):
    path, _ = assembled
    spec = json.loads(path.read_text())
    if mutation == 'duplicate':
        spec['endpoints'][1] = spec['endpoints'][0]
    elif mutation == 'budget':
        spec['budget']['total_ps'] = 0.
    elif mutation == 'missing':
        spec['endpoints'].pop()
    else:
        spec['endpoints'][0]['ladder_directory'] = '../unrelated'
    changed = tmp_path/'campaign.json'
    changed.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match=message):
        campaign.inspect(changed)


def test_transfer_rechecks_endpoint_identity_not_just_restart_hash(completed_campaign, tmp_path):
    _, report = completed_campaign
    original = Path(report['cells'][0]['ladder_path'])
    spec = json.loads(original.read_text())
    ladder = chain.inspect_plan(original)
    spec['series']['replica'] = 99
    with pytest.raises(ValueError, match='series identities'):
        campaign.validate_transfer(original, spec, ladder['series']['windows'])
    spec['endpoint_origin']['plan'] = str(original)
    spec['endpoint_origin']['sha256'] = cp.fingerprint(original)
    with pytest.raises(ValueError, match='not another ladder'):
        campaign.validate_transfer(original, spec, ladder['series']['windows'])


@pytest.mark.parametrize('timestep,steps', [(1., 6880000), (.5, 13760000)])
def test_production_profile_budget_is_explicit_without_running_it(timestep, steps):
    budget = campaign.protocol('feasibility_pilot', timestep)
    assert budget['total_ps'] == 6880.
    assert budget['total_steps'] == steps
    assert budget['canonical_weights'] == [i/10 for i in range(11)]
    assert budget['discard_ps_per_window'] == 10.
