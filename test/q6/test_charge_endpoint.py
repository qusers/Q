"""Fixed-endpoint planning and tiny native continuation; no equilibrium claim."""
import json
from pathlib import Path
import struct
import subprocess
import sys

import pytest

from QligFEP import charge_analysis, charge_build, charge_chain as chain
from QligFEP import charge_endpoint as endpoint, charge_probe as probe, charge_protocol as cp
from test_charge_chain import isolated_build


@pytest.fixture
def endpoint_factory(isolated_build, tmp_path):
    build = charge_build.validate(isolated_build)
    qdyn = isolated_build.parent/build['binaries']['qdyn']['path']
    qprep = isolated_build.parent/build['binaries']['qprep']['path']
    prepared = tmp_path/'prepared'
    probe.prepare(prepared, qprep, 10., 758971)
    def make(sign=1, direction='forward', total_ps=.1, segment_ps=.03, timestep=1.):
        seed = tmp_path/f'seed-{sign}-{direction}'
        probe.timing(prepared, seed, qdyn, sign=sign,
                     weight=0. if direction == 'forward' else 1., steps=20, timestep=timestep)
        report = endpoint.generate(prepared/'prepared.json', seed/'timing.json', isolated_build,
                                   tmp_path/f'endpoint-{sign}-{direction}', sign=sign,
                                   direction=direction, replica=1, total_ps=total_ps, segment_ps=segment_ps)
        return Path(report['plan_path']), report
    return make


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('direction', ['forward', 'reverse'])
def test_native_endpoint_schedule_is_not_a_free_energy_ladder(endpoint_factory, sign, direction):
    path, report = endpoint_factory(sign, direction)
    origin = report['preparation_origin']
    assert origin['segment_steps'] == [30, 30, 20]
    assert origin['total_steps'] == 100 and report['total_steps'] == 80
    assert origin['total_ps'] == .1 and origin['equilibrated'] is False
    assert origin['independent_replica_established'] is False
    assert not (path.parent/'p000/charge-started.json').exists()
    with pytest.raises(ValueError, match='not a free-energy ladder'):
        charge_analysis.analyze_chain(path, discard_frames=0, bootstrap=50)
    for index in range(3):
        result = chain.run_next(path, max_steps=30, timeout=30)
        assert result['completed_windows'] == index+1
        receipt = json.loads((path.parent/f'p{index:03d}/charge-completed.json').read_text())
        assert receipt['result']['final_restart']['offsets_radians'] == [0.]*len(report['initial_offsets']['offsets_radians'])
    assert chain.run_next(path)['gate'] == 'chain_outputs_consistent'
    assert chain.inspect_plan(path)['preparation_origin']['equilibrated'] is False


@pytest.mark.parametrize('timestep', [.5, 1.])
def test_full_endpoint_budget_includes_grid_seed_without_launch(endpoint_factory, timestep):
    path, report = endpoint_factory(total_ps=100., segment_ps=20., timestep=timestep)
    segment = int(20000/timestep)
    assert report['preparation_origin']['segment_steps'] == [segment]*4+[segment-20]
    assert report['total_steps'] == int(100000/timestep)-20
    assert report['preparation_origin']['total_ps'] == 100.
    with pytest.raises(ValueError, match='step budget'):
        chain.run_next(path)
    assert not (path.parent/'p000/charge-started.json').exists()


def test_endpoint_cli_only_writes_plan_and_refuses_overwrite(endpoint_factory):
    path, report = endpoint_factory()
    origin = report['preparation_origin']
    destination = path.parent/'cli-plan'
    command = [sys.executable, '-m', 'QligFEP.charge_endpoint', origin['prepared_report'],
               origin['seed_report'], report['build_report'], str(destination),
               '--sign', '1', '--direction', 'forward', '--replica', '1',
               '--total-ps', '.1', '--segment-ps', '.03']
    result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=30)
    output = json.loads(result.stdout)
    assert output['purpose'] == 'endpoint_preparation' and output['total_steps'] == 80
    assert not (destination/'p000/charge-started.json').exists()
    before = cp.fingerprint(destination/'plan.json')
    repeat = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert repeat.returncode == 2 and 'File exists' in repeat.stderr
    assert cp.fingerprint(destination/'plan.json') == before


@pytest.mark.parametrize('total,segment', [(101, 20), (.02, .02), (.09, .03), (.025, .02), (float('nan'), 20)])
def test_invalid_preparation_budget_never_launches(endpoint_factory, total, segment):
    with pytest.raises(ValueError, match='cap|remainder|Durations'):
        endpoint_factory(total_ps=total, segment_ps=segment)


@pytest.mark.parametrize('mutation,message', [
    ('budget', 'preparation budget'), ('endpoint', 'fixed starting endpoint'),
    ('velocity', 'retain restart velocities'), ('seed', 'asset changed'),
    ('restraint', 'frozen continuation protocol'),
])
def test_endpoint_identity_and_origin_are_not_self_certified(endpoint_factory, mutation, message):
    path, report = endpoint_factory()
    spec = json.loads(path.read_text())
    if mutation == 'budget':
        spec['preparation']['total_steps'] = 200
    elif mutation == 'seed':
        seed_log = Path(report['preparation_origin']['seed_report']).parent/'native.log'
        seed_log.write_text(seed_log.read_text()+'changed\n')
    else:
        for item in spec['series']['windows']:
            inp = path.parent/item['input']
            old, new = {'endpoint': ('1.00000000 0.00000000', '0.50000000 0.50000000'),
                        'velocity': ('random_seed 0', 'random_seed 112'),
                        'restraint': ('1 0 0 0 10 10 10 0', '1 0 0 0 20 20 20 0')}[mutation]
            assert old in inp.read_text()
            inp.write_text(inp.read_text().replace(old, new))
            item['sha256'] = cp.fingerprint(inp)
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match=message):
        chain.run_next(path)
    assert not (path.parent/'p000/charge-started.json').exists()


@pytest.mark.parametrize('retain', [False, True])
def test_native_seed_zero_uses_restart_velocities(endpoint_factory, retain):
    # Fixed coordinates and offsets; change only restart velocities in a test
    # copy. Positive seeds erase that difference; zero must retain its effect.
    path, report = endpoint_factory()
    original = Path(report['initial_restart']).read_bytes()
    coordinate_record_size = struct.unpack('<i', original[:4])[0]+8
    size = struct.unpack('<i', original[coordinate_record_size:coordinate_record_size+4])[0]
    start = coordinate_record_size+8  # skip record marker and coordinate count
    stop = coordinate_record_size+4+size
    velocities = struct.unpack('<'+'d'*((stop-start)//8), original[start:stop])
    changed = original[:start]+struct.pack('<'+'d'*len(velocities), *(1.2*v for v in velocities))+original[stop:]
    outputs = []
    qdyn = report['engine']['binary']
    for index, restart in enumerate((original, changed)):
        run = path.parent/f'velocity-control-{index}'
        run.mkdir()
        first = path.parent/'p000'
        for name in ('system.top', 'charge.fep'):
            (run/name).write_bytes((first/name).read_bytes())
        (run/'start.re').write_bytes(restart)
        inp = probe.md_input(report['preparation_origin']['effective_radius_angstrom'],
                             20, 1., 0 if retain else 112, 0.)
        (run/'run.inp').write_text(inp.replace('[files]\n', '[files]\nrestart start.re\n'))
        result = subprocess.run([qdyn, 'run.inp'], cwd=run, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0 and 'terminated normally.' in result.stdout
        outputs.append((run/'final.re').read_bytes())
    assert (outputs[0] != outputs[1]) == retain
