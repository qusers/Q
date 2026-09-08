"""Adversarial input-contract checks; synthetic assets do not validate physics."""
import json
import struct
import subprocess
import sys

import pytest

from QligFEP import charge_protocol as cp


INPUT = '''[MD]
steps 100
stepsize 1
temperature 298
bath_coupling 10
initial_temperature 298
random_seed 112
shake_solvent on
shake_hydrogens on
shake_solute off
lrf off
[cut-offs]
solute_solvent 99
solute_solute 99
solvent_solvent 99
q_atom 99
lrf 99
[sphere]
shell_force 10
shell_radius 20
[solvent]
radius 20.1
radial_force 60
polarization on
polarization_force 20
charge_correction on
perstate_polarization on
polarization_adaptation off
perstate_born_correction on
born_dielectric 80
[intervals]
output 10
non_bond 10
energy 10
[files]
topology system.top
fep charge.fep
restart start.re
final final.re
energy energy.en
[lambdas]
{weight} {complement}
'''


def _record(data):
    marker = struct.pack('<i', len(data))
    return marker+data+marker


@pytest.fixture
def staged(tmp_path):
    # Explicitly synthetic: tests the gate, not a runnable Q system.
    binary = tmp_path/'qdyn'
    binary.write_bytes(b'not an executable')
    series = []
    for sign in (-1, 1):
        for direction in ('forward', 'reverse'):
            windows = []
            for index, weight in enumerate((1, .5, 0) if direction == 'forward' else (0, .5, 1)):
                directory = tmp_path/f'{sign}-{direction}-{index}'
                directory.mkdir()
                (directory/'system.top').write_text('synthetic topology; native properties unverified')
                payload = _record(struct.pack('<i6d', 6, *([0.]*6)))
                (directory/'start.re').write_bytes(payload*2+_record(struct.pack('<i3f', 3, .01, -.02, .03)))
                (directory/'charge.fep').write_text(f'[FEP]\nstates 2\n[atoms]\n1 1\n[change_charges]\n1 0 {sign}\n')
                path = directory/'run.inp'
                path.write_text(INPUT.format(weight=weight, complement=1-weight))
                windows.append({'input': str(path), 'sha256': cp.fingerprint(path), 'assets_sha256': {
                    key: cp.fingerprint(directory/name) for key, name in
                    [('topology', 'system.top'), ('fep', 'charge.fep'), ('restart', 'start.re')]}})
            series.append({'id': f'{sign}-{direction}', 'system': 'water-r20', 'sign': sign,
                           'direction': direction, 'replica': 1, 'born_mode': 'integrated',
                           'apply_born_posthoc': False, 'windows': windows})
    spec = {'schema_version': 1, 'engine': {'binary': str(binary), 'sha256': cp.fingerprint(binary),
                                         'source_commit': 'a'*40}, 'series': series}
    path = tmp_path/'manifest.json'
    path.write_text(json.dumps(spec))
    return path, spec


def rewrite(path, spec):
    path.write_text(json.dumps(spec))


def change_input(spec, old, new, *, series=0, window=1, refresh=True):
    from pathlib import Path
    item = spec['series'][series]['windows'][window]
    path = Path(item['input'])
    assert old in path.read_text()
    path.write_text(path.read_text().replace(old, new))
    if refresh:
        item['sha256'] = cp.fingerprint(path)


def test_staged_consistency_is_not_production_qualification(staged):
    path, _ = staged
    report = cp.validate(path)
    assert report['gate'] == 'staged_input_consistency_passed'
    assert report['production_ready'] is False
    assert len(report['series']) == 4
    assert 'native nonzero unchanged LJ parameters and interaction coverage' in report['unverified']


def test_restart_velocity_mode_is_explicit_and_consistent(staged):
    path, spec = staged
    for s in range(4):
        for w in range(3):
            change_input(spec, 'random_seed 112', 'random_seed 0', series=s, window=w)
    rewrite(path, spec)
    report = cp.validate(path)
    assert report['series'][0]['windows'][0]['signature']['velocity_initialization'] == 'restart'
    change_input(spec, 'random_seed 0', 'random_seed 112')
    rewrite(path, spec)
    with pytest.raises(ValueError, match='settings differ'):
        cp.validate(path)


@pytest.mark.parametrize('value', ['-1', '100000000', '0.5'])
def test_unsupported_restart_seed_rejected(staged, value):
    path, spec = staged
    change_input(spec, 'random_seed 112', 'random_seed '+value)
    rewrite(path, spec)
    with pytest.raises(ValueError, match='random_seed'):
        cp.validate(path)


def test_shared_position_is_part_of_system_identity(staged):
    path, spec = staged
    for s in range(4):
        for w in range(3):
            change_input(spec, '[lambdas]', '[atom_restraints]\n1 0.1 -0.2 0.3 10 20 30 0\n[lambdas]', series=s, window=w)
    rewrite(path, spec)
    report = cp.validate(path)
    assert report['series'][0]['windows'][0]['signature']['atom_restraints'] == [[1, '0.1', '-0.2', '0.3', '10', '20', '30', 0]]
    change_input(spec, '0.1 -0.2 0.3', '0.1 -0.2 0.4')
    rewrite(path, spec)
    with pytest.raises(ValueError, match='settings differ'):
        cp.validate(path)


@pytest.mark.parametrize('row,message', [
    ('1 0 0 0 10 10 10 1', 'shared state 0'),
    ('1 0 0 0 10 10 10 -1', 'shared state 0'),
    ('1 0 0 0 10 10 10 256', 'shared state 0'),
    ('1 0 0 0 10 0 10 0', 'positive Cartesian'),
    ('1 0 0 0 10 -1 10 0', 'positive Cartesian'),
    ('2 0 0 0 10 10 10 0', 'unique mapped Q'),
    ('1.5 0 0 0 10 10 10 0', 'unique mapped Q'),
    ('1 NaN 0 0 10 10 10 0', 'Nonfinite'),
    ('1 1e999 0 0 10 10 10 0', 'finite doubles'),
    ('1 0 0 0 10 10 10', 'needs atom'),
    ('1 0 0 0 10 10 10 0\n1 0 0 0 10 10 10 0', 'unique mapped Q'),
])
def test_invalid_or_state_dependent_positions_are_rejected(staged, row, message):
    path, spec = staged
    change_input(spec, '[lambdas]', '[atom_restraints]\n'+row+'\n[lambdas]')
    rewrite(path, spec)
    with pytest.raises(ValueError, match=message):
        cp.validate(path)


@pytest.mark.parametrize('valid', [True, False])
def test_documented_command_line_entrypoint(staged, valid):
    path, spec = staged
    if not valid:
        spec['series'][0]['apply_born_posthoc'] = True
        rewrite(path, spec)
    result = subprocess.run([sys.executable, '-m', 'QligFEP.charge_protocol', str(path)],
                            capture_output=True, text=True, timeout=15)
    if valid:
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['production_ready'] is False
    else:
        assert result.returncode == 2
        assert 'double-count' in result.stderr
        assert not result.stdout


@pytest.mark.parametrize('old,new,message', [
    ('polarization_adaptation off', 'polarization_adaptation on', 'polarization_adaptation'),
    ('perstate_born_correction on', 'perstate_born_correction off', 'perstate_born_correction'),
    ('born_dielectric 80', 'born_dielectric 40', 'dielectric 80'),
    ('radius 20.1', 'radius 21.1', 'settings differ'),
    ('temperature 298', 'temperature 300', 'settings differ'),
    ('lrf off', 'lrf on', 'lrf off'),
    ('energy 10', 'energy 100', 'retain production'),
    ('steps 100', 'steps 100.5', 'integer'),
    ('0.5 0.5', '0.4 0.5', 'normalized'),
    ('[MD]', '[MD]\nsteps 10', 'Duplicate'),
    ('[MD]', '[MD]\nunknown_option yes', 'Unsupported'),
    ('[MD]', '[MD]\n[MD]', 'Duplicate section'),
])
def test_unsafe_inputs_rejected(staged, old, new, message):
    path, spec = staged
    change_input(spec, old, new)
    rewrite(path, spec)
    with pytest.raises(ValueError, match=message):
        cp.validate(path)


def test_modified_input_hash_rejected(staged):
    path, spec = staged
    change_input(spec, 'temperature 298', 'temperature 299', refresh=False)
    with pytest.raises(ValueError, match='input hash mismatch'):
        cp.validate(path)


@pytest.mark.parametrize('mode,posthoc,enabled', [('posthoc', True, 'off'), ('control', False, 'off')])
def test_explicit_nonintegrated_accounting_modes(staged, mode, posthoc, enabled):
    path, spec = staged
    for s, series in enumerate(spec['series']):
        series.update(born_mode=mode, apply_born_posthoc=posthoc)
        for w in range(3):
            change_input(spec, 'perstate_born_correction on', f'perstate_born_correction {enabled}', series=s, window=w)
    rewrite(path, spec)
    assert cp.validate(path)['production_ready'] is False


def test_double_born_rejected(staged):
    path, spec = staged
    spec['series'][0]['apply_born_posthoc'] = True
    rewrite(path, spec)
    with pytest.raises(ValueError, match='double-count'):
        cp.validate(path)


@pytest.mark.parametrize('mutation,message', [
    ('[softcore]\n1 0 0\n', 'Charge-only FEP'),
    ('[change_atoms]\n1 REAL DUMMY\n', 'Charge-only FEP'),
    ('[change_bonds]\n1 2 1 2\n', 'Charge-only FEP'),
])
def test_non_charge_fep_sections_rejected(staged, mutation, message):
    from pathlib import Path
    path, spec = staged
    fep = Path(spec['series'][0]['windows'][0]['input']).with_name('charge.fep')
    fep.write_text(fep.read_text()+mutation)
    with pytest.raises(ValueError, match=message):
        cp.validate(path)


def test_different_frozen_offsets_rejected(staged):
    from pathlib import Path
    path, spec = staged
    restart = Path(spec['series'][1]['windows'][0]['input']).with_name('start.re')
    restart.write_bytes(restart.read_bytes()[:-24]+_record(struct.pack('<i3f', 3, .01, -.02, .04)))
    spec['series'][1]['windows'][0]['assets_sha256']['restart'] = cp.fingerprint(restart)
    rewrite(path, spec)
    with pytest.raises(ValueError, match='offsets differ'):
        cp.validate(path)


def test_existing_output_rejected(staged):
    from pathlib import Path
    path, spec = staged
    Path(spec['series'][0]['windows'][0]['input']).with_name('energy.en').touch()
    with pytest.raises(ValueError, match='existing output'):
        cp.validate(path)


def test_modified_asset_hash_rejected(staged):
    from pathlib import Path
    path, spec = staged
    Path(spec['series'][0]['windows'][0]['input']).with_name('system.top').write_text('changed')
    with pytest.raises(ValueError, match='asset hash mismatch'):
        cp.validate(path)


def test_modified_engine_hash_rejected(staged):
    from pathlib import Path
    path, spec = staged
    Path(spec['engine']['binary']).write_bytes(b'changed engine')
    with pytest.raises(ValueError, match='Engine hash mismatch'):
        cp.validate(path)


def test_shared_output_path_rejected(staged):
    from pathlib import Path
    path, spec = staged
    first_output = Path(spec['series'][0]['windows'][0]['input']).parent/'energy.en'
    # Use a short relative path, within Q's native filename field.
    change_input(spec, 'energy energy.en', f'energy ../{first_output.parent.name}/energy.en')
    rewrite(path, spec)
    with pytest.raises(ValueError, match='Output collision'):
        cp.validate(path)


@pytest.mark.parametrize('kind', ['missing', 'nan_offset', 'bad_count', 'nan_coordinate', 'short_marker', 'bad_trailer'])
def test_bad_restart_is_rejected(tmp_path, kind):
    coordinates = _record(struct.pack('<i6d', 6, *([0.]*6)))
    offset = _record(struct.pack('<i3f', 3, .1, .2, .3))
    payload = coordinates*2+offset
    if kind == 'missing':
        payload = coordinates*2
    elif kind == 'nan_offset':
        payload = coordinates*2+_record(struct.pack('<i3f', 3, .1, float('nan'), .3))
    elif kind == 'bad_count':
        payload = coordinates*2+_record(struct.pack('<i3f', 2, .1, .2, .3))
    elif kind == 'nan_coordinate':
        payload = _record(struct.pack('<i6d', 6, float('nan'), *([0.]*5)))+coordinates+offset
    elif kind == 'short_marker':
        payload += b'x'
    elif kind == 'bad_trailer':
        payload = payload[:-1]
    path = tmp_path/'bad.re'
    path.write_bytes(payload)
    with pytest.raises(ValueError):
        cp.restart_offsets(path)


def test_exact_endpoints_required_for_this_charge_only_gate(staged):
    path, spec = staged
    spec['series'][0]['windows'] = spec['series'][0]['windows'][1:]
    rewrite(path, spec)
    with pytest.raises(ValueError, match='including both endpoints'):
        cp.validate(path)
