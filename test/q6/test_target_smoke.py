"""Target adapter and archive-transfer gates, independent of private input data."""
from pathlib import Path
import struct

import pytest

from QligFEP import charge_build, charge_target_smoke as smoke, charge_protocol as cp
from test_charge_chain import isolated_build


def test_zero_offsets_preserves_coordinate_and_velocity_bytes(tmp_path):
    def record(payload):
        marker = struct.pack('<i', len(payload))
        return marker+payload+marker
    first = record(struct.pack('<i9d', 9, *range(9)))
    second = record(struct.pack('<i9d', 9, *range(9, 18)))
    source, destination = tmp_path/'old.re', tmp_path/'new.re'
    original = first+second+record(struct.pack('<i2f', 2, .25, -.5))
    source.write_bytes(original)
    before = smoke.zero_offsets(source, destination)
    assert before['offsets_radians'] == [.25, -.5]
    assert source.read_bytes() == original
    assert destination.read_bytes()[:len(first+second)] == first+second
    assert cp.restart_offsets(destination)['offsets_radians'] == [0., 0.]
    with pytest.raises(FileExistsError):
        smoke.zero_offsets(source, destination)


@pytest.fixture
def template(tmp_path):
    text = """[MD]
steps 5000
stepsize 2.0
temperature T_VAR
bath_coupling 10
shake_solvent on
shake_hydrogens on
shake_solute off
lrf on
separate_scaling on
[cut-offs]
q_atom 99
[sphere]
shell_radius 20
shell_force 10
[solvent]
radial_force 60
polarisation on
polarisation_force 20
[intervals]
energy 10
[files]
fep FEP_VAR
[lambdas]
0.500 0.500
[distance_restraints]
1 2 0 .1 .5 0
"""
    path = tmp_path/'original.inp'
    path.write_text(text)
    return path


def test_target_adaptation_is_explicit_and_keeps_shared_restraints(template, tmp_path):
    original = template.read_bytes()
    rendered = tmp_path/'new.inp'
    rendered.write_text(smoke.render(template, 18.54, True))
    sections = cp.sections(rendered)
    md = cp.keyed(sections['md'])
    assert md['steps'] == '20' and md['stepsize'] == '1.0'
    assert md['random_seed'] == '0' and md['constraint_algorithm'] == 'shake shake'
    assert md['separate_scaling'] == 'on'
    assert md['lrf'] == 'off'
    assert set(cp.keyed(sections['cut-offs']).values()) == {'99'}
    assert sections['distance_restraints'] == cp.sections(template)['distance_restraints']
    assert template.read_bytes() == original


@pytest.mark.parametrize('weight', [.0001, .5, .9999])
def test_interior_checks_render_requested_steps_and_weights(template, tmp_path, weight):
    output = tmp_path/'run.inp'
    output.write_text(smoke.render(template, 18.54, True, steps=2000, weight=weight))
    raw = cp.sections(output)
    assert cp.keyed(raw['md'])['steps'] == '2000'
    assert list(map(float, raw['lambdas'][0])) == pytest.approx([1-weight, weight])


@pytest.mark.parametrize('weight', [0., 1., -1., float('nan')])
def test_exact_endpoints_and_invalid_weights_rejected(template, weight):
    with pytest.raises(ValueError, match='exclude exact endpoints'):
        smoke.render(template, 20., True, weight=weight)


@pytest.mark.parametrize('steps', [0, 19, 21, 2001, 2000.])
def test_target_budget_cannot_silently_expand(template, steps):
    with pytest.raises(ValueError, match='20..2000'):
        smoke.render(template, 20., True, steps=steps)


def test_remove_softcore_preserves_other_tables(tmp_path):
    source = tmp_path/'old.fep'
    source.write_text('''[FEP]
states 2
softcore_use_max_potential on
softcore_method gapsys
[atoms]
1 10
[change_charges]
1 1.0 0.0
[atom_types]
X 1 2 3 4 5 6
[softcore]
1 20 20
[change_atoms]
1 X X
''')
    before = source.read_bytes()
    output = tmp_path/'new.fep'
    output.write_text(smoke.without_softcore(source))
    original, converted = cp.sections(source), cp.sections(output)
    assert 'softcore' not in converted
    assert cp.keyed(converted['fep']) == {
        'states': '2', 'softcore_use_max_potential': 'off', 'softcore_method': 'standard'}
    for key in original.keys()-{'fep', 'softcore'}:
        assert converted[key] == original[key]
    assert source.read_bytes() == before
    with pytest.raises(ValueError, match='Unsupported target FEP dialect'):
        smoke.without_softcore(output)


def test_no_softcore_matrix_covers_both_targets_legs_directions_and_weights(tmp_path, monkeypatch):
    import json
    calls = []
    monkeypatch.setattr(smoke, 'stage_case', lambda *args, **kwargs: calls.append((args, kwargs)))
    destination = tmp_path/'matrix'
    smoke.stage(tmp_path/'references', destination, no_softcore=True,
                weights=(.0001, .5, .9999), steps=2000)
    matrix = json.loads((destination/'matrix.json').read_text())
    assert len(calls) == len(set(matrix['cases'])) == 24
    assert matrix['total_steps'] == 96000 and matrix['aggregate_ps'] == 96.
    identities = {(args[3]['target'], args[3]['leg'], args[3]['direction'], kwargs['weight'])
                  for args, kwargs in calls}
    assert identities == {(target, leg, direction, weight)
                          for target in ('cmet', 'eg5') for leg in ('1.water', '2.protein')
                          for direction in ('fwd', 'rev') for weight in (.0001, .5, .9999)}
    assert all(kwargs['no_softcore'] is True and kwargs['steps'] == 2000 for _, kwargs in calls)


@pytest.mark.parametrize('extra', ['[unknown]\nkey value\n', '[atom_restraints]\n1 0 0 0 1 1 1 0\n'])
def test_unreviewed_target_options_are_rejected(template, extra):
    template.write_text(template.read_text()+extra)
    with pytest.raises(ValueError, match='Unsupported target MD'):
        smoke.render(template, 20., True)


def test_transferred_archive_is_validated_and_builds(isolated_build, tmp_path):
    old = charge_build.validate(isolated_build)
    archive = isolated_build.parent/'native-source.tar'
    with pytest.raises(ValueError, match='archive hash'):
        charge_build.build_archive(archive, tmp_path/'bad-hash', 'gfortran-11',
                                   expected_commit=old['source_commit'], expected_sha256='0'*64)
    with pytest.raises(ValueError, match='commit differs'):
        charge_build.build_archive(archive, tmp_path/'bad-commit', 'gfortran-11',
                                   expected_commit='0'*40, expected_sha256=cp.fingerprint(archive))
    assert not (tmp_path/'bad-hash').exists() and not (tmp_path/'bad-commit').exists()
    report = charge_build.build_archive(archive, tmp_path/'build', old['compiler']['binary'],
                                        expected_commit=old['source_commit'], expected_sha256=cp.fingerprint(archive))
    assert charge_build.validate(tmp_path/'build/build.json') == report
    assert report['source_files_sha256'] == old['source_files_sha256']
