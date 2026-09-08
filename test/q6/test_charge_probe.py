"""Fresh Qprep probe inputs and bounded existing-MD endpoint checks."""
import json
from pathlib import Path
import struct

import pytest

from QligFEP import charge_probe as probe
from QligFEP.charge_protocol import fingerprint
from test_born_serialization import PROJECT_ROOT, QDYN

QPREP = PROJECT_ROOT/'src/q6/qprep'


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    assert QPREP.is_file(), 'Build serial Qprep before running native probe tests'
    assert QDYN.is_file(), 'Build serial Qdyn before running native probe tests'
    assert QDYN.stat().st_mtime >= max(p.stat().st_mtime for p in (PROJECT_ROOT/'src/q6').glob('*.f90')), 'Rebuild stale Qdyn'
    directory = tmp_path_factory.mktemp('probe')/'r10'
    report = probe.prepare(directory, QPREP, 10., 758971)
    assert report['production_ready'] is False
    assert (report['atoms'], report['waters']) == (439, 146)
    assert report['effective_radius_angstrom'] == 10.15
    assert report['effective_radius_angstrom'] != report['grid_radius_angstrom']
    return directory


def test_preparation_preserves_forcefield_and_rejects_overwrite(prepared):
    for suffix in ('lib', 'prm'):
        assert fingerprint(prepared/f'water.{suffix}') == fingerprint(PROJECT_ROOT/f'src/QligFEP/FF/AMBER14sb.{suffix}')
    before = fingerprint(prepared/'prepared.json')
    with pytest.raises(FileExistsError):
        probe.prepare(prepared, QPREP, 10., 758971)
    assert fingerprint(prepared/'prepared.json') == before


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('weight', [0., .5, 1.])
def test_native_probe_has_finite_real_lj_endpoints(prepared, tmp_path, sign, weight):
    result = probe.timing(prepared, tmp_path/'run', QDYN, sign=sign, weight=weight, steps=100)
    assert result['production_ready'] is False
    assert result['saved_frames'] == 9
    assert result['native']['meta'][-2:] == [2, 0]  # arithmetic LJ, optimized water compatible with zero H LJ
    assert result['native']['parameters'][2:4] == [332., 80.]
    qatom = result['native']['qatom'][0]
    assert qatom[4] == 12.01
    assert qatom[5:8] == [1.908, 0, 1.908]
    assert qatom[8] == pytest.approx(.1094**.5)
    assert qatom[-2:] == [0, sign]
    assert all(site[4:] == [0]*6 for site in result['native']['water_atom'][1:])
    with pytest.raises(FileExistsError):
        probe.timing(prepared, tmp_path/'run', QDYN, sign=sign, weight=weight, steps=100)


def test_prepared_input_tampering_is_rejected(prepared, tmp_path):
    import shutil
    copy = tmp_path/'copy'
    shutil.copytree(prepared, copy)
    with (copy/'negative.fep').open('a') as stream:
        stream.write('[softcore]\n1 0 0\n')
    with pytest.raises(ValueError, match='Prepared asset changed'):
        probe.timing(copy, tmp_path/'run', QDYN, sign=-1, weight=1., steps=100)
    assert not (tmp_path/'run').exists()


@pytest.mark.parametrize('steps', [0, 19, 21, 2001, 1_000_000])
def test_timing_cannot_be_used_for_long_calculations(tmp_path, steps):
    with pytest.raises(ValueError, match='20..2000'):
        probe.timing(tmp_path, tmp_path/'run', QDYN, sign=1, weight=1., steps=steps)
    assert not (tmp_path/'run').exists()


def test_orientation_seed_reproducibility(prepared, tmp_path):
    same = tmp_path/'same'
    different = tmp_path/'different'
    probe.prepare(same, QPREP, 10., 758971)
    probe.prepare(different, QPREP, 10., 758972)
    assert (same/'system.pdb').read_bytes() == (prepared/'system.pdb').read_bytes()
    assert (different/'system.pdb').read_bytes() != (prepared/'system.pdb').read_bytes()
    # This only changes generated hydrogen directions, not an equilibrated ensemble.
    assert json.loads((different/'prepared.json').read_text())['waters'] == 146


@pytest.mark.parametrize('kind', ['truncated', 'wrong-state', 'nonfinite', 'wrong-weight', 'offdiagonal'])
def test_energy_record_checker_rejects_corruption(tmp_path, kind):
    def record(payload):
        marker = struct.pack('<i', len(payload))
        return marker+payload+marker
    states = [struct.pack('<i15d', i, .5, *([0.]*14)) for i in (1, 2)]
    if kind == 'wrong-state':
        states[1] = states[0]
    elif kind == 'nonfinite':
        states[0] = struct.pack('<i15d', 1, .5, float('nan'), *([0.]*13))
    elif kind == 'wrong-weight':
        states[0] = struct.pack('<i15d', 1, .9, *([0.]*14))
    data = b''.join(map(record, states))+record(b'BAD!' if kind == 'offdiagonal' else b'')
    (tmp_path/'states.en').write_bytes(data[:-1] if kind == 'truncated' else data)
    with pytest.raises(ValueError):
        probe._check_energies(tmp_path/'states.en', [.5, .5], 1)
