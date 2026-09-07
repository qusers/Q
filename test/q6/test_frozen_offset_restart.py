"""Native MD restart safety: frozen offsets must not silently reset to zero."""
from pathlib import Path
import shutil
import struct
import subprocess

import pytest

from test_born_serialization import DATA, PROJECT_ROOT, QDYN, _md_input, _read_fortran_records


def _record(payload):
    marker = struct.pack('=i', len(payload))
    return marker+payload+marker


def _run(directory, *, restart=None, adapt=False):
    directory.mkdir(exist_ok=True)
    shutil.copyfile(DATA/'topology/Na-benzene-water.top', directory/'system.top')
    (directory/'charge.fep').write_text('[FEP]\nstates 2\n[atoms]\n1 1\n'
                                       '[change_charges]\n1 0.0 1.0\n')
    inp = _md_input(Path('system.top'), Path('charge.fep'), Path('final.re'), .5)
    inp = inp.replace('steps 1\n', 'steps 2\n').replace('charge_correction off',
                      'charge_correction on\nperstate_polarization on\n'
                      f'polarization_adaptation {"on" if adapt else "off"}')
    if restart is not None:
        (directory/'start.re').write_bytes(restart)
        inp = inp.replace('[files]\n', '[files]\nrestart start.re\n')
    (directory/'run.inp').write_text(inp)
    return subprocess.run([str(QDYN), 'run.inp'], cwd=directory,
                          capture_output=True, text=True, timeout=30)


@pytest.fixture(scope='module')
def initial_restart(tmp_path_factory):
    assert QDYN.is_file(), 'Build serial Qdyn before running native restart tests'
    assert QDYN.stat().st_mtime >= max(p.stat().st_mtime for p in (PROJECT_ROOT/'src/q6').glob('*.f90')), 'Rebuild stale Qdyn'
    directory = tmp_path_factory.mktemp('offset-restart')
    result = _run(directory)
    assert result.returncode == 0, result.stdout+result.stderr
    records = _read_fortran_records(directory/'final.re')
    assert len(records) == 3
    assert struct.unpack('=i3f', records[2]) == (3, 0., 0., 0.)
    return records


def test_valid_nonzero_frozen_offsets_survive_native_restart(initial_restart, tmp_path):
    offsets = struct.pack('=i3f', 3, .004204615484923124, .0362846776843071, .0221804678440094)
    restart = b''.join(map(_record, initial_restart[:2]+[offsets]))
    result = _run(tmp_path, restart=restart)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'Loaded polarization restraint data from restart file.' in result.stdout
    assert _read_fortran_records(tmp_path/'final.re')[2] == offsets


@pytest.mark.parametrize('offsets,message', [
    (None, 'requires restart offsets with matching shell count'),
    (struct.pack('=i2f', 2, .1, .2), 'requires restart offsets with matching shell count'),
    (struct.pack('=i2f', 3, .1, .2), 'Incomplete polarization restart offsets'),
    (struct.pack('=i3f', 3, float('nan'), .1, .2), 'Nonfinite polarization restart offsets'),
    (struct.pack('=i3f', 3, .1, float('inf'), .2), 'Nonfinite polarization restart offsets'),
])
def test_bad_frozen_restart_stops_before_dynamics(initial_restart, tmp_path, offsets, message):
    records = initial_restart[:2]+([] if offsets is None else [offsets])
    result = _run(tmp_path, restart=b''.join(map(_record, records)))
    assert result.returncode != 0
    assert message in result.stdout+result.stderr
    assert (tmp_path/'audit.en').stat().st_size == 0


@pytest.mark.parametrize('offsets', [None, struct.pack('=i2f', 2, .1, .2)])
def test_adaptive_legacy_fallback_is_preserved(initial_restart, tmp_path, offsets):
    records = initial_restart[:2]+([] if offsets is None else [offsets])
    result = _run(tmp_path, restart=b''.join(map(_record, records)), adapt=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'WARNING: Failed to read polarization restraint data' in result.stdout
    assert struct.unpack('=i3f', _read_fortran_records(tmp_path/'final.re')[2]) == (3, 0., 0., 0.)
