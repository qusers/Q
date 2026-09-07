"""Full-Q fixed-coordinate state bookkeeping; no dynamics or new sampler."""
from pathlib import Path
import shutil
import struct
import subprocess

import numpy as np
import pytest

from test_born_serialization import DATA, _md_input, _energy_state_records, _read_fortran_records

ROOT = Path(__file__).resolve().parents[2]
OBJECTS = ('md', 'boundary_corrections', 'mpiglob', 'qatom', 'sizes', 'trj',
           'topo', 'misc', 'nrgy', 'prmfile', 'index', 'mask')


@pytest.fixture(scope='module')
def executable(tmp_path_factory):
    compiler = shutil.which('gfortran-11') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Native state audit requires a Fortran compiler')
    build = tmp_path_factory.mktemp('state-audit-build')
    objects = [ROOT/'src/q6'/f'{name}.o' for name in OBJECTS]
    for name, path in zip(OBJECTS, objects):
        assert path.exists(), 'Build serial Qdyn before the native state audit'
        assert path.stat().st_mtime >= (ROOT/'src/q6'/f'{name}.f90').stat().st_mtime, 'Stale native object: '+name
    exe = build/'state-audit'
    subprocess.run([compiler, '-O2', '-fcheck=all', '-ffree-line-length-none', '-I', str(ROOT/'src/q6'),
                    str(Path(__file__).with_name('state_energy_audit.f90')), *map(str, objects), '-o', str(exe)],
                   cwd=build, capture_output=True, text=True, check=True, timeout=60)
    return exe


@pytest.fixture(scope='module', params=[-1, 1])
def audit(request, executable, tmp_path_factory):
    run = tmp_path_factory.mktemp('state-audit')
    shutil.copyfile(DATA/'topology/Na-benzene-water.top', run/'system.top')
    # Atom 1 is a real solute atom. No type/mass/LJ changes; a nonzero non-Q
    # environment tests the Born cross term rather than only an even q^2 case.
    (run/'audit.fep').write_text(f'[FEP]\nstates 2\n[atoms]\n1 1\n[change_charges]\n1 0.0 {request.param}.0\n')
    text = _md_input(Path('system.top'), Path('audit.fep'), Path('audit.re'), .5)
    text = text.replace('charge_correction off', 'charge_correction on\nperstate_polarization on\n'
                        'polarization_adaptation off\nperstate_born_correction on\nborn_dielectric 80')
    text = text.replace('\ntemperature 1\n', '\ntemperature 298\n')
    (run/'audit.inp').write_text(text)
    result = subprocess.run([str(executable), 'audit.inp'], cwd=run,
                            capture_output=True, text=True, check=True, timeout=45)
    metadata = [line for line in result.stdout.splitlines() if line.startswith('AUDIT_META ')]
    assert len(metadata) == 1
    meta = np.array([float(v) for v in metadata[0].split()[1:]])
    rows = np.array([[float(v) for v in line.split()[1:]] for line in result.stdout.splitlines()
                     if line.startswith('STATE_AUDIT ')])
    assert rows.shape == (28, 13) and np.isfinite(rows).all()
    assert [(int(r[0]), int(r[1])) for r in rows] == [(m, s) for m in range(1, 5) for s in range(1, 8)]
    records = _read_fortran_records(run/'gradients.bin')
    natom = struct.unpack('=i', records[0])[0]
    assert len(records) == 29
    gradients = []
    for row, raw in zip(rows, records[1:]):
        assert struct.unpack('=2i', raw[:8]) == (int(row[0]), int(row[1]))
        gradient = np.frombuffer(raw[8:], dtype='=f8')
        assert len(gradient) == 3*natom and np.isfinite(gradient).all()
        gradients.append(gradient)
    return request.param, meta, rows.reshape(4, 7, 13), np.array(gradients).reshape(4, 7, -1), run


def test_pure_states_do_not_change_with_lambda(audit):
    _, _, rows, _, _ = audit
    for mode in rows:
        for row in mode[1:]:
            # EQ totals, restraints, angular energies, Born constants (not mixture bucket).
            np.testing.assert_allclose(row[4:12], mode[0, 4:12], atol=1e-10, rtol=0)


def test_full_energy_and_force_follow_same_lambda_mixture(audit):
    _, _, rows, gradients, _ = audit
    for mode, forces in zip(rows, gradients):
        for row, gradient in zip(mode, forces):
            weight = row[2]
            assert row[3] == pytest.approx(weight*mode[0, 3]+(1-weight)*mode[1, 3], abs=1e-8, rel=0)
            np.testing.assert_allclose(gradient, weight*forces[0]+(1-weight)*forces[1], atol=1e-9, rtol=1e-11)
            assert mode[0, 3]-mode[1, 3] == pytest.approx(row[4]-row[5], abs=1e-8, rel=0)


def test_boundary_terms_appear_once_in_pure_and_total_energies(audit):
    sign, meta, rows, gradients, _ = audit
    radius, ke, eps, coefficient, environment, q0, q1 = meta
    assert q0 == 0 and q1 == sign
    assert environment == pytest.approx(1.115, abs=1e-6)
    assert coefficient == pytest.approx(ke*(1-1/eps)/(2*radius), rel=0, abs=1e-12)
    expected_born = -coefficient*(environment+np.array([q0, q1]))**2
    for setting in range(7):
        both, pol, born, neither = rows[:, setting]
        weight = both[2]
        np.testing.assert_allclose(both[10:12], expected_born, atol=1e-12, rtol=0)
        np.testing.assert_allclose(both[4:6]-pol[4:6], expected_born, atol=1e-10, rtol=0)
        np.testing.assert_allclose(born[4:6]-neither[4:6], expected_born, atol=1e-10, rtol=0)
        np.testing.assert_allclose(both[4:6]-born[4:6], both[8:10], atol=1e-10, rtol=0)
        np.testing.assert_allclose(pol[4:6]-neither[4:6], pol[8:10], atol=1e-10, rtol=0)
        mixture = np.array([weight, 1-weight])
        assert both[3]-pol[3] == pytest.approx(mixture@expected_born, abs=1e-8, rel=0)
        assert both[3]-born[3] == pytest.approx(mixture@both[8:10], abs=1e-8, rel=0)
        assert both[12] == pytest.approx(mixture@both[8:10], abs=1e-10, rel=0)
        np.testing.assert_array_equal(gradients[0, setting], gradients[1, setting])
        np.testing.assert_array_equal(gradients[2, setting], gradients[3, setting])
    assert np.max(abs(gradients[0]-gradients[2])) > 1e-4  # Polarization is not a force-free Born term.


def test_saved_energy_records_match_audited_pure_states(audit):
    _, _, rows, _, run = audit
    serialized = _energy_state_records(run/'audit.en')
    assert len(serialized) == 28
    for row, states in zip(rows.reshape(28, 13), serialized):
        assert states[1][0] == row[2] and states[2][0] == 1-row[2]
        np.testing.assert_array_equal([states[1][1], states[2][1]], row[4:6])
        np.testing.assert_array_equal([states[1][-1], states[2][-1]], row[6:8])
