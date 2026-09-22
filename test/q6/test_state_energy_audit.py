"""Full-Q fixed-coordinate state bookkeeping; no dynamics or new sampler."""
from pathlib import Path
import os
import shutil
import struct
import subprocess

import numpy as np
import pytest

from test_born_serialization import DATA, _md_input, _energy_state_records, _read_fortran_records

ROOT = Path(__file__).resolve().parents[2]
OBJECTS = ('md', 'boundary_corrections', 'lincs', 'settle', 'mpiglob', 'qatom', 'sizes', 'trj',
           'topo', 'misc', 'nrgy', 'prmfile', 'index', 'mask')


@pytest.fixture(scope='module')
def executable(tmp_path_factory):
    compiler = shutil.which('gfortran-11') or shutil.which('gfortran')
    if not compiler:
        pytest.skip('Native state audit requires a Fortran compiler')
    build = tmp_path_factory.mktemp('state-audit-build')
    native = Path(os.environ.get('Q_NATIVE_BUILD', ROOT/'src/q6')).resolve()
    objects = [native/f'{name}.o' for name in OBJECTS]
    for name, path in zip(OBJECTS, objects):
        assert path.exists(), 'Build serial Qdyn before the native state audit'
        assert path.stat().st_mtime >= (ROOT/'src/q6'/f'{name}.f90').stat().st_mtime, 'Stale native object: '+name
        assert (native/f'{name}.f90').read_bytes() == (ROOT/'src/q6'/f'{name}.f90').read_bytes(), 'Different native source: '+name
    exe = build/'state-audit'
    subprocess.run([compiler, '-O2', '-fcheck=all', '-ffree-line-length-none', '-I', str(native),
                    str(Path(__file__).with_name('state_energy_audit.f90')), *map(str, objects), '-o', str(exe)],
                   cwd=build, capture_output=True, text=True, check=True, timeout=60)
    return exe


@pytest.fixture(scope='module', params=[-1, 1])
def audit(request, executable, tmp_path_factory):
    run = tmp_path_factory.mktemp('state-audit')
    return run_audit(request.param, executable, run)


def run_audit(sign, executable, run, *, fixed_charge_as_q=False, general_water=False, position=False, smooth=False):
    run.mkdir(exist_ok=True)
    shutil.copyfile(DATA/'topology/Na-benzene-water.top', run/'system.top')
    if general_water:
        path = run/'system.top'
        marker = '       0 = solvent type (0=SPC,1=3-atom,2=general)'
        text = path.read_text()
        assert text.count(marker) == 1
        path.write_text(text.replace(marker, marker.replace('       0', '       1', 1)))
    # Atom 1 is a real solute atom. No type/mass/LJ changes; a nonzero non-Q
    # environment tests the Born cross term rather than only an even q^2 case.
    extra_atom = '2 13\n' if fixed_charge_as_q else ''
    extra_charge = '2 1.0 1.0\n' if fixed_charge_as_q else ''
    (run/'audit.fep').write_text(f'[FEP]\nstates 2\n[atoms]\n1 1\n{extra_atom}'
                               f'[change_charges]\n1 0.0 {sign}.0\n{extra_charge}')
    text = _md_input(Path('system.top'), Path('audit.fep'), Path('audit.re'), .5)
    text = text.replace('charge_correction off', 'charge_correction on\nperstate_polarization on\n'
                        'polarization_adaptation off\nperstate_born_correction on\nborn_dielectric 80')
    text = text.replace('\ntemperature 1\n', '\ntemperature 298\n')
    if smooth:
        text = text.replace('[solvent]', '[solvent]\nsmooth_polarization on')
    if position:
        text += '\n[atom_restraints]\n1 0.13 -0.27 0.41 10 20 30 0\n'
    (run/'audit.inp').write_text(text)
    result = subprocess.run([str(executable), 'audit.inp'], cwd=run,
                            capture_output=True, text=True, check=True, timeout=45)
    (run/'audit.log').write_text(result.stdout)
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
    return sign, meta, rows.reshape(4, 7, 13), np.array(gradients).reshape(4, 7, -1), run


def test_shared_position_energy_force_and_serialization(audit, executable, tmp_path):
    sign, _, baseline, gradients, base_run = audit
    _, _, restrained, restrained_gradients, run = run_audit(sign, executable, tmp_path, position=True)
    lines = (run/'system.top').read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if '= Total no. of atoms, no. of solute atoms.' in line)
    xyz = np.array(list(map(float, lines[start+1].split()[:3])))
    displacement = xyz-np.array([.13, -.27, .41])
    expected_gradient = np.array([10, 20, 30])*displacement
    expected_energy = .5*np.dot(displacement, expected_gradient)
    # Every control and lambda: whole-system U and both pure-state U/restraint
    # buckets gain exactly one identical term, not a lambda-scaled duplicate.
    for column in (3, 4, 5, 6, 7):
        np.testing.assert_allclose(restrained[:, :, column]-baseline[:, :, column], expected_energy, atol=1e-9, rtol=0)
    np.testing.assert_allclose(restrained[:, :, 8:]-baseline[:, :, 8:], 0, atol=1e-10, rtol=0)
    difference = restrained_gradients-gradients
    np.testing.assert_allclose(difference[:, :, :3], np.broadcast_to(expected_gradient, (4, 7, 3)), atol=1e-9, rtol=0)
    np.testing.assert_allclose(difference[:, :, 3:], 0, atol=1e-10, rtol=0)
    base_frames = _energy_state_records(base_run/'audit.en')
    frames = _energy_state_records(run/'audit.en')
    assert len(base_frames) == len(frames) == 28
    for before, after in zip(base_frames, frames):
        for state in (1, 2):
            # q_energies order: lambda, total, ... restraint (last).
            assert after[state][1]-before[state][1] == pytest.approx(expected_energy, abs=1e-9)
            assert after[state][-1]-before[state][-1] == pytest.approx(expected_energy, abs=1e-9)


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
