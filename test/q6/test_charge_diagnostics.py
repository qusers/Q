"""Native observation coverage, independent geometry checks and no dynamics change."""
import copy
import shutil
import struct

import numpy as np
import pytest

from QligFEP import charge_build, charge_diagnostics as diag, charge_probe as probe
from test_charge_probe import prepared, QDYN, PROJECT_ROOT


@pytest.fixture(scope='module', params=[-1, 1])
def observed(request, prepared, tmp_path_factory):
    run = tmp_path_factory.mktemp('charge-observations')/'run'
    report = probe.timing(prepared, run, QDYN, sign=request.param, weight=1., steps=100)
    window = {'signature': {'md': {'steps': '100'}, 'intervals': {'output': '10'}}}
    return run, report, window


def test_native_trace_covers_every_force_geometry(observed):
    run, report, window = observed
    result = diag.assess(window, report['native'], run/'native.log')
    assert result['force_geometries_observed'] == 101
    assert result['snapshot_steps'] == list(range(0, 101, 10))
    assert result['cutoff_margin_angstrom'] > 70
    assert result['water_distance_extrema_drift_angstrom'] < .005
    assert result['production_ready'] is False
    assert result['total_temperature_all_evaluation_range_kelvin'][0] > 0


def test_last_snapshot_geometry_matches_independent_restart_calculation(observed):
    run, report, _ = observed
    native = report['native']
    data = (run/'final.re').read_bytes()
    count = struct.unpack('<i', data[4:8])[0]
    xyz = np.frombuffer(data[8:8+8*count], dtype='<f8').reshape(-1, 3)
    center = np.array(native['center'])
    oxygen = xyz[1::3]-center
    h1, h2 = xyz[2::3]-xyz[1::3], xyz[3::3]-xyz[1::3]
    radius = np.linalg.norm(oxygen, axis=1)
    dipole = h1+h2
    cosine = np.einsum('ij,ij->i', dipole, oxygen)/(np.linalg.norm(dipole, axis=1)*radius)
    final = diag.trace(run/'native.log', steps=100, interval=10, waters=146, shells=len(native['shell']))[-1]
    edges = np.linspace(0, native['parameters'][0], 21)
    counts = np.histogram(radius[radius < edges[-1]], edges)[0].tolist()+[int(sum(radius >= edges[-1]))]
    assert final['density'] == counts
    assert final['values'][1] == pytest.approx(np.linalg.norm(xyz[0]-center), abs=1e-12)
    assert final['values'][0] >= np.linalg.norm(xyz-center, axis=1).max()-1e-12
    for s, shell in enumerate(native['shell']):
        # Shells are indexed outermost first. The outer shell includes escaped
        # waters; inner boundaries have the engine's strict lower inequality.
        lower = float(np.float32(shell[1])-np.float32(shell[2])) if s == len(native['shell'])-1 else native['shell'][s+1][1]
        mask = radius > lower
        if s:
            mask &= radius <= shell[1]
        assert final['shells'][s][0] == int(mask.sum())
        assert final['shells'][s][1] == pytest.approx(cosine[mask].sum(), abs=1e-10)
        assert final['shells'][s][2] == pytest.approx((cosine[mask]**2).sum(), abs=1e-10)


@pytest.mark.parametrize('mutation,message', [
    ('missing', 'snapshots'), ('duplicate', 'snapshots'), ('density', 'every water'),
    ('invalid', 'Invalid/nonfinite'), ('evaluation', 'Invalid/nonfinite'),
    ('shell', 'orientation moments'), ('version', 'Unsupported'), ('geometry', '0.005'),
    ('hot_atom', 'velocity reset'),
])
def test_damaged_or_failed_trace_is_not_success(observed, tmp_path, mutation, message):
    run, report, window = observed
    lines = (run/'native.log').read_text().splitlines()
    indices = [i for i, line in enumerate(lines) if line.startswith('Q_CHARGE_TRACE_V1')]
    last = indices[-1]
    if mutation == 'hot_atom':
        lines.append('>>> WARNING: hot atom, i = 1 Temp(i)= 100000.0')
    elif mutation == 'missing':
        lines = lines[:last]  # preserve an incomplete log, not a recovered pass
    elif mutation == 'duplicate':
        lines.append(lines[last])
    elif mutation == 'version':
        lines[last] = lines[last].replace('V1', 'V2')
    elif mutation in ('invalid', 'evaluation', 'geometry'):
        tokens = lines[last].split()
        if mutation == 'invalid':
            tokens[3] = '1'
        elif mutation == 'evaluation':
            tokens[2] = '100'
        else:
            tokens[8] = str(float(tokens[8])+.01)  # running OH maximum
        lines[last] = ' '.join(tokens)
    else:
        index = next(i for i in range(last+1, len(lines)) if lines[i].startswith('QCT_'+('DENSITY' if mutation == 'density' else 'SHELL')))
        tokens = lines[index].split()
        tokens[2 if mutation == 'density' else 4] = '10000'
        lines[index] = ' '.join(tokens)
    log = tmp_path/'damaged.log'
    log.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError, match=message):
        diag.assess(window, report['native'], log)


def test_global_radius_bound_must_fit_every_active_cutoff(observed):
    run, report, window = observed
    native = copy.deepcopy(report['native'])
    native['cutoffs'][2] = 10.
    with pytest.raises(ValueError, match='untruncated pair coverage'):
        diag.assess(window, native, run/'native.log')


@pytest.fixture(scope='module')
def pre_observation_binary(tmp_path_factory):
    compiler = shutil.which('gfortran-11') or shutil.which('gfortran')
    directory = tmp_path_factory.mktemp('prior-native')/'build'
    charge_build.build(PROJECT_ROOT, directory, compiler, commit='4ad905d1')
    report = charge_build.validate(directory/'build.json')
    return directory/report['binaries']['qdyn']['path']


def test_observation_does_not_change_energies_or_restart(observed, prepared, pre_observation_binary, tmp_path):
    run, report, _ = observed
    old = tmp_path/'old'
    probe.timing(prepared, old, pre_observation_binary, sign=report['sign'], weight=1., steps=100)
    assert 'Q_CHARGE_TRACE' not in (old/'native.log').read_text()
    for name in ('states.en', 'final.re'):
        assert (run/name).read_bytes() == (old/name).read_bytes()
