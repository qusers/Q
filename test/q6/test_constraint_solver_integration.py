"""Short, matched-restart solver checks; not an equilibrium comparison."""
import json
import shutil
import struct
import subprocess

import numpy as np
import pytest

from QligFEP import boundary_native as bn, charge_completion as completion
from QligFEP import charge_probe as probe, charge_protocol as cp
from test_charge_probe import QDYN, QPREP


def test_solver_selector_is_the_only_two_value_setting():
    assert cp.keyed([['constraint_algorithm', 'shake', 'settle']]) == {
        'constraint_algorithm': 'shake settle'}
    for rows in ([['constraint_algorithm', 'shake']],
                 [['constraint_algorithm', 'shake', 'settle', 'extra']],
                 [['steps', '10', '20']],
                 [['constraint_algorithm', 'shake', 'shake']]*2):
        with pytest.raises(ValueError, match='Duplicate or non-scalar'):
            cp.keyed(rows)


@pytest.fixture(scope='module')
def common_start(tmp_path_factory):
    root = tmp_path_factory.mktemp('constraint-comparison')
    prepared = root/'prepared'
    probe.prepare(prepared, QPREP, 14., 758982)
    probe.timing(prepared, root/'seed', QDYN, sign=1, weight=0., steps=20, seed=123)
    return root


@pytest.fixture(scope='module', params=[
    (sign, solver) for sign in (-1, 1)
    for solver in ('shake shake', 'shake settle', 'lincs lincs')
])
def completed_solver(request, common_start):
    sign, solver = request.param
    run = common_start/f'{sign}-{solver.replace(" ", "-")}'
    run.mkdir()
    shutil.copyfile(common_start/'prepared/system.top', run/'system.top')
    shutil.copyfile(common_start/'prepared'/('positive.fep' if sign == 1 else 'negative.fep'), run/'charge.fep')
    shutil.copyfile(common_start/'seed/final.re', run/'start.re')
    radius = json.loads((common_start/'prepared/prepared.json').read_text())['effective_radius_angstrom']
    text = probe.md_input(radius, 100, 1., 0, 1.).replace(
        'constraint_algorithm shake shake', 'constraint_algorithm '+solver)
    text = text.replace('[files]\n', '[files]\nrestart start.re\n')
    path = run/'run.inp'
    path.write_text(text)
    window = cp.inspect_window(path, 'integrated')
    result = subprocess.run([str(QDYN), 'run.inp'], cwd=run,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr
    log = run/'native.log'
    log.write_text(result.stdout)
    report = completion.validate_window(window, 'integrated', log)
    return window, log, report, solver


def test_solvers_preserve_geometry_and_charge_accounting(completed_solver):
    window, log, report, solver = completed_solver
    assert report['native_initialization']['native']['constraint_algorithms'] == solver.split()
    assert report['saved_frames'] == 9
    assert report['maximum_accounting_residual_kcal_mol'] < 1e-8
    assert report['trajectory_diagnostics']['water_distance_extrema_drift_angstrom'] < .005
    assert report['production_ready'] is False
    data = (log.parent/'final.re').read_bytes()
    count = struct.unpack('<i', data[4:8])[0]
    xyz = np.frombuffer(data[8:8+8*count], dtype='<f8').reshape(-1, 3)
    water = xyz[1:].reshape(-1, 3, 3)
    residuals = []
    for left, right, target in ((0, 1, .9572), (0, 2, .9572), (1, 2, 1.5136)):
        squared = np.sum((water[:, left]-water[:, right])**2, axis=1)
        residuals.extend(abs(squared-target**2)/target**2)
    tolerance = 1e-10 if solver.endswith('settle') else 1e-4
    assert max(residuals) <= tolerance
    assert window['signature']['velocity_initialization'] == 'restart'


@pytest.mark.parametrize('mutation', ['missing', 'different', 'duplicate'])
def test_solver_provenance_cannot_be_removed_or_relabelled(completed_solver, tmp_path, mutation):
    window, log, _, solver = completed_solver
    line = 'Constraint algorithms = '+solver.upper()
    text = log.read_text()
    assert text.count(line) == 1
    if mutation == 'missing':
        text = text.replace(line, '')
    elif mutation == 'duplicate':
        text = text.replace(line, line+'\n'+line)
    else:
        replacement = 'SHAKE SETTLE' if solver == 'shake shake' else 'SHAKE SHAKE'
        text = text.replace(line, 'Constraint algorithms = '+replacement)
    changed = tmp_path/'native.log'
    changed.write_text(text)
    with pytest.raises(ValueError, match='constraint algorithms'):
        bn.validate_window(window, 'integrated', changed)
