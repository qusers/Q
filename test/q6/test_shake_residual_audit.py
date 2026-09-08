"""Isolate the pre-existing shared-constraint convergence defect without dynamics."""
from pathlib import Path
import shutil
import subprocess

import pytest

from test_state_energy_audit import ROOT, OBJECTS


def test_native_shake_reports_all_ready_only_when_all_final_constraints_converge(tmp_path):
    compiler = shutil.which('gfortran-11') or shutil.which('gfortran')
    assert compiler, 'Constraint audit needs a Fortran compiler'
    objects = [ROOT/'src/q6'/f'{name}.o' for name in OBJECTS]
    for name, path in zip(OBJECTS, objects):
        assert path.is_file() and path.stat().st_mtime >= (ROOT/'src/q6'/f'{name}.f90').stat().st_mtime
    executable = tmp_path/'constraint-audit'
    subprocess.run([compiler, '-O2', '-fcheck=all', '-ffree-line-length-none', '-I', str(ROOT/'src/q6'),
                    str(Path(__file__).with_name('shake_residual_audit.f90')), *map(str, objects), '-o', str(executable)],
                   capture_output=True, text=True, check=True, timeout=60)
    result = subprocess.run([str(executable)], capture_output=True, text=True, check=True, timeout=10)
    fields = result.stdout.split()
    assert fields[0] == 'SHAKE_MAX' and len(fields) == 6
    _, sample, bond, ready, residual, tolerance = fields
    assert 1 <= int(sample) <= 64 and 1 <= int(bond) <= 3 and int(ready) == 1
    assert float(residual) <= float(tolerance), result.stdout


test_native_shake_reports_all_ready_only_when_all_final_constraints_converge = pytest.mark.xfail(
    strict=True, raises=AssertionError,
    reason='Native SHAKE keeps ready flags after later coupled constraints move shared atoms; repair requires scope approval'
)(test_native_shake_reports_all_ready_only_when_all_final_constraints_converge)
