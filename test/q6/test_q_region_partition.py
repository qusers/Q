"""Reproduce the current Q-label dependence without changing any physical charges."""
import json

import numpy as np
import pytest

from test_state_energy_audit import executable, run_audit


@pytest.fixture(scope='module', params=[-1, 1])
def comparison(request, executable, tmp_path_factory):
    directory = tmp_path_factory.mktemp('partition-audit')
    first = run_audit(request.param, executable, directory/'original')
    second = run_audit(request.param, executable, directory/'promoted', fixed_charge_as_q=True)
    return first, second


def test_same_enclosed_charges_and_born_constants(comparison):
    first, second = comparison
    _, meta0, rows0, _, _ = first
    _, meta1, rows1, _, _ = second
    np.testing.assert_array_equal(meta0[:4], meta1[:4])
    assert meta0[4]-meta1[4] == pytest.approx(1., abs=1e-12, rel=0)
    np.testing.assert_allclose(meta1[5:]-meta0[5:], 1., atol=1e-12, rtol=0)
    np.testing.assert_allclose(meta0[4]+meta0[5:], meta1[4]+meta1[5:], atol=1e-12, rtol=0)
    np.testing.assert_allclose(rows0[:, :, 10:12], rows1[:, :, 10:12], atol=1e-12, rtol=0)


def _components(audit):
    values = np.array([[float(v) for v in line.split()[1:]]
                       for line in (audit[4]/'audit.log').read_text().splitlines()
                       if line.startswith('AUDIT_COMPONENT ')])
    assert values.shape == (28, 10) and np.isfinite(values).all()
    return values.reshape(4, 7, 10)[:, :, 2:]


def test_nonangular_difference_is_accounted_for(comparison):
    first, second = comparison
    c0, c1 = map(_components, comparison)
    for audit, components in zip(comparison, (c0, c1)):
        np.testing.assert_allclose(components.sum(axis=2), audit[2][:, :, 3], atol=1e-8, rtol=0)
    delta = c1[2]-c0[2]  # Born-only control, across all weights.
    np.testing.assert_allclose(delta[:, [0, 3, 4, 5, 6, 7]], 0., atol=1e-10, rtol=0)
    np.testing.assert_allclose(second[2][2, :, 3]-first[2][2, :, 3],
                               delta[:, 1]+delta[:, 2], atol=1e-8, rtol=0)
    print('NONANGULAR_COMPONENT_DIFFERENCE', delta[0].tolist())


@pytest.mark.xfail(strict=True, reason='Known Q/non-Q nonbonded energy mismatch; see PARTITION_AUDIT.md')
def test_full_nonangular_potential_invariant(comparison):
    first, second = comparison
    for mode in (2, 3):  # Born-only and neither; shared physical configurations/interactions.
        np.testing.assert_allclose(first[2][mode, :, 3], second[2][mode, :, 3], atol=1e-8, rtol=0)


def test_reproduces_current_angular_partition_dependence(comparison):
    first, second = comparison
    sign, _, rows0, forces0, _ = first
    _, _, rows1, forces1, _ = second
    angular_change = rows1[0, :, 8:10]-rows0[0, :, 8:10]
    # Difference of within-partition controls isolates the angular contribution
    # without assuming the Q/non-Q nonbonded kernels are identical.
    isolated_energy = (rows1[0, :, 3]-rows1[2, :, 3])-(rows0[0, :, 3]-rows0[2, :, 3])
    isolated_gradient = (forces1[0]-forces1[2])-(forces0[0]-forces0[2])
    assert np.max(np.abs(angular_change)) > 1e-4
    assert np.max(np.abs(isolated_gradient)) > 1e-4
    for index in range(7):
        weight = rows0[0, index, 2]
        expected = np.array([weight, 1-weight])@angular_change[index]
        assert isolated_energy[index] == pytest.approx(expected, abs=1e-8, rel=0)
        np.testing.assert_allclose(isolated_gradient[index],
                                   weight*isolated_gradient[0]+(1-weight)*isolated_gradient[1],
                                   atol=1e-9, rtol=1e-11)
    gap_change = angular_change[0, 1]-angular_change[0, 0]
    assert abs(gap_change) > 1e-4
    print(json.dumps(dict(charge=sign, angular_state_changes=angular_change[0].tolist(),
                          state2_minus_state1_gap_change=float(gap_change),
                          maximum_angular_gradient_change=float(np.max(np.abs(isolated_gradient))),
                          maximum_nonangular_gradient_change=float(np.max(np.abs(forces1[2]-forces0[2])))), sort_keys=True))
    # Passing reproduces an observed defect/ambiguity. It is NOT a physical
    # partition-invariance pass or permission to change the production target.
