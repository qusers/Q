"""Reproduce the current Q-label dependence without changing any physical charges."""
import json

import numpy as np
import pytest

from test_state_energy_audit import executable, run_audit
from test_state_energy_audit import (
    test_pure_states_do_not_change_with_lambda as assert_pure_states,
    test_full_energy_and_force_follow_same_lambda_mixture as assert_mixture,
    test_boundary_terms_appear_once_in_pure_and_total_energies as assert_boundary_terms,
    test_saved_energy_records_match_audited_pure_states as assert_saved_states,
)


@pytest.fixture(scope='module', params=[-1, 1])
def comparison(request, executable, tmp_path_factory):
    directory = tmp_path_factory.mktemp('partition-audit')
    first = run_audit(request.param, executable, directory/'original')
    second = run_audit(request.param, executable, directory/'promoted', fixed_charge_as_q=True)
    return first, second


@pytest.fixture(scope='module', params=[-1, 1])
def general_comparison(request, executable, tmp_path_factory):
    directory = tmp_path_factory.mktemp('general-water-partition')
    first = run_audit(request.param, executable, directory/'original', general_water=True)
    second = run_audit(request.param, executable, directory/'promoted',
                       fixed_charge_as_q=True, general_water=True)
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


def test_nonangular_difference_matches_omitted_h_lj_and_coulomb_rounding(comparison):
    first, second = comparison
    c0, c1 = map(_components, comparison)
    logs = [(run[4]/'audit.log').read_text().splitlines() for run in comparison]
    omitted = [float(line.split()[1]) for line in logs[1] if line.startswith('PARTITION_H_LJ ')]
    assert len(omitted) == 1 and omitted[0] < 0
    rounding = [np.array([float(line.split()[2]) for line in log
                         if line.startswith('PARTITION_COULOMB_ROUNDING ')]) for log in logs]
    assert all(row.shape == (2,) for row in rounding)
    for mode in (2, 3):
        delta = c1[mode]-c0[mode]
        np.testing.assert_allclose(delta[:, 2], -omitted[0], atol=1e-11, rtol=0)
        for index, row in enumerate(first[2][mode]):
            weight = row[2]
            predicted = np.array([weight, 1-weight])@(rounding[1]-rounding[0])
            assert delta[index, 1] == pytest.approx(predicted, abs=1e-11, rel=0)


@pytest.mark.xfail(strict=True, reason='Archived water type omits H LJ in Q-water; also known charge-product rounding')
def test_full_nonangular_potential_invariant(comparison):
    first, second = comparison
    for mode in (2, 3):  # Born-only and neither; shared physical configurations/interactions.
        np.testing.assert_allclose(first[2][mode, :, 3], second[2][mode, :, 3], atol=1e-8, rtol=0)


def test_general_water_partition_has_only_identified_charge_rounding(general_comparison):
    first, second = general_comparison
    c0, c1 = map(_components, general_comparison)
    rounding = [np.array([float(line.split()[2]) for line in (run[4]/'audit.log').read_text().splitlines()
                         if line.startswith('PARTITION_COULOMB_ROUNDING ')]) for run in general_comparison]
    for mode in (2, 3):
        delta = c1[mode]-c0[mode]
        np.testing.assert_allclose(delta[:, 2], 0., atol=1e-11, rtol=0)
        for index, row in enumerate(first[2][mode]):
            predicted = np.array([row[2], 1-row[2]])@(rounding[1]-rounding[0])
            assert delta[index, 1] == pytest.approx(predicted, abs=1e-11, rel=0)
            assert second[2][mode, index, 3]-row[3] == pytest.approx(predicted, abs=1e-8, rel=0)
    # Removing the non-angular model inconsistency does not fix the distinct
    # angular Q-label dependence. Keep it visible, without changing the target.
    angular = second[2][0, 0, 8:10]-first[2][0, 0, 8:10]
    assert abs(angular[1]-angular[0]) > 1e-4


def test_general_water_native_gradient_matches_energy_derivative(general_comparison):
    for run in general_comparison:
        rows = np.array([[float(v) for v in line.split()[1:]] for line in (run[4]/'audit.log').read_text().splitlines()
                         if line.startswith('PARTITION_FORCE_FD ')])
        assert rows.shape == (3, 3) and np.isfinite(rows).all()
        np.testing.assert_allclose(rows[:, 1], rows[:, 2], atol=1e-7, rtol=0)


def test_general_water_retains_state_energy_bookkeeping(general_comparison):
    for run in general_comparison:
        assert_pure_states(run)
        assert_mixture(run)
        assert_saved_states(run)
    assert_boundary_terms(general_comparison[0])


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
