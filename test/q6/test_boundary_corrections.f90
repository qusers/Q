program test_boundary_corrections
  use boundary_corrections
  implicit none

  real(8), parameter :: pi = 3.1415926535897932384626433832795_8
  real(8), parameter :: ke = 332.0637_8
  real(8), parameter :: eps = 80.0_8
  real(8), parameter :: radius = 25.0_8
  real(8), parameter :: force_constant = 20.0_8
  real(8), parameter :: theta = 1.2_8
  real(8), parameter :: offset = 0.03_8
  real(8) :: c, t1, t2, e1, e2, mixed_gradient, finite_difference
  real(8) :: shift, shifted1, shifted2, angle_sample, target1, target2
  real(8) :: unit_strength, base_angle, required_shift(101), analytic_shift(101)
  real(8), parameter :: lambda1 = 0.35_8, lambda2 = 0.65_8, delta = 1.0e-6_8
  integer :: failures, sample, field

  failures = 0
  c = born_coefficient(ke, eps, radius)
  call expect_close('Born coefficient', c, 6.558258075_8, 1.0e-12_8)
  call expect_close('Born state difference', &
    born_self_energy(4.0_8, c) - born_self_energy(3.0_8, c), -7.0_8*c, 1.0e-12_8)
  call expect_close('Born neutral edge', &
    born_self_energy(-2.0_8, c) - born_self_energy(-2.0_8, c), 0.0_8, 0.0_8)

  call expect_close('zero-field target', polarization_target(pi/2.0_8, 0.0_8), pi/2.0_8, 1.0e-12_8)
  call expect_close('target lower clamp', polarization_target(pi/2.0_8, 2.0_8), 0.0_8, 0.0_8)
  call expect_close('target upper clamp', polarization_target(pi/2.0_8, -2.0_8), pi, 0.0_8)
  call expect_close('safe acos upper clamp', safe_acos(1.1_8), 0.0_8, 0.0_8)
  call expect_close('safe acos lower clamp', safe_acos(-1.1_8), pi, 0.0_8)

  t1 = polarization_target(pi/2.0_8, 0.10_8)
  t2 = polarization_target(pi/2.0_8, -0.08_8)
  e1 = polarization_energy(theta, t1, offset, force_constant)
  e2 = polarization_energy(theta, t2, offset, force_constant)
  mixed_gradient = lambda1 * polarization_gradient(theta, t1, offset, force_constant) + &
                   lambda2 * polarization_gradient(theta, t2, offset, force_constant)
  finite_difference = (&
    lambda1 * polarization_energy(theta + delta, t1, offset, force_constant) + &
    lambda2 * polarization_energy(theta + delta, t2, offset, force_constant) - &
    lambda1 * polarization_energy(theta - delta, t1, offset, force_constant) - &
    lambda2 * polarization_energy(theta - delta, t2, offset, force_constant)) / (2.0_8 * delta)
  call expect_close('mixed force is derivative of mixed endpoint energies', &
    mixed_gradient, finite_difference, 1.0e-8_8)
  call expect_close('endpoint energy remains distinct', e2 - e1, &
    polarization_energy(theta, t2, offset, force_constant) - &
    polarization_energy(theta, t1, offset, force_constant), 0.0_8)

  ! A common offset shift changes both states' energies and forces. The gap
  ! shift is independent of observed angle, but depends on the rank targets.
  ! Include target clamps and the sign of a historical protein-shell mismatch.
  shift = -0.0037077669985592365_8
  do field = -2, 2
    target1 = polarization_target(pi/2, real(field,8))
    target2 = polarization_target(pi/3, -real(field,8))
    do sample = 1, 9
      angle_sample = pi*sample/10
      e1 = polarization_energy(angle_sample, target1, offset, force_constant)
      e2 = polarization_energy(angle_sample, target2, offset, force_constant)
      shifted1 = polarization_energy(angle_sample, target1, offset+shift, force_constant)
      shifted2 = polarization_energy(angle_sample, target2, offset+shift, force_constant)
      call expect_close('offset changes single-state energy', shifted1-e1, &
        force_constant*shift*(angle_sample-target1+offset)+0.5_8*force_constant*shift**2, 1.e-12_8)
      call expect_close('offset changes pure-state gap', (shifted2-shifted1)-(e2-e1), &
        force_constant*shift*(target1-target2), 1.e-12_8)
      call expect_close('offset changes angular gradient', &
        polarization_gradient(angle_sample,target1,offset+shift,force_constant)- &
        polarization_gradient(angle_sample,target1,offset,force_constant), force_constant*shift, 1.e-12_8)
    end do
  end do

  ! A scalar shell offset cannot generally absorb a missing background charge:
  ! matching effective target t-a requires a rank-dependent offset change.
  ! This is an algebraic counterfactual, not a parameter fit or production target.
  unit_strength = polarization_strength(1._8, .98750_8, .0335_8, .489_8, 23.3_8)
  do field = -3, 3, 3
    do sample = 1, size(required_shift)
      base_angle = acos(1._8+(1._8-2._8*sample)/size(required_shift))
      target1 = polarization_target(base_angle, unit_strength)
      target2 = polarization_target(base_angle, (1+field)*unit_strength)
      required_shift(sample) = target1-target2
      analytic_shift(sample) = 1.5_8*field*unit_strength*sin(base_angle)
      call expect_close('background target shift is rank dependent', &
        required_shift(sample), analytic_shift(sample), 1.e-12_8)
    end do
    if (field == 0) then
      call expect_close('neutral background conventions agree', maxval(abs(required_shift)), 0._8, 0._8)
    else
      if (maxval(required_shift)-minval(required_shift) <= 1.e-3_8) then
        print '(a)', 'FAIL: nonzero background unexpectedly reducible to one scalar offset'
        failures = failures+1
      end if
    end if
  end do

  if (failures /= 0) stop 1
  print '(a)', 'PASS: spherical-boundary correction helpers'

contains

  subroutine expect_close(label, got, expected, tolerance)
    character(*), intent(in) :: label
    real(8), intent(in) :: got, expected, tolerance

    if (abs(got - expected) > tolerance) then
      print '(a,a,2(a,es16.8))', 'FAIL: ', label, ' got=', got, ' expected=', expected
      failures = failures + 1
    end if
  end subroutine expect_close

end program test_boundary_corrections
