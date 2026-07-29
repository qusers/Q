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
  real(8), parameter :: lambda1 = 0.35_8, lambda2 = 0.65_8, delta = 1.0e-6_8
  integer :: failures

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
