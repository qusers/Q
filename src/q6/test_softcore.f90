program test_softcore
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qsoftcore
  implicit none

  integer :: failures
  real(8) :: sigma, r_lj, r_q, r2, energy, coeff, energy_p, energy_m, coeff_dummy
  real(8) :: c12, c6, qq, h, numeric_dvdr, analytic_dvdr, delta_energy
  logical :: overlap, overlap_dummy

  failures = 0

  ! Gapsys paper carbon/TIP3P oxygen example, converted from nm to angstrom.
  sigma = 0.5_8*(3.0_8 + 3.15061_8)
  r_lj = gapsys_lj_radius(GAPSYS_ALPHA_LJ_DEFAULT, sigma, 0.5_8)
  call assert_close('paper LJ radius (A)', r_lj, 2.898108709891854_8, 1.0e-12_8, failures)

  r_q = gapsys_coulomb_radius(GAPSYS_ALPHA_Q_ANGSTROM_DEFAULT, &
                              GAPSYS_SIGMA_Q_DEFAULT, 0.25_8, -0.834_8, 0.5_8)
  call assert_close('paper Coulomb radius (A)', r_q, 3.229953302617800_8, 1.0e-12_8, failures)

  ! A representative LJ pair: epsilon=0.5 kcal/mol and sigma=3 A.
  c6 = 4.0_8*0.5_8*3.0_8**6
  c12 = 4.0_8*0.5_8*3.0_8**12
  r_lj = gapsys_lj_radius(GAPSYS_ALPHA_LJ_DEFAULT, 3.0_8, 0.5_8)

  call eval_gapsys_lj(c12, c6, r_lj, 0.0_8, energy, coeff, overlap)
  call assert_true('LJ exact overlap detected', overlap, failures)
  call assert_true('LJ exact-overlap energy finite', ieee_is_finite(energy), failures)
  call assert_close('LJ exact-overlap Cartesian coefficient', coeff, 0.0_8, 0.0_8, failures)

  ! Continuity at the switching point.
  h = 1.0e-7_8
  call eval_gapsys_lj(c12, c6, r_lj, (r_lj-h)**2, energy_m, coeff_dummy, overlap_dummy)
  call eval_gapsys_lj(c12, c6, r_lj, (r_lj+h)**2, energy_p, coeff_dummy, overlap_dummy)
  call assert_close('LJ potential continuity', energy_m, energy_p, 2.0e-5_8, failures)

  ! Force/energy consistency inside the soft region.
  r2 = (0.6_8*r_lj)**2
  call eval_gapsys_lj(c12, c6, r_lj, r2, energy, coeff, overlap)
  h = 1.0e-5_8
  call eval_gapsys_lj(c12, c6, r_lj, (sqrt(r2)+h)**2, energy_p, coeff_dummy, overlap_dummy)
  call eval_gapsys_lj(c12, c6, r_lj, (sqrt(r2)-h)**2, energy_m, coeff_dummy, overlap_dummy)
  numeric_dvdr = (energy_p-energy_m)/(2.0_8*h)
  analytic_dvdr = coeff*sqrt(r2)
  call assert_close('LJ finite-difference force', analytic_dvdr, numeric_dvdr, 1.0e-5_8, failures)

  ! Coulomb uses energy-scaled qq but a radius selected independently above.
  qq = -20.0_8
  call eval_gapsys_coulomb(qq, r_q, 0.0_8, energy, coeff, overlap)
  call assert_true('Coulomb exact overlap detected', overlap, failures)
  call assert_close('Coulomb exact-overlap energy', energy, 3.0_8*qq/r_q, 1.0e-12_8, failures)
  call assert_close('Coulomb exact-overlap Cartesian coefficient', coeff, 0.0_8, 0.0_8, failures)

  r2 = (0.6_8*r_q)**2
  call eval_gapsys_coulomb(qq, r_q, r2, energy, coeff, overlap)
  h = 1.0e-5_8
  call eval_gapsys_coulomb(qq, r_q, (sqrt(r2)+h)**2, energy_p, coeff_dummy, overlap_dummy)
  call eval_gapsys_coulomb(qq, r_q, (sqrt(r2)-h)**2, energy_m, coeff_dummy, overlap_dummy)
  numeric_dvdr = (energy_p-energy_m)/(2.0_8*h)
  analytic_dvdr = coeff*sqrt(r2)
  call assert_close('Coulomb finite-difference force', analytic_dvdr, numeric_dvdr, 1.0e-8_8, failures)

  ! Coupled endpoint is exactly hard-core away from overlap.
  call eval_gapsys_coulomb(qq, 0.0_8, 4.0_8, energy, coeff, overlap)
  call assert_close('Coulomb coupled endpoint energy', energy, qq/2.0_8, 1.0e-12_8, failures)
  call eval_gapsys_lj(c12, c6, 0.0_8, 4.0_8, energy, coeff, overlap)
  call assert_close('LJ coupled endpoint energy', energy, c12/2.0_8**12-c6/2.0_8**6, &
                    1.0e-8_8, failures)

  ! Energy-only shape corrections agree with two complete pair evaluations.
  r2 = (0.6_8*r_lj)**2
  call eval_gapsys_lj(c12,c6,0.8_8*r_lj,r2,energy_m,coeff_dummy,overlap_dummy)
  call eval_gapsys_lj(c12,c6,r_lj,r2,energy_p,coeff_dummy,overlap_dummy)
  delta_energy = gapsys_lj_energy_difference(c12,c6,0.8_8*r_lj,r_lj,r2)
  call assert_close('LJ shape energy difference', delta_energy, energy_p-energy_m, &
                    1.0e-10_8, failures)

  call eval_gapsys_coulomb(qq,0.8_8*r_q,r2,energy_m,coeff_dummy,overlap_dummy)
  call eval_gapsys_coulomb(qq,r_q,r2,energy_p,coeff_dummy,overlap_dummy)
  delta_energy = gapsys_coulomb_energy_difference(qq,0.8_8*r_q,r_q,r2)
  call assert_close('Coulomb shape energy difference', delta_energy, energy_p-energy_m, &
                    1.0e-12_8, failures)

  r2 = (1.1_8*r_lj)**2
  call assert_close('LJ hard/hard correction is exactly zero', &
                    gapsys_lj_energy_difference(c12,c6,0.8_8*r_lj,r_lj,r2), &
                    0.0_8, 0.0_8, failures)
  call assert_close('Coulomb hard/hard correction is exactly zero', &
                    gapsys_coulomb_energy_difference(qq,0.8_8*r_lj,r_lj,r2), &
                    0.0_8, 0.0_8, failures)

  ! An almost fully decoupled atom can approach a solvent atom arbitrarily
  ! closely. Work for the final weight 0.0001 -> 0 must stay bounded.
  r_lj = gapsys_lj_radius(GAPSYS_ALPHA_LJ_DEFAULT, 3.0_8, 0.9999_8)
  call eval_gapsys_lj(c12,c6,r_lj,1.0e-12_8,energy,coeff,overlap)
  call assert_true('near-overlap final-window work is finite', &
                    ieee_is_finite(0.0001_8*energy), failures)
  call assert_true('near-overlap final-window work stays bounded', &
                    abs(0.0001_8*energy) < 1.0_8, failures)

  if(failures /= 0) then
    write(*,'(a,i0)') 'soft-core tests failed: ', failures
    error stop 1
  end if
  write(*,'(a)') 'soft-core pair tests passed'

contains

subroutine assert_close(name, actual, expected, tolerance, nfail)
  character(*), intent(in) :: name
  real(8), intent(in) :: actual, expected, tolerance
  integer, intent(inout) :: nfail
  if(.not. ieee_is_finite(actual) .or. abs(actual-expected) > tolerance) then
    write(*,'(a,2(1x,es24.16),a,es12.4)') 'FAIL '//trim(name)//':', actual, expected, &
                                         ' tolerance=', tolerance
    nfail = nfail+1
  end if
end subroutine assert_close

subroutine assert_true(name, condition, nfail)
  character(*), intent(in) :: name
  logical, intent(in) :: condition
  integer, intent(inout) :: nfail
  if(.not. condition) then
    write(*,'(a)') 'FAIL '//trim(name)
    nfail = nfail+1
  end if
end subroutine assert_true

end program test_softcore
