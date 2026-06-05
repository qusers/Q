program test_perstate_born
  ! Regression test for the per-state finite-sphere Born self-energy helpers in
  ! src/q6/qatom.f90 (born_coefficient, born_self_energy), the basis of the
  ! perstate_born_correction fix for net-charge-changing FEP edges.
  !
  ! Asserts:
  !   - born_coefficient(k_e, eps, R) = k_e*(1-1/eps)/(2R)   (continuum Born C).
  !   - born_self_energy(Q, C) = -C*Q^2.
  !   - the per-state difference E(Q2)-E(Q1) = -C*(Q2^2 - Q1^2) (what enters ddG),
  !     with the correct sign for both positive and negative environment charge:
  !       * cmet-like 0->+1 edge, Q_env=+3: Q1=3,Q2=4 -> ddG_Born = -7C ~ -45.9
  !       * eg5-like  0->+1 edge, Q_env=-5: Q1=-5,Q2=-4 -> ddG_Born = +9C ~ +59.0
  !   - a neutral edge (Q1=Q2) gives EXACTLY zero (the correction is invisible).
  !
  ! Exit code: 0 = all pass, 1 = some failure.

  use qatom, only: born_coefficient, born_self_energy
  implicit none

  real(8), parameter :: ke  = 332.0637_8     ! kcal*A/mol/e^2
  real(8), parameter :: eps = 80.0_8         ! bulk water dielectric
  real(8), parameter :: R   = 25.0_8         ! sphere radius, A
  real(8), parameter :: C_expected = ke * (1.0_8 - 1.0_8/eps) / (2.0_8 * R)  ! ~6.5583

  integer :: failures = 0
  real(8) :: C

  C = born_coefficient(ke, eps, R)
  call expect_close('born_coefficient = k_e(1-1/eps)/(2R)', C, C_expected, 1.0e-9_8)
  call expect_close('continuum C ~ 6.5583', C, 6.55825815_8, 1.0e-5_8)

  ! born_self_energy(Q,C) = -C*Q^2
  call expect_close('E_born(Q=3)  = -C*9',  born_self_energy(3.0_8, C), -C*9.0_8,  1.0e-9_8)
  call expect_close('E_born(Q=4)  = -C*16', born_self_energy(4.0_8, C), -C*16.0_8, 1.0e-9_8)
  call expect_close('E_born(Q=-5) = -C*25', born_self_energy(-5.0_8, C), -C*25.0_8, 1.0e-9_8)

  ! per-state difference = what enters ddG = -C*(Q2^2 - Q1^2)
  ! cmet-like 0->+1 edge, Q_env=+3: should be ~ -45.9 (matches observed ~ -45/-46)
  call expect_close('ddG_Born cmet (Q1=3,Q2=4) = -7C ~ -45.9', &
       born_self_energy(4.0_8, C) - born_self_energy(3.0_8, C), -C*7.0_8, 1.0e-9_8)
  call expect_close('ddG_Born cmet magnitude ~ -45.9', &
       born_self_energy(4.0_8, C) - born_self_energy(3.0_8, C), -45.907807_8, 1.0e-4_8)

  ! eg5-like 0->+1 edge, Q_env=-5: sign FLIPS to positive (matches observed +70)
  call expect_close('ddG_Born eg5 (Q1=-5,Q2=-4) = +9C ~ +59.0', &
       born_self_energy(-4.0_8, C) - born_self_energy(-5.0_8, C), C*9.0_8, 1.0e-9_8)
  if (born_self_energy(-4.0_8, C) - born_self_energy(-5.0_8, C) <= 0.0_8) then
    print '(a)', 'FAIL: eg5-like edge ddG_Born should be POSITIVE (sign flip vs cmet)'
    failures = failures + 1
  end if

  ! neutral edge (Q1=Q2): EXACTLY zero, correction invisible
  call expect_exact('neutral edge (Q1=Q2) ddG_Born = 0', &
       born_self_energy(2.5_8, C) - born_self_energy(2.5_8, C), 0.0_8)

  ! general identity: E(Q2)-E(Q1) = -C*(Q2^2 - Q1^2) for arbitrary charges
  call expect_close('identity -C*(Q2^2-Q1^2)', &
       born_self_energy(1.7_8, C) - born_self_energy(-2.3_8, C), &
       -C*(1.7_8**2 - (-2.3_8)**2), 1.0e-9_8)

  if (failures == 0) then
    print '(a)', 'PASS: all per-state Born self-energy checks passed.'
    stop 0
  else
    print '(a,i0,a)', 'FAIL: ', failures, ' per-state Born check(s) failed.'
    stop 1
  end if

contains

  subroutine expect_close(label, got, expected, tol)
    character(*), intent(in) :: label
    real(8),      intent(in) :: got, expected, tol
    real(8) :: diff
    diff = abs(got - expected)
    if (diff > tol) then
      print '(a,a,a,es15.7,a,es15.7,a,es15.7)', 'FAIL: ', label, &
            '  got=', got, '  expected=', expected, '  diff=', diff
      failures = failures + 1
    else
      print '(a,a,a,es15.7)', 'PASS: ', label, '  value=', got
    end if
  end subroutine expect_close

  subroutine expect_exact(label, got, expected)
    character(*), intent(in) :: label
    real(8),      intent(in) :: got, expected
    if (got /= expected) then
      print '(a,a,a,es15.7,a,es15.7)', 'FAIL: ', label, '  got=', got, '  expected=', expected
      failures = failures + 1
    else
      print '(a,a)', 'PASS: ', label
    end if
  end subroutine expect_exact

end program test_perstate_born
