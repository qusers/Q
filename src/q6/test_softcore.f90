program test_softcore
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qsoftcore
  implicit none

  integer :: failures
  real(8) :: sigma, r_lj, r_q, r2, energy, coeff, energy_p, energy_m, coeff_dummy
  real(8) :: c12, c6, qq, h, numeric_dvdr, analytic_dvdr, delta_energy
  logical :: overlap, overlap_dummy

  failures = 0
  call test_manuscript_equations(failures)

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

subroutine test_manuscript_equations(nfail)
  integer, intent(inout) :: nfail
  real(8), parameter :: weights(5) = [0.0_8,0.0001_8,0.5_8,0.9999_8,1.0_8]
  real(8), parameter :: distances(4) = [0.03_8,0.7_8,3.0_8,10.0_8]
  real(8), parameter :: sig=3.0_8, eps=0.5_8, alpha=0.3_8, alpha_c=100.0_8, charge=-20.0_8
  real(8) :: w,r,rlj,rq,u,dv,ref,den,up,um,dummy,step,a,b,c,hard_dv
  real(8) :: a12,a6,rs,delta
  integer :: iw,ir
  logical :: overlap
  a12 = 4*eps*sig**12
  a6 = 4*eps*sig**6
  do iw=1,size(weights)
    w = weights(iw)
    rlj = beutler_lj_radius(alpha,sig,1-w)
    rq = beutler_coulomb_radius(alpha_c,1-w)
    call assert_close('Beutler LJ offset units', rlj**6,alpha*sig**6*(1-w)**2,1.0e-10_8,nfail)
    call assert_close('Beutler Coulomb offset units', rq**2,alpha_c*(1-w)**2,1.0e-12_8,nfail)
    do ir=1,size(distances)
      r = distances(ir)
      den = alpha*(1-w)**2+(r/sig)**6
      ref = 4*eps*(1/den**2-1/den)
      call eval_beutler_lj(a12,a6,rlj,r*r,u,dv,overlap)
      call assert_close('Beutler Eq 7 energy',u,ref,1.0e-12_8*max(1.0_8,abs(ref)),nfail)
      ref = 4*eps*(-2/den**3+1/den**2)*6*r**5/sig**6
      call assert_close('Beutler Eq 7 derivative',dv*r,ref,1.0e-12_8*max(1.0_8,abs(ref)),nfail)
      step = r*1.0e-5_8
      call eval_beutler_lj(a12,a6,rlj,(r+step)**2,up,dummy,overlap)
      call eval_beutler_lj(a12,a6,rlj,(r-step)**2,um,dummy,overlap)
      ref = (up-um)/(2*step)
      ! Include subtraction roundoff near the flat soft-core origin.
      call assert_close('Beutler LJ force',dv*r,ref, &
                        max(1.0e-7_8*max(1.0_8,abs(ref)),10*epsilon(u)*max(abs(up),abs(um))/step),nfail)
      ref = charge/sqrt(alpha_c*(1-w)**2+r*r)
      call eval_beutler_coulomb(charge,rq,r*r,u,dv,overlap)
      call assert_close('Beutler Eq 8 energy',u,ref,1.0e-12_8*max(1.0_8,abs(ref)),nfail)
      call eval_beutler_coulomb(charge,rq,(r+step)**2,up,dummy,overlap)
      call eval_beutler_coulomb(charge,rq,(r-step)**2,um,dummy,overlap)
      ref = (up-um)/(2*step)
      call assert_close('Beutler Coulomb force',dv*r,ref,1.0e-7_8*max(1.0_8,abs(ref)),nfail)
    end do
  end do
  call eval_beutler_lj(a12,a6,beutler_lj_radius(alpha,sig,0.9999_8),0.0_8,u,dv,overlap)
  call assert_true('Beutler softened overlap regular',.not. overlap,nfail)
  call assert_close('Beutler LJ overlap gradient',dv,0.0_8,0.0_8,nfail)
  call assert_true('Beutler last-window work bounded',abs(0.0001_8*u)<1.0_8,nfail)
  call eval_beutler_lj(0.0_8,0.0_8,0.0_8,0.0_8,u,dv,overlap)
  call assert_close('zero LJ atom contributes zero',u,0.0_8,0.0_8,nfail)
  call eval_beutler_coulomb(charge,10.0_8,0.0_8,u,dv,overlap)
  call assert_close('charged zero-LJ atom still softened',u,charge/10.0_8,1.0e-12_8,nfail)
  ! Unlike Gapsys, Beutler shape corrections must extend outside any soft radius.
  call eval_beutler_lj(a12,a6,2.0_8,100.0_8,um,dummy,overlap)
  call eval_beutler_lj(a12,a6,3.0_8,100.0_8,up,dummy,overlap)
  delta = beutler_lj_energy_difference(a12,a6,2.0_8,3.0_8,100.0_8)
  call assert_true('Beutler long-range LJ correction nonzero',delta /= 0.0_8,nfail)
  call assert_close('Beutler LJ work',delta,up-um,1.0e-12_8,nfail)
  call eval_beutler_coulomb(charge,2.0_8,100.0_8,um,dummy,overlap)
  call eval_beutler_coulomb(charge,3.0_8,100.0_8,up,dummy,overlap)
  delta = beutler_coulomb_energy_difference(charge,2.0_8,3.0_8,100.0_8)
  call assert_true('Beutler long-range Coulomb correction nonzero',delta /= 0.0_8,nfail)
  call assert_close('Beutler Coulomb work',delta,up-um,1.0e-12_8,nfail)

  ! Gapsys SI integrated force polynomial, independently transcribed.
  rs = gapsys_lj_radius(0.85_8,sig,0.5_8)
  a = 78*a12/rs**14-21*a6/rs**8
  b = 168*a12/rs**13-48*a6/rs**7
  c = 91*a12/rs**12-28*a6/rs**6
  r = 0.4_8*rs
  call eval_gapsys_lj(a12,a6,rs,r*r,u,dv,overlap)
  call assert_close('Gapsys SI LJ polynomial',u,a*r*r-b*r+c,1.0e-10_8,nfail)
  call assert_close('Gapsys SI LJ force',dv*r,2*a*r-b,1.0e-10_8,nfail)
  call eval_gapsys_lj(a12,a6,0.0_8,rs*rs,u,hard_dv,overlap)
  call assert_close('Gapsys LJ switch force continuity',2*a*rs-b,hard_dv*rs,1.0e-10_8,nfail)
  rs = gapsys_coulomb_radius(3.0_8,1.0_8,0.25_8,-0.834_8,0.5_8)
  r = 0.4_8*rs
  call eval_gapsys_coulomb(charge,rs,r*r,u,dv,overlap)
  call assert_close('Gapsys SI Coulomb polynomial',u,charge*(r*r/rs**3-3*r/rs**2+3/rs),1.0e-12_8,nfail)
  call assert_close('Gapsys SI Coulomb force',dv*r,charge*(2*r/rs**3-3/rs**2),1.0e-12_8,nfail)
  call assert_close('Gapsys Coulomb switch force continuity',charge*(2/rs**2-3/rs**2), &
                    -charge/rs**2,1.0e-12_8,nfail)
end subroutine

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
