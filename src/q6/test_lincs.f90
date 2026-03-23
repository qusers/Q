program test_lincs
!!-------------------------------------------------------------------------------
!!  Unit tests for the LINCS constraint algorithm.
!!
!!  Tests:
!!  1. Single C-H bond constraint
!!  2. Coupled constraints (CH3 group: 3 C-H bonds sharing C)
!!  3. Zero-displacement (positions already at equilibrium)
!!  4. Momentum conservation
!!-------------------------------------------------------------------------------
  use lincs
  implicit none

  integer :: nfailed, ntests
  nfailed = 0
  ntests = 0

  call test_single_bond(nfailed, ntests)
  call test_coupled_constraints(nfailed, ntests)
  call test_zero_displacement(nfailed, ntests)
  call test_momentum_conservation(nfailed, ntests)

  write(*,'(a)') '=================================='
  write(*,'(a,i3,a,i3,a)') 'LINCS tests: ', ntests - nfailed, '/', ntests, ' passed'
  write(*,'(a)') '=================================='

  if (nfailed > 0) then
    write(*,'(a)') 'FAIL'
    stop 1
  else
    write(*,'(a)') 'OK'
  end if

contains

subroutine assert_close(name, actual, expected, tol, nfailed, ntests)
  character(*), intent(in) :: name
  real(8), intent(in) :: actual, expected, tol
  integer, intent(inout) :: nfailed, ntests
  real(8) :: diff

  ntests = ntests + 1
  diff = abs(actual - expected)
  if (diff > tol) then
    write(*,'(a,a,a,es12.5,a,es12.5,a,es12.5)') &
      '  FAIL: ', name, '  got=', actual, '  expected=', expected, '  diff=', diff
    nfailed = nfailed + 1
  end if
end subroutine assert_close


subroutine test_single_bond(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! single C-H bond constraint
  real(8), parameter :: dCH = 1.09d0  ! Angstroms
  real(8), parameter :: massC = 12.0d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-6

  integer :: ai(1), aj(1), fep_idx(1)
  real(8) :: target_len(1), inv_mi(1), inv_mj(1)
  real(8) :: xx(6), x(6), vec(3), d

  write(*,'(a)') 'Test: single C-H bond'

  ai(1) = 1;  aj(1) = 2
  target_len(1) = dCH
  inv_mi(1) = 1.0d0 / massC
  inv_mj(1) = 1.0d0 / massH

  call setup_lincs(1, ai, aj, target_len, inv_mi, inv_mj, 0, fep_idx, 4)

  ! place atoms: C at origin, H along x-axis at dCH
  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0
  xx(4) = dCH;    xx(5) = 0.0d0;  xx(6) = 0.0d0

  ! displace H by 0.05 A in y (stretching the bond)
  x = xx
  x(2) = x(2) + 0.01d0   ! C moves a bit in y
  x(4) = x(4) + 0.03d0   ! H moves in x
  x(5) = x(5) + 0.05d0   ! H moves in y

  call lincs_positions(xx, x)

  ! check distance
  vec = x(1:3) - x(4:6)
  d = sqrt(dot_product(vec, vec))
  call assert_close('C-H distance', d, dCH, tol, nfailed, ntests)

  ! cleanup
  if (allocated(lincs_data%atom_i)) deallocate(lincs_data%atom_i)
  if (allocated(lincs_data%atom_j)) deallocate(lincs_data%atom_j)
  if (allocated(lincs_data%target_len)) deallocate(lincs_data%target_len)
  if (allocated(lincs_data%inv_mass_i)) deallocate(lincs_data%inv_mass_i)
  if (allocated(lincs_data%inv_mass_j)) deallocate(lincs_data%inv_mass_j)
  if (allocated(lincs_data%sdiag)) deallocate(lincs_data%sdiag)
  if (allocated(lincs_data%coupled_start)) deallocate(lincs_data%coupled_start)
  if (allocated(lincs_data%coupled_index)) deallocate(lincs_data%coupled_index)
  if (allocated(lincs_data%coupled_coeff)) deallocate(lincs_data%coupled_coeff)
  if (allocated(lincs_data%fep_idx)) deallocate(lincs_data%fep_idx)
  lincs_data%nconstraints = 0
  lincs_data%nfep = 0

end subroutine test_single_bond


subroutine test_coupled_constraints(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! CH3 group: C at center, 3 H atoms bonded to it
  ! 3 coupled constraints (all share the C atom)
  real(8), parameter :: dCH = 1.09d0
  real(8), parameter :: massC = 12.0d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-5

  integer :: ai(3), aj(3), fep_idx(1)
  real(8) :: target_len(3), inv_mi(3), inv_mj(3)
  real(8) :: xx(12), x(12), vec(3), d
  integer :: k

  write(*,'(a)') 'Test: coupled constraints (CH3)'

  ! atom 1 = C, atoms 2,3,4 = H
  do k = 1, 3
    ai(k) = 1
    aj(k) = k + 1
    target_len(k) = dCH
    inv_mi(k) = 1.0d0 / massC
    inv_mj(k) = 1.0d0 / massH
  end do

  call setup_lincs(3, ai, aj, target_len, inv_mi, inv_mj, 0, fep_idx, 4)

  ! C at origin, H atoms arranged tetrahedrally
  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0  ! C
  xx(4) = dCH;    xx(5) = 0.0d0;  xx(6) = 0.0d0   ! H1 along x
  xx(7) = -dCH*0.333d0;  xx(8) = dCH*0.943d0;  xx(9) = 0.0d0   ! H2
  xx(10) = -dCH*0.333d0; xx(11) = -dCH*0.471d0; xx(12) = dCH*0.816d0  ! H3

  ! apply displacements
  x = xx
  x(1) = x(1) + 0.02d0;  x(2) = x(2) - 0.01d0  ! C
  x(4) = x(4) + 0.05d0;  x(5) = x(5) + 0.03d0   ! H1
  x(8) = x(8) - 0.04d0;  x(9) = x(9) + 0.02d0   ! H2
  x(10) = x(10) - 0.01d0; x(12) = x(12) + 0.03d0 ! H3

  call lincs_positions(xx, x)

  ! check all 3 C-H distances
  do k = 1, 3
    vec = x(1:3) - x(3*k+1 : 3*k+3)
    d = sqrt(dot_product(vec, vec))
    call assert_close('C-H distance (coupled)', d, dCH, tol, nfailed, ntests)
  end do

  ! cleanup
  if (allocated(lincs_data%atom_i)) deallocate(lincs_data%atom_i)
  if (allocated(lincs_data%atom_j)) deallocate(lincs_data%atom_j)
  if (allocated(lincs_data%target_len)) deallocate(lincs_data%target_len)
  if (allocated(lincs_data%inv_mass_i)) deallocate(lincs_data%inv_mass_i)
  if (allocated(lincs_data%inv_mass_j)) deallocate(lincs_data%inv_mass_j)
  if (allocated(lincs_data%sdiag)) deallocate(lincs_data%sdiag)
  if (allocated(lincs_data%coupled_start)) deallocate(lincs_data%coupled_start)
  if (allocated(lincs_data%coupled_index)) deallocate(lincs_data%coupled_index)
  if (allocated(lincs_data%coupled_coeff)) deallocate(lincs_data%coupled_coeff)
  if (allocated(lincs_data%fep_idx)) deallocate(lincs_data%fep_idx)
  lincs_data%nconstraints = 0
  lincs_data%nfep = 0

end subroutine test_coupled_constraints


subroutine test_zero_displacement(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! positions already satisfy constraints — LINCS should make no change
  real(8), parameter :: dCH = 1.09d0
  real(8), parameter :: massC = 12.0d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-12

  integer :: ai(1), aj(1), fep_idx(1)
  real(8) :: target_len(1), inv_mi(1), inv_mj(1)
  real(8) :: xx(6), x(6), x_before(6)
  integer :: k

  write(*,'(a)') 'Test: zero displacement'

  ai(1) = 1;  aj(1) = 2
  target_len(1) = dCH
  inv_mi(1) = 1.0d0 / massC
  inv_mj(1) = 1.0d0 / massH

  call setup_lincs(1, ai, aj, target_len, inv_mi, inv_mj, 0, fep_idx, 4)

  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0
  xx(4) = dCH;    xx(5) = 0.0d0;  xx(6) = 0.0d0

  x = xx  ! no displacement
  x_before = x

  call lincs_positions(xx, x)

  do k = 1, 6
    call assert_close('no change', x(k), x_before(k), tol, nfailed, ntests)
  end do

  ! cleanup
  if (allocated(lincs_data%atom_i)) deallocate(lincs_data%atom_i)
  if (allocated(lincs_data%atom_j)) deallocate(lincs_data%atom_j)
  if (allocated(lincs_data%target_len)) deallocate(lincs_data%target_len)
  if (allocated(lincs_data%inv_mass_i)) deallocate(lincs_data%inv_mass_i)
  if (allocated(lincs_data%inv_mass_j)) deallocate(lincs_data%inv_mass_j)
  if (allocated(lincs_data%sdiag)) deallocate(lincs_data%sdiag)
  if (allocated(lincs_data%coupled_start)) deallocate(lincs_data%coupled_start)
  if (allocated(lincs_data%coupled_index)) deallocate(lincs_data%coupled_index)
  if (allocated(lincs_data%coupled_coeff)) deallocate(lincs_data%coupled_coeff)
  if (allocated(lincs_data%fep_idx)) deallocate(lincs_data%fep_idx)
  lincs_data%nconstraints = 0
  lincs_data%nfep = 0

end subroutine test_zero_displacement


subroutine test_momentum_conservation(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! LINCS corrections should conserve total momentum
  ! (corrections are equal and opposite, weighted by inverse mass)
  real(8), parameter :: dCH = 1.09d0
  real(8), parameter :: massC = 12.0d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-10

  integer :: ai(1), aj(1), fep_idx(1)
  real(8) :: target_len(1), inv_mi(1), inv_mj(1)
  real(8) :: xx(6), x(6), x_before(6)
  real(8) :: dp_before(3), dp_after(3)
  real(8) :: dt
  integer :: k

  write(*,'(a)') 'Test: momentum conservation'

  ai(1) = 1;  aj(1) = 2
  target_len(1) = dCH
  inv_mi(1) = 1.0d0 / massC
  inv_mj(1) = 1.0d0 / massH

  call setup_lincs(1, ai, aj, target_len, inv_mi, inv_mj, 0, fep_idx, 4)

  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0
  xx(4) = dCH;    xx(5) = 0.0d0;  xx(6) = 0.0d0

  x = xx
  x(1) = x(1) + 0.01d0
  x(4) = x(4) + 0.05d0;  x(5) = x(5) + 0.03d0

  x_before = x

  ! total momentum before (proportional to mass * displacement from xx)
  dp_before(1) = massC * (x_before(1) - xx(1)) + massH * (x_before(4) - xx(4))
  dp_before(2) = massC * (x_before(2) - xx(2)) + massH * (x_before(5) - xx(5))
  dp_before(3) = massC * (x_before(3) - xx(3)) + massH * (x_before(6) - xx(6))

  call lincs_positions(xx, x)

  ! total momentum after
  dp_after(1) = massC * (x(1) - xx(1)) + massH * (x(4) - xx(4))
  dp_after(2) = massC * (x(2) - xx(2)) + massH * (x(5) - xx(5))
  dp_after(3) = massC * (x(3) - xx(3)) + massH * (x(6) - xx(6))

  call assert_close('momentum x', dp_after(1), dp_before(1), tol, nfailed, ntests)
  call assert_close('momentum y', dp_after(2), dp_before(2), tol, nfailed, ntests)
  call assert_close('momentum z', dp_after(3), dp_before(3), tol, nfailed, ntests)

  ! cleanup
  if (allocated(lincs_data%atom_i)) deallocate(lincs_data%atom_i)
  if (allocated(lincs_data%atom_j)) deallocate(lincs_data%atom_j)
  if (allocated(lincs_data%target_len)) deallocate(lincs_data%target_len)
  if (allocated(lincs_data%inv_mass_i)) deallocate(lincs_data%inv_mass_i)
  if (allocated(lincs_data%inv_mass_j)) deallocate(lincs_data%inv_mass_j)
  if (allocated(lincs_data%sdiag)) deallocate(lincs_data%sdiag)
  if (allocated(lincs_data%coupled_start)) deallocate(lincs_data%coupled_start)
  if (allocated(lincs_data%coupled_index)) deallocate(lincs_data%coupled_index)
  if (allocated(lincs_data%coupled_coeff)) deallocate(lincs_data%coupled_coeff)
  if (allocated(lincs_data%fep_idx)) deallocate(lincs_data%fep_idx)
  lincs_data%nconstraints = 0
  lincs_data%nfep = 0

end subroutine test_momentum_conservation

end program test_lincs
