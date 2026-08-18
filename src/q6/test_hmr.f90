program test_hmr
  use hmr
  implicit none

  integer :: passed, total

  passed = 0
  total = 0
  call test_methyl_and_methylene
  call test_hydrogen_hydrogen_bond_is_ignored
  call test_ambiguous_parent_fails_without_mutation
  call test_parent_too_light_fails_without_mutation
  call test_invalid_hydrogen_mass_fails

  write(*,'(a,i0,a,i0,a)') 'HMR tests: ', passed, '/', total, ' passed'
  if (passed /= total) stop 1

contains

subroutine check(condition, description)
  logical, intent(in) :: condition
  character(*), intent(in) :: description
  total = total + 1
  if (condition) then
    passed = passed + 1
  else
    write(*,'(a)') 'FAILED: '//description
  end if
end subroutine check

subroutine test_methyl_and_methylene
  real(8) :: mass(7), before, minimum
  logical :: hydrogen(7)
  integer :: bi(5), bj(5), status, failed, nh, np

  ! CH3 and CH2 parents deliberately start with the same nominal carbon mass.
  mass = [12.010d0, 1.008d0, 1.008d0, 1.008d0, 12.010d0, 1.008d0, 1.008d0]
  hydrogen = [.false., .true., .true., .true., .false., .true., .true.]
  bi = [1, 1, 1, 5, 5]
  bj = [2, 3, 4, 6, 7]
  before = sum(mass)

  call repartition_hydrogen_masses(mass, hydrogen, bi, bj, 3.024d0, status, failed, nh, np, minimum)

  call check(status == HMR_OK, 'methyl/methylene repartition succeeds')
  call check(all(abs(mass([2,3,4,6,7])-3.024d0) < 1.0d-12), 'all hydrogen masses reach target')
  call check(abs(mass(1)-5.962d0) < 1.0d-12, 'methyl parent loses three increments')
  call check(abs(mass(5)-7.978d0) < 1.0d-12, 'methylene parent loses two increments')
  call check(abs(sum(mass)-before) < 1.0d-12, 'total mass is conserved')
  call check(nh == 5 .and. np == 2, 'summary counts hydrogens and unique parents')
  call check(abs(minimum-5.962d0) < 1.0d-12, 'summary reports minimum parent mass')
end subroutine test_methyl_and_methylene

subroutine test_hydrogen_hydrogen_bond_is_ignored
  real(8) :: mass(3), minimum
  logical :: hydrogen(3)
  integer :: bi(3), bj(3), status, failed, nh, np

  ! Three-distance rigid-water-style graph: each H still has one heavy parent.
  mass = [16.0d0, 1.008d0, 1.008d0]
  hydrogen = [.false., .true., .true.]
  bi = [1, 1, 2]
  bj = [2, 3, 3]
  call repartition_hydrogen_masses(mass, hydrogen, bi, bj, 3.024d0, status, failed, nh, np, minimum)

  call check(status == HMR_OK, 'H-H geometry bond does not create a second heavy parent')
  call check(abs(mass(1)-11.968d0) < 1.0d-12, 'water-style oxygen receives both decrements')
end subroutine test_hydrogen_hydrogen_bond_is_ignored

subroutine test_ambiguous_parent_fails_without_mutation
  real(8) :: mass(3), original(3), minimum
  logical :: hydrogen(3)
  integer :: bi(2), bj(2), status, failed, nh, np

  mass = [12.0d0, 1.008d0, 14.0d0]
  original = mass
  hydrogen = [.false., .true., .false.]
  bi = [1, 2]
  bj = [2, 3]
  call repartition_hydrogen_masses(mass, hydrogen, bi, bj, 3.024d0, status, failed, nh, np, minimum)

  call check(status == HMR_PARENT_NOT_UNIQUE .and. failed == 2, 'bridging H is rejected')
  call check(all(abs(mass-original) < 1.0d-15), 'ambiguous-parent failure is atomic')
end subroutine test_ambiguous_parent_fails_without_mutation

subroutine test_parent_too_light_fails_without_mutation
  real(8) :: mass(2), original(2), minimum
  logical :: hydrogen(2)
  integer :: bi(1), bj(1), status, failed, nh, np

  mass = [4.0d0, 1.008d0]
  original = mass
  hydrogen = [.false., .true.]
  bi = [1]
  bj = [2]
  call repartition_hydrogen_masses(mass, hydrogen, bi, bj, 3.024d0, status, failed, nh, np, minimum)

  call check(status == HMR_PARENT_TOO_LIGHT .and. failed == 1, 'too-light parent is rejected')
  call check(all(abs(mass-original) < 1.0d-15), 'too-light-parent failure is atomic')
end subroutine test_parent_too_light_fails_without_mutation

subroutine test_invalid_hydrogen_mass_fails
  real(8) :: mass(2), minimum
  logical :: hydrogen(2)
  integer :: bi(1), bj(1), status, failed, nh, np

  mass = [12.0d0, 3.024d0]
  hydrogen = [.false., .true.]
  bi = [1]
  bj = [2]
  call repartition_hydrogen_masses(mass, hydrogen, bi, bj, 3.024d0, status, failed, nh, np, minimum)
  call check(status == HMR_INVALID_HYDROGEN_MASS .and. failed == 2, 'already-heavy H is rejected')
end subroutine test_invalid_hydrogen_mass_fails

end program test_hmr
