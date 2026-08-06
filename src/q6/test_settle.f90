program test_settle
  use settle
  implicit none

  integer :: failures = 0
  integer :: assertions = 0

  call test_geometry_and_center_of_mass
  call test_matches_converged_shake
  call test_selected_water_list
  call test_impossible_displacement

  write(*,'(a,i0,a,i0,a)') 'SETTLE tests: ', assertions-failures, '/', assertions, ' passed'
  if (failures > 0) error stop 'SETTLE tests failed'

contains

subroutine test_geometry_and_center_of_mass
  real(8), parameter :: d_oh = 0.9572d0
  real(8), parameter :: angle = 104.52d0*acos(-1.0d0)/180.0d0
  real(8), parameter :: d_hh = 2.0d0*d_oh*sin(0.5d0*angle)
  real(8), parameter :: m_o = 15.9994d0, m_h = 1.008d0
  real(8) :: old(9), unconstrained(9), center_before(3), center_after(3)
  logical :: success

  call canonical_water(old, d_oh, d_hh)
  old = old + (/ 4.0d0, -2.0d0, 7.0d0, &
                 4.0d0, -2.0d0, 7.0d0, &
                 4.0d0, -2.0d0, 7.0d0 /)
  unconstrained = old + (/ 0.03d0, -0.02d0, 0.01d0, &
                           -0.01d0, 0.02d0, 0.03d0, &
                           0.01d0, -0.01d0, -0.02d0 /)
  center_before = center_of_mass(unconstrained, m_o, m_h)

  call setup_settle(d_oh, d_hh, m_o, m_h, (/ 1 /))
  success = settle_positions(old, unconstrained)
  center_after = center_of_mass(unconstrained, m_o, m_h)

  call assert_true('ordinary SETTLE step succeeds', success)
  call assert_close('O-H1 distance', bond_length(unconstrained, 1, 2), d_oh, 2.0d-14)
  call assert_close('O-H2 distance', bond_length(unconstrained, 1, 3), d_oh, 2.0d-14)
  call assert_close('H1-H2 distance', bond_length(unconstrained, 2, 3), d_hh, 2.0d-14)
  call assert_vector_close('center of mass', center_after, center_before, 2.0d-14)
end subroutine test_geometry_and_center_of_mass


subroutine test_matches_converged_shake
!!  SETTLE and converged SHAKE solve the same finite-step constraint equations.
!!  This catches rigid-body projections that preserve geometry but do not apply
!!  corrections along the old bond directions.
  real(8), parameter :: d_oh = 0.9572d0
  real(8), parameter :: angle = 104.52d0*acos(-1.0d0)/180.0d0
  real(8), parameter :: d_hh = 2.0d0*d_oh*sin(0.5d0*angle)
  real(8), parameter :: m_o = 15.9994d0, m_h = 1.008d0
  real(8) :: canonical(9), old(9), by_settle(9), by_shake(9)
  real(8) :: rotation(3,3), random_values(15), displacement_scale
  integer :: trial, seed_size
  integer, allocatable :: seed(:)
  logical :: success

  call canonical_water(canonical, d_oh, d_hh)
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 20260806
  call random_seed(put=seed)

  call setup_settle(d_oh, d_hh, m_o, m_h, (/ 1 /))
  do trial = 1, 250
    call random_number(random_values)
    call rotation_matrix(2.0d0*acos(-1.0d0)*random_values(1), &
                         acos(2.0d0*random_values(2)-1.0d0), &
                         2.0d0*acos(-1.0d0)*random_values(3), rotation)
    call rotate_water(canonical, rotation, old)
    old(1:3) = old(1:3) + random_values(4:6)*20.0d0 - 10.0d0
    old(4:6) = old(4:6) + random_values(4:6)*20.0d0 - 10.0d0
    old(7:9) = old(7:9) + random_values(4:6)*20.0d0 - 10.0d0

    displacement_scale = 0.02d0
    by_settle = old + displacement_scale*(2.0d0*random_values(7:15)-1.0d0)
    by_shake = by_settle
    success = settle_positions(old, by_settle)
    call converged_water_shake(old, by_shake, d_oh, d_hh, m_o, m_h)

    call assert_true('random SETTLE step succeeds', success)
    call assert_vector_close('SETTLE agrees with converged SHAKE', &
                             by_settle, by_shake, 2.0d-12)
  end do
end subroutine test_matches_converged_shake


subroutine test_selected_water_list
  real(8), parameter :: d_oh = 1.0d0, d_hh = 1.6d0
  real(8), parameter :: m_o = 16.0d0, m_h = 1.0d0
  real(8) :: old(18), unconstrained(18), untouched(9)
  logical :: success

  call canonical_water(old(1:9), d_oh, d_hh)
  call canonical_water(old(10:18), d_oh, d_hh)
  old(10:18) = old(10:18) + (/ 5.0d0, 4.0d0, 3.0d0, &
                                     5.0d0, 4.0d0, 3.0d0, &
                                     5.0d0, 4.0d0, 3.0d0 /)
  unconstrained = old
  unconstrained(1:9) = unconstrained(1:9) + 0.01d0
  unconstrained(10:18) = unconstrained(10:18) + &
                         (/ 0.02d0, 0.00d0, 0.00d0, 0.00d0, 0.03d0, 0.00d0, &
                            0.00d0, 0.00d0, -0.01d0 /)
  untouched = unconstrained(1:9)

  call setup_settle(d_oh, d_hh, m_o, m_h, (/ 4 /))
  success = settle_positions(old, unconstrained)

  call assert_true('selected-water SETTLE step succeeds', success)
  call assert_vector_close('unselected triplet is untouched', unconstrained(1:9), untouched, 0.0d0)
  call assert_close('selected O-H1', bond_length(unconstrained, 4, 5), d_oh, 2.0d-14)
  call assert_close('selected O-H2', bond_length(unconstrained, 4, 6), d_oh, 2.0d-14)
  call assert_close('selected H1-H2', bond_length(unconstrained, 5, 6), d_hh, 2.0d-14)
end subroutine test_selected_water_list


subroutine test_impossible_displacement
  real(8), parameter :: d_oh = 1.0d0, d_hh = 1.6d0
  real(8) :: old(9), unconstrained(9)
  integer :: failed_oxygen
  logical :: success

  call canonical_water(old, d_oh, d_hh)
  unconstrained = old
  unconstrained(3) = unconstrained(3) + 2.0d0
  call setup_settle(d_oh, d_hh, 16.0d0, 1.0d0, (/ 1 /))
  success = settle_positions(old, unconstrained, failed_oxygen)

  call assert_true('impossible SETTLE step reports failure', .not. success)
  call assert_true('failed oxygen is identified', failed_oxygen == 1)
end subroutine test_impossible_displacement


subroutine canonical_water(coordinates, d_oh, d_hh)
  real(8), intent(out) :: coordinates(9)
  real(8), intent(in) :: d_oh, d_hh
  real(8) :: height

  height = sqrt(d_oh*d_oh - 0.25d0*d_hh*d_hh)
  coordinates = (/ 0.0d0, 0.0d0, 0.0d0, &
                   -0.5d0*d_hh, -height, 0.0d0, &
                    0.5d0*d_hh, -height, 0.0d0 /)
end subroutine canonical_water


subroutine rotation_matrix(phi, theta, psi, rotation)
  real(8), intent(in) :: phi, theta, psi
  real(8), intent(out) :: rotation(3,3)

  rotation(1,1) = cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi)
  rotation(1,2) = -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi)
  rotation(1,3) = cos(phi)*sin(theta)
  rotation(2,1) = sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi)
  rotation(2,2) = -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi)
  rotation(2,3) = sin(phi)*sin(theta)
  rotation(3,1) = -sin(theta)*cos(psi)
  rotation(3,2) = sin(theta)*sin(psi)
  rotation(3,3) = cos(theta)
end subroutine rotation_matrix


subroutine rotate_water(input, rotation, output)
  real(8), intent(in) :: input(9), rotation(3,3)
  real(8), intent(out) :: output(9)
  integer :: atom

  do atom = 1, 3
    output(3*atom-2:3*atom) = matmul(rotation, input(3*atom-2:3*atom))
  end do
end subroutine rotate_water


subroutine converged_water_shake(old, new, d_oh, d_hh, m_o, m_h)
  real(8), intent(in) :: old(9), d_oh, d_hh, m_o, m_h
  real(8), intent(inout) :: new(9)
  integer, parameter :: first_atom(3) = (/ 1, 1, 2 /)
  integer, parameter :: second_atom(3) = (/ 2, 3, 3 /)
  real(8) :: inverse_mass(3), target2(3), current_bond(3), old_bond(3)
  real(8) :: difference, scalar_product, correction
  integer :: iteration, constraint, atom_i, atom_j

  inverse_mass = (/ 1.0d0/m_o, 1.0d0/m_h, 1.0d0/m_h /)
  target2 = (/ d_oh*d_oh, d_oh*d_oh, d_hh*d_hh /)

  do iteration = 1, 10000
    do constraint = 1, 3
      atom_i = first_atom(constraint)
      atom_j = second_atom(constraint)
      current_bond = new(3*atom_i-2:3*atom_i) - new(3*atom_j-2:3*atom_j)
      difference = target2(constraint) - dot_product(current_bond, current_bond)
      old_bond = old(3*atom_i-2:3*atom_i) - old(3*atom_j-2:3*atom_j)
      scalar_product = dot_product(current_bond, old_bond)
      correction = difference / &
                   (2.0d0*scalar_product*(inverse_mass(atom_i)+inverse_mass(atom_j)))
      new(3*atom_i-2:3*atom_i) = new(3*atom_i-2:3*atom_i) + &
                                  old_bond*correction*inverse_mass(atom_i)
      new(3*atom_j-2:3*atom_j) = new(3*atom_j-2:3*atom_j) - &
                                  old_bond*correction*inverse_mass(atom_j)
    end do
    if (max(abs(bond_length(new, 1, 2)-d_oh), &
            abs(bond_length(new, 1, 3)-d_oh), &
            abs(bond_length(new, 2, 3)-d_hh)) < 2.0d-14) return
  end do
  error stop 'reference SHAKE failed to converge'
end subroutine converged_water_shake


pure real(8) function bond_length(coordinates, atom_i, atom_j)
  real(8), intent(in) :: coordinates(:)
  integer, intent(in) :: atom_i, atom_j

  bond_length = sqrt(sum((coordinates(3*atom_i-2:3*atom_i) - &
                          coordinates(3*atom_j-2:3*atom_j))**2))
end function bond_length


pure function center_of_mass(coordinates, m_o, m_h) result(center)
  real(8), intent(in) :: coordinates(9), m_o, m_h
  real(8) :: center(3)

  center = (m_o*coordinates(1:3) + m_h*coordinates(4:6) + &
            m_h*coordinates(7:9)) / (m_o + 2.0d0*m_h)
end function center_of_mass


subroutine assert_close(name, actual, expected, tolerance)
  character(*), intent(in) :: name
  real(8), intent(in) :: actual, expected, tolerance

  assertions = assertions + 1
  if (abs(actual-expected) > tolerance) then
    failures = failures + 1
    write(*,'(a,a,a,es14.6,a,es14.6,a,es14.6)') &
      'FAIL: ', name, ': actual=', actual, ', expected=', expected, &
      ', difference=', abs(actual-expected)
  end if
end subroutine assert_close


subroutine assert_vector_close(name, actual, expected, tolerance)
  character(*), intent(in) :: name
  real(8), intent(in) :: actual(:), expected(:), tolerance

  assertions = assertions + 1
  if (maxval(abs(actual-expected)) > tolerance) then
    failures = failures + 1
    write(*,'(a,a,a,es14.6)') 'FAIL: ', name, ': max difference=', &
                              maxval(abs(actual-expected))
  end if
end subroutine assert_vector_close


subroutine assert_true(name, condition)
  character(*), intent(in) :: name
  logical, intent(in) :: condition

  assertions = assertions + 1
  if (.not. condition) then
    failures = failures + 1
    write(*,'(a,a)') 'FAIL: ', name
  end if
end subroutine assert_true

end program test_settle
