program test_settle
!!-------------------------------------------------------------------------------
!!  Unit tests for the SETTLE constraint algorithm.
!!
!!  Tests:
!!  1. Single water molecule: verify O-H and H-H distances after constraint
!!  2. Center-of-mass preservation
!!  3. Multiple water molecules
!!  4. Large displacement (stress test)
!!-------------------------------------------------------------------------------
  use settle
  implicit none

  integer :: nfailed, ntests
  nfailed = 0
  ntests = 0

  call test_single_water(nfailed, ntests)
  call test_com_preservation(nfailed, ntests)
  call test_multiple_waters(nfailed, ntests)
  call test_large_displacement(nfailed, ntests)
  call test_random_orientations(nfailed, ntests)
  call test_rotated_water(nfailed, ntests)
  call test_extreme_displacement(nfailed, ntests)

  write(*,'(a)') '=================================='
  write(*,'(a,i3,a,i3,a)') 'SETTLE tests: ', ntests - nfailed, '/', ntests, ' passed'
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


subroutine test_single_water(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! SPC water: dOH = 1.0 A, angle = 109.47 degrees
  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: angle_deg = 109.47d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = angle_deg * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-10

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH_target, d1, d2, d3
  real(8) :: vec(3)

  write(*,'(a)') 'Test: single water molecule'

  ! setup SETTLE parameters
  call setup_settle(dOH, angle_rad, massO, massH)

  dHH_target = 2.0d0 * dOH * sin(angle_rad / 2.0d0)

  ! place water at a known geometry
  ! O at origin, H1 and H2 symmetric about y-axis
  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0     ! O
  xx(4) = -dHH_target/2.0d0; xx(5) = -dOH*cos(angle_rad/2.0d0); xx(6) = 0.0d0  ! H1
  xx(7) =  dHH_target/2.0d0; xx(8) = -dOH*cos(angle_rad/2.0d0); xx(9) = 0.0d0  ! H2

  ! apply a small displacement (simulating unconstrained Verlet step)
  x = xx
  x(1) = x(1) + 0.01d0   ! displace O in x
  x(2) = x(2) - 0.005d0  ! displace O in y
  x(4) = x(4) + 0.02d0   ! displace H1 in x
  x(5) = x(5) + 0.01d0   ! displace H1 in y
  x(6) = x(6) + 0.005d0  ! displace H1 in z
  x(7) = x(7) - 0.01d0   ! displace H2 in x
  x(8) = x(8) + 0.015d0  ! displace H2 in y

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  ! apply SETTLE (nat_solute=0, so water starts at atom 1)
  call settle_positions(1, 0, xx, x, winv)

  ! check O-H1 distance
  vec = x(1:3) - x(4:6)
  d1 = sqrt(dot_product(vec, vec))
  call assert_close('O-H1 distance', d1, dOH, tol, nfailed, ntests)

  ! check O-H2 distance
  vec = x(1:3) - x(7:9)
  d2 = sqrt(dot_product(vec, vec))
  call assert_close('O-H2 distance', d2, dOH, tol, nfailed, ntests)

  ! check H1-H2 distance
  vec = x(4:6) - x(7:9)
  d3 = sqrt(dot_product(vec, vec))
  call assert_close('H1-H2 distance', d3, dHH_target, tol, nfailed, ntests)

end subroutine test_single_water


subroutine test_com_preservation(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-10

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH, com_before(3), com_after(3), total_mass

  write(*,'(a)') 'Test: center-of-mass preservation'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)

  ! place water at known geometry, offset from origin
  xx(1) = 5.0d0;  xx(2) = 3.0d0;  xx(3) = 7.0d0     ! O
  xx(4) = 5.0d0 - dHH/2.0d0; xx(5) = 3.0d0 - dOH*cos(angle_rad/2.0d0); xx(6) = 7.0d0  ! H1
  xx(7) = 5.0d0 + dHH/2.0d0; xx(8) = 3.0d0 - dOH*cos(angle_rad/2.0d0); xx(9) = 7.0d0  ! H2

  ! apply displacement
  x = xx
  x(1) = x(1) + 0.03d0
  x(2) = x(2) - 0.02d0
  x(3) = x(3) + 0.01d0
  x(4) = x(4) - 0.01d0
  x(5) = x(5) + 0.02d0
  x(6) = x(6) + 0.03d0
  x(7) = x(7) + 0.01d0
  x(8) = x(8) - 0.01d0
  x(9) = x(9) - 0.02d0

  ! compute COM of unconstrained positions (this is what SETTLE should preserve)
  total_mass = massO + 2.0d0 * massH
  com_before(1) = (massO*x(1) + massH*x(4) + massH*x(7)) / total_mass
  com_before(2) = (massO*x(2) + massH*x(5) + massH*x(8)) / total_mass
  com_before(3) = (massO*x(3) + massH*x(6) + massH*x(9)) / total_mass

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  call settle_positions(1, 0, xx, x, winv)

  ! compute COM after SETTLE
  com_after(1) = (massO*x(1) + massH*x(4) + massH*x(7)) / total_mass
  com_after(2) = (massO*x(2) + massH*x(5) + massH*x(8)) / total_mass
  com_after(3) = (massO*x(3) + massH*x(6) + massH*x(9)) / total_mass

  call assert_close('COM x preserved', com_after(1), com_before(1), tol, nfailed, ntests)
  call assert_close('COM y preserved', com_after(2), com_before(2), tol, nfailed, ntests)
  call assert_close('COM z preserved', com_after(3), com_before(3), tol, nfailed, ntests)

end subroutine test_com_preservation


subroutine test_multiple_waters(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-10

  integer, parameter :: nw = 3
  real(8) :: xx(9*nw), x(9*nw)
  real :: winv(3*nw)
  real(8) :: dHH, d1, d2, d3, vec(3)
  integer :: w, i3O

  write(*,'(a)') 'Test: multiple water molecules'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)

  ! place 3 waters at different locations
  do w = 1, nw
    i3O = 9*(w-1)
    xx(i3O+1) = 10.0d0 * w;  xx(i3O+2) = 0.0d0;  xx(i3O+3) = 0.0d0  ! O
    xx(i3O+4) = 10.0d0*w - dHH/2.0d0;  xx(i3O+5) = -dOH*cos(angle_rad/2.0d0);  xx(i3O+6) = 0.0d0  ! H1
    xx(i3O+7) = 10.0d0*w + dHH/2.0d0;  xx(i3O+8) = -dOH*cos(angle_rad/2.0d0);  xx(i3O+9) = 0.0d0  ! H2

    winv(3*(w-1)+1) = real(1.0d0 / massO)
    winv(3*(w-1)+2) = real(1.0d0 / massH)
    winv(3*(w-1)+3) = real(1.0d0 / massH)
  end do

  ! apply different displacements to each water
  x = xx
  x(1) = x(1) + 0.01d0; x(5) = x(5) + 0.02d0; x(9) = x(9) - 0.01d0
  x(10) = x(10) - 0.02d0; x(14) = x(14) + 0.03d0; x(18) = x(18) + 0.01d0
  x(19) = x(19) + 0.005d0; x(23) = x(23) - 0.01d0; x(27) = x(27) + 0.02d0

  call settle_positions(nw, 0, xx, x, winv)

  ! check all 3 waters
  do w = 1, nw
    i3O = 9*(w-1)

    vec = x(i3O+1:i3O+3) - x(i3O+4:i3O+6)
    d1 = sqrt(dot_product(vec, vec))

    vec = x(i3O+1:i3O+3) - x(i3O+7:i3O+9)
    d2 = sqrt(dot_product(vec, vec))

    vec = x(i3O+4:i3O+6) - x(i3O+7:i3O+9)
    d3 = sqrt(dot_product(vec, vec))

    call assert_close('O-H1 dist (multi)', d1, dOH, tol, nfailed, ntests)
    call assert_close('O-H2 dist (multi)', d2, dOH, tol, nfailed, ntests)
    call assert_close('H1-H2 dist (multi)', d3, dHH, tol, nfailed, ntests)
  end do

end subroutine test_multiple_waters


subroutine test_large_displacement(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-8  ! relaxed tolerance for large displacement

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH, d1, d2, d3, vec(3)

  write(*,'(a)') 'Test: large displacement'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)

  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0
  xx(4) = -dHH/2.0d0; xx(5) = -dOH*cos(angle_rad/2.0d0); xx(6) = 0.0d0
  xx(7) =  dHH/2.0d0; xx(8) = -dOH*cos(angle_rad/2.0d0); xx(9) = 0.0d0

  ! large displacement (0.1 A, typical of a few fs with fast H motion)
  x = xx
  x(1) = x(1) + 0.1d0
  x(2) = x(2) - 0.05d0
  x(3) = x(3) + 0.08d0
  x(4) = x(4) - 0.12d0
  x(5) = x(5) + 0.1d0
  x(6) = x(6) + 0.07d0
  x(7) = x(7) + 0.08d0
  x(8) = x(8) - 0.09d0
  x(9) = x(9) - 0.06d0

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  call settle_positions(1, 0, xx, x, winv)

  vec = x(1:3) - x(4:6)
  d1 = sqrt(dot_product(vec, vec))
  vec = x(1:3) - x(7:9)
  d2 = sqrt(dot_product(vec, vec))
  vec = x(4:6) - x(7:9)
  d3 = sqrt(dot_product(vec, vec))

  call assert_close('O-H1 large disp', d1, dOH, tol, nfailed, ntests)
  call assert_close('O-H2 large disp', d2, dOH, tol, nfailed, ntests)
  call assert_close('H1-H2 large disp', d3, dHH, tol, nfailed, ntests)

end subroutine test_large_displacement


subroutine test_random_orientations(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-8

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH, d1, d2, d3, vec(3)
  real(8) :: ra, phi, theta, psi
  real(8) :: R(3,3)  ! rotation matrix
  real(8) :: canonical(3,3)  ! canonical water: O, H1, H2
  real(8) :: com(3), disp(3)
  integer :: trial, nfailed_here, seed
  real(8) :: rng

  write(*,'(a)') 'Test: random orientations (100 trials)'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)
  ra = dOH * cos(angle_rad / 2.0d0)

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  ! canonical water geometry
  canonical(1,:) = (/ 0.0d0, 0.0d0, 0.0d0 /)       ! O at origin
  canonical(2,:) = (/ -dHH/2.0d0, -ra, 0.0d0 /)     ! H1
  canonical(3,:) = (/  dHH/2.0d0, -ra, 0.0d0 /)     ! H2

  nfailed_here = 0
  seed = 42

  do trial = 1, 100
    ! simple LCG random number generator
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    phi = dble(mod(abs(seed), 10000)) / 10000.0d0 * 2.0d0 * pi
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    theta = dble(mod(abs(seed), 10000)) / 10000.0d0 * pi
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    psi = dble(mod(abs(seed), 10000)) / 10000.0d0 * 2.0d0 * pi

    ! build rotation matrix (Euler angles ZYZ convention)
    R(1,1) = cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi)
    R(1,2) = -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi)
    R(1,3) = cos(phi)*sin(theta)
    R(2,1) = sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi)
    R(2,2) = -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi)
    R(2,3) = sin(phi)*sin(theta)
    R(3,1) = -sin(theta)*cos(psi)
    R(3,2) = sin(theta)*sin(psi)
    R(3,3) = cos(theta)

    ! random translation
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    com(1) = dble(mod(abs(seed), 10000)) / 10000.0d0 * 20.0d0 - 10.0d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    com(2) = dble(mod(abs(seed), 10000)) / 10000.0d0 * 20.0d0 - 10.0d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    com(3) = dble(mod(abs(seed), 10000)) / 10000.0d0 * 20.0d0 - 10.0d0

    ! rotate and translate water
    xx(1:3) = matmul(R, canonical(1,:)) + com
    xx(4:6) = matmul(R, canonical(2,:)) + com
    xx(7:9) = matmul(R, canonical(3,:)) + com

    ! random small displacements (0.01-0.05 A)
    x = xx
    do seed = seed, seed + 8
    end do
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(1) = x(1) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(2) = x(2) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(3) = x(3) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(4) = x(4) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(5) = x(5) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(6) = x(6) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(7) = x(7) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(8) = x(8) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0
    seed = mod(seed * 1103515245 + 12345, 2147483647)
    x(9) = x(9) + dble(mod(abs(seed), 1000)) / 1000.0d0 * 0.05d0 - 0.025d0

    call settle_positions(1, 0, xx, x, winv)

    vec = x(1:3) - x(4:6)
    d1 = sqrt(dot_product(vec, vec))
    vec = x(1:3) - x(7:9)
    d2 = sqrt(dot_product(vec, vec))
    vec = x(4:6) - x(7:9)
    d3 = sqrt(dot_product(vec, vec))

    ntests = ntests + 3
    if (abs(d1 - dOH) > tol) then
      write(*,'(a,i4,a,es12.5)') '  FAIL trial ', trial, ' O-H1 diff=', abs(d1-dOH)
      nfailed = nfailed + 1
      nfailed_here = nfailed_here + 1
    end if
    if (abs(d2 - dOH) > tol) then
      write(*,'(a,i4,a,es12.5)') '  FAIL trial ', trial, ' O-H2 diff=', abs(d2-dOH)
      nfailed = nfailed + 1
      nfailed_here = nfailed_here + 1
    end if
    if (abs(d3 - dHH) > tol) then
      write(*,'(a,i4,a,es12.5)') '  FAIL trial ', trial, ' H1-H2 diff=', abs(d3-dHH)
      nfailed = nfailed + 1
      nfailed_here = nfailed_here + 1
    end if
  end do

  if (nfailed_here == 0) then
    write(*,'(a)') '  All 100 random orientations passed'
  else
    write(*,'(a,i4,a)') '  ', nfailed_here, ' failures in random orientations'
  end if

end subroutine test_random_orientations


subroutine test_rotated_water(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! water rotated 90 degrees around z-axis — known to trigger sign bug
  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-10

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH, ra, d1, d2, d3, vec(3)

  write(*,'(a)') 'Test: 90-degree rotated water'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)
  ra = dOH * cos(angle_rad / 2.0d0)

  ! O at origin, H1-H2 along Y axis (rotated 90 degrees from canonical)
  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0     ! O
  xx(4) = -ra;    xx(5) = -dHH/2.0d0; xx(6) = 0.0d0  ! H1
  xx(7) = -ra;    xx(8) =  dHH/2.0d0; xx(9) = 0.0d0  ! H2

  x = xx
  x(1) = x(1) + 0.02d0;   x(2) = x(2) - 0.01d0
  x(4) = x(4) - 0.015d0;  x(5) = x(5) + 0.02d0;  x(6) = x(6) + 0.01d0
  x(7) = x(7) + 0.01d0;   x(8) = x(8) - 0.015d0; x(9) = x(9) - 0.005d0

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  call settle_positions(1, 0, xx, x, winv)

  vec = x(1:3) - x(4:6)
  d1 = sqrt(dot_product(vec, vec))
  vec = x(1:3) - x(7:9)
  d2 = sqrt(dot_product(vec, vec))
  vec = x(4:6) - x(7:9)
  d3 = sqrt(dot_product(vec, vec))

  call assert_close('O-H1 rotated', d1, dOH, tol, nfailed, ntests)
  call assert_close('O-H2 rotated', d2, dOH, tol, nfailed, ntests)
  call assert_close('H1-H2 rotated', d3, dHH, tol, nfailed, ntests)

end subroutine test_rotated_water


subroutine test_extreme_displacement(nfailed, ntests)
  integer, intent(inout) :: nfailed, ntests

  ! large asymmetric H displacements — stresses the alpha determination
  real(8), parameter :: dOH = 1.0d0
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), parameter :: angle_rad = 109.47d0 * pi / 180.0d0
  real(8), parameter :: massO = 15.9994d0
  real(8), parameter :: massH = 1.008d0
  real(8), parameter :: tol = 1.0d-6

  real(8) :: xx(9), x(9)
  real :: winv(3)
  real(8) :: dHH, ra, d1, d2, d3, vec(3)
  real(8) :: com_before(3), com_after(3), total_mass

  write(*,'(a)') 'Test: extreme asymmetric displacement'

  call setup_settle(dOH, angle_rad, massO, massH)
  dHH = 2.0d0 * dOH * sin(angle_rad / 2.0d0)
  ra = dOH * cos(angle_rad / 2.0d0)
  total_mass = massO + 2.0d0 * massH

  winv(1) = real(1.0d0 / massO)
  winv(2) = real(1.0d0 / massH)
  winv(3) = real(1.0d0 / massH)

  ! canonical water
  xx(1) = 0.0d0;  xx(2) = 0.0d0;  xx(3) = 0.0d0
  xx(4) = -dHH/2.0d0; xx(5) = -ra; xx(6) = 0.0d0
  xx(7) =  dHH/2.0d0; xx(8) = -ra; xx(9) = 0.0d0

  ! large opposing H displacements + O moved beyond COM
  x = xx
  x(2) = x(2) - 0.3d0    ! O pushed far down (past COM)
  x(3) = x(3) + 0.15d0   ! O pushed out of plane
  x(4) = x(4) + 0.2d0    ! H1 pushed right
  x(5) = x(5) + 0.15d0   ! H1 pushed up
  x(6) = x(6) - 0.1d0    ! H1 out of plane
  x(7) = x(7) - 0.15d0   ! H2 pushed left
  x(8) = x(8) + 0.2d0    ! H2 pushed up
  x(9) = x(9) + 0.1d0    ! H2 out of plane

  com_before(1) = (massO*x(1) + massH*x(4) + massH*x(7)) / total_mass
  com_before(2) = (massO*x(2) + massH*x(5) + massH*x(8)) / total_mass
  com_before(3) = (massO*x(3) + massH*x(6) + massH*x(9)) / total_mass

  call settle_positions(1, 0, xx, x, winv)

  vec = x(1:3) - x(4:6)
  d1 = sqrt(dot_product(vec, vec))
  vec = x(1:3) - x(7:9)
  d2 = sqrt(dot_product(vec, vec))
  vec = x(4:6) - x(7:9)
  d3 = sqrt(dot_product(vec, vec))

  call assert_close('O-H1 extreme', d1, dOH, tol, nfailed, ntests)
  call assert_close('O-H2 extreme', d2, dOH, tol, nfailed, ntests)
  call assert_close('H1-H2 extreme', d3, dHH, tol, nfailed, ntests)

  ! COM must be preserved
  com_after(1) = (massO*x(1) + massH*x(4) + massH*x(7)) / total_mass
  com_after(2) = (massO*x(2) + massH*x(5) + massH*x(8)) / total_mass
  com_after(3) = (massO*x(3) + massH*x(6) + massH*x(9)) / total_mass

  call assert_close('COM x extreme', com_after(1), com_before(1), 1.0d-10, nfailed, ntests)
  call assert_close('COM y extreme', com_after(2), com_before(2), 1.0d-10, nfailed, ntests)
  call assert_close('COM z extreme', com_after(3), com_before(3), 1.0d-10, nfailed, ntests)

end subroutine test_extreme_displacement


end program test_settle
