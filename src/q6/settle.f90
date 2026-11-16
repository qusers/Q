module settle
!!-------------------------------------------------------------------------------
!!  Analytical position constraints for rigid three-site water.
!!
!!  Implements Miyamoto and Kollman's SETTLE construction for an O-H-H
!!  molecule with equal hydrogen masses:
!!
!!    S. Miyamoto and P. A. Kollman, J. Comput. Chem. 13, 952-962 (1992).
!!
!!  Q's leap-frog integrator supplies the positions before the unconstrained
!!  step (x_old) and the unconstrained positions after it (x_new).  SETTLE
!!  replaces each selected HOH triplet in x_new with the exact rigid geometry
!!  while preserving the unconstrained center of mass.
!!-------------------------------------------------------------------------------
  implicit none
  private

  real(8), parameter :: numerical_tolerance = 1.0d-12

  real(8) :: mass_o = 0.0d0
  real(8) :: mass_h = 0.0d0
  real(8) :: distance_oh = 0.0d0
  real(8) :: distance_hh = 0.0d0
  real(8) :: radius_a = 0.0d0
  real(8) :: radius_b = 0.0d0
  real(8) :: radius_c = 0.0d0
  real(8) :: inverse_total_mass = 0.0d0
  integer, allocatable :: oxygen_atoms(:)

  public :: setup_settle, settle_positions, settle_is_active

contains

subroutine setup_settle(d_oh, d_hh, m_o, m_h, waters)
!!-------------------------------------------------------------------------------
!!  Validate and cache the common geometry, masses, and oxygen atom indices for
!!  the HOH molecules assigned to SETTLE.  Each oxygen must be followed by H1
!!  and H2 in Q's topology.
!!-------------------------------------------------------------------------------
  real(8), intent(in) :: d_oh, d_hh, m_o, m_h
  integer, intent(in) :: waters(:)
  real(8) :: height

  if (d_oh <= 0.0d0 .or. d_hh <= 0.0d0 .or. d_hh >= 2.0d0*d_oh) then
    error stop 'SETTLE requires a valid triangular O-H-H geometry'
  end if
  if (m_o <= 0.0d0 .or. m_h <= 0.0d0) then
    error stop 'SETTLE requires positive oxygen and hydrogen masses'
  end if

  distance_oh = d_oh
  distance_hh = d_hh
  mass_o = m_o
  mass_h = m_h
  inverse_total_mass = 1.0d0 / (mass_o + 2.0d0*mass_h)

  ! Paper Appendix A canonical coordinates:
  !   O  = (0,  r_a, 0)
  !   H1 = (-r_c, -r_b, 0)
  !   H2 = ( r_c, -r_b, 0)
  ! with the origin at the molecular center of mass.
  radius_c = 0.5d0 * distance_hh
  height = sqrt(distance_oh*distance_oh - radius_c*radius_c)
  radius_a = 2.0d0 * mass_h * height * inverse_total_mass
  radius_b = height - radius_a

  if (allocated(oxygen_atoms)) deallocate(oxygen_atoms)
  allocate(oxygen_atoms(size(waters)))
  oxygen_atoms = waters
end subroutine setup_settle


logical function settle_is_active()
  if (allocated(oxygen_atoms)) then
    settle_is_active = size(oxygen_atoms) > 0
  else
    settle_is_active = .false.
  end if
end function settle_is_active


logical function settle_positions(x_old, x_new, failed_oxygen) result(success)
!!-------------------------------------------------------------------------------
!!  Apply SETTLE to all configured HOH molecules.
!!
!!  Returns .false. when the unconstrained displacement lies outside SETTLE's
!!  geometrical solution domain.  The caller should terminate the MD step rather
!!  than silently producing NaNs.  failed_oxygen identifies the affected oxygen.
!!-------------------------------------------------------------------------------
  real(8), intent(in) :: x_old(:)
  real(8), intent(inout) :: x_new(:)
  integer, intent(out), optional :: failed_oxygen

  integer :: water, atom_o, atom_h1, atom_h2
  integer :: coordinate_o, coordinate_h1, coordinate_h2
  real(8) :: old_o(3), old_h1(3), old_h2(3)
  real(8) :: displacement_o(3), displacement_h1(3), displacement_h2(3)
  real(8) :: old_oh1(3), old_oh2(3), center(3)
  real(8) :: a1(3), b1(3), c1(3)
  real(8) :: axis_x(3), axis_y(3), axis_z(3)
  real(8) :: axis_x_length, axis_y_length, axis_z_length
  real(8) :: xb0, yb0, xc0, yc0
  real(8) :: za1, xb1, yb1, zb1, xc1, yc1, zc1
  real(8) :: sin_phi, cos_phi, sin_psi, cos_psi
  real(8) :: ya2, xb2, yb2, yc2, hh2, delta_x
  real(8) :: alpha, beta, gamma, alpha2_beta2, radicand
  real(8) :: sin_theta, cos_theta
  real(8) :: a3_local(3), b3_local(3), c3_local(3)
  real(8) :: a3(3), b3(3), c3(3)

  success = .true.
  if (present(failed_oxygen)) failed_oxygen = 0
  if (.not. settle_is_active()) return

  do water = 1, size(oxygen_atoms)
    atom_o = oxygen_atoms(water)
    atom_h1 = atom_o + 1
    atom_h2 = atom_o + 2
    coordinate_o = 3*atom_o - 3
    coordinate_h1 = 3*atom_h1 - 3
    coordinate_h2 = 3*atom_h2 - 3

    old_o = x_old(coordinate_o+1:coordinate_o+3)
    old_h1 = x_old(coordinate_h1+1:coordinate_h1+3)
    old_h2 = x_old(coordinate_h2+1:coordinate_h2+3)
    displacement_o = x_new(coordinate_o+1:coordinate_o+3) - old_o
    displacement_h1 = x_new(coordinate_h1+1:coordinate_h1+3) - old_h1
    displacement_h2 = x_new(coordinate_h2+1:coordinate_h2+3) - old_h2

    ! Work relative to the old oxygen.  center is the unconstrained new center
    ! of mass in that frame and must be preserved by the constraint forces.
    old_oh1 = old_h1 - old_o
    old_oh2 = old_h2 - old_o
    center = (mass_o*displacement_o + &
              mass_h*(old_oh1 + displacement_h1) + &
              mass_h*(old_oh2 + displacement_h2)) * inverse_total_mass
    a1 = displacement_o - center
    b1 = old_oh1 + displacement_h1 - center
    c1 = old_oh2 + displacement_h2 - center

    ! Appendix A transformed frame: z is normal to the old water plane, x is
    ! perpendicular to both z and the unconstrained O-COM vector, and y=z*x.
    axis_z = cross_product(old_oh1, old_oh2)
    axis_x = cross_product(a1, axis_z)
    axis_y = cross_product(axis_z, axis_x)
    axis_x_length = sqrt(dot_product(axis_x, axis_x))
    axis_y_length = sqrt(dot_product(axis_y, axis_y))
    axis_z_length = sqrt(dot_product(axis_z, axis_z))
    if (min(axis_x_length, axis_y_length, axis_z_length) <= tiny(1.0d0)) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    axis_x = axis_x / axis_x_length
    axis_y = axis_y / axis_y_length
    axis_z = axis_z / axis_z_length

    xb0 = dot_product(axis_x, old_oh1)
    yb0 = dot_product(axis_y, old_oh1)
    xc0 = dot_product(axis_x, old_oh2)
    yc0 = dot_product(axis_y, old_oh2)
    za1 = dot_product(axis_z, a1)
    xb1 = dot_product(axis_x, b1)
    yb1 = dot_product(axis_y, b1)
    zb1 = dot_product(axis_z, b1)
    xc1 = dot_product(axis_x, c1)
    yc1 = dot_product(axis_y, c1)
    zc1 = dot_product(axis_z, c1)

    ! Appendix A, equations A8-A10.
    sin_phi = za1 / radius_a
    if (.not. unit_interval(sin_phi)) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    sin_phi = clamp_unit(sin_phi)
    cos_phi = sqrt(max(0.0d0, 1.0d0 - sin_phi*sin_phi))
    if (cos_phi <= numerical_tolerance) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if

    sin_psi = (zb1 - zc1) / (2.0d0*radius_c*cos_phi)
    if (.not. unit_interval(sin_psi)) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    sin_psi = clamp_unit(sin_psi)
    cos_psi = sqrt(max(0.0d0, 1.0d0 - sin_psi*sin_psi))

    ya2 = radius_a*cos_phi
    xb2 = -radius_c*cos_psi
    yb2 = -radius_b*cos_phi - radius_c*sin_psi*sin_phi
    yc2 = -radius_b*cos_phi + radius_c*sin_psi*sin_phi

    ! Correct roundoff in the reconstructed H-H distance before solving theta.
    hh2 = 4.0d0*xb2*xb2 + (yb2-yc2)**2 + (zb1-zc1)**2
    radicand = 4.0d0*xb2*xb2 - hh2 + distance_hh*distance_hh
    if (radicand < -numerical_tolerance*distance_hh*distance_hh) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    delta_x = 2.0d0*xb2 + sqrt(max(0.0d0, radicand))
    xb2 = xb2 - 0.5d0*delta_x

    ! Appendix A, equations A14-A17.
    alpha = xb2*(xb0-xc0) + yb0*yb2 + yc0*yc2
    beta = xb2*(yc0-yb0) + xb0*yb2 + xc0*yc2
    gamma = xb0*yb1 - xb1*yb0 + xc0*yc1 - xc1*yc0
    alpha2_beta2 = alpha*alpha + beta*beta
    if (alpha2_beta2 <= tiny(1.0d0)) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    radicand = alpha2_beta2 - gamma*gamma
    if (radicand < -numerical_tolerance*alpha2_beta2) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    sin_theta = (alpha*gamma - beta*sqrt(max(0.0d0, radicand))) / alpha2_beta2
    if (.not. unit_interval(sin_theta)) then
      call report_failure(atom_o, failed_oxygen)
      success = .false.
      return
    end if
    sin_theta = clamp_unit(sin_theta)
    cos_theta = sqrt(max(0.0d0, 1.0d0 - sin_theta*sin_theta))

    a3_local = (/ -ya2*sin_theta, ya2*cos_theta, za1 /)
    b3_local = (/ xb2*cos_theta-yb2*sin_theta, &
                  xb2*sin_theta+yb2*cos_theta, zb1 /)
    c3_local = (/ -xb2*cos_theta-yc2*sin_theta, &
                  -xb2*sin_theta+yc2*cos_theta, zc1 /)

    a3 = axis_x*a3_local(1) + axis_y*a3_local(2) + axis_z*a3_local(3)
    b3 = axis_x*b3_local(1) + axis_y*b3_local(2) + axis_z*b3_local(3)
    c3 = axis_x*c3_local(1) + axis_y*c3_local(2) + axis_z*c3_local(3)

    x_new(coordinate_o+1:coordinate_o+3) = old_o + center + a3
    x_new(coordinate_h1+1:coordinate_h1+3) = old_o + center + b3
    x_new(coordinate_h2+1:coordinate_h2+3) = old_o + center + c3
  end do
end function settle_positions


pure function cross_product(a, b) result(c)
  real(8), intent(in) :: a(3), b(3)
  real(8) :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product


pure logical function unit_interval(value)
  real(8), intent(in) :: value

  unit_interval = abs(value) <= 1.0d0 + numerical_tolerance
end function unit_interval


pure real(8) function clamp_unit(value)
  real(8), intent(in) :: value

  clamp_unit = max(-1.0d0, min(1.0d0, value))
end function clamp_unit


subroutine report_failure(atom_o, failed_oxygen)
  integer, intent(in) :: atom_o
  integer, intent(out), optional :: failed_oxygen

  if (present(failed_oxygen)) failed_oxygen = atom_o
end subroutine report_failure

end module settle
