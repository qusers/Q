module boundary_corrections
!!-----------------------------------------------------------------------------
!! Pure helpers for spherical-boundary corrections used by qdyn.
!!-----------------------------------------------------------------------------
  implicit none

contains

  pure function born_coefficient(ke, eps, radius) result(c)
    real(8), intent(in) :: ke, eps, radius
    real(8)             :: c

    c = ke * (1.0_8 - 1.0_8/eps) / (2.0_8 * radius)
  end function born_coefficient

  pure function born_self_energy(q_sphere, c) result(energy)
    real(8), intent(in) :: q_sphere, c
    real(8)             :: energy

    energy = -c * q_sphere * q_sphere
  end function born_self_energy

  pure function polarization_strength(q_region, dielectric_factor, density, dipole, radius) result(c)
    real(8), intent(in) :: q_region, dielectric_factor, density, dipole, radius
    real(8)             :: c
    real(8), parameter  :: pi = 3.1415926535897932384626433832795_8

    c = q_region * dielectric_factor / (density * dipole * 4.0_8 * pi * radius**2)
  end function polarization_strength

  pure function polarization_target(theta_field_free, strength) result(theta_target)
    real(8), intent(in) :: theta_field_free, strength
    real(8)             :: theta_target
    real(8), parameter  :: pi = 3.1415926535897932384626433832795_8

    theta_target = theta_field_free - 1.5_8 * strength * sin(theta_field_free)
    theta_target = max(0.0_8, min(pi, theta_target))
  end function polarization_target

  pure function polarization_energy(theta, theta_target, theta_offset, force_constant) result(energy)
    real(8), intent(in) :: theta, theta_target, theta_offset, force_constant
    real(8)             :: energy
    real(8)             :: displacement

    displacement = theta - theta_target + theta_offset
    energy = 0.5_8 * force_constant * displacement * displacement
  end function polarization_energy

  pure function polarization_gradient(theta, theta_target, theta_offset, force_constant) result(gradient)
    real(8), intent(in) :: theta, theta_target, theta_offset, force_constant
    real(8)             :: gradient

    gradient = force_constant * (theta - theta_target + theta_offset)
  end function polarization_gradient

  pure function safe_acos(value) result(angle)
    real(8), intent(in) :: value
    real(8)             :: angle

    angle = acos(max(-1.0_8, min(1.0_8, value)))
  end function safe_acos

end module boundary_corrections
