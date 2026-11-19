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

end module boundary_corrections
