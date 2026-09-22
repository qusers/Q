module boundary_corrections
!!-----------------------------------------------------------------------------
!! Pure helpers for spherical-boundary corrections used by qdyn.
!!-----------------------------------------------------------------------------
  implicit none

contains

  pure subroutine polarization_switch(radius, boundary, halfwidth, weight, derivative)
    ! Quintic transition: value and first two derivatives join at both ends.
    real(8), intent(in) :: radius, boundary, halfwidth
    real(8), intent(out) :: weight, derivative
    real(8) :: t
    t=(radius-boundary+halfwidth)/(2*halfwidth)
    weight=0; derivative=0
    if (t <= 0) return
    if (t >= 1) then
      weight=1
      return
    end if
    ! Complement form near the upper end avoids loss of tiny shell weights.
    if (t <= .5_8) then
      weight=t**3*(10+t*(-15+6*t))
    else
      weight=1-(1-t)**3*(10+(1-t)*(-15+6*(1-t)))
    end if
    derivative=30*t*t*(1-t)**2/(2*halfwidth)
  end subroutine polarization_switch

  pure subroutine smooth_rank_energy(weight, angle, strength, offset, force_constant, bandwidth, energy, dw, da)
    ! Smooth weighted empirical angular ranks. No sorting or integer population.
    ! Inputs must have nonnegative weights, positive bandwidth and |1.5*strength|<1.
    ! dw and da differentiate the COMPLETE shell energy, including all targets.
    real(8), intent(in) :: weight(:), angle(:), strength(:), offset, force_constant, bandwidth
    real(8), intent(out) :: energy(size(strength)), dw(size(weight),size(strength)), da(size(weight),size(strength))
    real(8) :: population, rank(size(weight)), base(size(weight)), z, h, hp, displacement, target, adjoint, term
    real(8) :: kernel(size(weight),size(weight)), kernel_prime(size(weight),size(weight))
    real(8), parameter :: pi=3.1415926535897932384626433832795_8
    integer :: i,j,s,n
    energy=0; dw=0; da=0
    population=sum(weight); n=size(weight)
    if (population <= 0) return
    do i=1,n
      do j=1,n
        z=(angle(i)-angle(j))/bandwidth
        ! atan2 avoids cancellation of .5+atan(z)/pi for negative large z.
        h=atan2(1.0_8,-z)/pi
        hp=1/(pi*bandwidth*(1+z*z))
        kernel(i,j)=h; kernel_prime(i,j)=hp
      end do
      rank(i)=dot_product(weight,kernel(i,:))/population
      base(i)=acos(1-2*rank(i))
    end do
    do s=1,size(strength)
      do i=1,n
        target=base(i)-1.5_8*strength(s)*sin(base(i))
        displacement=angle(i)-target+offset
        energy(s)=energy(s)+.5_8*force_constant*weight(i)*displacement**2
        dw(i,s)=dw(i,s)+.5_8*force_constant*displacement**2
        da(i,s)=da(i,s)+force_constant*weight(i)*displacement
        adjoint=-force_constant*weight(i)*displacement*(1-1.5_8*strength(s)*cos(base(i))) &
          /sqrt(rank(i)*(1-rank(i)))
        do j=1,n
          dw(j,s)=dw(j,s)+adjoint*(kernel(i,j)-rank(i))/population
          term=adjoint*weight(j)*kernel_prime(i,j)/population
          da(i,s)=da(i,s)+term
          da(j,s)=da(j,s)-term
        end do
      end do
    end do
  end subroutine smooth_rank_energy

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
