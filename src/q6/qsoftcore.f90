!------------------------------------------------------------------------------!
! Shared soft-core pair potentials for Q.                                      !
!------------------------------------------------------------------------------!
module qsoftcore
  implicit none
  private

  ! Gapsys, Seeliger & de Groot (JCTC 2012) defaults converted from nm to
  ! Q's internal distance unit (angstrom). alpha_LJ and sigma_Q are dimensionless.
  real(8), parameter, public :: GAPSYS_ALPHA_LJ_DEFAULT = 0.85_8
  real(8), parameter, public :: GAPSYS_ALPHA_Q_ANGSTROM_DEFAULT = 3.0_8
  real(8), parameter, public :: GAPSYS_SIGMA_LJ_ANGSTROM_DEFAULT = 3.0_8
  real(8), parameter, public :: GAPSYS_SIGMA_Q_DEFAULT = 1.0_8

  public :: gapsys_lj_radius
  public :: gapsys_coulomb_radius
  public :: eval_gapsys_lj
  public :: eval_gapsys_coulomb
  public :: gapsys_lj_energy_difference
  public :: gapsys_coulomb_energy_difference

contains

pure real(8) function gapsys_lj_radius(alpha_lj, sigma_lj, decoupling)
  real(8), intent(in) :: alpha_lj, sigma_lj, decoupling
  real(8) :: omw

  omw = max(0.0_8, decoupling)
  if(alpha_lj <= 0.0_8 .or. sigma_lj <= 0.0_8 .or. omw <= 0.0_8) then
    gapsys_lj_radius = 0.0_8
  else
    gapsys_lj_radius = alpha_lj * ((26.0_8/7.0_8) * sigma_lj**6 * omw)**(1.0_8/6.0_8)
  end if
end function gapsys_lj_radius


pure real(8) function gapsys_coulomb_radius(alpha_q, sigma_q, qi_e, qj_e, decoupling)
  ! qi_e and qj_e are physical partial charges in elementary-charge units.
  ! They must not include Q's sqrt(coulomb_constant) energy scaling.
  real(8), intent(in) :: alpha_q, sigma_q, qi_e, qj_e, decoupling
  real(8) :: omw

  omw = max(0.0_8, decoupling)
  if(alpha_q <= 0.0_8 .or. omw <= 0.0_8) then
    gapsys_coulomb_radius = 0.0_8
  else
    gapsys_coulomb_radius = alpha_q * (1.0_8 + sigma_q*abs(qi_e*qj_e)) * &
                             omw**(1.0_8/6.0_8)
  end if
end function gapsys_coulomb_radius


pure subroutine eval_gapsys_lj(c12, c6, r_sc, r2, energy, dVdr_over_r, exact_overlap)
  ! Evaluate LJ energy and the Cartesian-force coefficient dV/dr/r.
  ! A caller using displacement x_j-x_i can retain Q's existing updates:
  !   d_i = d_i - (dV/dr/r) * displacement
  !   d_j = d_j + (dV/dr/r) * displacement
  real(8), intent(in)  :: c12, c6, r_sc, r2
  real(8), intent(out) :: energy, dVdr_over_r
  logical, intent(out) :: exact_overlap
  real(8) :: r, invr2, invr6, va, vb
  real(8) :: inv, inv2, inv6, inv7, inv8, inv12, inv13, inv14
  real(8) :: a_lj, b_lj, c_lj

  exact_overlap = (r2 <= 0.0_8)
  if(c12 == 0.0_8 .and. c6 == 0.0_8) then
    energy = 0.0_8
    dVdr_over_r = 0.0_8
    return
  end if

  if(exact_overlap) then
    if(r_sc > 0.0_8) then
      inv   = 1.0_8/r_sc
      inv2  = inv*inv
      inv6  = inv2*inv2*inv2
      inv12 = inv6*inv6
      energy = 91.0_8*c12*inv12 - 28.0_8*c6*inv6
    else if(c12 == 0.0_8 .and. c6 == 0.0_8) then
      energy = 0.0_8
    else
      energy = huge(energy)
    end if
    ! Gapsys has a finite radial force but no defined direction at exact overlap.
    dVdr_over_r = 0.0_8
    return
  end if

  r = sqrt(r2)
  if(r_sc > 0.0_8 .and. r < r_sc) then
    inv   = 1.0_8/r_sc
    inv2  = inv*inv
    inv6  = inv2*inv2*inv2
    inv7  = inv6*inv
    inv8  = inv6*inv2
    inv12 = inv6*inv6
    inv13 = inv12*inv
    inv14 = inv12*inv2
    a_lj  = 78.0_8*c12*inv14 - 21.0_8*c6*inv8
    b_lj  = 168.0_8*c12*inv13 - 48.0_8*c6*inv7
    c_lj  = 91.0_8*c12*inv12 - 28.0_8*c6*inv6
    energy = a_lj*r2 - b_lj*r + c_lj
    dVdr_over_r = (2.0_8*a_lj*r - b_lj)/r
  else
    invr2 = 1.0_8/r2
    invr6 = invr2*invr2*invr2
    va = c12*invr6*invr6
    vb = c6*invr6
    energy = va-vb
    dVdr_over_r = (-12.0_8*va + 6.0_8*vb)*invr2
  end if
end subroutine eval_gapsys_lj


pure real(8) function gapsys_lj_energy(c12, c6, r_sc, r2)
  ! Energy-only Gapsys LJ evaluator used by foreign-lambda work corrections.
  ! Callers screen the common hard/hard case before entering this routine.
  real(8), intent(in) :: c12, c6, r_sc, r2
  real(8) :: r, invr2, invr6, inv, inv2, inv6, inv7, inv8, inv12, inv13, inv14
  real(8) :: a_lj, b_lj, c_lj

  if(c12 == 0.0_8 .and. c6 == 0.0_8) then
    gapsys_lj_energy = 0.0_8
  else if(r2 <= 0.0_8) then
    if(r_sc <= 0.0_8) then
      gapsys_lj_energy = huge(gapsys_lj_energy)
    else
      inv = 1.0_8/r_sc
      inv2 = inv*inv
      inv6 = inv2*inv2*inv2
      inv12 = inv6*inv6
      gapsys_lj_energy = 91.0_8*c12*inv12 - 28.0_8*c6*inv6
    end if
  else
    r = sqrt(r2)
    if(r_sc > 0.0_8 .and. r < r_sc) then
      inv = 1.0_8/r_sc
      inv2 = inv*inv
      inv6 = inv2*inv2*inv2
      inv7 = inv6*inv
      inv8 = inv6*inv2
      inv12 = inv6*inv6
      inv13 = inv12*inv
      inv14 = inv12*inv2
      a_lj = 78.0_8*c12*inv14 - 21.0_8*c6*inv8
      b_lj = 168.0_8*c12*inv13 - 48.0_8*c6*inv7
      c_lj = 91.0_8*c12*inv12 - 28.0_8*c6*inv6
      gapsys_lj_energy = a_lj*r2 - b_lj*r + c_lj
    else
      invr2 = 1.0_8/r2
      invr6 = invr2*invr2*invr2
      gapsys_lj_energy = c12*invr6*invr6 - c6*invr6
    end if
  end if
end function gapsys_lj_energy


pure real(8) function gapsys_coulomb_energy(qq_energy, r_sc, r2)
  ! Energy-only Gapsys Coulomb evaluator used by work corrections.
  real(8), intent(in) :: qq_energy, r_sc, r2
  real(8) :: r, invr, inv, inv2, inv3, a_q, b_q, c_q

  if(qq_energy == 0.0_8) then
    gapsys_coulomb_energy = 0.0_8
  else if(r2 <= 0.0_8) then
    if(r_sc <= 0.0_8) then
      gapsys_coulomb_energy = sign(huge(gapsys_coulomb_energy), qq_energy)
    else
      gapsys_coulomb_energy = 3.0_8*qq_energy/r_sc
    end if
  else
    r = sqrt(r2)
    if(r_sc > 0.0_8 .and. r < r_sc) then
      inv = 1.0_8/r_sc
      inv2 = inv*inv
      inv3 = inv2*inv
      a_q = qq_energy*inv3
      b_q = 3.0_8*qq_energy*inv2
      c_q = 3.0_8*qq_energy*inv
      gapsys_coulomb_energy = a_q*r2 - b_q*r + c_q
    else
      invr = 1.0_8/r
      gapsys_coulomb_energy = qq_energy*invr
    end if
  end if
end function gapsys_coulomb_energy


pure real(8) function gapsys_lj_energy_difference(c12, c6, r_sc_source, r_sc_target, r2)
  ! Exact target-shape minus source-shape energy. Above both switching radii
  ! both Hamiltonians use the same hard potential, so the correction is zero.
  real(8), intent(in) :: c12, c6, r_sc_source, r_sc_target, r2

  if((c12 == 0.0_8 .and. c6 == 0.0_8) .or. &
     r2 >= max(r_sc_source, r_sc_target)**2) then
    gapsys_lj_energy_difference = 0.0_8
  else
    gapsys_lj_energy_difference = gapsys_lj_energy(c12,c6,r_sc_target,r2) - &
                                  gapsys_lj_energy(c12,c6,r_sc_source,r2)
  end if
end function gapsys_lj_energy_difference


pure real(8) function gapsys_coulomb_energy_difference(qq_energy, r_sc_source, r_sc_target, r2)
  ! Exact target-shape minus source-shape Coulomb energy with hard/hard sparsity.
  real(8), intent(in) :: qq_energy, r_sc_source, r_sc_target, r2

  if(qq_energy == 0.0_8 .or. r2 >= max(r_sc_source, r_sc_target)**2) then
    gapsys_coulomb_energy_difference = 0.0_8
  else
    gapsys_coulomb_energy_difference = &
      gapsys_coulomb_energy(qq_energy,r_sc_target,r2) - &
      gapsys_coulomb_energy(qq_energy,r_sc_source,r2)
  end if
end function gapsys_coulomb_energy_difference


pure subroutine eval_gapsys_coulomb(qq_energy, r_sc, r2, energy, dVdr_over_r, exact_overlap)
  ! qq_energy is the fully energy-scaled charge product, including any requested
  ! 1-4 factor. Radius selection is deliberately separate and uses physical charges.
  real(8), intent(in)  :: qq_energy, r_sc, r2
  real(8), intent(out) :: energy, dVdr_over_r
  logical, intent(out) :: exact_overlap
  real(8) :: r, invr, invr2, inv, inv2, inv3
  real(8) :: a_q, b_q, c_q

  exact_overlap = (r2 <= 0.0_8)
  if(qq_energy == 0.0_8) then
    energy = 0.0_8
    dVdr_over_r = 0.0_8
    return
  end if

  if(exact_overlap) then
    if(r_sc > 0.0_8) then
      energy = 3.0_8*qq_energy/r_sc
    else if(qq_energy == 0.0_8) then
      energy = 0.0_8
    else
      energy = sign(huge(energy), qq_energy)
    end if
    dVdr_over_r = 0.0_8
    return
  end if

  r = sqrt(r2)
  if(r_sc > 0.0_8 .and. r < r_sc) then
    inv  = 1.0_8/r_sc
    inv2 = inv*inv
    inv3 = inv2*inv
    a_q  = qq_energy*inv3
    b_q  = 3.0_8*qq_energy*inv2
    c_q  = 3.0_8*qq_energy*inv
    energy = a_q*r2 - b_q*r + c_q
    dVdr_over_r = (2.0_8*a_q*r - b_q)/r
  else
    invr = 1.0_8/r
    invr2 = 1.0_8/r2
    energy = qq_energy*invr
    dVdr_over_r = -energy*invr2
  end if
end subroutine eval_gapsys_coulomb

end module qsoftcore
