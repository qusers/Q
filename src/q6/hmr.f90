module hmr
!!-------------------------------------------------------------------------------
!!  Hydrogen mass repartitioning on per-atom topology masses.
!!
!!  The caller selects the atoms and bonds to transform. Qprep passes only the
!!  solute slice, which deliberately leaves rigid solvent water unchanged.
!!-------------------------------------------------------------------------------
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: HMR_OK = 0
  integer, parameter :: HMR_INVALID_TARGET = 1
  integer, parameter :: HMR_INVALID_HYDROGEN_MASS = 2
  integer, parameter :: HMR_PARENT_NOT_UNIQUE = 3
  integer, parameter :: HMR_PARENT_TOO_LIGHT = 4
  integer, parameter :: HMR_MASS_NOT_CONSERVED = 5

contains

subroutine repartition_hydrogen_masses(mass, is_hydrogen, bond_i, bond_j, target_mass, &
                                       status, failed_atom, hydrogen_count, parent_count, &
                                       minimum_parent_mass)
  real(8), intent(inout) :: mass(:)
  logical, intent(in) :: is_hydrogen(:)
  integer, intent(in) :: bond_i(:), bond_j(:)
  real(8), intent(in) :: target_mass
  integer, intent(out) :: status, failed_atom, hydrogen_count, parent_count
  real(8), intent(out) :: minimum_parent_mass

  real(8), allocatable :: trial_mass(:)
  logical, allocatable :: parent_changed(:)
  real(8) :: delta, mass_before, mass_after, tolerance
  integer :: atom, bond, parent, candidates

  status = HMR_OK
  failed_atom = 0
  hydrogen_count = 0
  parent_count = 0
  minimum_parent_mass = huge(1.0d0)

  if (size(mass) /= size(is_hydrogen) .or. size(bond_i) /= size(bond_j) .or. &
      .not. ieee_is_finite(target_mass) .or. target_mass <= 0.0d0) then
    status = HMR_INVALID_TARGET
    return
  end if

  allocate(trial_mass(size(mass)), parent_changed(size(mass)))
  trial_mass = mass
  parent_changed = .false.
  mass_before = sum(mass)

  do atom = 1, size(mass)
    if (.not. is_hydrogen(atom)) cycle

    if (.not. ieee_is_finite(mass(atom)) .or. mass(atom) <= 0.0d0 .or. &
        mass(atom) >= target_mass) then
      status = HMR_INVALID_HYDROGEN_MASS
      failed_atom = atom
      return
    end if

    parent = 0
    candidates = 0
    do bond = 1, size(bond_i)
      if (bond_i(bond) == atom) then
        if (bond_j(bond) >= 1 .and. bond_j(bond) <= size(mass)) then
          if (.not. is_hydrogen(bond_j(bond))) then
            parent = bond_j(bond)
            candidates = candidates + 1
          end if
        end if
      else if (bond_j(bond) == atom) then
        if (bond_i(bond) >= 1 .and. bond_i(bond) <= size(mass)) then
          if (.not. is_hydrogen(bond_i(bond))) then
            parent = bond_i(bond)
            candidates = candidates + 1
          end if
        end if
      end if
    end do

    if (candidates /= 1) then
      status = HMR_PARENT_NOT_UNIQUE
      failed_atom = atom
      return
    end if

    delta = target_mass - mass(atom)
    trial_mass(atom) = target_mass
    trial_mass(parent) = trial_mass(parent) - delta
    parent_changed(parent) = .true.
    hydrogen_count = hydrogen_count + 1
  end do

  do atom = 1, size(mass)
    if (.not. parent_changed(atom)) cycle
    parent_count = parent_count + 1
    minimum_parent_mass = min(minimum_parent_mass, trial_mass(atom))
    if (.not. ieee_is_finite(trial_mass(atom)) .or. trial_mass(atom) <= target_mass) then
      status = HMR_PARENT_TOO_LIGHT
      failed_atom = atom
      return
    end if
  end do

  if (parent_count == 0) minimum_parent_mass = 0.0d0
  mass_after = sum(trial_mass)
  tolerance = 1.0d-10*max(1.0d0, abs(mass_before))
  if (.not. ieee_is_finite(mass_after) .or. abs(mass_after-mass_before) > tolerance) then
    status = HMR_MASS_NOT_CONSERVED
    return
  end if

  mass = trial_mass
end subroutine repartition_hydrogen_masses

end module hmr
