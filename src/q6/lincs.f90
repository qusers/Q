module lincs
!!-------------------------------------------------------------------------------
!!  LINCS constraint algorithm for solute bonds
!!
!!  Hess B, Bekker H, Berendsen HJC, Fraaije JGEM (1997). "LINCS: A linear
!!  constraint solver for molecular simulations." J Comput Chem 18(12):1463-1472.
!!
!!  P-LINCS extension for better handling of coupled constraints:
!!  Hess B (2008). "P-LINCS: A parallel linear constraint solver for molecular
!!  simulation." J Chem Theory Comput 4(1):116-122.
!!-------------------------------------------------------------------------------
  implicit none

  type lincs_data_type
    integer :: nconstraints = 0           ! total solute constraints
    integer, allocatable :: atom_i(:)     ! first atom of each constraint
    integer, allocatable :: atom_j(:)     ! second atom of each constraint
    real(8), allocatable :: target_len(:) ! target bond lengths
    real(8), allocatable :: inv_mass_i(:) ! cached inverse mass of atom i
    real(8), allocatable :: inv_mass_j(:) ! cached inverse mass of atom j
    real(8), allocatable :: sdiag(:)      ! 1/sqrt(mi) + 1/sqrt(mj) diagonal factor

    ! coupling matrix in CSR format
    integer, allocatable :: coupled_start(:) ! (nconstraints+1) row pointers
    integer, allocatable :: coupled_index(:) ! column indices of coupled constraints
    real(8), allocatable :: coupled_coeff(:) ! coupling coefficients

    ! FEP constraints with lambda-dependent distances
    integer :: nfep = 0
    integer, allocatable :: fep_idx(:)    ! indices into constraint arrays

    integer :: lincs_order = 4            ! expansion order (default 4)
  end type lincs_data_type

  type(lincs_data_type) :: lincs_data
  logical :: use_lincs = .false.

contains

subroutine setup_lincs(nc, ai, aj, target_len, inv_mi, inv_mj, nfep, fep_idx, order)
!!-------------------------------------------------------------------------------
!!  Initialize LINCS constraint data.
!!
!!  nc         : number of constraints
!!  ai, aj     : atom index pairs (1-based)
!!  target_len : equilibrium bond lengths
!!  inv_mi/mj  : inverse masses
!!  nfep       : number of FEP constraints
!!  fep_idx    : indices of FEP constraints in the arrays
!!  order      : LINCS expansion order (default 4)
!!-------------------------------------------------------------------------------
  integer, intent(in) :: nc
  integer, intent(in) :: ai(:), aj(:)
  real(8), intent(in) :: target_len(:)
  real(8), intent(in) :: inv_mi(:), inv_mj(:)
  integer, intent(in) :: nfep
  integer, intent(in) :: fep_idx(:)
  integer, intent(in) :: order

  integer :: k, m, ik, jk, im, jm, ncoupled, idx
  real(8) :: total_inv_mass

  lincs_data%nconstraints = nc
  lincs_data%lincs_order = order
  lincs_data%nfep = nfep

  if (nc == 0) return

  ! allocate arrays
  allocate(lincs_data%atom_i(nc))
  allocate(lincs_data%atom_j(nc))
  allocate(lincs_data%target_len(nc))
  allocate(lincs_data%inv_mass_i(nc))
  allocate(lincs_data%inv_mass_j(nc))
  allocate(lincs_data%sdiag(nc))
  allocate(lincs_data%coupled_start(nc + 1))

  ! copy constraint data
  lincs_data%atom_i(1:nc) = ai(1:nc)
  lincs_data%atom_j(1:nc) = aj(1:nc)
  lincs_data%target_len(1:nc) = target_len(1:nc)
  lincs_data%inv_mass_i(1:nc) = inv_mi(1:nc)
  lincs_data%inv_mass_j(1:nc) = inv_mj(1:nc)

  ! compute diagonal factors: sdiag(k) = 1 / (inv_mi(k) + inv_mj(k))
  do k = 1, nc
    total_inv_mass = inv_mi(k) + inv_mj(k)
    if (total_inv_mass > 0.0d0) then
      lincs_data%sdiag(k) = 1.0d0 / total_inv_mass
    else
      lincs_data%sdiag(k) = 0.0d0
    end if
  end do

  ! FEP indices
  if (nfep > 0) then
    allocate(lincs_data%fep_idx(nfep))
    lincs_data%fep_idx(1:nfep) = fep_idx(1:nfep)
  end if

  ! build coupling matrix in CSR format
  ! two constraints are coupled if they share an atom
  ! first pass: count coupled pairs per constraint
  lincs_data%coupled_start(1) = 1
  do k = 1, nc
    ncoupled = 0
    ik = lincs_data%atom_i(k)
    jk = lincs_data%atom_j(k)
    do m = 1, nc
      if (m == k) cycle
      im = lincs_data%atom_i(m)
      jm = lincs_data%atom_j(m)
      if (ik == im .or. ik == jm .or. jk == im .or. jk == jm) then
        ncoupled = ncoupled + 1
      end if
    end do
    lincs_data%coupled_start(k+1) = lincs_data%coupled_start(k) + ncoupled
  end do

  ! allocate CSR arrays
  ncoupled = lincs_data%coupled_start(nc+1) - 1
  allocate(lincs_data%coupled_index(ncoupled))
  allocate(lincs_data%coupled_coeff(ncoupled))

  ! second pass: fill in indices and coefficients
  do k = 1, nc
    idx = lincs_data%coupled_start(k)
    ik = lincs_data%atom_i(k)
    jk = lincs_data%atom_j(k)
    do m = 1, nc
      if (m == k) cycle
      im = lincs_data%atom_i(m)
      jm = lincs_data%atom_j(m)
      ! coupling coefficient: sign depends on which atoms are shared
      ! A_km = -w_shared / (w_ik + w_jk) where w = inverse mass
      if (ik == im) then
        ! constraints k and m share atom ik=im
        lincs_data%coupled_index(idx) = m
        lincs_data%coupled_coeff(idx) = -inv_mi(k) * lincs_data%sdiag(k)
        idx = idx + 1
      else if (ik == jm) then
        lincs_data%coupled_index(idx) = m
        lincs_data%coupled_coeff(idx) = -inv_mi(k) * lincs_data%sdiag(k)
        idx = idx + 1
      else if (jk == im) then
        lincs_data%coupled_index(idx) = m
        lincs_data%coupled_coeff(idx) = -inv_mj(k) * lincs_data%sdiag(k)
        idx = idx + 1
      else if (jk == jm) then
        lincs_data%coupled_index(idx) = m
        lincs_data%coupled_coeff(idx) = -inv_mj(k) * lincs_data%sdiag(k)
        idx = idx + 1
      end if
    end do
  end do

  use_lincs = .true.

end subroutine setup_lincs


subroutine lincs_positions(xx, x)
!!-------------------------------------------------------------------------------
!!  Apply LINCS constraints to positions.
!!
!!  xx : old (pre-step) positions, dimension(3*natom)
!!  x  : unconstrained positions to be corrected, dimension(3*natom)
!!
!!  Algorithm (Hess et al. 1997):
!!  1. Compute unconstrained bond vectors and project onto old directions
!!  2. Apply coupling matrix correction (order N expansion)
!!  3. Correct positions
!!  4. Apply rotation correction (second pass)
!!-------------------------------------------------------------------------------
  real(8), intent(in) :: xx(:)
  real(8), intent(inout) :: x(:)

  integer :: nc, k, m, iter, i3, j3, idx
  real(8) :: rhs(lincs_data%nconstraints)
  real(8) :: sol(lincs_data%nconstraints)
  real(8) :: tmp(lincs_data%nconstraints)
  real(8) :: dir(3, lincs_data%nconstraints)  ! unit bond directions from old positions
  real(8) :: bond(3)
  real(8) :: len2, invlen, proj, p
  real(8) :: correction

  nc = lincs_data%nconstraints
  if (nc == 0) return

  ! Step 1: compute bond unit vectors from old positions
  do k = 1, nc
    i3 = lincs_data%atom_i(k) * 3 - 3
    j3 = lincs_data%atom_j(k) * 3 - 3

    bond(1) = xx(i3+1) - xx(j3+1)
    bond(2) = xx(i3+2) - xx(j3+2)
    bond(3) = xx(i3+3) - xx(j3+3)

    len2 = bond(1)**2 + bond(2)**2 + bond(3)**2
    invlen = 1.0d0 / sqrt(len2)

    dir(1,k) = bond(1) * invlen
    dir(2,k) = bond(2) * invlen
    dir(3,k) = bond(3) * invlen
  end do

  ! Step 2: compute RHS = B * (x_unconstrained - x_old) projected onto constraint directions
  ! rhs(k) = sdiag(k) * (d_k . (r_ij_new - r_ij_old))
  ! Actually, rhs is based on the difference between current and target length:
  ! rhs(k) = sdiag(k) * (projection of new bond onto old direction - target_length)
  do k = 1, nc
    i3 = lincs_data%atom_i(k) * 3 - 3
    j3 = lincs_data%atom_j(k) * 3 - 3

    bond(1) = x(i3+1) - x(j3+1)
    bond(2) = x(i3+2) - x(j3+2)
    bond(3) = x(i3+3) - x(j3+3)

    ! project new bond vector onto old direction
    proj = bond(1)*dir(1,k) + bond(2)*dir(2,k) + bond(3)*dir(3,k)

    rhs(k) = lincs_data%sdiag(k) * (proj - lincs_data%target_len(k))
  end do

  ! Step 3: solve (I - A) * sol = rhs using iterative expansion
  ! sol = rhs + A*rhs + A^2*rhs + ... (lincs_order terms)
  sol = rhs
  tmp = rhs
  do iter = 1, lincs_data%lincs_order
    ! tmp = A * tmp
    rhs = tmp  ! save for matrix-vector multiply
    tmp = 0.0d0
    do k = 1, nc
      do idx = lincs_data%coupled_start(k), lincs_data%coupled_start(k+1) - 1
        m = lincs_data%coupled_index(idx)
        ! A_km coefficient times direction dot product
        p = lincs_data%coupled_coeff(idx) * &
            (dir(1,k)*dir(1,m) + dir(2,k)*dir(2,m) + dir(3,k)*dir(3,m))
        tmp(k) = tmp(k) + p * rhs(m)
      end do
    end do
    sol = sol + tmp
  end do

  ! Step 4: apply position corrections
  do k = 1, nc
    i3 = lincs_data%atom_i(k) * 3 - 3
    j3 = lincs_data%atom_j(k) * 3 - 3

    correction = sol(k)

    x(i3+1) = x(i3+1) - correction * lincs_data%inv_mass_i(k) * dir(1,k)
    x(i3+2) = x(i3+2) - correction * lincs_data%inv_mass_i(k) * dir(2,k)
    x(i3+3) = x(i3+3) - correction * lincs_data%inv_mass_i(k) * dir(3,k)

    x(j3+1) = x(j3+1) + correction * lincs_data%inv_mass_j(k) * dir(1,k)
    x(j3+2) = x(j3+2) + correction * lincs_data%inv_mass_j(k) * dir(2,k)
    x(j3+3) = x(j3+3) + correction * lincs_data%inv_mass_j(k) * dir(3,k)
  end do

  ! Step 5: rotation correction (iterated for coupled constraints)
  ! Scale each bond to exact target length, distributed by inverse mass.
  ! For coupled constraints, one correction affects neighbors, so iterate.
  do iter = 1, 3
    do k = 1, nc
      i3 = lincs_data%atom_i(k) * 3 - 3
      j3 = lincs_data%atom_j(k) * 3 - 3

      bond(1) = x(i3+1) - x(j3+1)
      bond(2) = x(i3+2) - x(j3+2)
      bond(3) = x(i3+3) - x(j3+3)

      len2 = bond(1)**2 + bond(2)**2 + bond(3)**2
      invlen = 1.0d0 / sqrt(len2)
      p = lincs_data%sdiag(k) * (1.0d0/invlen - lincs_data%target_len(k))

      bond = bond * invlen  ! normalize

      x(i3+1) = x(i3+1) - p * lincs_data%inv_mass_i(k) * bond(1)
      x(i3+2) = x(i3+2) - p * lincs_data%inv_mass_i(k) * bond(2)
      x(i3+3) = x(i3+3) - p * lincs_data%inv_mass_i(k) * bond(3)

      x(j3+1) = x(j3+1) + p * lincs_data%inv_mass_j(k) * bond(1)
      x(j3+2) = x(j3+2) + p * lincs_data%inv_mass_j(k) * bond(2)
      x(j3+3) = x(j3+3) + p * lincs_data%inv_mass_j(k) * bond(3)
    end do
  end do

end subroutine lincs_positions

end module lincs
