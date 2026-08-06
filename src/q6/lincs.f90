module lincs
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!!-------------------------------------------------------------------------------
!!  LINear Constraint Solver (LINCS) for bond constraints.
!!
!!  B. Hess, H. Bekker, H. J. C. Berendsen, and J. G. E. M. Fraaije,
!!  J. Comput. Chem. 18, 1463-1472 (1997).
!!
!!  The normalized sparse coupling matrix is inverted through the truncated
!!  Neumann series described in equations 20-22 of the paper.  Nonlinear bond
!!  rotation is handled by the paper's projected-length correction.
!!-------------------------------------------------------------------------------
  implicit none
  private

  type lincs_data_type
    integer :: constraint_count = 0
    integer :: atom_count = 0
    integer :: expansion_order = 4
    integer :: rotation_iterations = 1
    integer, allocatable :: atom_i(:)
    integer, allocatable :: atom_j(:)
    real(8), allocatable :: target_length(:)
    real(8), allocatable :: inverse_mass(:)
    real(8), allocatable :: scale(:)
    integer, allocatable :: coupled_start(:)
    integer, allocatable :: coupled_index(:)
    real(8), allocatable :: coupled_prefactor(:)
    real(8), allocatable :: coupling(:)
    real(8), allocatable :: direction(:,:)
    real(8), allocatable :: rhs_a(:)
    real(8), allocatable :: rhs_b(:)
    real(8), allocatable :: solution(:)
  end type lincs_data_type

  type(lincs_data_type) :: data

  public :: setup_lincs, initialize_lincs_positions, lincs_positions, lincs_is_active

contains

subroutine setup_lincs(atom_i, atom_j, target_length, inverse_mass, &
                       expansion_order, rotation_iterations)
!!-------------------------------------------------------------------------------
!!  Cache the constraint topology and build its sparse coupling graph in O(K+C)
!!  storage, where K is the number of constraints and C the number of directly
!!  coupled constraint pairs.
!!-------------------------------------------------------------------------------
  integer, intent(in) :: atom_i(:), atom_j(:)
  real(8), intent(in) :: target_length(:), inverse_mass(:)
  integer, intent(in), optional :: expansion_order, rotation_iterations

  integer :: constraints, atoms, k, m, atom, other_atom, total_couplings
  integer :: position
  integer, allocatable :: degree(:), atom_start(:), atom_cursor(:)
  integer, allocatable :: atom_constraints(:)
  real(8) :: inverse_mass_sum

  call reset_lincs

  constraints = size(atom_i)
  atoms = size(inverse_mass)
  if (size(atom_j) /= constraints .or. size(target_length) /= constraints) then
    error stop 'LINCS constraint arrays have inconsistent sizes'
  end if
  if (present(expansion_order)) then
    if (expansion_order < 0) error stop 'LINCS expansion order must be nonnegative'
    data%expansion_order = expansion_order
  end if
  if (present(rotation_iterations)) then
    if (rotation_iterations < 0) error stop 'LINCS rotation iterations must be nonnegative'
    data%rotation_iterations = rotation_iterations
  end if

  data%constraint_count = constraints
  data%atom_count = atoms
  if (constraints == 0) return

  allocate(data%atom_i(constraints), data%atom_j(constraints), &
           data%target_length(constraints), data%inverse_mass(atoms), &
           data%scale(constraints), data%coupled_start(constraints+1), &
           data%direction(3,constraints), data%rhs_a(constraints), &
           data%rhs_b(constraints), data%solution(constraints))
  data%atom_i = atom_i
  data%atom_j = atom_j
  data%target_length = target_length
  data%inverse_mass = inverse_mass

  allocate(degree(atoms), atom_start(atoms+1), atom_cursor(atoms))
  degree = 0
  do k = 1, constraints
    if (atom_i(k) < 1 .or. atom_i(k) > atoms .or. &
        atom_j(k) < 1 .or. atom_j(k) > atoms .or. atom_i(k) == atom_j(k)) then
      error stop 'LINCS constraint contains invalid atom indices'
    end if
    if (target_length(k) <= 0.0d0) error stop 'LINCS target lengths must be positive'
    inverse_mass_sum = inverse_mass(atom_i(k)) + inverse_mass(atom_j(k))
    if (inverse_mass_sum <= 0.0d0) then
      error stop 'LINCS constraint must contain at least one mobile atom'
    end if
    data%scale(k) = 1.0d0/sqrt(inverse_mass_sum)
    degree(atom_i(k)) = degree(atom_i(k)) + 1
    degree(atom_j(k)) = degree(atom_j(k)) + 1
  end do

  atom_start(1) = 1
  do atom = 1, atoms
    atom_start(atom+1) = atom_start(atom) + degree(atom)
  end do
  allocate(atom_constraints(2*constraints))
  atom_cursor = atom_start(1:atoms)
  do k = 1, constraints
    atom = atom_i(k)
    atom_constraints(atom_cursor(atom)) = k
    atom_cursor(atom) = atom_cursor(atom) + 1
    atom = atom_j(k)
    atom_constraints(atom_cursor(atom)) = k
    atom_cursor(atom) = atom_cursor(atom) + 1
  end do

  ! Reject duplicate atom pairs.  Without this check they would appear twice in
  ! each other's coupling rows and make the matrix representation ambiguous.
  do k = 1, constraints
    atom = atom_i(k)
    other_atom = atom_j(k)
    do position = atom_start(atom), atom_start(atom+1)-1
      m = atom_constraints(position)
      if (m >= k) cycle
      if ((atom_i(m) == atom .and. atom_j(m) == other_atom) .or. &
          (atom_j(m) == atom .and. atom_i(m) == other_atom)) then
        error stop 'LINCS received a duplicate constraint atom pair'
      end if
    end do
  end do

  data%coupled_start(1) = 1
  do k = 1, constraints
    total_couplings = degree(atom_i(k)) + degree(atom_j(k)) - 2
    data%coupled_start(k+1) = data%coupled_start(k) + total_couplings
  end do
  total_couplings = data%coupled_start(constraints+1)-1
  allocate(data%coupled_index(total_couplings), &
           data%coupled_prefactor(total_couplings), data%coupling(total_couplings))

  do k = 1, constraints
    position = data%coupled_start(k)
    call append_atom_couplings(k, atom_i(k), atom_start, atom_constraints, position)
    call append_atom_couplings(k, atom_j(k), atom_start, atom_constraints, position)
    if (position /= data%coupled_start(k+1)) then
      error stop 'LINCS internal coupling-list size mismatch'
    end if
  end do

  deallocate(degree, atom_start, atom_cursor, atom_constraints)

contains

subroutine append_atom_couplings(row, shared_atom, starts, incident, output_position)
  integer, intent(in) :: row, shared_atom, starts(:), incident(:)
  integer, intent(inout) :: output_position
  integer :: entry, column, sign_factor
  logical :: row_first, column_first

  do entry = starts(shared_atom), starts(shared_atom+1)-1
    column = incident(entry)
    if (column == row) cycle
    data%coupled_index(output_position) = column
    row_first = data%atom_i(row) == shared_atom
    column_first = data%atom_i(column) == shared_atom
    if (row_first .eqv. column_first) then
      sign_factor = -1
    else
      sign_factor = 1
    end if
    data%coupled_prefactor(output_position) = &
      dble(sign_factor)*data%inverse_mass(shared_atom)* &
      data%scale(row)*data%scale(column)
    output_position = output_position + 1
  end do
end subroutine append_atom_couplings

end subroutine setup_lincs


logical function lincs_is_active()
  lincs_is_active = data%constraint_count > 0
end function lincs_is_active


logical function initialize_lincs_positions(coordinates, failed_constraint, &
                                            max_relative_error) result(success)
!!-------------------------------------------------------------------------------
!!  Project an arbitrary starting structure onto the LINCS constraint manifold.
!!  LINCS updates assume the reference coordinates are already constrained, so
!!  this one-time tightly converged SHAKE projection establishes that invariant.
!!-------------------------------------------------------------------------------
  real(8), intent(inout) :: coordinates(:)
  integer, intent(out), optional :: failed_constraint
  real(8), intent(out), optional :: max_relative_error

  integer, parameter :: maximum_iterations = 1000
  real(8), parameter :: convergence_tolerance = 1.0d-12
  integer :: iteration, k, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: reference(3), current(3), reference_dot_current
  real(8) :: difference, correction, relative_error, maximum_error

  success = .true.
  if (present(failed_constraint)) failed_constraint = 0
  if (present(max_relative_error)) max_relative_error = 0.0d0
  if (.not. lincs_is_active()) return

  ! Keep the initial bond directions fixed, as in SHAKE and finite-step LINCS.
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    reference = coordinates(coordinate_i+1:coordinate_i+3) - &
                coordinates(coordinate_j+1:coordinate_j+3)
    if (dot_product(reference,reference) <= tiny(1.0d0)) then
      call report_failure(k,failed_constraint)
      success = .false.
      return
    end if
    data%direction(:,k) = reference
  end do

  do iteration = 1, maximum_iterations
    do k = 1, data%constraint_count
      atom_i = data%atom_i(k)
      atom_j = data%atom_j(k)
      coordinate_i = 3*atom_i-3
      coordinate_j = 3*atom_j-3
      current = coordinates(coordinate_i+1:coordinate_i+3) - &
                coordinates(coordinate_j+1:coordinate_j+3)
      reference = data%direction(:,k)
      reference_dot_current = dot_product(reference,current)
      if (abs(reference_dot_current) <= tiny(1.0d0)) then
        call report_failure(k,failed_constraint)
        success = .false.
        return
      end if
      difference = data%target_length(k)**2-dot_product(current,current)
      correction = difference/(2.0d0*reference_dot_current* &
                   (data%inverse_mass(atom_i)+data%inverse_mass(atom_j)))
      coordinates(coordinate_i+1:coordinate_i+3) = &
        coordinates(coordinate_i+1:coordinate_i+3) + &
        data%inverse_mass(atom_i)*reference*correction
      coordinates(coordinate_j+1:coordinate_j+3) = &
        coordinates(coordinate_j+1:coordinate_j+3) - &
        data%inverse_mass(atom_j)*reference*correction
    end do

    maximum_error = 0.0d0
    do k = 1, data%constraint_count
      atom_i = data%atom_i(k)
      atom_j = data%atom_j(k)
      coordinate_i = 3*atom_i-3
      coordinate_j = 3*atom_j-3
      current = coordinates(coordinate_i+1:coordinate_i+3) - &
                coordinates(coordinate_j+1:coordinate_j+3)
      relative_error = abs(sqrt(dot_product(current,current))/ &
                           data%target_length(k)-1.0d0)
      maximum_error = max(maximum_error,relative_error)
    end do
    if (maximum_error <= convergence_tolerance) then
      if (present(max_relative_error)) max_relative_error = maximum_error
      return
    end if
  end do

  success = .false.
  if (present(max_relative_error)) max_relative_error = maximum_error
end function initialize_lincs_positions


logical function lincs_positions(x_old, x_new, failed_constraint, &
                                 max_relative_error) result(success)
!!-------------------------------------------------------------------------------
!!  Constrain one unconstrained position update.  Returns .false. if an old bond
!!  is degenerate or rotational lengthening exceeds LINCS' solution domain.
!!-------------------------------------------------------------------------------
  real(8), intent(in) :: x_old(:)
  real(8), intent(inout) :: x_new(:)
  integer, intent(out), optional :: failed_constraint
  real(8), intent(out), optional :: max_relative_error

  integer :: k, iteration, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: bond(3), length2, projected_length, radicand
  real(8) :: relative_error, maximum_error
  real(8), parameter :: domain_tolerance = 1.0d-12

  success = .true.
  if (present(failed_constraint)) failed_constraint = 0
  if (present(max_relative_error)) max_relative_error = 0.0d0
  if (.not. lincs_is_active()) return

  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = x_old(coordinate_i+1:coordinate_i+3) - &
           x_old(coordinate_j+1:coordinate_j+3)
    length2 = dot_product(bond, bond)
    if (length2 <= tiny(1.0d0)) then
      call report_failure(k, failed_constraint)
      success = .false.
      return
    end if
    data%direction(:,k) = bond/sqrt(length2)
  end do
  call update_coupling_matrix

  ! Linear projection, equations 19-22.
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = x_new(coordinate_i+1:coordinate_i+3) - &
           x_new(coordinate_j+1:coordinate_j+3)
    projected_length = dot_product(data%direction(:,k), bond)
    data%rhs_a(k) = data%scale(k)*(projected_length-data%target_length(k))
  end do
  call solve_expansion
  call apply_solution(x_new)

  ! Nonlinear correction for rotational lengthening, equations 17-18.
  do iteration = 1, data%rotation_iterations
    do k = 1, data%constraint_count
      atom_i = data%atom_i(k)
      atom_j = data%atom_j(k)
      coordinate_i = 3*atom_i-3
      coordinate_j = 3*atom_j-3
      bond = x_new(coordinate_i+1:coordinate_i+3) - &
             x_new(coordinate_j+1:coordinate_j+3)
      length2 = dot_product(bond, bond)
      radicand = 2.0d0*data%target_length(k)**2-length2
      if (radicand < -domain_tolerance*data%target_length(k)**2 .or. &
          .not. ieee_is_finite(radicand)) then
        call report_failure(k, failed_constraint)
        success = .false.
        return
      end if
      projected_length = sqrt(max(0.0d0, radicand))
      data%rhs_a(k) = data%scale(k)*(data%target_length(k)-projected_length)
    end do
    call solve_expansion
    call apply_solution(x_new)
  end do

  maximum_error = 0.0d0
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = x_new(coordinate_i+1:coordinate_i+3) - &
           x_new(coordinate_j+1:coordinate_j+3)
    relative_error = abs(sqrt(dot_product(bond,bond))/data%target_length(k)-1.0d0)
    maximum_error = max(maximum_error, relative_error)
  end do
  if (present(max_relative_error)) max_relative_error = maximum_error
end function lincs_positions


subroutine update_coupling_matrix
  integer :: row, position, column

  do row = 1, data%constraint_count
    do position = data%coupled_start(row), data%coupled_start(row+1)-1
      column = data%coupled_index(position)
      data%coupling(position) = data%coupled_prefactor(position)* &
        dot_product(data%direction(:,row), data%direction(:,column))
    end do
  end do
end subroutine update_coupling_matrix


subroutine solve_expansion
  integer :: order, row, position, column

  data%solution = data%rhs_a
  do order = 1, data%expansion_order
    data%rhs_b = 0.0d0
    do row = 1, data%constraint_count
      do position = data%coupled_start(row), data%coupled_start(row+1)-1
        column = data%coupled_index(position)
        data%rhs_b(row) = data%rhs_b(row) + &
                         data%coupling(position)*data%rhs_a(column)
      end do
    end do
    data%solution = data%solution + data%rhs_b
    data%rhs_a = data%rhs_b
  end do
end subroutine solve_expansion


subroutine apply_solution(coordinates)
  real(8), intent(inout) :: coordinates(:)
  integer :: k, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: correction(3)

  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    correction = data%direction(:,k)*data%scale(k)*data%solution(k)
    coordinates(coordinate_i+1:coordinate_i+3) = &
      coordinates(coordinate_i+1:coordinate_i+3) - &
      data%inverse_mass(atom_i)*correction
    coordinates(coordinate_j+1:coordinate_j+3) = &
      coordinates(coordinate_j+1:coordinate_j+3) + &
      data%inverse_mass(atom_j)*correction
  end do
end subroutine apply_solution


subroutine report_failure(constraint, failed_constraint)
  integer, intent(in) :: constraint
  integer, intent(out), optional :: failed_constraint

  if (present(failed_constraint)) failed_constraint = constraint
end subroutine report_failure


subroutine reset_lincs
  if (allocated(data%atom_i)) deallocate(data%atom_i)
  if (allocated(data%atom_j)) deallocate(data%atom_j)
  if (allocated(data%target_length)) deallocate(data%target_length)
  if (allocated(data%inverse_mass)) deallocate(data%inverse_mass)
  if (allocated(data%scale)) deallocate(data%scale)
  if (allocated(data%coupled_start)) deallocate(data%coupled_start)
  if (allocated(data%coupled_index)) deallocate(data%coupled_index)
  if (allocated(data%coupled_prefactor)) deallocate(data%coupled_prefactor)
  if (allocated(data%coupling)) deallocate(data%coupling)
  if (allocated(data%direction)) deallocate(data%direction)
  if (allocated(data%rhs_a)) deallocate(data%rhs_a)
  if (allocated(data%rhs_b)) deallocate(data%rhs_b)
  if (allocated(data%solution)) deallocate(data%solution)
  data%constraint_count = 0
  data%atom_count = 0
  data%expansion_order = 4
  data%rotation_iterations = 1
end subroutine reset_lincs

end module lincs
