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
    integer :: maximum_rotation_iterations = 1
    real(8) :: accuracy_tolerance = 1.0d-4
    integer, allocatable :: atom_i(:)
    integer, allocatable :: atom_j(:)
    real(8), allocatable :: target_length(:)
    real(8), allocatable :: inverse_mass(:)
    real(8), allocatable :: scale(:)
    integer, allocatable :: coupled_start(:)
    integer, allocatable :: coupled_index(:)
    integer, allocatable :: complex_rows(:), pair_i(:), pair_j(:)
    real(8), allocatable :: pair_even(:), pair_odd(:)
    integer :: triple_count = 0
    integer, allocatable :: triple_rows(:,:)
    real(8), allocatable :: triple_polynomial(:,:,:)
    real(8), allocatable :: coupled_prefactor(:)
    real(8), allocatable :: coupling(:)
    real(8), allocatable :: direction(:,:)
    real(8), allocatable :: rhs_a(:)
    real(8), allocatable :: rhs_b(:)
    real(8), allocatable :: solution(:)
  end type lincs_data_type

  type(lincs_data_type), target :: default_data

  public :: copy_lincs_data, lincs_data_type, lincs_linear_step, lincs_rotation_step, measure_constraint_error
  public :: setup_lincs, initialize_lincs_positions, lincs_positions, lincs_is_active

contains

subroutine copy_lincs_data(context)
  type(lincs_data_type), intent(out) :: context
  context = default_data
end subroutine copy_lincs_data


subroutine setup_lincs(atom_i, atom_j, target_length, inverse_mass, &
                       expansion_order, rotation_iterations, accuracy_tolerance, &
                       maximum_rotation_iterations, context)
!!-------------------------------------------------------------------------------
!!  Cache the constraint topology and build its sparse coupling graph in O(K+C)
!!  storage, where K is the number of constraints and C the number of directly
!!  coupled constraint pairs.
!!-------------------------------------------------------------------------------
  type(lincs_data_type), intent(inout), optional, target :: context
  type(lincs_data_type), pointer :: data
  integer, intent(in) :: atom_i(:), atom_j(:)
  real(8), intent(in) :: target_length(:), inverse_mass(:)
  integer, intent(in), optional :: expansion_order, rotation_iterations
  integer, intent(in), optional :: maximum_rotation_iterations
  real(8), intent(in), optional :: accuracy_tolerance

  integer :: constraints, atoms, k, m, atom, other_atom, total_couplings
  integer :: position, pairs, complex_count
  integer, allocatable :: degree(:), atom_start(:), atom_cursor(:)
  integer, allocatable :: atom_constraints(:)
  real(8) :: inverse_mass_sum

  data => default_data
  if (present(context)) data => context

  data = lincs_data_type()

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
    data%maximum_rotation_iterations = rotation_iterations
  end if
  if (present(accuracy_tolerance)) then
    if (.not. ieee_is_finite(accuracy_tolerance) .or. accuracy_tolerance <= 0.0d0) then
      error stop 'LINCS accuracy tolerance must be finite and positive'
    end if
    data%accuracy_tolerance = accuracy_tolerance
  end if
  if (present(maximum_rotation_iterations)) then
    if (maximum_rotation_iterations < data%rotation_iterations) then
      error stop 'LINCS maximum rotations must be at least the requested rotations'
    end if
    data%maximum_rotation_iterations = maximum_rotation_iterations
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
    if (.not. ieee_is_finite(target_length(k)) .or. target_length(k) <= 0.0d0) then
      error stop 'LINCS target lengths must be finite and positive'
    end if
    if (.not. ieee_is_finite(inverse_mass(atom_i(k))) .or. &
        .not. ieee_is_finite(inverse_mass(atom_j(k))) .or. &
        inverse_mass(atom_i(k)) < 0.0d0 .or. inverse_mass(atom_j(k)) < 0.0d0) then
      error stop 'LINCS inverse masses must be finite and nonnegative'
    end if
    inverse_mass_sum = inverse_mass(atom_i(k)) + inverse_mass(atom_j(k))
    if (.not. ieee_is_finite(inverse_mass_sum) .or. inverse_mass_sum <= 0.0d0) then
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

  ! Isolated two-bond components are common in hydrogen-only constraints.
  ! Their Neumann polynomial reduces to even/odd scalar powers of A^2.
  ! Retain the configured finite expansion, including odd and zero orders.
  pairs = 0
  complex_count = 0
  do k = 1, constraints
    position = data%coupled_start(k)
    if (data%coupled_start(k+1)-position == 0) cycle
    if (data%coupled_start(k+1)-position == 1) then
      m = data%coupled_index(position)
      if (data%coupled_start(m+1)-data%coupled_start(m) == 1) then
        if (k < m) pairs = pairs+1
        cycle
      end if
    end if
    complex_count = complex_count+1
  end do
  allocate(data%complex_rows(complex_count), data%pair_i(pairs), data%pair_j(pairs), &
           data%pair_even(pairs), data%pair_odd(pairs))
  pairs = 0
  complex_count = 0
  do k = 1, constraints
    position = data%coupled_start(k)
    if (data%coupled_start(k+1)-position == 0) cycle
    if (data%coupled_start(k+1)-position == 1) then
      m = data%coupled_index(position)
      if (data%coupled_start(m+1)-data%coupled_start(m) == 1) then
        if (k < m) then
          pairs = pairs+1
          data%pair_i(pairs) = k
          data%pair_j(pairs) = m
        end if
        cycle
      end if
    end if
    complex_count = complex_count+1
    data%complex_rows(complex_count) = k
  end do

  call find_three_bond_components(data)

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


subroutine find_three_bond_components(data)
  type(lincs_data_type), intent(inout) :: data
  integer, allocatable :: groups(:,:), membership(:)
  integer :: k, row, position, neighbor, queue(4), found, next, groups_found
  groups_found = 0
  do k = 1, data%constraint_count
    if (data%coupled_start(k+1)-data%coupled_start(k) > 2) cycle
    queue(1) = k; found = 1; next = 1
    ! Bounded search: stop as soon as a fourth constraint is encountered.
    do while (next <= found .and. found < 4)
      row = queue(next)
      do position = data%coupled_start(row), data%coupled_start(row+1)-1
        neighbor = data%coupled_index(position)
        if (any(queue(1:found) == neighbor)) cycle
        found = found+1
        queue(found) = neighbor
        if (found == 4) exit
      end do
      next = next+1
    end do
    if (found /= 3) cycle
    if (minval(queue(1:3)) /= k) cycle
    if (groups_found == 0) then
      allocate(groups(3,data%constraint_count/3),membership(data%constraint_count))
      membership = 0
    end if
    groups_found = groups_found+1
    groups(:,groups_found) = queue(1:3)
    membership(queue(1:3)) = groups_found
  end do
  ! Leave allocations unchanged for graphs that cannot use this shortcut.
  ! Besides saving storage, this avoids perturbing subsequent MD allocations.
  data%triple_count = groups_found
  if (groups_found == 0) return
  data%triple_rows = groups(:,1:groups_found)
  allocate(data%triple_polynomial(3,3,groups_found))
  data%complex_rows = pack(data%complex_rows,membership(data%complex_rows) == 0)
end subroutine find_three_bond_components

subroutine update_three_bond_polynomials(data)
  type(lincs_data_type), intent(inout) :: data
  integer :: group, i, j, position, order, rows(3)
  real(8) :: a(3,3), square(3,3), power(3), coefficients(3), second, determinant
  do group = 1, data%triple_count
    rows = data%triple_rows(:,group)
    a = 0.0d0
    do i = 1, 3
      do position = data%coupled_start(rows(i)),data%coupled_start(rows(i)+1)-1
        do j = 1, 3
          if (data%coupled_index(position) == rows(j)) a(i,j) = data%coupling(position)
        end do
      end do
    end do
    ! For any zero-diagonal 3x3 matrix, Cayley-Hamilton gives
    ! A^3 = second*A + det(A)*I. Reduce the finite Neumann polynomial
    ! to c0*I+c1*A+c2*A^2.
    second = a(1,2)*a(2,1)+a(1,3)*a(3,1)+a(2,3)*a(3,2)
    determinant = a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(3,2)*a(2,1)
    power = [1.0d0,0.0d0,0.0d0]
    coefficients = power
    do order = 1, data%expansion_order
      power = [determinant*power(3),power(1)+second*power(3),power(2)]
      coefficients = coefficients+power
    end do
    square = matmul(a,a)
    data%triple_polynomial(:,:,group) = coefficients(2)*a+coefficients(3)*square
    do i = 1, 3
      data%triple_polynomial(i,i,group) = data%triple_polynomial(i,i,group)+coefficients(1)
    end do
  end do
end subroutine update_three_bond_polynomials


logical function lincs_is_active(context)
  type(lincs_data_type), intent(in), optional :: context

  if (present(context)) then
    lincs_is_active = context%constraint_count > 0
  else
    lincs_is_active = default_data%constraint_count > 0
  end if
end function lincs_is_active


logical function initialize_lincs_positions(coordinates, failed_constraint, &
                                            max_relative_error, context) result(success)
!!-------------------------------------------------------------------------------
!!  Project an arbitrary starting structure onto the LINCS constraint manifold.
!!  LINCS updates assume the reference coordinates are already constrained, so
!!  this one-time tightly converged SHAKE projection establishes that invariant.
!!-------------------------------------------------------------------------------
  type(lincs_data_type), intent(inout), optional, target :: context
  type(lincs_data_type), pointer :: data
  real(8), intent(inout) :: coordinates(:)
  integer, intent(out), optional :: failed_constraint
  real(8), intent(out), optional :: max_relative_error

  integer, parameter :: maximum_iterations = 1000
  real(8), parameter :: convergence_tolerance = 1.0d-12
  integer :: iteration, k, atom_i, atom_j, coordinate_i, coordinate_j
  integer :: worst_constraint
  real(8) :: reference(3), current(3), reference_dot_current, length2
  real(8) :: difference, correction, maximum_error, inverse_mass_sum
  logical :: finite_geometry

  data => default_data
  if (present(context)) data => context

  success = .true.
  if (present(failed_constraint)) failed_constraint = 0
  if (present(max_relative_error)) max_relative_error = 0.0d0
  if (data%constraint_count == 0) return

  ! Keep the initial bond directions fixed, as in SHAKE and finite-step LINCS.
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    reference = coordinates(coordinate_i+1:coordinate_i+3) - &
                coordinates(coordinate_j+1:coordinate_j+3)
    length2 = dot_product(reference,reference)
    if (.not. ieee_is_finite(length2) .or. length2 <= tiny(1.0d0)) then
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
      inverse_mass_sum = data%inverse_mass(atom_i)+data%inverse_mass(atom_j)
      if (.not. ieee_is_finite(reference_dot_current) .or. &
          reference_dot_current <= tiny(1.0d0) .or. &
          .not. ieee_is_finite(inverse_mass_sum) .or. inverse_mass_sum <= 0.0d0) then
        call report_failure(k,failed_constraint)
        success = .false.
        return
      end if
      length2 = dot_product(current,current)
      difference = data%target_length(k)**2-length2
      correction = difference/(2.0d0*reference_dot_current*inverse_mass_sum)
      if (.not. ieee_is_finite(length2) .or. .not. ieee_is_finite(difference) .or. &
          .not. ieee_is_finite(correction)) then
        call report_failure(k,failed_constraint)
        success = .false.
        return
      end if
      coordinates(coordinate_i+1:coordinate_i+3) = &
        coordinates(coordinate_i+1:coordinate_i+3) + &
        data%inverse_mass(atom_i)*reference*correction
      coordinates(coordinate_j+1:coordinate_j+3) = &
        coordinates(coordinate_j+1:coordinate_j+3) - &
        data%inverse_mass(atom_j)*reference*correction
    end do

    call measure_constraint_error(coordinates, maximum_error, worst_constraint, &
                                  finite_geometry, data)
    if (.not. finite_geometry) then
      call report_failure(worst_constraint,failed_constraint)
      success = .false.
      if (present(max_relative_error)) max_relative_error = maximum_error
      return
    end if
    if (maximum_error <= convergence_tolerance) then
      if (present(max_relative_error)) max_relative_error = maximum_error
      return
    end if
  end do

  success = .false.
  call report_failure(worst_constraint,failed_constraint)
  if (present(max_relative_error)) max_relative_error = maximum_error
end function initialize_lincs_positions


logical function lincs_positions(x_old, x_new, failed_constraint, &
                                 max_relative_error, context) result(success)
!!-------------------------------------------------------------------------------
!!  Constrain one unconstrained position update.  The requested number of
!!  nonlinear corrections is a minimum; difficult steps receive bounded
!!  iterative refinement until the configured accuracy is reached.
!!-------------------------------------------------------------------------------
  type(lincs_data_type), intent(inout), optional, target :: context
  type(lincs_data_type), pointer :: data
  real(8), intent(in) :: x_old(:)
  real(8), intent(inout) :: x_new(:)
  integer, intent(out), optional :: failed_constraint
  real(8), intent(out), optional :: max_relative_error

  integer :: iteration
  integer :: worst_constraint
  real(8) :: maximum_error
  logical :: finite_geometry

  data => default_data
  if (present(context)) data => context

  success = .true.
  if (present(failed_constraint)) failed_constraint = 0
  if (present(max_relative_error)) max_relative_error = 0.0d0
  if (data%constraint_count == 0) return
  if (present(max_relative_error)) max_relative_error = -1.0d0

  call lincs_linear_step(x_old, x_new, worst_constraint, data)
  if (worst_constraint /= 0) then
    call report_failure(worst_constraint, failed_constraint)
    success = .false.
    return
  end if

  ! Nonlinear correction for rotational lengthening, equations 17-18.
  maximum_error = huge(1.0d0)
  worst_constraint = 1
  do iteration = 1, data%maximum_rotation_iterations
    call lincs_rotation_step(x_new, worst_constraint, data)
    if (worst_constraint /= 0) then
      call report_failure(worst_constraint, failed_constraint)
      success = .false.
      return
    end if

    if (iteration >= data%rotation_iterations) then
      call measure_constraint_error(x_new, maximum_error, worst_constraint, &
                                    finite_geometry, data)
      if (.not. finite_geometry) then
        call report_failure(worst_constraint, failed_constraint)
        success = .false.
        return
      end if
      if (maximum_error <= data%accuracy_tolerance) exit
    end if
  end do

  ! Also handle the explicitly supported zero-rotation configuration.
  if (data%maximum_rotation_iterations == 0) then
    call measure_constraint_error(x_new, maximum_error, worst_constraint, &
                                  finite_geometry, data)
  end if
  if (.not. finite_geometry .or. maximum_error > data%accuracy_tolerance) then
    call report_failure(worst_constraint, failed_constraint)
    success = .false.
  end if
  if (present(max_relative_error)) max_relative_error = maximum_error
end function lincs_positions


subroutine lincs_linear_step(x_old, x_new, failed_constraint, data)
  type(lincs_data_type), intent(inout) :: data
  real(8), intent(in) :: x_old(:)
  real(8), intent(inout) :: x_new(:)
  integer, intent(out) :: failed_constraint
  integer :: k, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: bond(3), length2, projected_length
  failed_constraint = 0
  if (data%constraint_count == 0) return
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = x_old(coordinate_i+1:coordinate_i+3) - &
           x_old(coordinate_j+1:coordinate_j+3)
    length2 = dot_product(bond, bond)
    if (.not. ieee_is_finite(length2) .or. length2 <= tiny(1.0d0)) then
      failed_constraint = k
      return
    end if
    data%direction(:,k) = bond/sqrt(length2)
  end do
  call update_coupling_matrix(data)
  if (.not. all(ieee_is_finite(data%coupling))) then
    failed_constraint = 1
    return
  end if

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
    if (.not. ieee_is_finite(data%rhs_a(k))) then
      failed_constraint = k
      return
    end if
  end do
  call solve_expansion(data)
  if (.not. all(ieee_is_finite(data%solution))) then
    failed_constraint = 1
    return
  end if
  call apply_solution(x_new, data)

end subroutine lincs_linear_step

subroutine lincs_rotation_step(x_new, failed_constraint, data)
  type(lincs_data_type), intent(inout) :: data
  real(8), intent(inout) :: x_new(:)
  integer, intent(out) :: failed_constraint
  integer :: k, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: bond(3), length2, projected_length, radicand
  real(8), parameter :: domain_tolerance = 1.0d-12
  failed_constraint = 0
  if (data%constraint_count == 0) return
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = x_new(coordinate_i+1:coordinate_i+3) - &
           x_new(coordinate_j+1:coordinate_j+3)
    length2 = dot_product(bond, bond)
    radicand = 2.0d0*data%target_length(k)**2-length2
    if (.not. ieee_is_finite(length2) .or. .not. ieee_is_finite(radicand) .or. &
        radicand < -domain_tolerance*data%target_length(k)**2) then
      failed_constraint = k
      return
    end if
    projected_length = sqrt(max(0.0d0, radicand))
    data%rhs_a(k) = data%scale(k)*(data%target_length(k)-projected_length)
  end do
  call solve_expansion(data)
  if (.not. all(ieee_is_finite(data%solution))) then
    failed_constraint = 1
    return
  end if
  call apply_solution(x_new, data)

end subroutine lincs_rotation_step

subroutine measure_constraint_error(coordinates, maximum_error, worst_constraint, &
                                    finite_geometry, context)
  type(lincs_data_type), intent(in), optional, target :: context
  type(lincs_data_type), pointer :: data
  real(8), intent(in) :: coordinates(:)
  real(8), intent(out) :: maximum_error
  integer, intent(out) :: worst_constraint
  logical, intent(out) :: finite_geometry

  integer :: k, atom_i, atom_j, coordinate_i, coordinate_j
  real(8) :: bond(3), length2, relative_error

  data => default_data
  if (present(context)) data => context

  maximum_error = 0.0d0
  worst_constraint = 1
  finite_geometry = .true.
  do k = 1, data%constraint_count
    atom_i = data%atom_i(k)
    atom_j = data%atom_j(k)
    coordinate_i = 3*atom_i-3
    coordinate_j = 3*atom_j-3
    bond = coordinates(coordinate_i+1:coordinate_i+3) - &
           coordinates(coordinate_j+1:coordinate_j+3)
    length2 = dot_product(bond,bond)
    if (.not. ieee_is_finite(length2) .or. length2 < 0.0d0) then
      finite_geometry = .false.
      worst_constraint = k
      maximum_error = huge(1.0d0)
      return
    end if
    relative_error = abs(sqrt(length2)/data%target_length(k)-1.0d0)
    if (relative_error > maximum_error) then
      maximum_error = relative_error
      worst_constraint = k
    end if
  end do
end subroutine measure_constraint_error


subroutine update_coupling_matrix(data)
  type(lincs_data_type), intent(inout) :: data
  integer :: row, position, column, pair, order
  real(8) :: power, square

  do row = 1, data%constraint_count
    do position = data%coupled_start(row), data%coupled_start(row+1)-1
      column = data%coupled_index(position)
      data%coupling(position) = data%coupled_prefactor(position)* &
        dot_product(data%direction(:,row), data%direction(:,column))
    end do
  end do
  if (data%triple_count > 0) call update_three_bond_polynomials(data)
  do pair = 1, size(data%pair_i)
    row = data%pair_i(pair)
    column = data%pair_j(pair)
    square = data%coupling(data%coupled_start(row))*data%coupling(data%coupled_start(column))
    power = 1.0d0
    data%pair_even(pair) = 1.0d0
    data%pair_odd(pair) = 0.0d0
    do order = 1, data%expansion_order, 2
      data%pair_odd(pair) = data%pair_odd(pair)+power
      if (order < data%expansion_order) then
        power = power*square
        data%pair_even(pair) = data%pair_even(pair)+power
      end if
    end do
  end do
end subroutine update_coupling_matrix


subroutine solve_expansion(data)
  type(lincs_data_type), intent(inout) :: data
  integer :: order, pair, i, j, k

  data%solution = data%rhs_a
  do pair = 1, size(data%pair_i)
    i = data%pair_i(pair); j = data%pair_j(pair)
    data%solution(i) = data%pair_even(pair)*data%rhs_a(i) + &
      data%pair_odd(pair)*data%coupling(data%coupled_start(i))*data%rhs_a(j)
    data%solution(j) = data%pair_even(pair)*data%rhs_a(j) + &
      data%pair_odd(pair)*data%coupling(data%coupled_start(j))*data%rhs_a(i)
  end do
  do pair = 1, data%triple_count
    i = data%triple_rows(1,pair); j = data%triple_rows(2,pair); k = data%triple_rows(3,pair)
    ! Scalar expressions avoid allocating MATMUL temporaries for each group.
    data%solution(i) = data%triple_polynomial(1,1,pair)*data%rhs_a(i) + &
      data%triple_polynomial(1,2,pair)*data%rhs_a(j) + data%triple_polynomial(1,3,pair)*data%rhs_a(k)
    data%solution(j) = data%triple_polynomial(2,1,pair)*data%rhs_a(i) + &
      data%triple_polynomial(2,2,pair)*data%rhs_a(j) + data%triple_polynomial(2,3,pair)*data%rhs_a(k)
    data%solution(k) = data%triple_polynomial(3,1,pair)*data%rhs_a(i) + &
      data%triple_polynomial(3,2,pair)*data%rhs_a(j) + data%triple_polynomial(3,3,pair)*data%rhs_a(k)
  end do
  if (size(data%complex_rows) == 0) return
  do order = 1, data%expansion_order
    if (mod(order,2) == 1) then
      call expansion_product(data%rhs_a, data%rhs_b)
    else
      call expansion_product(data%rhs_b, data%rhs_a)
    end if
  end do
contains
  subroutine expansion_product(input, output)
    real(8), intent(in) :: input(:)
    real(8), intent(out) :: output(:)
    real(8) :: value
    integer :: row, position, entry
    do entry = 1, size(data%complex_rows)
      row = entry
      if (size(data%complex_rows) /= data%constraint_count) row = data%complex_rows(entry)
      value = 0.0d0
      do position = data%coupled_start(row), data%coupled_start(row+1)-1
        value = value + data%coupling(position)*input(data%coupled_index(position))
      end do
      output(row) = value
      data%solution(row) = data%solution(row) + value
    end do
  end subroutine expansion_product
end subroutine solve_expansion


subroutine apply_solution(coordinates, data)
  type(lincs_data_type), intent(in) :: data
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


end module lincs
