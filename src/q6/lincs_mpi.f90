module lincs_mpi
  ! Overlapping graph partitions: expansion_order hops suffice for each linear
  ! solve. Exchange corrected owned atoms before every nonlinear correction.
  ! No communication occurs inside the Neumann expansion, even for one large
  ! connected solute. All indices in diagnostics refer to the original topology.
  use mpi
  use lincs
  implicit none
  private
  public :: setup_lincs_mpi, lincs_positions_mpi, lincs_mpi_active, lincs_mpi_work
  type(lincs_data_type), save :: local_data
  integer, save :: comm, rank, ranks
  integer, allocatable, save :: active_atoms(:), local_active(:), owned_local(:), constraint_ids(:)
  integer, allocatable, save :: counts(:), offsets(:)
  real(8), allocatable, save :: gathered_coordinates(:,:), old_local(:), new_local(:), send_coords(:)
  integer, allocatable, save :: send_local(:), receive_local(:), send_counts(:), receive_counts(:)
  integer, allocatable, save :: send_offsets(:), receive_offsets(:)
  real(8), allocatable, save :: send_ghost(:), receive_ghost(:)
  integer, allocatable, save :: input_atoms(:), input_counts(:), input_offsets(:)
  real(8), allocatable, save :: input_global(:,:), input_local(:,:)
  logical, save :: exchange_ghosts = .false.
  logical, save :: enabled = .false.
contains

subroutine check_mpi(status)
  integer, intent(in) :: status
  integer :: ierr
  if (status /= MPI_SUCCESS) call MPI_Abort(comm, status, ierr)
end subroutine check_mpi

logical function lincs_mpi_active()
  lincs_mpi_active = enabled
end function lincs_mpi_active

subroutine lincs_mpi_work(owned_atoms, local_constraints)
  integer, intent(out) :: owned_atoms, local_constraints
  owned_atoms = 0
  local_constraints = 0
  if (.not. enabled) return
  owned_atoms = size(owned_local)
  local_constraints = local_data%constraint_count
end subroutine lincs_mpi_work

subroutine setup_lincs_mpi(communicator, root_data, minimum_constraints)
  integer, intent(in) :: communicator
  integer, intent(in), optional :: minimum_constraints
  type(lincs_data_type), intent(in), optional :: root_data
  type(lincs_data_type) :: topology
  integer :: options(5), ierr, atoms, constraints, minimum, nactive, nlocal
  integer :: owned_begin, owned_end, k, a, p, layer, j
  integer, allocatable :: ai(:), aj(:), degree(:), active_index(:), distance(:), map(:)
  real(8), allocatable :: target(:), mass(:)
  real(8) :: tolerance

  comm = communicator
  call MPI_Comm_rank(comm, rank, ierr)
  call check_mpi(ierr)
  call MPI_Comm_size(comm, ranks, ierr)
  call check_mpi(ierr)
  if (allocated(active_atoms)) then
    deallocate(input_atoms, input_counts, input_offsets, input_global, input_local)
    deallocate(send_local, receive_local, send_counts, receive_counts, send_offsets, receive_offsets, &
               send_ghost, receive_ghost)
    deallocate(active_atoms, local_active, owned_local, constraint_ids, counts, offsets, &
               gathered_coordinates, old_local, new_local, send_coords)
  end if
  enabled = .false.
  options = 0
  if (rank == 0) then
    if (.not. present(root_data)) error stop 'Root LINCS MPI setup requires topology'
    options = [root_data%atom_count, root_data%constraint_count, root_data%expansion_order, &
               root_data%rotation_iterations, root_data%maximum_rotation_iterations]
    tolerance = root_data%accuracy_tolerance
  end if
  call MPI_Bcast(options, 5, MPI_INTEGER, 0, comm, ierr)
  call check_mpi(ierr)
  atoms = options(1); constraints = options(2)
  minimum = 0
  if (present(minimum_constraints)) minimum = minimum_constraints
  if (ranks == 1 .or. constraints == 0 .or. constraints < minimum) then
    if (rank == 0) local_data = root_data
    return
  end if
  call MPI_Bcast(tolerance, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call check_mpi(ierr)
  allocate(ai(constraints), aj(constraints), target(constraints), mass(atoms))
  if (rank == 0) then
    ai = root_data%atom_i; aj = root_data%atom_j
    target = root_data%target_length; mass = root_data%inverse_mass
  end if
  call MPI_Bcast(ai, constraints, MPI_INTEGER, 0, comm, ierr)
  call check_mpi(ierr)
  call MPI_Bcast(aj, constraints, MPI_INTEGER, 0, comm, ierr)
  call check_mpi(ierr)
  call MPI_Bcast(target, constraints, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call check_mpi(ierr)
  call MPI_Bcast(mass, atoms, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call check_mpi(ierr)
  call setup_lincs(ai, aj, target, mass, options(3), options(4), tolerance, options(5), topology)

  allocate(degree(atoms), active_index(atoms), distance(constraints), map(atoms))
  degree = 0
  do k = 1, constraints
    degree(ai(k)) = degree(ai(k))+1
    degree(aj(k)) = degree(aj(k))+1
  end do
  nactive = count(degree > 0)
  allocate(active_atoms(nactive), counts(ranks), offsets(ranks))
  active_index = 0
  j = 0
  do a = 1, atoms
    if (degree(a) == 0) cycle
    j = j+1
    active_atoms(j) = a
    active_index(a) = j
  end do
  do p = 0, ranks-1
    offsets(p+1) = 3*((nactive*p)/ranks)
    counts(p+1) = 3*((nactive*(p+1))/ranks)-offsets(p+1)
  end do
  owned_begin = offsets(rank+1)/3+1
  owned_end = (offsets(rank+1)+counts(rank+1))/3
  distance = -1
  do k = 1, constraints
    a = active_index(ai(k)); j = active_index(aj(k))
    if ((a >= owned_begin .and. a <= owned_end) .or. &
        (j >= owned_begin .and. j <= owned_end)) distance(k) = 0
  end do
  do layer = 0, options(3)-1
    do k = 1, constraints
      if (distance(k) /= layer) cycle
      do p = topology%coupled_start(k), topology%coupled_start(k+1)-1
        j = topology%coupled_index(p)
        if (distance(j) < 0) distance(j) = layer+1
      end do
    end do
  end do
  allocate(constraint_ids(count(distance >= 0)))
  degree = 0
  j = 0
  do k = 1, constraints
    if (distance(k) < 0) cycle
    j = j+1
    constraint_ids(j) = k
    degree(ai(k)) = 1; degree(aj(k)) = 1
  end do
  nlocal = count(degree > 0)
  allocate(local_active(nlocal), owned_local(owned_end-owned_begin+1))
  map = 0
  j = 0
  do a = 1, atoms
    if (degree(a) == 0) cycle
    j = j+1
    map(a) = j
    local_active(j) = active_index(a)
  end do
  do a = owned_begin, owned_end
    owned_local(a-owned_begin+1) = map(active_atoms(a))
  end do
  call setup_lincs(map(ai(constraint_ids)), map(aj(constraint_ids)), target(constraint_ids), &
                   mass(active_atoms(local_active)), options(3), options(4), tolerance, options(5), local_data)
  allocate(gathered_coordinates(3,nactive), old_local(3*nlocal), new_local(3*nlocal), send_coords(counts(rank+1)))
  call setup_input_scatter()
  call setup_ghost_exchange(map)
  enabled = .true.
end subroutine setup_lincs_mpi

subroutine setup_input_scatter()
  integer :: nlocal, p, total, ierr
  nlocal = size(local_active)
  allocate(input_counts(ranks), input_offsets(ranks))
  input_counts = 0; input_offsets = 0
  call MPI_Gather(nlocal,1,MPI_INTEGER,input_counts,1,MPI_INTEGER,0,comm,ierr)
  call check_mpi(ierr)
  do p = 2, ranks
    input_offsets(p) = input_offsets(p-1)+input_counts(p-1)
  end do
  total = sum(input_counts)
  allocate(input_atoms(total),input_global(6,total),input_local(6,nlocal))
  call MPI_Gatherv(local_active,nlocal,MPI_INTEGER,input_atoms,input_counts,input_offsets,MPI_INTEGER,0,comm,ierr)
  call check_mpi(ierr)
  if (rank == 0) input_atoms = active_atoms(input_atoms)
  input_counts = 6*input_counts; input_offsets = 6*input_offsets
end subroutine setup_input_scatter

subroutine setup_ghost_exchange(map)
  integer, intent(in) :: map(:)
  integer, allocatable :: owner(:), requested(:), incoming(:), cursor(:)
  integer :: p, j, a, ierr, ghosts, all_ghosts
  allocate(owner(size(active_atoms)), cursor(ranks), send_counts(ranks), receive_counts(ranks), &
           send_offsets(ranks), receive_offsets(ranks))
  do p = 1, ranks
    owner(offsets(p)/3+1:(offsets(p)+counts(p))/3) = p
  end do
  receive_counts = 0
  do j = 1, size(local_active)
    p = owner(local_active(j))
    if (p /= rank+1) receive_counts(p) = receive_counts(p)+1
  end do
  call MPI_Alltoall(receive_counts,1,MPI_INTEGER,send_counts,1,MPI_INTEGER,comm,ierr)
  call check_mpi(ierr)
  send_offsets(1) = 0; receive_offsets(1) = 0
  do p = 2, ranks
    send_offsets(p) = send_offsets(p-1)+send_counts(p-1)
    receive_offsets(p) = receive_offsets(p-1)+receive_counts(p-1)
  end do
  ghosts = sum(receive_counts)
  allocate(requested(ghosts), receive_local(ghosts), incoming(sum(send_counts)), send_local(sum(send_counts)))
  cursor = receive_offsets
  do j = 1, size(local_active)
    a = local_active(j); p = owner(a)
    if (p == rank+1) cycle
    cursor(p) = cursor(p)+1
    requested(cursor(p)) = a
    receive_local(cursor(p)) = j
  end do
  call MPI_Alltoallv(requested,receive_counts,receive_offsets,MPI_INTEGER, &
                     incoming,send_counts,send_offsets,MPI_INTEGER,comm,ierr)
  call check_mpi(ierr)
  do j = 1, size(incoming)
    send_local(j) = map(active_atoms(incoming(j)))
  end do
  send_counts = 3*send_counts; receive_counts = 3*receive_counts
  send_offsets = 3*send_offsets; receive_offsets = 3*receive_offsets
  allocate(send_ghost(sum(send_counts)), receive_ghost(sum(receive_counts)))
  call MPI_Allreduce(ghosts,all_ghosts,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call check_mpi(ierr)
  exchange_ghosts = all_ghosts > 0
end subroutine setup_ghost_exchange

subroutine exchange_coordinates()
  integer :: j, a, ierr
  if (.not. exchange_ghosts) return
  do j = 1, size(send_local)
    a = send_local(j)
    send_ghost(3*j-2:3*j) = new_local(3*a-2:3*a)
  end do
  call MPI_Alltoallv(send_ghost,send_counts,send_offsets,MPI_DOUBLE_PRECISION, &
                     receive_ghost,receive_counts,receive_offsets,MPI_DOUBLE_PRECISION,comm,ierr)
  call check_mpi(ierr)
  do j = 1, size(receive_local)
    a = receive_local(j)
    new_local(3*a-2:3*a) = receive_ghost(3*j-2:3*j)
  end do
end subroutine exchange_coordinates

logical function lincs_positions_mpi(x_old, x_new, failed_constraint, max_relative_error) result(success)
  real(8), intent(in) :: x_old(:)
  real(8), intent(inout) :: x_new(:)
  integer, intent(out) :: failed_constraint
  real(8), intent(out) :: max_relative_error
  integer :: a, j, ierr, iteration, failed, global_failed, worst
  real(8) :: error
  logical :: finite

  success = .true.
  failed_constraint = 0
  max_relative_error = 0.0d0
  if (.not. enabled) then
    if (rank == 0) success = lincs_positions(x_old, x_new, failed_constraint, max_relative_error, local_data)
    call MPI_Bcast(success, 1, MPI_LOGICAL, 0, comm, ierr)
    call check_mpi(ierr)
    call MPI_Bcast(failed_constraint, 1, MPI_INTEGER, 0, comm, ierr)
    call check_mpi(ierr)
    call MPI_Bcast(max_relative_error, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call check_mpi(ierr)
    return
  end if
  if (rank == 0) then
    do j = 1, size(input_atoms)
      a = input_atoms(j)
      input_global(1:3,j) = x_old(3*a-2:3*a)
      input_global(4:6,j) = x_new(3*a-2:3*a)
    end do
  end if
  call MPI_Scatterv(input_global,input_counts,input_offsets,MPI_DOUBLE_PRECISION, &
                    input_local,size(input_local),MPI_DOUBLE_PRECISION,0,comm,ierr)
  call check_mpi(ierr)
  do j = 1, size(local_active)
    old_local(3*j-2:3*j) = input_local(1:3,j)
    new_local(3*j-2:3*j) = input_local(4:6,j)
  end do
  max_relative_error = -1.0d0
  do iteration = 0, local_data%maximum_rotation_iterations
    if (iteration == 0) then
      call lincs_linear_step(old_local, new_local, failed, local_data)
    else
      call lincs_rotation_step(new_local, failed, local_data)
    end if
    global_failed = huge(1)
    if (failed > 0) global_failed = constraint_ids(failed)
    call MPI_Allreduce(MPI_IN_PLACE, global_failed, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    call check_mpi(ierr)
    if (global_failed /= huge(1)) then
      failed_constraint = global_failed
      success = .false.
      return
    end if
    call exchange_coordinates()
    if (local_data%maximum_rotation_iterations > 0 .and. &
        iteration < max(1,local_data%rotation_iterations)) cycle
    call measure_constraint_error(new_local, error, worst, finite, local_data)
    call MPI_Allreduce(error, max_relative_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    call check_mpi(ierr)
    if (max_relative_error <= local_data%accuracy_tolerance) exit
  end do
  if (max_relative_error > local_data%accuracy_tolerance) then
    global_failed = huge(1)
    if (local_data%constraint_count > 0 .and. error == max_relative_error) global_failed = constraint_ids(worst)
    call MPI_Allreduce(global_failed, failed_constraint, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    call check_mpi(ierr)
    success = .false.
  end if
  do j = 1, size(owned_local)
    a = owned_local(j)
    send_coords(3*j-2:3*j) = new_local(3*a-2:3*a)
  end do
  call MPI_Gatherv(send_coords,size(send_coords),MPI_DOUBLE_PRECISION,gathered_coordinates, &
                   counts,offsets,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call check_mpi(ierr)
  ! Only rank zero owns the global integrator coordinates in Q's MPI model.
  if (rank /= 0) return
  do j = 1, size(active_atoms)
    a = active_atoms(j)
    x_new(3*a-2:3*a) = gathered_coordinates(:,j)
  end do
end function lincs_positions_mpi
end module lincs_mpi
