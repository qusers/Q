program test_lincs_mpi
  use mpi
  use lincs
  use lincs_mpi
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(lincs_data_type) :: reference_data
  integer :: rank, ranks, ierr, n, mode, nc, na, k, j, b, failed, case_id, assertions
  integer :: configuration, expansion, minimum
  integer, allocatable :: ai(:), aj(:)
  real(8), allocatable :: old(:), candidate(:), reference(:), target(:), mass(:)
  real(8) :: error, difference, global_difference
  logical :: success
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,ranks,ierr)
  assertions = 0
  do configuration = 1, 3
    expansion = 8
    minimum = 1
    if (configuration == 1) minimum = 2
    if (configuration == 2) expansion = 12
    do case_id = 1, 5
      n = 257
      if (case_id == 1) n = 1
      if (case_id == 4) n = 4096
      if (case_id == 5) n = 8
      do mode = 1, 2
        na = 3*n
        nc = 2*n
        if (mode == 2) nc = nc+n-1
        if (case_id == 5 .and. mode == 2) nc = nc+1
        allocate(ai(nc),aj(nc),target(nc),mass(na),old(3*na),candidate(3*na),reference(3*na))
        old = 0.0d0
        k = 0
        do j = 1, n
          b = 3*j-2
          old(3*b-2:3*b) = [1.2d0*(j-1), 0.9d0*mod(j,2), 0.0d0]
          old(3*b+1:3*b+3) = old(3*b-2:3*b)+[0.0d0,0.8d0,0.6d0]
          old(3*b+4:3*b+6) = old(3*b-2:3*b)+[0.0d0,-0.8d0,0.6d0]
          mass(b:b+2) = [1.0d0/12.0d0,1.0d0,1.0d0]
          k = k+1; ai(k) = b; aj(k) = b+1; target(k) = 1.0d0
          k = k+1; ai(k) = b; aj(k) = b+2; target(k) = 1.0d0
          if (mode == 2 .and. j > 1) then
            k = k+1; ai(k) = b-3; aj(k) = b; target(k) = 1.5d0
          end if
        end do
        if (case_id == 5) then
          do j=1,n
            b=3*j-2
            old(3*b-2:3*b)=[cos(2.0d0*acos(-1.0d0)*j/n),sin(2.0d0*acos(-1.0d0)*j/n),0.0d0] &
                             *1.5d0/(2.0d0*sin(acos(-1.0d0)/n))
            old(3*b+1:3*b+3)=old(3*b-2:3*b)+[0.0d0,0.8d0,0.6d0]
            old(3*b+4:3*b+6)=old(3*b-2:3*b)+[0.0d0,-0.8d0,0.6d0]
          end do
          if (mode == 2) then
            ai(nc)=1;aj(nc)=3*n-2;target(nc)=1.5d0
          end if
        end if
        if (case_id == 3) then
          mass(1) = 0.0d0
          ai = ai(nc:1:-1); aj = aj(nc:1:-1); target = target(nc:1:-1)
        end if
        do k = 1, size(old)
          candidate(k) = old(k)+merge(0.075d0,0.003d0,case_id == 3)*sin(1.7d0*k)
        end do
        reference = candidate
        call setup_lincs(ai,aj,target,mass,expansion,minimum,1.0d-6,8,reference_data)
        call setup_lincs_mpi(MPI_COMM_WORLD,reference_data)
        success = lincs_positions(old,reference,failed,error,reference_data)
        call require(success,'serial reference succeeds')
        success = lincs_positions_mpi(old,candidate,failed,error)
        call require(success,'parallel solver succeeds')
        difference = 0.0d0
        if (rank == 0) difference = maxval(abs(candidate-reference))
        ! The wrapper only modifies rank zero for the one-rank fallback.
        call MPI_Allreduce(difference,global_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        call require(global_difference < 2.0d-12,'serial and MPI coordinates agree')
        call require(error <= 1.0d-6,'production accuracy retained')
        if (rank == 0) write(*,'(a,3i8,es14.5)') 'case/mode/constraints/difference ',case_id,mode,nc,global_difference
        if (case_id == 1) then
          do j=0,1
            call setup_lincs(ai,aj,target,mass,1,0,huge(1.0d0),j,reference_data)
            call setup_lincs_mpi(MPI_COMM_WORLD,reference_data)
            candidate=old+0.02d0
            candidate(2)=candidate(2)+0.02d0
            reference=candidate
            success=lincs_positions(old,reference,failed,error,reference_data)
            call require(success,'zero-minimum-rotation serial succeeds')
            success=lincs_positions_mpi(old,candidate,failed,error)
            call require(success,'zero-minimum-rotation MPI succeeds')
            if (rank == 0) call require(maxval(abs(candidate-reference)) < 2.0d-14, &
                                       'zero-minimum-rotation semantics match')
          end do
        end if
        if (ranks > 1) then
          candidate = old
          candidate(size(candidate)) = ieee_value(0.0d0,ieee_quiet_nan)
          success = lincs_positions_mpi(old,candidate,failed,error)
          call require(.not.success .and. failed > 0,'nonfinite geometry fails collectively')
        end if
        deallocate(ai,aj,target,mass,old,candidate,reference)
      end do
    end do
  end do
  ! A three-bond component split across owned-atom ranges must still match
  ! the serial solve after halo construction and ghost exchange.
  block
    real(8) :: origin(12), trial(12), expected(12), lengths(3), weights(4)
    weights=[1.0d0/12.0d0,1.0d0,1.0d0,1.0d0]
    lengths=1.0d0
    origin=[0.0d0,0.0d0,0.0d0, 1.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0]
    do configuration=1,3
      expansion=8;minimum=1
      if (configuration==1) minimum=2
      if (configuration==2) expansion=12
      call setup_lincs([1,1,1],[2,3,4],lengths,weights,expansion,minimum,1.0d-6,8,reference_data)
      call setup_lincs_mpi(MPI_COMM_WORLD,reference_data)
      do k=1,12
        trial(k)=origin(k)+0.03d0*sin(1.7d0*k)
      end do
      expected=trial
      success=lincs_positions(origin,expected,failed,error,reference_data)
      call require(success,'three-bond serial solve succeeds')
      success=lincs_positions_mpi(origin,trial,failed,error)
      call require(success .and. error<=1.0d-6,'three-bond MPI solve succeeds')
      if (rank==0) call require(maxval(abs(trial-expected))<2.0d-14,'three-bond MPI coordinates agree')
    end do
  end block
  ! The rank-zero fallback must report the same failure on every rank.
  block
    real(8) :: origin(6), trial(6)
    origin = 0.0d0
    trial = 0.0d0
    call setup_lincs([1],[2],[1.0d0],[1.0d0,1.0d0],context=reference_data)
    call setup_lincs_mpi(MPI_COMM_WORLD,reference_data,huge(1))
    success=lincs_positions_mpi(origin,trial,failed,error)
    call require(.not.success .and. failed==1 .and. error==-1.0d0, &
                 'inactive fallback failure is collective')
  end block
  ! Reinitializing an active MPI solver with no constraints clears its state.
  call setup_lincs([integer::],[integer::],[real(8)::],[1.0d0],context=reference_data)
  call setup_lincs_mpi(MPI_COMM_WORLD,reference_data)
  block
    real(8) :: point(3)
    point=0.0d0
    success=lincs_positions_mpi([0.0d0,0.0d0,0.0d0],point,failed,error)
    call require(success .and. error == 0.0d0,'empty topology is inactive')
  end block
  if (rank == 0) write(*,'(a,i0)') 'MPI LINCS assertions on rank zero: ',assertions
  call MPI_Finalize(ierr)
contains
  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    assertions = assertions+1
    if (.not.ok) then
      write(*,*) 'FAIL rank ',rank,': ',message
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_lincs_mpi
