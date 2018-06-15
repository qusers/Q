!------------------------------------------------------------------------------!
!  Q version 5.7                                                               !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          !
!  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer,     !
!  Miha Purg, Irek Szeler, Mauricio Esguerra                                   !
!  latest update: August 29, 2017                                              !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!  Copyright (c) 2017 Johan Aqvist, John Marelius, Shina Caroline Lynn Kamerlin
!!  and Paul Bauer
!!  program: **qdyn**  
!!  by Johan Aqvist, John Marelius, Anders Kaplan & Martin Nervall  
!!  qdyn molecular dynamics main program  
!------------------------------------------------------------------------------!
program qdyn
  use iso_fortran_env

  use md
  use mpiglob

#if defined (_DF_VERSION_)
  use dfport ! portability lib for signals. Used in windows.
#endif

  implicit none
  ! version data
  character(*), parameter :: program_name = 'qdyn'
  character(*), parameter :: program_version = '5.7'
  character(*), parameter :: program_date = '2015-02-22'
!  character(*), parameter :: options = compiler_options()
  character(len=32)       :: arg
  integer                 :: k
 
  
#if defined (USE_MPI)
  character(10) :: program_suffix = '_parallel'
#else
  character(10) :: program_suffix = ''
#endif

#if defined (USE_MPI)
  ! MPI error code
  integer :: qdyn_ierr
#endif

  ! signal handler data and declarations
  integer(4) :: sigret

#if defined (_DF_VERSION_)
  ! nothing
#else
  integer(4), parameter :: SIGINT = 2  ! CTRL-C signal
  integer(4), parameter :: SIGKILL = 9 ! kill/CTRL-BREAK signal
  integer(4), parameter :: SIGABRT = 6 ! kill/CTRL-BREAK signal
#endif

  external sigint_handler
  external sigkill_handler
  external sigabrt_handler

#if defined (USE_MPI)
  ! initialize MPI
  call MPI_Init(qdyn_ierr)
  if (qdyn_ierr .ne. MPI_SUCCESS) call die('failure at MPI init')
  call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, qdyn_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numnodes, qdyn_ierr)
#else
  nodeid = 0
  numnodes = 1
#endif

  ! init signal handlers
  sigret = SIGNAL(SIGINT, sigint_handler, -1_4)
  sigret = SIGNAL(SIGKILL, sigkill_handler, -1_4)
  sigret = SIGNAL(SIGABRT, sigabrt_handler, -1_4)

  !  call commandlineoptions

  ! initialize static data, display banner etc
  call startup

  if (nodeid .eq. 0) then
    ! master node: read input and initialize
    if(.not. initialize()) call die('Invalid data in input file') ! read input data
    call open_files ! open necessary files
    call topology ! read topology
    call prep_coord ! read coords, solvates etc
    if ( nstates > 0 ) call get_fep ! read fep/evb strategy
    !remove things with code 0 and maybe excluded bonded interactions

    call prep_sim ! prepare for simulation (calc. inv. mass, total charge,...)
    call close_input_files ! close input files

    call init_shake
    call make_nbqqlist
    call shrink_topology
    call nbmonitorlist
    call init_trj

    ! generate Maxwellian velocities and shake initial x and v if necessary
    if ( iseed > 0 ) then
      call maxwell
      call initial_shaking
      call stop_cm_translation
    end if
  end if

#if defined (USE_MPI)
  ! initialize slave nodes
  if (numnodes .gt. 1) call init_nodes
#endif
  
  ! count non-bonded pairs to get the maximum number, then distribute them 
  ! among the nodes
  call distribute_nonbonds

  ! do the work!
  call md_run

  if (nodeid .eq. 0) then
    ! master node: close output files
    call close_output_files
  end if

  ! deallocate memory etc.
  call shutdown

#if defined (USE_MPI)
  ! shut down MPI
  call MPI_Finalize(qdyn_ierr)
#endif



contains

  !-----------------------------------------------------------------------------
  ! startup subroutine
  !-----------------------------------------------------------------------------
  subroutine startup
    integer :: i
    integer :: datum(8)

    if (nodeid .eq. 0) then
      ! start-of-header
      ! write (*,'(80a)') ('=',i=1,80)

#if defined (DUM)
      print '(a)',  '--------------------------------------------------------------------------------'
      print '(a)',  'Welcome to Q'
      print '(2a)', 'This version was compiled using: ', compiler_version()
      print '(a)',  'And the following compiler options: '
!      write (output_unit, *, delim='quote') options
      write ( *, '( A /)' ) trim ( compiler_options() )
      write(*,'(a,a,a)') 'qdum input checker version ', trim(program_version), ' initializing'
      print '(a)',  '--------------------------------------------------------------------------------'
#elif defined(EVAL)
      print '(a)',  '--------------------------------------------------------------------------------'
      print '(a)',  'Welcome to Q'
      print '(2a)', 'This version was compiled using: ', compiler_version()
      print '(a)',  'And the following compiler options: '
!      write (output_unit, *, delim='quote') options
      write ( *, '( A /)' ) trim ( compiler_options() )
      write(*,'(a,a,a)') 'qdyn evaluation version ', trim(program_version), ' initializing'
      write(*,'(a)') 'This version is for evaluation purposes only.'
      write(*,'(a)') 'Optimizations are disabled - runs at <20% of maximum speed.'
      print '(a)',  '--------------------------------------------------------------------------------'
#else
      print '(a)',  '--------------------------------------------------------------------------------'
      print '(a)',  'Welcome to Q'
      print '(2a)', 'This version was compiled using: ', compiler_version()
      print '(a)',  'And the following compiler options: '
!      write (output_unit, *, delim='quote') options
      write ( *, '( A /)' ) trim ( compiler_options() )
      write(*,'(a,a,a,a)') 'qdyn version ', trim(program_version), trim(program_suffix),' initializing'
      print '(a)',  '--------------------------------------------------------------------------------'
#endif

      call date_and_time(values=datum)
      write(*,130) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
130   format('Current date ',i4,'-',i2,'-',i2,' and time ',i2,':',i2,':',i2)
    end if

    ! initialize used modules
    call md_startup

  end subroutine startup

  !----------------------------------------------------------------------------!
  !!  subroutine: startup  
  !!  Startup  
  !----------------------------------------------------------------------------!
  !subroutine startup
  !if (nodeid .eq. 0) then
  !  print '(a)',  '--------------------------------------------------------------------------------'
  !  print '(4a)', 'Welcome to ', program_name, ' version: ', program_version
  !  print '(a)',  ' '
  !  print '(2a)', 'This version was compiled using: ', compiler_version()
  !  print '(a)',  ' '
  !  print '(a)',  'And using the following compiler options: '
  !  write (output_unit, *, delim='quote') options
  !!  write ( *, '( A /)' ) trim ( compiler_options() )
  !!  print '(a)',  trim(compiler_options())
  !  print '(a)',  ' '
  !  print '(a)',  'For command line options type qdyn --help  or qprep -h at the terminal.'
  !  print '(a)',  ' '
  !  print '(a)',  'If you are using the interactive mode you can type "help"'
  !  print '(a)',  'at the prompt. To quit type "quit".'
  !  print '(a)',  '--------------------------------------------------------------------------------'
  !end if
  !  ! initialize modules
  !  call md_startup
  !  write(*,*)
  !end subroutine startup



  !-----------------------------------------------------------------------------
  ! shutdown subroutine
  !-----------------------------------------------------------------------------
  subroutine shutdown
    integer :: i

    if (nodeid .eq. 0) then
#if defined (DUM)
      write(*,*) 'qdum input checker version ', program_version, ' terminated normally.'
#else
      write(*,*) 'qdyn version ', trim(program_version), trim(program_suffix), ' terminated normally.'
#endif
      write (*,'(80a)') ('=',i=1,80)
    end if

    ! call md's shutdown
    call md_shutdown

  end subroutine shutdown

#if ! defined(__PATHCC__)
  integer(4) function signal( signum, proc, flag )
    implicit none
    integer(4) :: signum, flag
    external proc
    signal = 1
  end function signal
#endif



  !----------------------------------------------------------------------------!
  !!  subroutine: commandlineoptions  
  !----------------------------------------------------------------------------!
  subroutine commandlineoptions
    do k = 2, command_argument_count()
      call get_command_argument(k, arg)
      select case (arg)
        case ('-v', '--version')
          print '(3a)', program_name, ' version ', program_version
          stop
        case ('-h', '--help')
          call print_help()
          stop
        case default
          print '(a,a,/)', 'Unrecognized command-line option: ', arg
          call print_help()
          stop
      end select
    end do
  end subroutine commandlineoptions

  !----------------------------------------------------------------------------!
  !!  subroutine: print_help  
  !----------------------------------------------------------------------------!
  subroutine print_help()
    print '(a)', 'usage:'
    print '(a)', 'qdyn [OPTION]'
    print '(a)', '  or'
    print '(a)', 'qdyn  inputfile.inp > outputfile.out'
    print '(a)', ''
    print '(a)', 'Without options, qdyn goes into interactive mode.'
    print '(a)', ''
    print '(a)', 'qdyn [OPTION]:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
  end subroutine print_help


end program qdyn



!-------------------------------------------------------------------------------
! signal handlers
!-------------------------------------------------------------------------------
integer(4) function sigint_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num

  call die('user request (control-C)')
  sigint_handler = 1
end function sigint_handler

integer(4) function sigkill_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num

  call die('kill signal')
  sigkill_handler = 1
end function sigkill_handler

integer(4) function sigabrt_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num

  call die('kill signal')
  sigabrt_handler = 1
end function sigabrt_handler

