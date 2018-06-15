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
!>  program: **qprep**
!!  by Johan Aqvist & John Marelius  
!!  qprep topology preparation main program  
!------------------------------------------------------------------------------!
program qprep
  use iso_fortran_env

  use prep
  use avetr

  implicit none
  character(*), parameter :: program_name = 'qprep'
  character(*), parameter :: program_version = '5.7'
  character(*), parameter :: program_date = '2015-04-01'
  character(*), parameter :: program_suffix  = ''
!  character(*), parameter :: options = compiler_options()
  character(len=200)      :: fileName = ''
  character(len=32)       :: arg  
  logical                 :: use_inputfile
  integer                 :: i

  call commandlineoptions

  call startup

! allocate memory for libraries
  call topo_set_max(0, max_lib, 0)

! reset residue count
  call clearlib

  use_inputfile = check_inputfile(fileName)
  if (use_inputfile) then
    call qprep_from_inputfile(fileName)
  else
    call qprep_from_commandline
  endif

  call shutdown

contains

  !----------------------------------------------------------------------------!
  !>  subroutine: **qprep_from_inputfile**
  !!  Read input from file and execute commands  
  !----------------------------------------------------------------------------!
  subroutine qprep_from_inputfile(filename)
  character(200)  :: filename
  character(200)  :: command
  logical         :: readable
  integer         :: INPF_U = 0
  integer         :: stat

!Fix so that file is read per line and parsed in the same way as in the command line
  if (.not. parse_open_file(filename)) then
    write( * , '(a,a)') 'Could not open ',filename
    return
  endif

  do
    call get_string_arg(command)
    select case (command)
      case('quit', 'q') 
        exit
      case default
        call parse_command(command)
    end select
  enddo

  end subroutine qprep_from_inputfile


  !----------------------------------------------------------------------------!
  !>  subroutine: qprep_from_commandline
  !!  Read input from command line and execute commands  
  !----------------------------------------------------------------------------!
  subroutine qprep_from_commandline
  character(200) :: command

  do
    call parse_reset
    call get_string_arg(command, 'qprep> ')
    call locase(command) !make all lower case
    select case (command)
    case('quit', 'q')
      exit
    case default
      call parse_command(command)
    end select
  enddo

  end subroutine qprep_from_commandline


  !----------------------------------------------------------------------------!
  !>  subroutine: parse_command(command)
  !!  Parse a command and call corresponding subroutine  
  !----------------------------------------------------------------------------!
  subroutine parse_command(command)
  character(*), intent(IN) :: command
  ! --- Command loop
  select case (command)
    case ('average', 'av')
          call avetr_calc
    case ('readlib', 'rl')
          call readlib
    case('clearlib', 'cl')
          call clearlib
    case('readpdb', 'rp')
          call readpdb
    case('readprm', 'readff', 'rff', 'rprm')
          call readff
    case('addbond', 'ab')
          call addbond
    case('clearbond', 'clearbonds')
          call clearbond
    case('maketop', 'mt')
          call maketop
    case('cleartop')
          call cleartop
    case('listseq', 'ls')
          call listseq
    case('listres', 'lr')
          call listres
    case('writetop', 'wt')
          call writetop
    case('checkbonds', 'cb')
          call checkbonds
    case('checkangs', 'ca')
          call checkangs
    case('checktors', 'ct')
          call checktors
    case('checkimps', 'ci')
          call checkimps
    case('changeimp')
          call changeimp
    case('readtop', 'rt')
          call readtop
    case('readx', 'rx')
          call readx
    case('makeshell', 'ms')
          call make_shell2
    case('mask', 'ma')
          call modify_mask
    case('trajectory', 'trj', 'tr')
          call trajectory
    case('readtrajectory', 'readframe','rf')
          call readframe
    case('readnext', 'rn')
          call readnext
    case('solvate', 'so')
          call solvate
    case('writepdb', 'wp')
          call writepdb 
    case('writemol2', 'wm')
          call writemol2
    case('xlink', 'crosslink', 'xl')
          call xlink
    case('prefs', 'preferences')
          call listprefs
    case('set', 's')
          call set
    case('help', '?', 'h')
          call help
    case('boundary', 'bc')
          call define_boundary_condition
    case default
          write( * , '(/,a,a,a,/)') 'unrecognized command "', trim(command), '" (type ? for help)'
  end select
  end subroutine parse_command


  !----------------------------------------------------------------------------!
  !>  subroutine: help
  !!  Give help on commands  
  !----------------------------------------------------------------------------!
  subroutine help
  write( *, * )
  write( * , '(a)') &
  'command      argument            description', &
  '----------- ------------------- -------------------------------------------------', &
  'addbond                          adds extra bonds(e.g. S-S)',&
  'average                          computes an average structure from a trajectory file',&
  'boundary    [boundary condition]  set boundary condition',&
  '            [center] ',&
  '            (box) [boxlengths]',&
  '            (sphere) [radius]',&
  '            (sphere) [inner radius]',&
  'changeimp                        redefine(specified) improper torsions',&
  'checkangs   [energy_threshold]   check angle energies',&
  'checkbonds  [energy_threshold]   check bonds energies',&
  'checkimps   [energy_threshold]   check improper torsion energies',&
  'checktors   [energy_threshold]   check torsion energies',&
  'clearbond                        clears extra bonds',&
  'clearlib                         unloads all libaries',&
  'cleartop                         clears topology & parameters',&
  'help                             shows command list',&
  'listseq                          lists the residue sequence',&
  'listres     [residue_number]     lists atoms in residue',&
  'makeshell                        fix the mask of the atoms in the restrained shell',&
  'maketop     [name]               generates the topology',&
  'mask        [mask_def|none]      add to or clear atom mask',&
  'preferences                      list preferences', &
  'quit                             quits the program',&
  'readframe   [frame]              reads coordinates from trajectory',&
  'readlib     [library_file]       reads library file',&
  'readnext                         reads next frame from trajectory', &
  'readpdb     [pdb_file]           reads pdb file',&
  'readprm     [param.file]         reads FF parameter file',&
  'readtop     [topology_file]      reads topology file',&
  'readx       [restart_file]       reads coord. file',&
  'set                              set preferences',&
  'solvate(boundry=sphere)          solvate sphere with specified options',&      
  '            [center]',&
  '            [radius]',&
  '            [grid|file|restart] ',&
  '            [solvent name] ',&
  '            [file name] ',&
  '',&
  'solvate(boundry=box)             solvate box with specified options',&
  '            [grid|file|restart]  ',&   
  '            [solvent name] ',&
  '            [file name] ',&
  '',&
  'trajectory  [trajectory_file]    open trajectory file',&
  'writetop    [topology_file]      writes topology file',&
  'writepdb    [pdb_file]           writes pdb file', &
  '            [hydrogen_flag [gap_flag [water_flag]]]',&
  '',&                                
  'writemol2   [mol_no]             writes molecule mol_no (0 for all) to',&
  '            [mol2_file]          SYBYL mol2 file', &
  '            [hydrogen_flag [water_flag]]',&
  '',&                                
  'xlink                            add crosslink bonds interactively', &
  '', &
  'short form meaning' ,&
  '---------- -------', &
  'ab          addbond', &
  'av          average', &
  'bc          boundary', &
  'cb,ca,ct,ci checkbonds, checkangs, checktors, checkimps', &
  'cl          clearlib', &
  'h, ?        help', &
  'lr          listres', &
  'ls          listseq', &
  'ma          mask', &
  'ma 0        clear mask', &
  'ms          makeshell', &
  'mt          maketop', &
  'prefs       preferences', &
  'q           quit', &
  'rf          readframe', &
  'rff         readprm', &
  'rl          readlib', &
  'rn          readnext', &
  'rp          readpdb', &
  'rt          readtop', &
  'rx          readx', &
  's           set',&
  'so          solvate', &
  'trj, tr     trajectory', &
  'wp          writepdb', &
  'wt          writetop', &
  'xl          xlink'
  write( *, * )

  end subroutine help


  !----------------------------------------------------------------------------!
  !>  subroutine: startup
  !!  Startup  
  !----------------------------------------------------------------------------!
subroutine startup
  print '(a)',  '--------------------------------------------------------------------------------'
  print '(4a)', 'Welcome to ', program_name, ' version: ', program_version
  print '(a)',  ' '
  print '(2a)', 'This version was compiled using: ', compiler_version()
  print '(a)',  ' '
  print '(a)',  'And using the following compiler options: '
!  write (output_unit, *, delim='quote') options
  write ( *, '( A /)' ) trim ( compiler_options() )
!  print '(a)',  trim(compiler_options())
  print '(a)',  ' '
  print '(a)',  'For command line options type qprep --help  or qprep -h at the terminal.'
  print '(a)',  ' '
  print '(a)',  'If you are using the interactive mode you can type "help"'
  print '(a)',  'at the prompt. To quit type "quit".'
  print '(a)',  '--------------------------------------------------------------------------------'
  call prep_startup
  write(*,*)
end subroutine startup


  !----------------------------------------------------------------------------!
  !>  subroutine: shutdown
  !!  Shutdown call
  !----------------------------------------------------------------------------!
  subroutine shutdown
    call prep_shutdown
    stop 'qprep ended normally'
  end subroutine shutdown

  
  !----------------------------------------------------------------------------!
  !>  subroutine: commandlineoptions
  !----------------------------------------------------------------------------!
  subroutine commandlineoptions
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    select case (arg)
    case ('-v', '--version')
      print '(3a)', PROGRAM_NAME, ' version ', PROGRAM_VERSION
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
  !>  subroutine: print_help
  !----------------------------------------------------------------------------!
  subroutine print_help()
    print '(a)', 'usage:'
    print '(a)', 'qprep [OPTION]'
    print '(a)', '  or'
    print '(a)', 'qprep < inputfile.inp > outputfile.out'
    print '(a)', ''
    print '(a)', 'Without options, qprep goes into interactive mode.'
    print '(a)', ''
    print '(a)', 'qprep [OPTION]:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
  end subroutine print_help


  !----------------------------------------------------------------------------!
  !> function: check_inputfile
  !! Determine if qprep is to be run from command line or from input file  
  !----------------------------------------------------------------------------!
  logical function check_inputfile(infilename)
    !local variables
    integer :: num_args
    character(200), intent(OUT) :: infilename
    character(300) :: text

    ! read name of input file from the command line
    num_args = command_argument_count()
    if (num_args .lt. 1) then
      check_inputfile = .false.
      return
    endif

    call getarg(num_args, infilename)
    text = 'Reading input from '//infilename
    call centered_heading(trim(text), '-')
    check_inputfile = .true.
  end function check_inputfile


end program qprep
