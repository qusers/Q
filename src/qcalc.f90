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
!!  program: **qcalc**  
!!  by John Marelius  
!!  qcalc trajectory analysis main program  
!------------------------------------------------------------------------------!
program qcalc
  use iso_fortran_env

  use trj
  use calc_base
  use calc_rms
  use calc_fit
  use calc_geom
  use calc_entropy
  use calc_nb
  use calc_chemscore
  use calc_xscore
  use calc_pmf
  use calc_rdf
  use calc_rmsf
  use calc_com_ke
  use calc_com
  !add more calculation kind modules here...

  implicit none

  ! version data
  character(*), parameter :: program_name = 'qcalc'
  character(*), parameter :: program_version = '5.7'
  character(*), parameter :: program_date = '2015-02-22'
!  character(*), parameter :: options = compiler_options()
  character(len=32)       :: arg 
  integer                 :: k

  !constants
  integer, parameter      :: max_calcs = 99
  integer, parameter      :: max_calc_kinds = 99
  
  !data types  
  type calc_type
    character(len=60)     :: desc
    integer               :: i
    integer               :: typ
  end type calc_type

  type calc_kind_type
    character(len=40)     :: desc
    character(len=14)     :: key
    logical               :: output
  end type calc_kind_type

  integer                 :: nkinds = 0, ncalcs = 0
  type(calc_kind_type)    :: cdef(max_calc_kinds)
  type(calc_type)         :: calcs(max_calcs)


  !register the calculation kinds
  call add_kind(desc='RMS coord. deviation',          key = 'rmsd', output=.true.)
  call add_kind(desc='Least squares fit',             key = 'fit', output=.false.)
  call add_kind(desc='Distance, bond energy',         key = 'dist', output=.true.)
  call add_kind(desc='Angle, angle energy',           key = 'angle', output=.true.)
  call add_kind(desc='Torsion, torsion energy',       key = 'torsion', output=.true.)
  call add_kind(desc='Entropy [Schlitters formula]',  key = 'entropy', output=.true.) 
  call add_kind(desc='Nonbonded monitor',             key = 'nonbond', output=.true.)
  call add_kind(desc='Residue nonbond monitor',       key = 'nb_prot_qatom', output=.true.)
  call add_kind(desc='ChemScore',                     key = 'chemscore', output=.true.)
  call add_kind(desc='X-Score',                       key = 'xscore', output=.true.)
  call add_kind(desc='PMF-Score',                     key = 'pmfscore', output=.true.)
  call add_kind(desc='RDF',                           key = 'rdf', output=.true.)
  call add_kind(desc='RMSF',                          key = 'rmsf', output=.true.)
  call add_kind(desc='Center of mass Kinetic Energy', key = 'com_ke', output=.true.)
  call add_kind(desc='Center of mass position',       key = 'com', output=.true.)
  !add more calc. kind registrations here...

  
  call commandlineoptions

  ! startup prints welcome msg and calls calc_help to print 
  ! helper info. mask_startup and topo_startup are dummy routines
  call startup

  if(get_topology()) then  ! attempt to load topology (req. user input) reads bnd(:) among other things
    call initialize        ! initialize modules, mostly allocate arrays and define constants
    call add_calcs         ! display calc menu and add calcs
    call do_pre_calcs
    call process_help
    call print_headings
    call process_data      ! processes data frame by frame
    call finalize_topology
    call finalize
  else
    write(*,900)
900 format('>>>>> ERROR: Failed to load topology. Exiting')
  end if

  call shutdown

contains

  !----------------------------------------------------------------------------!
  !!  subroutine: startup  
  !!  Startup  
  !----------------------------------------------------------------------------!
subroutine startup
  integer :: i

!  write(*,'(79a)')('#',i=1,79)
!  print '(a)',  '--------------------------------------------------------------------------------'
!  write(*,'(a,a)') 'Welcome to qcalc version ', program_version
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
  print '(a)',  'For command line options type qcalc --help  or qcalc -h at the terminal.'
  print '(a)',  ' '
  print '(a)',  'If you are using the interactive mode and want to quit type "quit" or Ctrl-C.' 
  print '(a)',  '--------------------------------------------------------------------------------'

  call mask_startup       ! mask_startup calls topo_startup, which is empty
  write(*,*)
  
  call calc_help          ! outputs help to screen
end subroutine startup


  !----------------------------------------------------------------------------!
  !!  subroutine: shutdown
  !!  Shutdown call
  !!  TODO: Implement fully.
  !----------------------------------------------------------------------------!
subroutine shutdown
  call trj_shutdown
end subroutine shutdown


  !----------------------------------------------------------------------------!
  !!  function: get_topology
  !!
  !----------------------------------------------------------------------------!
logical function get_topology()
  !locals
  character(len=200)                      :: topfile
  call getlin(topfile, '--> Topology file: ')
  write(*,100) trim(topfile)
100 format('Loading topology ',a)
  get_topology = topo_load(topfile, require_version=4.)
end function get_topology


subroutine initialize
  allocate(xin(3*nat_pro))   ! contains all coordinates, 3*no_atoms
  call rms_initialize
  call rmsf_initialize
  call fit_initialize
  call geom_initialize
  call entropy_initialize
  call nb_initialize
  call score_initialize
  call xscore_initialize
  call pmf_initialize
  call rdf_initialize
  call com_ke_initialize
  call com_initialize
  ! add more calls to module initialization routines if needed here
end subroutine initialize

subroutine do_pre_calcs
        integer                                         ::      c, i

        do c = 1, Ncalcs
                i = calcs(c)%i
                select case(cdef(calcs(c)%typ)%key)
                        case('chemscore')
                                call score_precalc
                        case('xscore')
                                call xscore_precalc
                        case('pmfscore')
                                call pmf_precalc

       !add more call to calculation pre-routines here...
                end select
         end do 
end subroutine

subroutine finalize
        integer                                         ::      c, i

        deallocate(xin)
        do c=1, Ncalcs
                if(cdef(calcs(c)%typ)%output) then
                        i = calcs(c)%i
                        select case(cdef(calcs(c)%typ)%key)
                        case('rmsd')
                                call rms_finalize(i)
                        case('fit')
                                call fit_finalize(i)
                        case('entropy')
                            call entropy_finalize(i)
                        case('nonbond')
                            call nb_finalize(i)
                        case('nb_prot_qatom')
                            call nb_qp_finalize()
                        case('dist', 'angle', 'torsion', 'improper')
                                !nothing needed
                                
                        case('chemscore')
                                call chemscore_finalize
                        case('xscore')
                                call xscore_finalize
                        case('pmfscore')
                                call pmf_finalize
                        case('rdf')
                                call rdf_finalize(i)
                        case('rmsf')
                                call RMSF_finalize(i)
                        case('com_ke')
                                call COM_KE_finalize(i)
                        case('com')
                                call COM_finalize(i)
                        !add more call to finalize routines here...
                    
                        end select
                end if
        end do
                
        call mask_shutdown
end subroutine finalize


subroutine calc_help
  WRITE( * , '(a)') &
  'qcalc works in three stages:', &
  '1. Load topology.', &
  '2. Set up a list of calculations to make.', &
  '3. Read coordinates from trajectory and/or restart files and calculate.' 
end subroutine calc_help


subroutine process_help
        call centered_heading('Help for reading coordinate files','-')
        WRITE( * , '(a)') &
        'Enter names of coordinate files (restart or trajectory) on separate lines', &
        'Terminate with EOF or period (.).'

end subroutine process_help

subroutine finishtrj_all
        ! this sub invokes corresponding end-of-trajectory functions in calc. modules
        ! when end of trj file is reached
        ! used by score to calculate mean
        
!       call score_calc_mean
        ! add more calls here
end subroutine finishtrj_all

subroutine process_data
        !locals
        character(len=256):: cf !coord file name
        integer                                         :: fn !LUN for cf
        real(8)                                         :: rms
        integer                                         :: frame,p1,p2,frame_start,frame_end,every_n
        character(len=80)       :: frames
    
        do while(get_coordfile(cf, fn,frames))
                every_n = 0
                if(frames.ne.'') then
                        p1 = index(frames,'=')
                        p2 = index(frames,'-')
                        
                        if(p2.ne.0) then
                                read(frames(p1+1:p2-1),*) frame_start
                                read(frames(p2+1:len_trim(frames)),*) frame_end
                        else
                                frames = frames(p1+1:80)
                                if(frames(1:5).eq.'every') then
                                        read(frames(6:80),*) every_n
                                end if
                        end if
                end if
                
                frame = 0
                if(fn /= 0) then
                        ! restart file
                        if(load_restart(fn)) then
                                call process_frame(cf, 0)               ! was process_frame(cf,1)
                        end if
                else
                        ! trajectory file
                        do while(trj_read(xin))
                                ! coordinates for frame "frame" are now in xin(:)
                                frame = frame+1
!                               write(*,*) 'frame= ', frame
                                if(frames.ne.'') then
                                        if(every_n.eq.0) then
                                                if(frame>=frame_start.and.frame<=frame_end) call process_frame(cf, frame)
                                        else
                                                if(modulo(frame,every_n).eq.0) call process_frame(cf,frame)
                                        end if
                                else
                                        call process_frame(cf, frame)
                                end if
                        end do
                        call trj_close
                end if
        end do
end subroutine process_data

subroutine process_frame(cf, frame)
        !arguments
        character(*)                    ::      cf
        integer                                         ::      frame
        character(len=32) :: filename

        if (len_trim(cf) .gt. 20) then
                filename(1:5)  = cf(1:5)
                filename(6:8)  = '...'
                filename(9:20) = cf(len_trim(cf)-11:len_trim(cf))
        else
                filename(:) = ' '
                filename(1:len_trim(cf)) = cf(1:len_trim(cf))
        end if
        
        
        
        write(*, 100, advance='no') filename, frame
100     format(a21, i4)
        call calc_all(frame)
    write(*,*) !new line                                            
end subroutine process_frame

subroutine make_mean_all                                ! calls mean routines in other modules
        integer :: c

        do c = 1, Ncalcs
                select case(cdef(calcs(c)%typ)%key)
                        case('chemscore')
                                call score_mean
                        case('xscore')
                                call xscore_mean
                        case('pmfscore')
                                call pmfscore_mean

       !add more call to mean routines here...
                end select
         end do 
end subroutine make_mean_all    

logical function get_coordfile(cf, fn, frames)
        !arguments
        character(len=256), intent(out) :: cf
        integer, intent(out)                                            :: fn
        character(len=80),optional                      :: frames

        !locals
        integer                                         :: filestat,p
        character(len=4)        :: trj_type

        get_coordfile = .false.
        
        
        do 
                read(*,'(a)', iostat=filestat) cf
                if(filestat /= 0) exit !EOF on input
                if(cf == 'end' .or. cf == 'END' .or. cf == '.') exit
                if(cf == 'mean' .or. cf == 'MEAN') then
                        call make_mean_all
                        cycle
                end if
                if(cf == 'ref' .or. cf == 'REF') then
                        call make_ref_all(xin)
                        write(*,'(a)') 'Updated reference coordinate set.'
                        cycle
                end if

                ! see if user specified frame interval
                p = index(cf,',')
                if(p.ne.0.and.present(frames)) then
                        frames = adjustl(trim(cf(p+1:256)))
                        cf = adjustl(trim(cf(1:p-1)))
                elseif(present(frames)) then
                        frames = ''
                end if

                !open file if filename was specified
                fn = freefile()
                open(unit=fn, file=cf, status='old', action='read', form='unformatted', &
                        iostat=filestat)
                if(filestat /= 0) then
                        write(*,900) trim(cf)
900                     format('>>>>> ERROR: failed to open ',a)
                        cycle
                else
                        !determine type
                        
                        read(fn) trj_type
                        if(trj_type == 'CORD') then !it's a trajectory
                                !it is read by the trajectory module
                                close(fn)
                                fn = 0
                                get_coordfile = trj_open(cf)
                                if(get_coordfile) then
                                        !write the reference coordinates into xin 
                                        !so that atoms not in the trajectory
                                        !get their reference coordinates
                                        xin(:) = xtop(:)
                                end if
                        else !is's a restrart file
                                backspace(fn)
                                get_coordfile = .true.
                        end if
                        exit
                end if
        end do
end function get_coordfile


logical function load_restart(fn)
        !arguments
        integer, intent(in)                     ::      fn
        
        !locals
        integer(4)                                      ::      nat3
        integer                                         ::      filestat

        load_restart = .false.
        rewind(fn)
        read(fn, iostat=filestat) nat3
        if(filestat > 0) then
                write(*,900) 
900             format('>>>>> ERROR: Coordinate file read failure.')
        elseif(filestat < 0) then
                ! EOF - do nothing
        elseif(nat3 == 3) then
                !We found a polarization restraint data record in a restart file
                !do nothing
        elseif(nat3 /= 3*nat_pro) then
                write(*,910) nat3/3
910             format('>>>>> ERROR: Wrong number of atoms in coordinate file:',i5)
        else
                backspace(fn)
                read(fn) nat3, xin(:)
                load_restart = .true.
        end if
        close(fn)
end function load_restart


subroutine finalize_topology
        call topo_deallocate
end subroutine finalize_topology

subroutine add_kind(desc, key, output)
        character(*)            ::      desc
        character(*)            ::      key
        logical                         ::      output
                
        If(Nkinds == MAX_CALC_KINDS) then
                stop 'Too many kinds of calculations.'
        end if
        Nkinds = Nkinds + 1
        cdef(Nkinds)%desc = desc
        cdef(Nkinds)%key = key
        cdef(Nkinds)%output = output
end subroutine add_kind

subroutine calc_menu
        integer                                         ::      i

        call centered_heading('List of calculation types','=')
        write(*,90) ' ', 'description', 'command'
90      format(a2,1x,a,t44,a,t20)
100     format(i2,1x,a,t44,a,t20)
        do i = 1, Nkinds
                write(*, 100) i, cdef(i)%desc, cdef(i)%key
        end do
        write(*,90) ' ', 'Show list of calculations+menu', 'list'
        write(*,90) ' ', 'Proceed to next phase', 'go'
end subroutine calc_menu


subroutine add_calcs
        
        character(len=200)                      ::      input
        integer                                         ::      kind, readstat

        call calc_menu                          ! calc_menu outputs menu and calc. descriptions

        do
                write(*,'(a)', advance='no') 'qcalc> '
                read(*,*) input
                read(input, '(i5)', iostat=readstat) kind
                if(readstat == 0) then !got numeric input
                        if(kind < 1 .or. kind > Nkinds) then
                                write(*,900) kind
900                             format('>>>>> ERROR: Command number ',i2,' is not available.')
                                cycle
                        end if
                else
                        input = adjustl(input)
                        call locase(input)
                        if(input == 'go') exit
                        if(input == 'list') then
                                call list_calcs         ! outputs list of calcs
                                call calc_menu
                                cycle
                        end if
                        do kind = 1, Nkinds
                                if(cdef(kind)%key == input) exit
                        end do
                        if(kind > Nkinds) then
                                write(*,910) trim(input)
910                             format('>>>>> ERROR: The command ',a,' is not available.')
                                cycle
                        end if
                end if
                !at this stage we do have a proper input
                call add_a_calc(kind)
        end do
end subroutine add_calcs

subroutine make_ref_all(x)
        real(8)                                         ::      x(:)

        integer                                         ::      i
        
        !copy current coordinates to xtop
        xtop(:) = xin(:)

        do i = 1, Ncalcs
                select case(cdef(calcs(i)%typ)%key)
                case('rmsd')
                        call rms_make_ref(calcs(i)%i)
                case('rmsf')
                        call RMSF_make_ref(calcs(i)%i)
                case('fit')
                        call fit_make_ref(calcs(i)%i)

                end select
        end do
end subroutine make_ref_all


subroutine list_calcs
        integer                                         ::      c

        call centered_heading('List of defined calculations','=')
        do c=1, Ncalcs
                write(*,100) c, calcs(c)%desc
100             format(i2, 1x, a)
        end do

end subroutine list_calcs

subroutine add_a_calc(kind)
        integer         ::      kind

        if(Ncalcs == MAX_CALCS) then
                write(*,'(a)') '>>>>> ERROR: Too many calculations.'
                return
        end if
        Ncalcs = NCalcs + 1
        calcs(Ncalcs)%typ = kind

        select case(cdef(kind)%key)
                case('rmsd')
                        calcs(Ncalcs)%i = RMS_add(calcs(Ncalcs)%desc)
                case('rmsf')
                        calcs(Ncalcs)%i = RMSF_add(calcs(Ncalcs)%desc)
                case('fit')
                        calcs(Ncalcs)%i = fit_add(calcs(Ncalcs)%desc)
                case('dist')
                        calcs(Ncalcs)%i = dist_add(calcs(Ncalcs)%desc)
                case('angle')
                        calcs(Ncalcs)%i = angle_add(calcs(Ncalcs)%desc)
                case('torsion')
                        calcs(Ncalcs)%i = torsion_add(calcs(Ncalcs)%desc)
                case('entropy')                                             
                  calcs(Ncalcs)%i = entropy_add(calcs(Ncalcs)%desc)
                case('nonbond')                                             
                        calcs(Ncalcs)%i = nb_add(calcs(Ncalcs)%desc)
                case('nb_prot_qatom')                                             
                  calcs(Ncalcs)%i = nb_qp_add(calcs(Ncalcs)%desc)
                case('chemscore')
                        calcs(Ncalcs)%i = score_add(calcs(Ncalcs)%desc)
                case('xscore')
                        calcs(Ncalcs)%i = xscore_add(calcs(Ncalcs)%desc)
                case('pmfscore')
                        calcs(Ncalcs)%i = pmf_add(calcs(Ncalcs)%desc)
                case('rdf')
                        calcs(Ncalcs)%i = rdf_add(calcs(Ncalcs)%desc)
                case('com_ke')
                        calcs(Ncalcs)%i = com_ke_add(calcs(Ncalcs)%desc)
                case('com')
                        calcs(Ncalcs)%i = com_add(calcs(Ncalcs)%desc)
                !add more calls to add routines here...
    end select

        if(calcs(Ncalcs)%i == 0) then !add failed
                !remove the empty entry
                Ncalcs = Ncalcs - 1
        else !successful
                write(*,100) Ncalcs, calcs(Ncalcs)%desc
100             format('Added calc. #',i2,': ',a)
        end if
end subroutine add_a_calc

subroutine print_headings
        integer                                         ::      c, i

        call list_calcs

        call centered_heading('Calculation results', '-')
        write(*,90, advance='no') 'file', 'frame'
90                      format(a, t22, a)

        do c=1, Ncalcs
                if(cdef(calcs(c)%typ)%output) then
                        i = calcs(c)%i
                        write(*, 100, advance='no') c
100                     format(1x,i2,':')
                        select case(cdef(calcs(c)%typ)%key)
                        case('rmsd')
                                call rms_heading(i)
                        case('rmsf')
                                call RMSF_heading(i)
                        case('fit')
                                !no output - no heading
                        case('dist')
                                call dist_heading(i)
                        case('angle')
                                call angle_heading(i)
                        case('torsion')
                                call torsion_heading(i)
                        case('nonbond')
                                call nb_heading(i)
                case('nb_prot_qatom')                             
!                               call nb_heading(i)
                        case('chemscore')
                                call score_heading(i)
                        case('xscore')
                                call xscore_heading(i)
                        case('pmfscore')
                                call pmf_heading(i)
                        case('rdf')
                                call rdf_heading(i)
                        case('com_ke')
                                call com_ke_heading(i)
                        case('com')
                                call com_heading(i)
                        !add more call to heading routines here...
                        end select
                end if
        end do
        write(*,*) !new line

end subroutine print_headings

subroutine calc_all(frame)
    integer             ::  frame
        integer                         ::      c, i

        do c = 1, Ncalcs
                i = calcs(c)%i
                if(cdef(calcs(c)%typ)%output) then
                        write(*,100, advance='no')
100                     format(1x)
                end if
                select case(cdef(calcs(c)%typ)%key)
                        case('rmsd')
                                call rms_calc(i)
                        case('rmsf')
                                call RMSF_calc(i)
                        case('fit')
                                call fit_calc(i)
                        case('dist')
                                call dist_calc(i)
                        case('angle')
                                call angle_calc(i)
                        case('torsion')
                                call torsion_calc(i)
                        case('entropy')                                                  
                                call entropy_calc(i)
                        case('nonbond')                                                         
                                call nb_calc(i)
                        case('nb_prot_qatom')                                                     
                                call nb_qp_calc()
                        case('chemscore')
                                call score_calc(i,frame)
                        case('xscore')
                                call xscore_calc(i,frame)
                        case('pmfscore')
                                call pmf_calc(i,frame)
                        case('rdf')
                                call rdf_calc(i)
                        case('com_ke')
                                call com_ke_calc(i)
                        case('com')
                                call com_calc(i)
       !add more call to calculation routines here...
                end select
         end do 
end subroutine calc_all

  !----------------------------------------------------------------------------!
  !!  subroutine: commandlineoptions  
  !----------------------------------------------------------------------------!
  subroutine commandlineoptions
  do k = 1, command_argument_count()
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
    print '(a)', 'qcalc [OPTION]'
    print '(a)', '  or'
    print '(a)', 'qcalc < inputfile.inp > outputfile.out'
    print '(a)', ''
    print '(a)', 'Without options, qcalc goes into interactive mode.'
    print '(a)', ''
    print '(a)', 'qcalc [OPTION]:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
  end subroutine print_help


end program qcalc

