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
!  parse.f90
!  by John Marelius
!  command parser 
!------------------------------------------------------------------------------!
module parse
  use misc
  implicit none

!constants
  character(*), private, parameter  :: MODULE_VERSION = '5.7'
  character(*), private, parameter  :: MODULE_DATE = '2015-02-22'
  integer, parameter                :: MAX_ARGS = 10

  type substring
   integer                          :: istart, iend
  end type substring

  character(200)                    :: inbuf
  type(substring)                   :: argv(0:MAX_ARGS)
  integer                           :: argc = 0, argp = 1
  
  logical,private                   :: read_from_file = .false.  !Enables input from file
  integer,private                   :: INFILE = 0


contains

!------------------------------------------------------------------------------!
!>  function: parse_startup
!>
!------------------------------------------------------------------------------!
subroutine parse_startup()
  argc = 0
  argp = 1
end subroutine parse_startup


!------------------------------------------------------------------------------!
!>  function: parse_open_file
!>
!------------------------------------------------------------------------------!
logical function parse_open_file(filename)
character(*)    :: filename
integer         :: stat


        parse_open_file = .false.
        read_from_file = .true.

        INFILE = freefile()
        open(unit=INFILE, file = filename, status = 'old', form='formatted', &
                action='read', iostat = stat, access='sequential')
        if (stat /= 0) then
                parse_open_file = .false.
                INFILE = 0
                return
        endif

        parse_open_file = .true.
end function parse_open_file


!------------------------------------------------------------------------------!
!>  function: openit
!!  Moved from Misc
!!  split text in inbuf into 
!------------------------------------------------------------------------------!
integer function openit(lun, fil, stat, frm, mode)
!arguments
        INTEGER lun
        CHARACTER( * ) fil, stat, frm
        character(*)                            ::      mode
!locals
        integer                                         ::      errcode
        
        if(trim(fil).eq.'') write(*,*) 'Invalid filename'

        errcode = -1

        do while(errcode /= 0)
                OPEN(unit = lun, file = fil, status = stat, form = frm, &
                        iostat = errcode, action=mode)
                if(errcode /= 0) then
                        WRITE( * , 10) trim(fil), lun, errcode
                        CALL get_string_arg(fil, ' >>>>> Enter a new file name (or nothing to cancel): ')
                        if(fil == '' .or. fil == 'nothing') exit
                end if
        end do
        openit = errcode

10      format('>>>>> Failed to open ',a,' (unit ',i2,' error code ',i4,')')
end function openit


!------------------------------------------------------------------------------!
!>  subroutine: split
!!  split text in inbuf into 
!------------------------------------------------------------------------------!
subroutine split
  integer                           :: p
  logical                           :: ws_flag 
  logical                           :: quote_flag 
  character, save                   :: TAB = achar(9)
  integer                           :: inlen
  integer                           :: trimlen(1)

        argc = 0
        argp = 1
        ws_flag = .true.
        quote_flag = .false.
        trimlen = len_trim(inbuf)
        inlen = trimlen(1)
        do p = 1, inlen
                if(inbuf(p:p) == '"') then
                        if(quote_flag) then !closing quotation mark
                                quote_flag = .false.
                                ws_flag = .true.
                        else !opening quotation mark
                                quote_flag = .true. !set to accept all characters (like ' ') in string
                                ws_flag = .true. !set to begin string at next positon
                        end if
                elseif((inbuf(p:p) == ' ' .or. inbuf(p:p) == TAB) .and. .not. quote_flag) then
                        if(.not. ws_flag) argv(argc)%iend = p-1
                        ws_flag = .true.
                elseif(ws_flag) then
                        argc = argc + 1
                        if(argc > MAX_ARGS) then
                                write(*,*) 'WARNING: Command argument overflow'
                                exit
                        end if
                        argv(argc)%istart = p
                        argv(argc)%iend = p !step up character by character
                        ws_flag = .false.
                else
                        argv(argc)%iend = p !step up character by character
                        ws_flag = .false.
                end if
        end do
end subroutine split


!-------------------------------------------------------------------------------!
!>  subroutine: **getline**
!!  goes line by line in input files.
!!  TODO: In case there is no quit card give a different error to:
!!    Fortran runtime error: End of file
!!  which is now the case.
!-------------------------------------------------------------------------------!
subroutine getline()
  do
     if (read_from_file) then
        read(INFILE, '(a200)') inbuf
     else
        read(*, '(a200)') inbuf
     endif
     inbuf = adjustl(inbuf)
     !only exit if not a comment line
     if(scan(inbuf(1:1), '!#*') == 0) exit
  end do
  call split
end subroutine getline


!-------------------------------------------------------------------------------!
!>  subroutine: **get_string_arg**
!! 
!-------------------------------------------------------------------------------!
subroutine get_string_arg(arg, prompt)
  !arguments
  character(*), intent(out)               :: arg
  character(*), optional, intent(in)      :: prompt

  do while(argp > argc)
     if(present(prompt))     write(unit=*, fmt='(a)', advance='no') prompt
     call getline
  end do
  arg = inbuf(argv(argp)%istart :argv(argp)%iend)
  argp = argp + 1                
end subroutine get_string_arg


!--------------------------------------------------------------------------
subroutine get_line_arg(arg, prompt)
!arguments
        character(*), intent(out)       ::      arg
        character(*), optional, intent(in)      ::      prompt
        
        do while(argp > argc)
                if(present(prompt))     write(unit=*, fmt='(a)', advance='no') prompt
                call getline
        end do
        !return the rest of the line    
        arg = inbuf(argv(argp)%istart : len_trim(inbuf))
        call parse_reset !reset argument counter - new line next time           
end subroutine get_line_arg

!--------------------------------------------------------------------------

logical function get_string_single_line(arg, prompt)
!arguments
        character(*), intent(out)       ::      arg
        character(*), optional, intent(in)      ::      prompt
        
        if(argc == 0) then
                if(present(prompt))     write(unit=*, fmt='(a)', advance='no') prompt
                call getline
        end if
        if(argp>argc) then
                arg = ''
                get_string_single_line = .false.
        else
                arg = inbuf(argv(argp)%istart :argv(argp)%iend)
                argp = argp + 1
                get_string_single_line = .true.
        end if          
                
end function get_string_single_line

!--------------------------------------------------------------------------

subroutine parse_reset
        argc = 0
        argp = 1
end subroutine parse_reset

!--------------------------------------------------------------------------

integer function get_int_arg(prompt)
!arguments
        character(*), optional, intent(in)      ::      prompt
!locals
        integer                                         ::      value

1       do while(argp > argc)
                if(present(prompt))     write(unit=*, fmt='(a)', advance='no') prompt
                call getline
        end do
        
        read(inbuf(argv(argp)%istart :argv(argp)%iend), fmt=*, err=100) value

        get_int_arg = value

        argp = argp + 1
        return

100     write(*,*) 'Please enter an integer value!'
        call parse_reset
        goto 1
end function get_int_arg

!--------------------------------------------------------------------------

real function get_real_arg(prompt)
!arguments
        character(*), optional, intent(in)      ::      prompt
!locals
        real                                            ::      value

1       do while(argp > argc)
                if(present(prompt))     write(unit=*, fmt='(a)', advance='no') prompt
                call getline
        end do
        
        read(inbuf(argv(argp)%istart :argv(argp)%iend), fmt=*, err=100) value
        get_real_arg = value

        argp = argp + 1
        return

100     write(*,*) 'Please enter a number!'
        call parse_reset
        goto 1
end function get_real_arg

!--------------------------------------------------------------------------

end module parse
