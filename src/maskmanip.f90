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
!  maskmanip.f90
!  by John Marelius
!  atom mask manipulation functions
!------------------------------------------------------------------------------!
module maskmanip
  use atom_mask
  use misc
  use parse
  implicit none

contains

  !*******************************************************

  integer function maskmanip_make(mask)
    !arguments
    type(MASK_TYPE)                         ::      mask

    !locals
    character(len=200)                      ::      line
    integer                                         ::      included, inthis
        
    call mask_clear(mask)

    do
      call getlin(line, 'Mask:')
      call upcase(line)
      if(line == 'END' .or. line == '.') exit
      if(line == 'CLEAR') then
        call mask_clear(mask)
        cycle
      elseif(line == 'HELP' .or. line == '?') then
        call  maskmanip_help
        cycle
      end if
      inthis=mask_add(mask, line)
      write(*,110) inthis, mask%included
110   format('Added',i5,' atoms to the mask which now contains',i5,' atoms.')
    end do

    maskmanip_make = mask%included
  end function maskmanip_make

!*******************************************************

subroutine maskmanip_help
        call centered_heading('Help for mask creation','-')
        WRITE( * , '(a)') &
        'The following syntax is used to specify the set of atoms to include:', &
        '[properties] [residue] [first [last]]', &
        '  properties are one or more of:',&
        '    solute: only solute atoms', &
        '    excluded: atoms outside the simulation sphere', &
        '    restrained: atoms outside the simulation sphere and in the restrained shell', &
        '    heavy: heavy (non-hydrogen) atoms', &
        '    sybyl syb.type: atoms of sybyl type "syb.type"', &
        '  Each property may be negated by prepending not.', &
        '  If no property is given, all atoms are selected.',&
        '  residue: first and last are residue numbers, not atom numbers', &
        '  first: first atom/residue of a sequence.', &
        '  last: last atom/residue of a sequence.',&
        '  If first and last are omitted, first=1 and last=number of atoms/residues', &
        '  If last is omitted, last=first.',&
        'Enter atom set specifications and terminate with end.'

end subroutine maskmanip_help

!*******************************************************

integer function maskmanip_make_pretop(mask)
        !arguments
        type(MASK_TYPE)                         ::      mask

        !locals
        character(len=200)                      ::      line
        integer                                         ::      included, inthis
        logical                     ::  pretop=.true.

        call mask_clear(mask)

        do
                call get_line_arg(line, 'Mask:')
                call upcase(line)
                if(line == 'END' .or. line == '.') exit
                if(line == 'CLEAR') then
                        call mask_clear(mask)
                        cycle
                elseif(line == 'HELP' .or. line == '?') then
                        call  maskmanip_help_pretop
                        cycle
                end if
                inthis=mask_add(mask, line, pretop)
                write(*,110) inthis, mask%included
110             format('Added',i5,' atoms to the mask which now contains',i5,' atoms.')
        end do

        maskmanip_make_pretop = mask%included
end function maskmanip_make_pretop

!*******************************************************

subroutine maskmanip_help_pretop
        call centered_heading('Help for mask creation','-')
        WRITE( * , '(a)') &
        'The following syntax is used to specify the set of atoms to include:', &
        '[properties] [residue] [first [last]]', &
        '  properties are one or more of:',&
        '    solute: only solute atoms', &
        '    heavy: heavy (non-hydrogen) atoms', &         
        '  Each property may be negated by prepending not.', &
        '  If no property is given, all atoms are selected.',&
        '  residue: first and last are residue numbers, not atom numbers', &
        '  first: first atom/residue of a sequence.', &
        '  last: last atom/residue of a sequence.',&
        '  If first and last are omitted, first=1 and last=number of atoms/residues', &
        '  If last is omitted, last=first.',&
        'Enter atom set specifications and terminate with end.'

end subroutine maskmanip_help_pretop

!*******************************************************

end module maskmanip
