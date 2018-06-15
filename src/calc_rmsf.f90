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
!  calc_rmsf.f90
!  by Martin Almlof
!  root mean square coordinate deviation
!------------------------------------------------------------------------------!
module calc_rmsf
  use calc_base
  use maskmanip
  implicit none

!constants
  integer, parameter                   :: max_masks = 10

!module variables
  integer, private                     :: frames(max_masks), apa
  real(8), allocatable                 :: rmsf(:)
  type(mask_type), private, target     :: masks(max_masks)
  integer, private                     :: Nmasks = 0
  type rmsf_coord_type
    real, pointer                      :: x(:), x0(:), x2(:)
  end type rmsf_coord_type
  type(rmsf_coord_type), private       :: coords(max_masks)
contains

subroutine rmsf_initialize
end subroutine rmsf_initialize

subroutine rmsf_finalize(i)
  integer                              :: i

  call mask_finalize(masks(i))
end subroutine rmsf_finalize

integer function rmsf_add(desc)
  !arguments
  character(*)                         :: desc
  integer                              :: ats
  if(Nmasks == MAX_MASKS) then
    write(*,10) MAX_MASKS
    return
  end if
10      format('Sorry, the maximum number of RMSF calculations is ',i2)
  !add a new RMSF mask
  Nmasks = Nmasks + 1
  call mask_initialize(masks(Nmasks))
  ats =  maskmanip_make(masks(Nmasks))
  !discard if no atoms in mask
  if(ats == 0) then
    call mask_finalize(masks(Nmasks))
    Nmasks = Nmasks - 1
    RMSF_add = 0
    return
  end if

  allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%x0(3*ats),coords(Nmasks)%x2(3*ats), rmsf(masks(Nmasks)%included))
  coords(Nmasks)%x2(:) = 0
  coords(Nmasks)%x0(:) = 0
  frames(Nmasks) = 0
  call RMSF_make_ref(Nmasks)
  RMSF_add = Nmasks
  write(desc, 20) masks(Nmasks)%included
20      format('RMSF coordinate deviation for ',i6,' atoms')
 end function rmsf_add


subroutine rmsf_calc(i)
  !arguments
  integer, intent(in)                  :: i
  !locals
  real(8)                              :: r
  if(i < 1 .or. i > Nmasks) return
  frames(i)=frames(i) + 1
  call mask_get(masks(i), xin, coords(i)%x)
  !calculate the sum of the squared coordinates up to this frame
  coords(i)%x2 = coords(i)%x2 + (coords(i)%x)**2
  !calculate the average coordinates up to this frame
  coords(i)%x0 = coords(i)%x0*(1._8 - 1._8/frames(i)) + coords(i)%x/frames(i)
  do apa = 1, masks(i)%included
    rmsf(apa) = sqrt(1._8/frames(i) * sum(coords(i)%x2(3*apa-2:3*apa)) - sum((coords(i)%x0(3*apa-2:3*apa))**2))
  end do
  
  r = 1._8/masks(i)%included * sum(rmsf)
  
  write(*,100, advance='no') r
100     format(f10.3)
end subroutine rmsf_calc


subroutine rmsf_make_ref(i)
  integer                              :: i,at

  if(i < 1 .or. i > Nmasks) return
  call mask_get(masks(i), xtop, coords(i)%x0)

end subroutine rmsf_make_ref

subroutine rmsf_heading(i)
  integer                              :: i

  write(*,'(a)', advance='no') 'RMSFd(A)'
end subroutine rmsf_heading

end module calc_rmsf
