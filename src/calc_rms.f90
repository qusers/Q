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
!  calc_rms.f90
!  by John Marelius
!  root mean square coordinate deviation
!------------------------------------------------------------------------------!
module calc_rms
  use calc_base
  use maskmanip
  implicit none

!constants
  integer, parameter              :: MAX_MASKS = 10

!module variables
  type(MASK_TYPE), private, target :: masks(MAX_MASKS)
  integer, private                 :: Nmasks = 0
  type RMS_COORD_TYPE
    real, pointer                  :: x(:), x0(:)
  end type RMS_COORD_TYPE
  type(RMS_COORD_TYPE), private    :: coords(MAX_MASKS)
contains

subroutine rms_initialize
end subroutine rms_initialize

subroutine rms_finalize(i)
  integer                         :: i
  call mask_finalize(masks(i))
end subroutine rms_finalize

integer function rms_add(desc)
  !arguments
  character(*)                    :: desc
  integer                         :: ats
  if(Nmasks == MAX_MASKS) then
    write(*,10) MAX_MASKS
    return
  end if
10      format('Sorry, the maximum number of RMSD calculations is ',i2)
  !add a new RMSD mask
  Nmasks = Nmasks + 1
  call mask_initialize(masks(Nmasks))
  ats =  maskmanip_make(masks(Nmasks))
  !discard if no atoms in mask
  if(ats == 0) then
    call mask_finalize(masks(Nmasks))
    Nmasks = Nmasks - 1
    RMS_add = 0
    return
  end if

  allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%x0(3*ats))

  call rms_make_ref(Nmasks)
  RMS_add = Nmasks
  write(desc, 20) masks(Nmasks)%included
20      format('RMS coordinate deviation for ',i6,' atoms')
 end function rms_add


subroutine rms_calc(i)
  !arguments
  integer, intent(in)             :: i

  !locals
  real(8)                         :: r

  open(10, file='rmsd.dat')

  if(i < 1 .or. i > Nmasks) return

  call mask_get(masks(i), xin, coords(i)%x)
  r = sqrt(  sum((coords(i)%x-coords(i)%x0)**2)/(masks(i)%included) )  ! removed 3 in front of mask(i)
  write(*,100, advance='no') r
  write(10, 100, advance='yes') r
100 format(f8.6)
end subroutine rms_calc


subroutine rms_make_ref(i)
  integer                         :: i,at
  if(i < 1 .or. i > Nmasks) return
  call mask_get(masks(i), xtop, coords(i)%x0)
end subroutine rms_make_ref


subroutine rms_heading(i)
  integer                         :: i
  write(*,'(a)', advance='no') ' RMSD(A)'
end subroutine rms_heading


end module calc_rms
