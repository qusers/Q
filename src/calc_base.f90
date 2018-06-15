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
!  (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!  calc_base.f90
!  by John Marelius
!  shared data for trajectory calculation modules
!------------------------------------------------------------------------------!
module calc_base
use topo
implicit none

! In Fortran 90 and later, intrinsic types such as real and integer have a
! kind attribute which guarantees a specific precision and/or range. real*8
! and counterparts should no longer be used (Chin, Worth, and Greenough, 2006, p. 5).
  real(kind=dp), allocatable :: xin(:) ! All functions share this vector which contains a frame (J)
! real(8), allocatable :: xin(:) ! All functions share this vector which contains a frame (J)
! real(8), allocatable :: xin2d(:,:) ! Schlitters formula needs all coordinates at the same time

end module calc_base
