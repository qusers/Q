!------------------------------------------------------------------------------!
!  Q version 5.7                                                               !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          !
!  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer,     !
!  Miha Purg, Irek Szeler,  Mauricio Esguerra                                  !
!  latest update: Augut 29, 2017                                               !
!------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!!  Copyright (c) 2017 Johan Aqvist, John Marelius, Shina Caroline Lynn Kamerlin
!!  and Paul Bauer
!  sizes.f90
!  by John Marelius
!  data storage specifications for all Q programs
!-------------------------------------------------------------------------------
module sizes
  implicit none

  ! In Fortran 90 and later, intrinsic types such as real and integer have a
  ! kind attribute which guarantees a specific precision and/or range. real*8
  ! and counterparts should no longer be used (Chin, Worth, and Greenough, 2006, p. 5).

  ! Martin Nervall 2002-11-11
  ! set a nice, portable standard for variables
  integer, parameter :: wp8 = selected_real_kind(15,307)

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)

  integer, parameter :: q_int2 = selected_int_kind(4)
  integer, parameter :: q_int4 = selected_int_kind(9)
  integer, parameter :: q_int8 = selected_int_kind(18) 
  integer, parameter :: int_byte = selected_int_kind(2)

  ! storage specifications for q
  ! change according to 
  ! 1) the hardware of your machine
  ! 2) maximum number of atoms desired

  ! max number of atoms is limited only by the size of the variables used 
  ! as indices to the atom arrays.
  ! here we use 16-bit signed integers for indexes to the atom arrays
  ! if more than 2**15-1=32767 atoms are required, change max_atom and ai
  ! integer, parameter :: ai = 2 !integer(ai) = integer(2)
  ! integer, parameter :: max_ai = 2**(8*ai-1)-1

  ! change to this setting to allow up to 2**31-1 = 2147483647 atoms
  integer, parameter :: ai = 4 !integer(ai) = integer(4)

  ! size of integer for flag (1 or 0) variables
  ! change this for machines do byte-wise memory access slowly
  integer, parameter :: flag = 1
  integer, parameter :: max_flag = 2**(8*flag-1)-1

  ! size of integer for small range variables 
  ! (atom types, qatom/qangle/qtorsion/qimproper types,...)
  ! change this for machines do byte-wise memory access slowly
  integer, parameter :: tiny = 2
  integer, parameter :: max_tiny = 2**(8*tiny-1)-1      

  ! size of integer for medium range variables (angle/torsion types, ...)
  ! change this for machines do doublebyte-wise memory access slowly
  integer, parameter :: shrt = 2
  integer, parameter :: max_shrt = 2**(8*shrt-1)-1

  ! note
  ! there are still a couple of hard-coded limitations left. 
  ! the most important of these are
  ! * maxbndtyp   number of entries in bond type library
  ! * maxangtyp   number of entries in angle type library
  ! * max_atyp    number of entries in atom type library
  ! * max_states  number of different states in perturbations
  ! * max_qat     number of q-atoms


  ! note 2  
  ! it is important that you change the types in subroutine 'init_nodes'
  ! in md.f90 when you change the types above. otherwise the datastructures 
  ! sent to the slave nodes will be incorrect. this can be dynamically 
  ! handled with later versions of mpi, but is not supported in the 
  ! current version.

end module sizes

