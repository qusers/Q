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
!  mpiglob.f90
!  by John Marelius & Anders Kaplan
!  global variables for MPI parallell qdyn
!------------------------------------------------------------------------------!
module  mpiglob

  use nrgy

  ! types and variables used for calculation assignment
  type PAIR_ASSIGNMENT_TYPE
     integer                                     :: start, end
     integer                                     :: max
  end type PAIR_ASSIGNMENT_TYPE


  type NODE_ASSIGNMENT_TYPE
     type(PAIR_ASSIGNMENT_TYPE) :: pp, pw, ww, qp, qw
  end type NODE_ASSIGNMENT_TYPE

  !Gathering
  type MPI_NB_ENERGIES
     sequence
     real(kind=dp)       :: lrf
     type(NB_ENERGIES)    :: pp,pw,ww
  end type MPI_NB_ENERGIES

  !Gathering
  type MPI_NBQ_ENERGIES
     sequence
     type(NB_ENERGIES)    :: qp,qw
  end type MPI_NBQ_ENERGIES

  ! global MPI data
  integer                  :: nodeid, numnodes,ierr
  integer,allocatable      :: mpi_status(:,:)  

  !Used for gathering d,E,EQ
  real(kind=dp),allocatable         :: d_recv(:,:)
  type(MPI_NB_ENERGIES),allocatable  :: E_recv(:),E_send(:)
  type(MPI_NBQ_ENERGIES),allocatable :: EQ_recv(:,:),EQ_send(:)
  integer,allocatable                :: request_recv(:,:)

  !Book keeping of nb-pairs
  integer  :: totnbpp,totnbpw,totnbww,totnbqp,totnbqw
  integer,dimension(5)   ::nbxx,nbxx_tot

  !Balancing nodes
  !Keep track of how many nb-pairs each chargegroup will generate
  integer,allocatable  :: nbpp_per_cgp(:)
  integer,allocatable  :: nbpw_per_cgp(:)
  integer,allocatable  :: nbww_per_cgp(:)
  integer,allocatable  :: nbqp_per_cgp(:)
  integer,allocatable  :: nbqw_per_cgp(:)

end module mpiglob

