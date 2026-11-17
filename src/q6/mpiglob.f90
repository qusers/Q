!------------------------------------------------------------------------------!
!  Q version 6.0.1                                                             !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          !
!  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer,     !
!  Miha Purg, Irek Szeler, Mauricio Esguerra                                   !
!  latest update: August 29, 2017                                              !
!------------------------------------------------------------------------------!


module  mpiglob
!!-------------------------------------------------------------------------------
!!  Copyright (c) 2017 Johan Aqvist, John Marelius, Shina Caroline Lynn Kamerlin
!!  and Paul Bauer
!!  **module mpiglob**
!!  by John Marelius & Anders Kaplan
!!  global variables for MPI parallell qdyn
!!-------------------------------------------------------------------------------
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

contains
! BEGIN GENERATED reset_state (tools/reset_state.py; do not edit by hand)
subroutine mpiglob_reset_state
!! Every variable of this module back to its value at program start (qdyn --serve).
  integer :: reset_status
  reset_status = 0
  ierr = 0
  if (allocated(mpi_status)) deallocate(mpi_status)
  if (allocated(d_recv)) deallocate(d_recv)
  if (allocated(E_recv)) deallocate(E_recv)
  if (allocated(E_send)) deallocate(E_send)
  if (allocated(EQ_recv)) deallocate(EQ_recv)
  if (allocated(EQ_send)) deallocate(EQ_send)
  if (allocated(request_recv)) deallocate(request_recv)
  totnbpp = 0
  totnbpw = 0
  totnbww = 0
  totnbqp = 0
  totnbqw = 0
  nbxx = 0
  nbxx_tot = 0
  if (allocated(nbpp_per_cgp)) deallocate(nbpp_per_cgp)
  if (allocated(nbpw_per_cgp)) deallocate(nbpw_per_cgp)
  if (allocated(nbww_per_cgp)) deallocate(nbww_per_cgp)
  if (allocated(nbqp_per_cgp)) deallocate(nbqp_per_cgp)
  if (allocated(nbqw_per_cgp)) deallocate(nbqw_per_cgp)
end subroutine mpiglob_reset_state
! END GENERATED reset_state
end module mpiglob

