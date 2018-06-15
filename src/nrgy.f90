!------------------------------------------------------------------------------!
!  Q version 5.7                                                               !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          !
!  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer,     !
!  Miha Purg, Irek Szeler,  Mauricio Esguerra                                  !
!  latest update: August 29, 2017                                              !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  Copyright (c) 2017 Johan Aqvist, John Marelius, Shina Caroline Lynn Kamerlin
!  and Paul Bauer
!>  **nrgy.f90**
!!  by John Marelius
!!  energy data and energy file I/O
!------------------------------------------------------------------------------!
module nrgy
  use sizes

  implicit none

  character(*), parameter     :: NRGY_VERSION = '5.7'
  character(*), parameter     :: NRGY_DATE = '2015-02-22'

  type bonded_energies
    sequence
    real(8)                   :: bond, angle, torsion, improper
  end type bonded_energies

  type nb_energies
     sequence
     real(8)                  :: el, vdw
  end type nb_energies

  type restraint_energies
     sequence
     real(8)                  :: total, fix, shell, protein
     real(8)                  :: solvent_radial, water_pol
  end type restraint_energies

  type energies
    sequence
    real(8)                   :: potential, kinetic, LRF
    type(bonded_energies)     :: p, w, q
    type(nb_energies)         :: pp, pw, ww, qx
    type(restraint_energies)  :: restraint
  end type energies

  type q_energies
     sequence
     real(8)                  :: lambda
     real(8)                  :: total
     type(bonded_energies)    :: q
     type(nb_energies)        :: qx, qq, qp, qw
     real(8)                  :: restraint
  end type q_energies

  type offdiag_save
     sequence        
     !integer(4) avoids unaligned access & is compatible with qdyn v2
     integer(4)               :: i, j
     real(8)                  :: Hij, rkl
  end type offdiag_save

  type offdiag_aux
     integer(4)               :: k, l
     real(8)                  :: A, mu
  end type offdiag_aux

  interface operator(+)
     module procedure add_ene
  end interface
!  end interface operator(+)

  interface operator(*)
     module procedure scale_ene
  end interface
!  end interface operator(*)


contains
  !----------------------------------------------------------------------
  subroutine nrgy_startup

  end subroutine nrgy_startup


!------------------------------------------------------------------------------!
!>  **subroutine: put_ene**
!!    defines what goes in the energy files.
!------------------------------------------------------------------------------!
  subroutine put_ene(unit, e2, offd)
    !arguments
    integer                       :: unit
    type(q_energies), intent(in)  :: e2(:)
    type(offdiag_save), intent(in):: offd(:)

    !local variables
    integer                       :: i, bound(1), first, last

    bound = lbound(e2)
    first = bound(1)
    bound = ubound(e2)
    last = bound(1)

    do i=first, last
       write (unit) i, e2(i)
    end do

    bound = lbound(offd)
    first = bound(1)
    bound = ubound(offd)
    last = bound(1)
    write(unit) offd(first:last)

    !       write (unit) (offd(i)%i, offd(i)%j, offd(i)%Hij, offd(i)%rkl, i=first, last)
  end subroutine put_ene


  !----------------------------------------------------------------------
  integer function get_ene(unit, e2, offd, nstates, noffd)
    !arguments
    integer                                 :: unit
    type(q_energies), intent(out)           :: e2(:)
    type(offdiag_save), intent(out)         :: offd(:)
    integer, optional                       :: nstates, noffd

    !local variables
    integer                                 :: i, bound(1), first, last, dummy
    if(present(nstates)) then
       first = 1
       last = nstates
    else
       bound = lbound(e2)
       first = bound(1)
       bound = ubound(e2)
       last = bound(1)
    end if

    do i=first, last
       read (unit, end=10) dummy, e2(i)
    end do

    if(present(noffd)) then 
       first = 1
       last = noffd
    else
       bound = lbound(offd)
       first = bound(1)
       bound = ubound(offd)
       last = bound(1)
    end if
    !       read(unit, end=20) (offd(i)%i, offd(i)%j, offd(i)%Hij, offd(i)%rkl, i=first, last)
    read(unit, end=20) offd(first:last)

    get_ene = 0 !it's OK
    return

10  get_ene = i !failed energies
    return  
20  get_ene = -1 !failed offd
    return  

  end function get_ene


  !----------------------------------------------------------------------
  function add_ene (e1, e2)
    type(q_energies), intent (in) :: e1 (:), e2 (size (e1))
    type(q_energies) :: add_ene (size (e1))

    add_ene(:)%total=e1(:)%total+e2(:)%total
    add_ene(:)%q%bond=e1(:)%q%bond+e2(:)%q%bond
    add_ene(:)%q%angle=e1(:)%q%angle+e2(:)%q%angle
    add_ene(:)%q%torsion=e1(:)%q%torsion+e2(:)%q%torsion
    add_ene(:)%q%improper=e1(:)%q%improper+e2(:)%q%improper
    add_ene(:)%qx%el=e1(:)%qx%el+e2(:)%qx%el
    add_ene(:)%qx%vdw=e1(:)%qx%vdw+e2(:)%qx%vdw
    add_ene(:)%qq%el=e1(:)%qq%el+e2(:)%qq%el
    add_ene(:)%qq%vdw=e1(:)%qq%vdw+e2(:)%qq%vdw
    add_ene(:)%qp%el=e1(:)%qp%el+e2(:)%qp%el
    add_ene(:)%qp%vdw=e1(:)%qp%vdw+e2(:)%qp%vdw
    add_ene(:)%qw%el=e1(:)%qw%el+e2(:)%qw%el
    add_ene(:)%qw%vdw=e1(:)%qw%vdw+e2(:)%qw%vdw
    add_ene(:)%restraint=e1(:)%restraint+e2(:)%restraint

  end function add_ene


  !----------------------------------------------------------------------
  function scale_ene (e1, k)
    type(q_energies), intent (in) :: e1 (:)
    real, intent(in)              :: k
    type(q_energies)              :: scale_ene (size (e1))

    scale_ene(:)%total=e1(:)%total*k
    scale_ene(:)%q%bond=e1(:)%q%bond*k
    scale_ene(:)%q%angle=e1(:)%q%angle*k
    scale_ene(:)%q%torsion=e1(:)%q%torsion*k
    scale_ene(:)%q%improper=e1(:)%q%improper*k
    scale_ene(:)%qx%el=e1(:)%qx%el*k
    scale_ene(:)%qx%vdw=e1(:)%qx%vdw*k
    scale_ene(:)%qq%el=e1(:)%qq%el*k
    scale_ene(:)%qq%vdw=e1(:)%qq%vdw*k
    scale_ene(:)%qp%el=e1(:)%qp%el*k
    scale_ene(:)%qp%vdw=e1(:)%qp%vdw*k
    scale_ene(:)%qw%el=e1(:)%qw%el*k
    scale_ene(:)%qw%vdw=e1(:)%qw%vdw*k
    scale_ene(:)%restraint=e1(:)%restraint*k

  end function scale_ene


  !----------------------------------------------------------------------
  real(8) function sum_bonded(r, eb)
    real(8), intent(in)                     :: r
    type(bonded_energies), intent(in)       :: eb

    sum_bonded = r + eb%bond + eb%angle + eb%torsion + eb%improper
  end function sum_bonded


  !----------------------------------------------------------------------
  real(8) function sum_non_bonded(r, enb)
    real(8), intent(in)                     :: r
    type(nb_energies), intent(in)           :: enb

    sum_non_bonded = r + enb%el + enb%vdw
  end function sum_non_bonded


end module nrgy
