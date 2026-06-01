module debug
  use md_test, only: dt, natom, nwat, restart, tau_T
  use qatom, only: nqat, nstates
  use topo, only: nat_solute, nmol
  use mpiglob, only: nodeid

  implicit none
  private
  public :: output_ctx_in_file

contains

  subroutine output_ctx_in_file()
    integer :: unit
    integer :: n_qatoms

    if (nodeid .ne. 0) return

    n_qatoms = 0
    if (nstates .gt. 0) n_qatoms = nqat

    open(newunit=unit, file='fortran_context_debug.txt', status='replace', action='write')

    ! write(unit,'(i0)') merge(1, 0, .not. restart)
    ! write(unit,'(i0)') natom
    ! write(unit,'(i0)') nat_solute
    ! write(unit,'(i0)') nat_solute - n_qatoms
    ! write(unit,'(i0)') n_qatoms
    ! write(unit,'(i0)') nwat
    ! write(unit,'(i0)') nmol
    write(unit,'(es24.16)') dt
    write(unit,'(es24.16)') tau_T

    close(unit)
  end subroutine output_ctx_in_file

end module debug
