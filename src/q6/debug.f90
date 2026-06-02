module debug
  use md, only: natom, nwat, restart, v
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
    call write_logical(unit, .not. restart)
    write(unit,'(i0)') natom
    write(unit,'(i0)') nat_solute
    write(unit,'(i0)') nat_solute - n_qatoms
    write(unit,'(i0)') n_qatoms
    write(unit,'(i0)') nwat
    write(unit,'(i0)') nmol
    call write_vectors(unit, v, natom)
    call write_zero_vectors(unit, natom)
    close(unit)
  end subroutine output_ctx_in_file

  subroutine write_vectors(unit, values, n_vectors)
    integer, intent(in) :: unit
    integer, intent(in) :: n_vectors
    real(8), intent(in) :: values(:)
    integer :: i, i3

    write(unit,'(i0)') n_vectors
    do i = 1, n_vectors
      i3 = 3 * i
      write(unit,'(f0.8,1x,f0.8,1x,f0.8)') values(i3 - 2), values(i3 - 1), values(i3)
    end do
  end subroutine write_vectors

  subroutine write_zero_vectors(unit, n_vectors)
    integer, intent(in) :: unit
    integer, intent(in) :: n_vectors
    integer :: i

    write(unit,'(i0)') n_vectors
    do i = 1, n_vectors
      write(unit,'(f0.8,1x,f0.8,1x,f0.8)') 0.0d0, 0.0d0, 0.0d0
    end do
  end subroutine write_zero_vectors

  subroutine write_logical(unit, value)
    integer, intent(in) :: unit
    logical, intent(in) :: value

    write(unit,'(i0)') merge(1, 0, value)
  end subroutine write_logical



end module debug
