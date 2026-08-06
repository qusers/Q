program test_lincs
  use lincs
  implicit none

  integer :: failures = 0
  integer :: assertions = 0

  call test_single_constraint
  call test_coupled_chain
  call test_expansion_convergence
  call test_orientation_invariance
  call test_rigid_triangle
  call test_momentum_conservation
  call test_large_sparse_setup
  call test_initial_projection
  call test_rotation_domain_failure

  write(*,'(a,i0,a,i0,a)') 'LINCS tests: ', assertions-failures, '/', assertions, ' passed'
  if (failures > 0) error stop 'LINCS tests failed'

contains

subroutine test_single_constraint
  integer :: atom_i(1), atom_j(1)
  real(8) :: target(1), inverse_mass(2), old(6), candidate(6), error
  logical :: success

  atom_i = 1; atom_j = 2; target = 1.09d0
  inverse_mass = (/ 1.0d0/12.0d0, 1.0d0/1.008d0 /)
  old = (/ 0.0d0, 0.0d0, 0.0d0, 1.09d0, 0.0d0, 0.0d0 /)
  candidate = old + (/ 0.01d0, -0.02d0, 0.00d0, 0.04d0, 0.05d0, 0.01d0 /)

  call setup_lincs(atom_i, atom_j, target, inverse_mass, 4, 1)
  success = lincs_positions(old, candidate, max_relative_error=error)

  call assert_true('single constraint succeeds', success)
  call assert_close('single constraint length', bond_length(candidate,1,2), target(1), 2.0d-14)
  call assert_close('single constraint reported error', error, 0.0d0, 2.0d-14)
end subroutine test_single_constraint


subroutine test_coupled_chain
  integer, parameter :: atom_count = 8, constraint_count = atom_count-1
  integer :: atom_i(constraint_count), atom_j(constraint_count), k
  real(8) :: target(constraint_count), inverse_mass(atom_count)
  real(8) :: old(3*atom_count), by_lincs(3*atom_count), by_shake(3*atom_count)
  real(8) :: error
  logical :: success

  call make_zigzag_chain(atom_count, old)
  do k = 1, constraint_count
    atom_i(k)=k; atom_j(k)=k+1; target(k)=1.0d0
  end do
  inverse_mass = 1.0d0/12.0d0
  inverse_mass(2:atom_count:2) = 1.0d0
  by_lincs = old
  do k = 1, 3*atom_count
    by_lincs(k) = by_lincs(k) + 0.015d0*sin(1.7d0*dble(k))
  end do
  by_shake = by_lincs

  call setup_lincs(atom_i, atom_j, target, inverse_mass, 8, 2)
  success = lincs_positions(old, by_lincs, max_relative_error=error)
  call converged_shake(old, by_shake, atom_i, atom_j, target, inverse_mass)

  call assert_true('coupled chain succeeds', success)
  call assert_true('coupled chain constraints are accurate', error < 5.0d-7)
  call assert_true('coupled LINCS agrees with converged SHAKE', &
                   maxval(abs(by_lincs-by_shake)) < 5.0d-7)
end subroutine test_coupled_chain


subroutine test_expansion_convergence
  integer, parameter :: atom_count=7, constraint_count=6
  integer :: atom_i(constraint_count), atom_j(constraint_count), k
  real(8) :: target(constraint_count), inverse_mass(atom_count), old(3*atom_count)
  real(8) :: low_order(3*atom_count), high_order(3*atom_count), low_error, high_error
  logical :: success

  call make_zigzag_chain(atom_count, old)
  do k=1,constraint_count
    atom_i(k)=k; atom_j(k)=k+1; target(k)=1.0d0
  end do
  inverse_mass=1.0d0
  low_order=old
  do k=1,3*atom_count
    low_order(k)=low_order(k)+0.03d0*cos(0.9d0*dble(k))
  end do
  high_order=low_order

  call setup_lincs(atom_i,atom_j,target,inverse_mass,0,1)
  success=lincs_positions(old,low_order,max_relative_error=low_error)
  call assert_true('zeroth-order chain succeeds',success)
  call setup_lincs(atom_i,atom_j,target,inverse_mass,8,1)
  success=lincs_positions(old,high_order,max_relative_error=high_error)

  call assert_true('high-order chain succeeds',success)
  call assert_true('higher expansion order reduces error',high_error < low_error)
  call assert_true('order-eight chain accuracy',high_error < 1.0d-5)
end subroutine test_expansion_convergence


subroutine test_orientation_invariance
  integer :: forward_i(3)=(/1,2,3/),forward_j(3)=(/2,3,4/)
  integer :: reversed_i(3)=(/2,2,4/),reversed_j(3)=(/1,3,3/)
  real(8) :: target(3)=(/1.0d0,1.0d0,1.0d0/),inverse_mass(4)=(/0.2d0,1.0d0,0.1d0,0.5d0/)
  real(8) :: old(12),forward(12),reversed(12)
  logical :: success

  call make_zigzag_chain(4,old)
  forward=old+(/0.01d0,0.00d0,-0.01d0, -0.02d0,0.03d0,0.01d0, &
                0.02d0,-0.01d0,0.00d0, 0.00d0,0.02d0,-0.02d0/)
  reversed=forward
  call setup_lincs(forward_i,forward_j,target,inverse_mass,8,2)
  success=lincs_positions(old,forward)
  call assert_true('forward-oriented constraints succeed',success)
  call setup_lincs(reversed_i,reversed_j,target,inverse_mass,8,2)
  success=lincs_positions(old,reversed)
  call assert_true('reversed constraints succeed',success)
  call assert_vector_close('constraint orientation does not change LINCS', &
                           reversed,forward,2.0d-14)
end subroutine test_orientation_invariance


subroutine test_rigid_triangle
  integer :: atom_i(3)=(/1,1,2/), atom_j(3)=(/2,3,3/)
  real(8) :: target(3)=(/0.9572d0,0.9572d0,1.5139d0/)
  real(8) :: inverse_mass(3)=(/1.0d0/15.9994d0,1.0d0/1.008d0,1.0d0/1.008d0/)
  real(8) :: old(9), candidate(9), height, error
  logical :: success

  height=sqrt(target(1)**2-0.25d0*target(3)**2)
  old=(/0.0d0,0.0d0,0.0d0,-0.5d0*target(3),-height,0.0d0, &
       0.5d0*target(3),-height,0.0d0/)
  candidate=old+(/0.01d0,-0.01d0,0.02d0,-0.02d0,0.03d0,-0.01d0, &
                  0.02d0,0.01d0,0.015d0/)

  call setup_lincs(atom_i,atom_j,target,inverse_mass,12,3)
  success=lincs_positions(old,candidate,max_relative_error=error)

  call assert_true('rigid triangle succeeds',success)
  call assert_true('rigid triangle accuracy',error < 1.0d-7)
end subroutine test_rigid_triangle


subroutine test_momentum_conservation
  integer :: atom_i(3)=(/1,1,1/), atom_j(3)=(/2,3,4/), atom
  real(8) :: target(3)=(/1.0d0,1.0d0,1.0d0/)
  real(8) :: inverse_mass(4)=(/1.0d0/12.0d0,1.0d0,1.0d0,1.0d0/)
  real(8) :: old(12), candidate(12), before(3), after(3), mass
  logical :: success

  old=(/0.0d0,0.0d0,0.0d0, 1.0d0,0.0d0,0.0d0, &
       -0.5d0,sqrt(0.75d0),0.0d0, -0.5d0,-sqrt(0.75d0),0.0d0/)
  candidate=old+(/0.01d0,0.02d0,-0.01d0, 0.03d0,-0.02d0,0.01d0, &
                  -0.02d0,0.01d0,0.03d0, 0.01d0,-0.03d0,-0.02d0/)
  before=0.0d0
  do atom=1,4
    mass=1.0d0/inverse_mass(atom)
    before=before+mass*(candidate(3*atom-2:3*atom)-old(3*atom-2:3*atom))
  end do
  call setup_lincs(atom_i,atom_j,target,inverse_mass,8,2)
  success=lincs_positions(old,candidate)
  after=0.0d0
  do atom=1,4
    mass=1.0d0/inverse_mass(atom)
    after=after+mass*(candidate(3*atom-2:3*atom)-old(3*atom-2:3*atom))
  end do

  call assert_true('momentum test succeeds',success)
  call assert_vector_close('constraint correction conserves momentum',after,before,2.0d-13)
end subroutine test_momentum_conservation


subroutine test_large_sparse_setup
  integer, parameter :: atom_count=10000, constraint_count=atom_count-1
  integer, allocatable :: atom_i(:), atom_j(:)
  real(8), allocatable :: target(:), inverse_mass(:), old(:), candidate(:)
  real(8) :: error
  integer :: k
  logical :: success

  allocate(atom_i(constraint_count),atom_j(constraint_count),target(constraint_count), &
           inverse_mass(atom_count),old(3*atom_count),candidate(3*atom_count))
  old=0.0d0
  do k=1,atom_count
    old(3*k-2)=dble(k-1)
  end do
  do k=1,constraint_count
    atom_i(k)=k;atom_j(k)=k+1;target(k)=1.0d0
  end do
  inverse_mass=1.0d0
  candidate=old

  call setup_lincs(atom_i,atom_j,target,inverse_mass,4,1)
  success=lincs_positions(old,candidate,max_relative_error=error)
  call assert_true('large sparse chain succeeds',success)
  call assert_close('large sparse chain remains exact',error,0.0d0,1.0d-15)
  deallocate(atom_i,atom_j,target,inverse_mass,old,candidate)
end subroutine test_large_sparse_setup


subroutine test_initial_projection
  integer :: atom_i(3)=(/1,2,3/),atom_j(3)=(/2,3,4/),failed
  real(8) :: target(3)=(/1.0d0,1.0d0,1.0d0/),inverse_mass(4)=(/1.0d0,1.0d0,1.0d0,1.0d0/)
  real(8) :: coordinates(12),error
  logical :: success

  coordinates=(/0.0d0,0.0d0,0.0d0, 1.25d0,0.1d0,0.0d0, &
                2.1d0,0.3d0,0.0d0, 3.3d0,0.2d0,0.1d0/)
  call setup_lincs(atom_i,atom_j,target,inverse_mass,8,2)
  success=initialize_lincs_positions(coordinates,failed,error)
  call assert_true('initial nonlinear projection succeeds',success)
  call assert_true('initial nonlinear projection converges tightly',error<1.0d-12)
end subroutine test_initial_projection


subroutine test_rotation_domain_failure
  integer :: atom_i(1)=(/1/),atom_j(1)=(/2/),failed
  real(8) :: target(1)=(/1.0d0/),inverse_mass(2)=(/1.0d0,1.0d0/)
  real(8) :: old(6)=(/0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0/),candidate(6)
  logical :: success

  candidate=(/0.0d0,0.0d0,0.0d0,1.0d0,2.0d0,0.0d0/)
  call setup_lincs(atom_i,atom_j,target,inverse_mass,4,1)
  success=lincs_positions(old,candidate,failed)
  call assert_true('impossible rotation reports failure',.not.success)
  call assert_true('failed LINCS constraint is identified',failed==1)
end subroutine test_rotation_domain_failure


subroutine make_zigzag_chain(atom_count,coordinates)
  integer,intent(in)::atom_count
  real(8),intent(out)::coordinates(3*atom_count)
  integer::atom

  coordinates=0.0d0
  do atom=2,atom_count
    coordinates(3*atom-2)=coordinates(3*(atom-1)-2)+0.8d0
    coordinates(3*atom-1)=coordinates(3*(atom-1)-1)+merge(0.6d0,-0.6d0,mod(atom,2)==0)
  end do
end subroutine make_zigzag_chain


subroutine converged_shake(old,new,atom_i,atom_j,target,inverse_mass)
  real(8),intent(in)::old(:),target(:),inverse_mass(:)
  real(8),intent(inout)::new(:)
  integer,intent(in)::atom_i(:),atom_j(:)
  integer::iteration,k,i,j
  real(8)::current(3),reference(3),difference,correction,max_error

  do iteration=1,10000
    do k=1,size(atom_i)
      i=atom_i(k);j=atom_j(k)
      current=new(3*i-2:3*i)-new(3*j-2:3*j)
      reference=old(3*i-2:3*i)-old(3*j-2:3*j)
      difference=target(k)**2-dot_product(current,current)
      correction=difference/(2.0d0*dot_product(current,reference)* &
                 (inverse_mass(i)+inverse_mass(j)))
      new(3*i-2:3*i)=new(3*i-2:3*i)+reference*correction*inverse_mass(i)
      new(3*j-2:3*j)=new(3*j-2:3*j)-reference*correction*inverse_mass(j)
    end do
    max_error=0.0d0
    do k=1,size(atom_i)
      max_error=max(max_error,abs(bond_length(new,atom_i(k),atom_j(k))/target(k)-1.0d0))
    end do
    if(max_error<2.0d-14)return
  end do
  error stop 'reference SHAKE failed to converge'
end subroutine converged_shake


pure real(8) function bond_length(coordinates,atom_i,atom_j)
  real(8),intent(in)::coordinates(:)
  integer,intent(in)::atom_i,atom_j
  real(8)::difference(3)
  difference=coordinates(3*atom_i-2:3*atom_i)-coordinates(3*atom_j-2:3*atom_j)
  bond_length=sqrt(dot_product(difference,difference))
end function bond_length


subroutine assert_close(name,actual,expected,tolerance)
  character(*),intent(in)::name
  real(8),intent(in)::actual,expected,tolerance
  assertions=assertions+1
  if(abs(actual-expected)>tolerance)then
    failures=failures+1
    write(*,'(a,a,a,es14.6)')'FAIL: ',name,': difference=',abs(actual-expected)
  end if
end subroutine assert_close


subroutine assert_vector_close(name,actual,expected,tolerance)
  character(*),intent(in)::name
  real(8),intent(in)::actual(:),expected(:),tolerance
  assertions=assertions+1
  if(maxval(abs(actual-expected))>tolerance)then
    failures=failures+1
    write(*,'(a,a,a,es14.6)')'FAIL: ',name,': max difference=',maxval(abs(actual-expected))
  end if
end subroutine assert_vector_close


subroutine assert_true(name,condition)
  character(*),intent(in)::name
  logical,intent(in)::condition
  assertions=assertions+1
  if(.not.condition)then
    failures=failures+1
    write(*,'(a,a)')'FAIL: ',name
  end if
end subroutine assert_true

end program test_lincs
