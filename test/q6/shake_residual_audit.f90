! Diagnostic only: call Q's unchanged constraint solver on one water.
program shake_residual_audit
  use md
  implicit none
  real(8) :: reference(9), trial(9), residual, oh, hh, a, b, worst
  integer :: iterations, ic, left, right, sample, worst_sample, worst_bond, ready
  natom=3; nmol=1; shake_molecules=1
  allocate(winv(3), shake_mol(1))
  winv=[1.0_8/15.9994_8,1.0_8/1.008_8,1.0_8/1.008_8]
  shake_mol(1)%nconstraints=3
  allocate(shake_mol(1)%bond(3))
  oh=.9572_8; hh=1.5136_8
  a=(2*oh**2-hh**2)/(2*oh); b=sqrt(oh**2-a**2)
  reference=[0.0_8,0.0_8,0.0_8,oh,0.0_8,0.0_8,a,b,0.0_8]
  do ic=1,3
    left=1; right=ic+1
    if (ic == 3) then
      left=2; right=3
    end if
    shake_mol(1)%bond(ic)%i=left
    shake_mol(1)%bond(ic)%j=right
    shake_mol(1)%bond(ic)%dist2=sum((reference(3*left-2:3*left)-reference(3*right-2:3*right))**2)
  end do
  worst=0; ready=0
  do sample=1,64
    trial=reference
    ! Both OH lengths start exact; perturb only the angle by up to 0.03 rad.
    ! Correcting HH moves their shared hydrogen atoms and must recheck OH.
    trial(7)=a*cos(.03_8*sin(real(sample,8)))-b*sin(.03_8*sin(real(sample,8)))
    trial(8)=b*cos(.03_8*sin(real(sample,8)))+a*sin(.03_8*sin(real(sample,8)))
    iterations=shake(reference,trial)
    do ic=1,3
      left=shake_mol(1)%bond(ic)%i; right=shake_mol(1)%bond(ic)%j
      residual=abs(sum((trial(3*left-2:3*left)-trial(3*right-2:3*right))**2)- &
                   shake_mol(1)%bond(ic)%dist2)/shake_mol(1)%bond(ic)%dist2
      if (residual > worst) then
        worst=residual; worst_sample=sample; worst_bond=ic
        ready=merge(1,0,all(shake_mol(1)%bond(:)%ready))
      end if
    end do
  end do
  write(*,'(a,3i8,2es26.17e3)') 'SHAKE_MAX ',worst_sample,worst_bond,ready,worst,shake_tol
end program shake_residual_audit
