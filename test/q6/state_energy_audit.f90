! Fixed-coordinate audit of actual Q state energies/gradients. No dynamics or sampler.
program state_energy_audit
  use md
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(8), parameter :: weights(7)=[1.0_8,0.0_8,.25_8,.5_8,.75_8,.0001_8,.9999_8]
  real(8) :: angular(2),born_terms(2)
  integer :: mode,setting,unit
  nodeid=0; numnodes=1
  call md_startup
  if (.not. initialize()) call die('Invalid state audit input')
  if (use_PBC .or. use_LRF .or. wpol_adapt) call die('Audit needs spherical direct interactions and frozen offsets')
  if (.not. perstate_wpol .or. .not. perstate_born) call die('Audit input must enable both state terms')
  call open_files
  call topology
  call prep_coord
  call get_fep
  if (nstates /= 2 .or. noffd /= 0) call die('Audit needs two uncoupled states')
  call prep_sim
  call close_input_files
  call init_shake
  call make_nbqqlist
  call shrink_topology
  call nbmonitorlist
  call distribute_nonbonds
  ! A common nonzero test offset checks bookkeeping away from the zero-offset special case.
  wshell(:)%theta_corr=0.03
  ! A fixed snapshot of the existing wall's temperature argument, not a wall modification.
  Tfree=Temp0; istep=1
  write(*,'(a,7es26.17e3)') 'AUDIT_META ',rwat,coulomb_constant,born_eps,born_C,born_crg_env,q_region_state
  open(newunit=unit,file='gradients.bin',status='new',form='unformatted')
  write(unit) natom
  do mode=1,4
    ! Explicit diagnostic subtraction controls, NOT proposed production targets.
    ! 1 both terms, 2 polarization only, 3 Born only, 4 neither term.
    perstate_born=(mode == 1 .or. mode == 3)
    perstate_wpol=(mode == 1 .or. mode == 2)
    wpol_restr=perstate_wpol
    do setting=1,size(weights)
      EQ(1)%lambda=weights(setting); EQ(2)%lambda=1-weights(setting)
      call make_pair_lists
      call pot_energy
      if (.not. ieee_is_finite(E%potential) .or. .not. all(ieee_is_finite(d))) call die('Nonfinite audit state')
      if (.not. all(ieee_is_finite(EQ(1:2)%total))) call die('Nonfinite pure-state energy')
      angular=0; born_terms=0
      if (perstate_wpol) angular=wpol_state_energy
      if (perstate_born) born_terms=born_self_state
      write(*,'(a,2i6,11es26.17e3)') 'STATE_AUDIT ',mode,setting,weights(setting),E%potential, &
        EQ(1:2)%total,EQ(1:2)%restraint,angular,born_terms,E%restraint%water_pol
      write(*,'(a,2i6,8es26.17e3)') 'AUDIT_COMPONENT ',mode,setting, &
        E%p%bond+E%w%bond+E%q%bond+E%p%angle+E%w%angle+E%q%angle+ &
        E%p%torsion+E%q%torsion+E%p%improper+E%q%improper, &
        E%pp%el+E%pw%el+E%ww%el+E%qx%el,E%pp%vdw+E%pw%vdw+E%ww%vdw+E%qx%vdw, &
        E%restraint%fix+E%restraint%shell+E%restraint%protein,E%restraint%solvent_radial, &
        E%restraint%water_pol,dot_product(EQ(1:2)%lambda,born_terms),E%LRF
      write(unit) mode,setting,d
      ! Actual Q serialization, in the same mode/setting order as the log.
      call put_ene(11,EQ,OFFD)
    end do
  end do
  call partition_trace
  close(unit)
  call close_output_files
contains
  subroutine partition_trace
    ! Fixture-only analytical trace of the unchanged sodium (topology atom 13).
    ! No dynamics or production interaction changes. Charges are already scaled.
    integer :: j, s, ip, iq, i, jq, axis
    real(8) :: r2, inv6, omitted_h_lj, qi, qj, rounding, saved_x, plus, minus, native_gradient(3)
    real(8), parameter :: step=1.e-5_8
    if (natom /= 3418 .or. nat_solute /= 13 .or. ivdw_rule /= VDW_GEOMETRIC) &
      call die('Partition trace requires the Na/benzene/water fixture')
    omitted_h_lj=0
    do j=nat_solute+1,natom
      if (mod(j-nat_solute,3) == 1 .or. excl(j)) cycle
      r2=sum((x(37:39)-x(3*j-2:3*j))**2)
      inv6=1/r2**3
      omitted_h_lj=omitted_h_lj+iaclib(iac(13))%avdw(1)*iaclib(iac(j))%avdw(1)*inv6**2 &
        -iaclib(iac(13))%bvdw(1)*iaclib(iac(j))%bvdw(1)*inv6
    end do
    write(*,'(a,es26.17e3)') 'PARTITION_H_LJ ',omitted_h_lj
    do s=1,nstates
      rounding=0
      do ip=1,nbqq_pair(s)
        iq=nbqq(ip,s)%iq; i=iqseq(iq); j=nbqq(ip,s)%j; jq=nbqq(ip,s)%jq
        if (i /= 13 .and. j /= 13) cycle
        qi=qcrg(iq,s); qj=crg(j)
        if (jq /= 0) qj=qcrg(jq,s)
        r2=sum((x(3*i-2:3*i)-x(3*j-2:3*j))**2)
        rounding=rounding+(qi*qj-real(real(qi,4)*real(qj,4),8))/sqrt(r2)
      end do
      write(*,'(a,i6,es26.17e3)') 'PARTITION_COULOMB_ROUNDING ',s,rounding
    end do
    ! Last audit setting has both boundary terms disabled. Probe the sodium's
    ! actual Cartesian gradient without changing pair membership (99 A cutoffs).
    native_gradient=d(37:39)
    do axis=1,3
      i=36+axis; saved_x=x(i)
      x(i)=saved_x+step
      call pot_energy
      plus=E%potential
      x(i)=saved_x-step
      call pot_energy
      minus=E%potential
      x(i)=saved_x
      write(*,'(a,i6,2es26.17e3)') 'PARTITION_FORCE_FD ',axis,native_gradient(axis),(plus-minus)/(2*step)
    end do
  end subroutine partition_trace
end program state_energy_audit
