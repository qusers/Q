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
!!  Copyright (c) 2017 Johan Aqvist, John Marelius, Shina Caroline Lynn Kamerlin
!!  and Paul Bauer
!!  program: **qfep**  
!!  by Johan Aqvist, Karin Kolmodin, John Marelius, Johan Sund  
!!  qfep free energy analysis program for FEP, EVB & Umbrella Sampling  
!------------------------------------------------------------------------------!
program qfep
  use iso_fortran_env, only : compiler_version, compiler_options
  
  use nrgy
  use parse

  implicit none
  character(*), parameter :: program_name = 'qfep'
  character(*), parameter :: program_version = '5.7'
  character(*), parameter :: program_date = '2015-04-01'
!  character(*), parameter :: options = compiler_options()
  character(len=32)       :: arg
  integer                 :: k
  
  integer,parameter :: mxpts = 20000000
  integer,parameter :: mxbin = 100
  integer,parameter :: mxstates = 4
  character(80)     :: filnam, line
  integer           :: i,j,ifile,ipt,istate,ibin,nfiles,nstates,ERR, &
                       nskip,nbins,nmin,idum,noffd,nnoffd,offel,ngroups, &
                       iexclg,igrp,mpts,iexclgn,astate,bstate

  type(OFFDIAG_SAVE), dimension(mxstates) :: offd

  real(8) :: rt,gapmin,gapmax,sum,dv,gaprange,xint,dvg,veff1,veff2, &
             dGa,dGb,dGg,alpha_B,scale_Hij,veff,min,dlam,sumf,sumb, &
             konst,fel,nfnr,nrnf
  real(8),dimension(mxbin) :: sumg,sumg2,avdvg,avc1,avc2,avc11,avc12,avc13, &
                              avc21,avc22,avc23,avc31,avc32,avc33,avr

  real(8),dimension(mxbin,4) :: binsum
  integer,dimension(mxbin)   :: nbinpts, ptsum

  type(Q_ENERGIES), dimension(mxstates) :: EQ
  type(Q_ENERGIES), dimension(mxstates) :: avEQ
  real(8),dimension(mxstates)           :: dvv,dGv,alfa,coeff,Vel,Vvdw

  real(8),dimension(3) :: u,y
  real(8),allocatable  :: dgf(:),dgr(:),dgfsum(:),dgrsum(:),dG(:),dgti(:), &
                          dgtisum(:),dglu(:),dglusum(:),dgbar(:),dgbarsum(:)
  real(8),dimension(mxstates,mxstates) :: A,mu,eta,rxy0
  real(8),allocatable                  :: Hij(:,:),d(:),e(:)

  type FEP_DATA_TYPE
     integer           :: npts
     real(8)           :: lambda(mxstates)
     real(8), pointer  :: v(:,:), r(:,:) !indices are state, point
     real(8), pointer  :: vg(:), gap(:), c1(:), c2(:) !index is point
  end type FEP_DATA_TYPE

  type(FEP_DATA_TYPE), allocatable :: FEP(:) !index is file
  type(FEP_DATA_TYPE) :: FEPtmp !temporary storage for one file

  real(8),dimension(mxstates,mxbin) :: avdvv,sumv,sumv2 

  integer        :: f, gas=0, error, dummyno ! masoud
  real           :: dummy !masoud
  character(100) :: iline !masoud

  
  call commandlineoptions
  call startup




  !------------------------------------------
  ! INPUT OF PARAMETERS
  !------------------------------------------
  call prompt ('--> Number of energy files: ')
  read (*,*) nfiles
  write (*,1) nfiles
1 format('# Number of files                 =',i6)
  call prompt ('--> No. of states, no. of predefined off-diag elements: ')
  read (*,*) nstates, noffd
  write (*,2) nstates, noffd

2 format('# Number of states                 =',i6,/, &
       '# Number of off-diagonal elements =',i6)
  !allocate size of secular determinant
  allocate(Hij(nstates,nstates),d(nstates),e(nstates),STAT=ERR)
  if(ERR /= 0) then
     write(*,*) 'ERROR: Out of memory when allocating Hij array.'
     stop 'qfep terminated abnormally: Out of memory.'
  end if

  ! Continue to read input
  call prompt ('--> Give kT & num. of pts to skip & calculation mode: ')
  read(*,'(a)') iline

  ! masoud reading an extra option for just calculating QQ energies
  !print*, iline
  do i=1,4
    read (iline,*,iostat=error)(dummy,j=1,i)
    !   print*, error,"***",i,dummy
    if (error .ne. 0) then
      dummyno=i-1
      exit
    endif
  enddo
  if (dummyno .eq. 2) then
    read (iline,*) rt,nskip
  elseif (dummyno .eq. 3) then
    read (iline,*) rt,nskip,gas
  else
    print*, "num. correct arguments for KT, data points to skip or calculation mode"
    stop
  end if
  ! masoud

! read (*,*) rt,nskip
  write (*,3) rt,nskip,gas
3 format('# kT                              =',f6.3,/, &
         '# Number of data points to skip   =',i6,/, '# Only QQ interactions will be considered   = ',i3)


  call prompt ('--> Give number of gap-bins: ')
  read (*,*) nbins
  write (*,4) nbins
4 format('# Number of gap-bins              =',i6)

  call prompt ('--> Give minimum # pts/bin: ')
  read (*,*) nmin
  write (*,5) nmin
5 format('# Minimum number of points per bin=',i6)

  do istate=2,nstates
    write(line,7) istate
    call prompt(line)
    read (*,*) alfa(istate)
    write (*,6) istate,alfa(istate)
  end do
6 format('# Alpha for state ',i2,'              =',f6.2)
7 format('--> Give alpha for state ',i2,':')

  scale_Hij=0.0
  if (noffd /=0) then
    call prompt ('--> Hij scaling:')
    read (*,*) scale_Hij
    write (*,8) scale_Hij
8   format('# Scale factor for Hij            =',f6.2)
  else
    call prompt ('--> Number of off-diagonal elements:')
    read (*,*) nnoffd
    write (*,9) nnoffd
    if(nnoffd >0) write(*,11) 
    do offel=1,nnoffd
      call prompt ('--> i, j, A_ij, mu_ij, eta_ij, r_xy0: ')
      read (*,*) i, j, A(i,j), mu(i,j), eta(i,j), rxy0(i,j)
      write(*,12) i, j, A(i,j), mu(i,j), eta(i,j), rxy0(i,j)
    end do
9   format('# Number of off-diagonal elements =',i6)
11  format('#   i   j   A(i,j)  mu(i,j) eta(i,j) rxy0(i,j)')
12  format('#',2i4,4f9.2) 
  end if

  call prompt ('--> linear combination of states defining reaction coord: ')
  read (*,*) (coeff(i),i=1,nstates)
  write(*,13) coeff(1:nstates)
13 format('# Linear combination coefficients=', 8f8.2)

  !allocate large arrays
  allocate(FEP(nfiles), FEPtmp%v(nstates,mxpts), FEPtmp%r(nstates,mxpts), &
       FEPtmp%vg(mxpts), FEPtmp%gap(mxpts), &
       FEPtmp%c1(mxpts), FEPtmp%c2(mxpts), &
       dgf(0:nfiles+1),dgr(0:nfiles+1),dgti(0:nfiles+1), &
       dglu(0:nfiles+1),dgbar(0:nfiles+1), &
       dgfsum(0:nfiles+1),dgrsum(0:nfiles+1),dG(0:nfiles+1),dgtisum(0:nfiles+1), &
       dglusum(0:nfiles+1),dgbarsum(0:nfiles+1), &
       STAT=ERR)
  if(ERR /= 0) then
    stop 'qfep terminated abnormally: Out of memory when allocating arrays.'
  end if

  !---------------------------------
  ! Energy files are opened and read
  !---------------------------------
  binsum(:,:)=0.
  EQ(:)%total=0.

  gapmin=999.
  gapmax=-999.
  f = freefile()
  FEPtmp%c1(:) = 0.
  FEPtmp%c2(:) = 0.

  write(*,*)
  write(*,*)
  write(*,15)
  write(*,16)
15 format('# Part 0: Average energies for all states in all files')
16 format('# file             state   pts   lambda    EQtot   EQbond',&
       '  EQang   EQtor   EQimp    EQel   EQvdW  Eel_qq  EvdW_qq Eel_qp  EvdW_qp Eel_qw EvdW_qw Eqrstr')

  do ifile=1,nfiles
    write(line,14) ifile
14  format('--> Name of file number',i4,':')
    call prompt(line)
    read(*,*) filnam
    write (*,*) ''
    if(openit(f,filnam,'old','unformatted','read') /= 0) then
      stop 'qfep terminated abnormally: Failed to open energy file.'
    end if

    !read 1st record to get lambdas
    idum = get_ene(f, EQ(:), offd, nstates,nnoffd)
    if(idum /= 0) then 
      !we have a problem
      write(*,'(a,a)') 'ERROR: Unexpected end-of-file in first record of ', &
           trim(filnam)
      if(idum > 0) then
        !the problem is in energies
        write(*,'(a,i1)') 'while reading energies of state ',idum
        write(*,'(a,i1,a)') 'Maybe there are only ',idum-1, ' states?'
        stop 'qfep terminated abnormally: Failed to read energies.'
      else
        !idum < 0 indicates problems with offdiags
        write(*,'(a)') 'while reading off-diagonal elements.'
        write(*,'(a,i1,a)') 'Maybe there are less than ',nnoffd, &
             ' off-diagonal elements?'
        stop 'qfep terminated abnormally: Failed to read off-diagonals.'
      end if
    end if

    FEP(ifile)%lambda(:) = EQ(:)%lambda

    rewind (f)

    ipt = 0
    avEQ(:) = avEQ(:) * 0. !set all fields to zero using multiplication operator
    FEPtmp%gap(:) = 0.
    do while(get_ene(f, EQ(:), offd, nstates, nnoffd) == 0) !keep reading till EOF
      ipt = ipt + 1
!   masoud
!if the gas flag is > 0 then total energy will be just q (bonded), qq(nonbonded) and restraint
      if (gas .gt. 0 ) then
        do i=1,nstates
          EQ(i)%total=EQ(i)%q%bond+EQ(i)%q%angle+EQ(i)%q%torsion+EQ(i)%q%improper+EQ(i)%qq%el+EQ(i)%qq%vdw+EQ(i)%restraint
!         print*, EQ(i)%total, EQ(i)%q,EQ(i)%qq,EQ(i)%restraint
!         print*, "''''''''''"
!         if (ipt > 10 ) stop
        end do
      end if
!   masoud
    if(ipt > nskip) then
      avEQ(:) = avEQ(:) + EQ(:) !use Qenergies + operator
    end if


      !-------------------------------------------
      ! Correct H_ii with alfa, and modify H_ij...
      !-------------------------------------------
      alfa(1)=0.
      EQ(:)%total=EQ(:)%total+alfa(:)
      FEPtmp%v(1:nstates, ipt) = EQ(1:nstates)%total
      if (nnoffd .ne. 0) then
        FEPtmp%r(:,ipt) = offd(:)%rkl
      end if

      do i=1, noffd
        Hij(offd(i)%i, offd(i)%j) = offd(i)%Hij
      end do

      if ( scale_Hij .gt. 0.0 ) then
        Hij(1,2) = scale_Hij*Hij(1,2)
      else
        do i=1,nstates
          do j=1,nstates
            if (i==j) then
              Hij(i,j)=EQ(i)%total
            else
              if (A(i,j)==0.0) then
                Hij(i,j) = A(j,i)*exp(-mu(j,i)*(offd(1)%rkl-rxy0(j,i)))* &
                exp(-eta(j,i)*(offd(1)%rkl-rxy0(j,i))**2)
              else 
                Hij(i,j) = A(i,j)*exp(-mu(i,j)*(offd(1)%rkl-rxy0(i,j)))* &
                exp(-eta(i,j)*(offd(1)%rkl-rxy0(i,j))**2)
              end if
            end if
          end do
        end do
      end if

      !-----------------------------------------------------------
      ! Ground state energy is calculated from secular determinant
      !-----------------------------------------------------------
      if (nstates==2) then
        FEPtmp%vg(ipt)=0.5*(EQ(1)%total+EQ(2)%total)-  &
        0.5*sqrt( (EQ(1)%total-EQ(2)%total)**2 + 4.*Hij(1,2)**2 )
        if(nnoffd > 0) then
          FEPtmp%c1(ipt)=1./(1.+((FEPtmp%vg(ipt)-EQ(1)%total)/Hij(1,2))**2)
          FEPtmp%c2(ipt)=1-FEPtmp%c1(ipt)
        end if
        else 
          call tred2(Hij,nstates,nstates,d,e)
          call tqli(d,e,nstates,nstates,Hij)
          FEPtmp%vg(ipt)=MINVAL(d)
      end if

      do istate=1,nstates
        FEPtmp%gap(ipt)=FEPtmp%gap(ipt)+FEPtmp%v(istate,ipt)*coeff(istate)
!        print*, "LINEAR COMBINATION"
!        print*, "FEPtmp%gap(ipt)        coeff(istate)                istate           nstates"        
!        print*, FEPtmp%gap(ipt) ,coeff(istate), istate, nstates
      end do

      if(ipt .gt. nskip) then
        if(FEPtmp%gap(ipt) .lt. gapmin) gapmin=FEPtmp%gap(ipt)
        if(FEPtmp%gap(ipt) .gt. gapmax) gapmax=FEPtmp%gap(ipt)
      end if
  end do  !(ipt)


    close(f)
    FEP(ifile)%npts = ipt



     !copy FEPtmp to D
     allocate(FEP(ifile)%v(nstates, FEP(ifile)%npts), &
          FEP(ifile)%vg(FEP(ifile)%npts), &
          FEP(ifile)%gap(FEP(ifile)%npts), &
          FEP(ifile)%c1(FEP(ifile)%npts), &
          FEP(ifile)%c2(FEP(ifile)%npts), &
          STAT=ERR)
     if(ERR /= 0) then
        stop 'qfep terminated abnormally: Out of memory when allocating arrays.'
     end if
     FEP(ifile)%v(:,:) = FEPtmp%v(1:nstates, 1:FEP(ifile)%npts)
     FEP(ifile)%vg(:) = FEPtmp%vg(1:FEP(ifile)%npts)
     FEP(ifile)%gap(:) = FEPtmp%gap(1:FEP(ifile)%npts)
     FEP(ifile)%c1(:) = FEPtmp%c1(1:FEP(ifile)%npts)
     FEP(ifile)%c2(:) = FEPtmp%c2(1:FEP(ifile)%npts)
     if(nnoffd > 0) then 
        allocate(FEP(ifile)%r(nstates, FEP(ifile)%npts), STAT=ERR)
        if(ERR /= 0) then
           stop 'qfep terminated abnormally: Out of memory when allocating arrays.'
        end if
        FEP(ifile)%r(:,:) = FEPtmp%r(1:nstates, 1:FEP(ifile)%npts)
     end if


     ! Average energies for each file are calculated
     if(ipt <= nskip) then !skipped all points in the file
        write(*,900) trim(filnam)
900     format('>>>>> ERROR: number of data sets in ',a,&
             ' is less than number of points to skip!')
        stop 'qfep terminated abnormally: Too few data points in file.'
     else
        avEQ(:) = avEQ(:) * (1./(FEP(ifile)%npts-nskip)) !use new * operator
     end if
     do istate=1,nstates
        write(*,17) trim(filnam), istate, FEP(ifile)%npts-nskip, FEP(ifile)%lambda(istate), &
             avEQ(istate)%total,avEQ(istate)%q%bond,avEQ(istate)%q%angle,&
             avEQ(istate)%q%torsion, &
             avEQ(istate)%q%improper,avEQ(istate)%qx%el,avEQ(istate)%qx%vdw, &
             avEQ(istate)%qq%el,avEQ(istate)%qq%vdw,avEQ(istate)%qp%el,&
             avEQ(istate)%qp%vdw, &
             avEQ(istate)%qw%el,avEQ(istate)%qw%vdw,avEQ(istate)%restraint

     end do
17   format(a,t23,i2,1x,i6,1x,f8.6,14f8.2)
  end do  !ifile

  if(nfiles > 1) then !the following is meaningless for a single file
     dgf=0.
     dgfsum=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     do ifile=1,nfiles-1
        !check that number of points > number to skip
        if(nskip >= FEP(ifile)%npts) then
           write(*,999) ifile, nskip, FEP(ifile)%npts
999        format('File',i5,' contains only',i5,' points. Can''t skip',i5)
        end if
        do ipt=nskip+1,FEP(ifile)%npts
           do istate=1,nstates
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
              veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do
           dv=veff2-veff1
           veff1=0.
           veff2=0.
           sum=sum+exp(-dv/rt)
        end do
        sum=sum/real(FEP(ifile)%npts-nskip)
        dgf(ifile)=-rt*dlog(sum)
        dgfsum(ifile+1)=dgfsum(ifile)+dgf(ifile)
        sum=0.
     end do
     dgrsum=0.
     dgr=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     do ifile=nfiles,2,-1
        do ipt=nskip+1,FEP(ifile)%npts
           do istate=1,nstates
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
              veff2=veff2+FEP(ifile-1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do
           dv=veff2-veff1
           veff1=0.
           veff2=0.
           sum=sum+exp(-dv/rt)
        end do
        sum=sum/real(FEP(ifile)%npts-nskip)
        dgr(ifile)=-rt*dlog(sum)
        dgrsum(ifile-1)=dgrsum(ifile)+dgr(ifile)
        sum=0.
     end do
     dgtisum=0.
     dgti=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     astate=1
     bstate=2 
     istate=1
     do ifile=1,nfiles                        
        do ipt=nskip+1,FEP(ifile)%npts
           veff1=FEP(ifile)%v(astate,ipt)
           veff2=FEP(ifile)%v(bstate,ipt)
           dv=veff2-veff1
           veff1=0.
           veff2=0.
           sum=sum+dv
        end do
        dgti(ifile)=sum/real(FEP(ifile)%npts-nskip)
        sum=0.
     end do
     do ifile=2,nfiles
        dlam=FEP(ifile-1)%lambda(istate)-FEP(ifile)%lambda(istate)
        dgtisum(ifile)=dgtisum(ifile-1)+(0.5*(dgti(ifile-1)+dgti(ifile))*dlam)
     end do
     dglu=0.
     dglusum=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     do ifile=1,nfiles-1
        do ipt=nskip+1,FEP(ifile)%npts
           do istate=1,nstates
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
              veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do
           dv=(veff2-veff1)/2
           veff1=0.
           veff2=0.
           sum=sum+exp(-dv/rt)
        end do
        sumf=sum/real(FEP(ifile)%npts-nskip)
        sum=0.
        do ipt=nskip+1,FEP(ifile+1)%npts
           do istate=1,nstates
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
              veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
           end do
           dv=(veff2-veff1)/2
           veff1=0.
           veff2=0.
           sum=sum+exp(dv/rt)
        end do
        sumb=sum/real(FEP(ifile+1)%npts-nskip)
        dglu(ifile)=-rt*dlog(sumf/sumb)
        dglusum(ifile+1)=dglusum(ifile)+dglu(ifile)
        sum=0.
     end do
     dgbar=0.
     dgbarsum=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     konst=0
     fel=1
     do ifile=1,nfiles-1
        konst=dglu(ifile)
        do while (fel>0.001)  
           nfnr=real(FEP(ifile)%npts-nskip)/real(FEP(ifile+1)%npts-nskip)
           nrnf=real(FEP(ifile+1)%npts-nskip)/real(FEP(ifile)%npts-nskip)
           do ipt=nskip+1,FEP(ifile)%npts
              do istate=1,nstates
                 veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
                 veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
              end do
              dv=(veff2-veff1)
              veff1=0.
              veff2=0.
              sum=sum+1/((1+(nfnr*exp((dv-konst)/rt))))
           end do
           sumf=sum/real(FEP(ifile)%npts-nskip)
           sum=0.
           do ipt=nskip+1,FEP(ifile+1)%npts
              do istate=1,nstates
                 veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
                 veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
              end do
              dv=(veff2-veff1)
              veff1=0.
              veff2=0.
              sum=sum+1/((1+(nrnf*exp((-dv+konst)/rt))))
           end do
           sumb=sum/real(FEP(ifile+1)%npts-nskip)
           dgbar(ifile)=-rt*dlog((sumf/sumb)*exp(-konst/rt)*nfnr)
           sum=0.
           fel=ABS(konst-dgbar(ifile))
           konst=dgbar(ifile)
        end do
        dgbarsum(ifile+1)=dgbarsum(ifile)+dgbar(ifile)
        fel=1
     end do


     write(*,*) 
     write(*,*) 
     write(*,21)
     write(*,22)
21   format('# Part 1: Free energy perturbation summary:')
22   format('# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>')
     dG(1)=0.0
     do ifile=2,nfiles
        dG(ifile)=dG(ifile-1)+0.5*(dgf(ifile-1)-dgr(ifile))
     end do
     do ifile=1,nfiles
        write (*,23) &
             FEP(ifile)%lambda(1),dgf(ifile-1),dgfsum(ifile), &
             dgr(ifile+1),dgrsum(ifile),dG(ifile)
     end do
23   format(3x,f8.6,5f9.3)


     write (*,*)
     write (*,'(a,f9.2)') '# Min energy-gap is: ',gapmin
     write (*,'(a,f9.2)') '# Max energy-gap is: ',gapmax

     !-----------------------------------
     ! Reaction free energy is calculated
     !-----------------------------------
     write (*,*)
     write (*,*)
     write(*,24)
     write(*,25)
24   format('# Part 2: Reaction free energy summary:')
25   format('# Lambda(1)  bin Energy gap      dGa     dGb     dGg    # pts    c1**2    c2**2')
26   format(3x,f8.6,i5,2x,4f9.2,2x,i5,2f9.3)
     gaprange=gapmax-gapmin      !Range of reaction coordinate
     xint=gaprange/real(nbins)   !Divide R.C. into bins
     do ifile=1,nfiles         
        avdvv=0.
        avdvg=0.
        sumv=0.
        sumg=0.
        avc1=0.
        avc2=0.
        avr=0.
        dvv=0.
        dvg=0.
        nbinpts=0

        do ipt=nskip+1,FEP(ifile)%npts
           ibin=int((FEP(ifile)%gap(ipt)-gapmin)/xint)+1
           veff=0.
           do istate=1,nstates
              veff=veff+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do  !states
           dvv(1:nstates)=FEP(ifile)%v(1:nstates,ipt)-veff
           dvg=FEP(ifile)%vg(ipt)-veff
           avdvv(:,ibin)=avdvv(:,ibin)+dvv(:)
           avdvg(ibin)=avdvg(ibin)+dvg
           avc1(ibin)=avc1(ibin)+FEP(ifile)%c1(ipt)
           avc2(ibin)=avc2(ibin)+FEP(ifile)%c2(ipt)
           !Only gives first r_xy distance
           if(nnoffd > 0)  avr(ibin)=avr(ibin)+FEP(ifile)%r(1,ipt)
           nbinpts(ibin)=nbinpts(ibin)+1
        end do           !ipt
        do ibin=1,nbins
           if ( nbinpts(ibin) .ne. 0 ) then
              avc1(ibin)=avc1(ibin)/real(nbinpts(ibin))
              avc2(ibin)=avc2(ibin)/real(nbinpts(ibin))
              avr(ibin)=avr(ibin)/real(nbinpts(ibin))
              avdvv(:,ibin)=avdvv(:,ibin)/nbinpts(ibin)
              avdvg(ibin)=avdvg(ibin)/nbinpts(ibin)

           end if
        end do !ibin
        do ipt=nskip+1,FEP(ifile)%npts
           ibin=int((FEP(ifile)%gap(ipt)-gapmin)/xint)+1
           veff=0.
           do istate=1,nstates
              veff=veff+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do  !istate

           do istate=1,nstates
              dvv(istate)=FEP(ifile)%v(istate,ipt)-veff-avdvv(istate,ibin)
           end do
           dvg=FEP(ifile)%vg(ipt)-veff-avdvg(ibin)
           sumv(:,ibin)=sumv(:,ibin)+exp(-dvv(:)/rt)
           sumg(ibin)=sumg(ibin)+exp(-dvg/rt)
        end do   !ipt

        do ibin=1,nbins
           if (nbinpts(ibin).ge.nmin) then
              binsum(ibin,2)=binsum(ibin,2)+avc1(ibin)*nbinpts(ibin)
              binsum(ibin,3)=binsum(ibin,3)+avc2(ibin)*nbinpts(ibin)
              binsum(ibin,4)=binsum(ibin,4)+avr(ibin)*nbinpts(ibin)  !Bin-averaged r_xy
              sumv(:,ibin)=sumv(:,ibin)/real(nbinpts(ibin)) 
              sumg(ibin)=sumg(ibin)/real(nbinpts(ibin))

              ptsum(ibin)=ptsum(ibin)+nbinpts(ibin)

              do istate=1,nstates
                 sumv2(istate,ibin)=-rt*dlog(sumv(istate,ibin))+avdvv(istate,ibin)
              end do
              sumg2(ibin)=-rt*dlog(sumg(ibin))+avdvg(ibin) 
              ! These are the diabatic free energy curves
              dGv(:)=dG(ifile)+sumv2(:,ibin) 
              ! This is the reaction free energy
              dGg=dG(ifile)+sumg2(ibin)

              binsum(ibin,1)=binsum(ibin,1)+dGg*int(nbinpts(ibin))

              write (*,26) FEP(ifile)%lambda(1),ibin, &
                   gapmin+real(ibin)*xint-xint/2., dGv(1),dGv(2),dGg, &
                   int(nbinpts(ibin)),avc1(ibin),avc2(ibin) 
           end if
        end do  !ibin
     end do      !ifile
     write(*,*)
     write(*,*)
     write(*,27)
     write(*,28)

27   format('# Part 3: Bin-averaged summary:')
28   format('# bin  energy gap  <dGg> <dGg norm> pts  <c1**2> <c2**2> <r_xy>')

     do ibin=1,nbins
        if (ptsum(ibin).ge.nmin) then
           binsum(ibin,1)=binsum(ibin,1)/real(ptsum(ibin)) ! Bin-averaged reaction free energy
           binsum(ibin,2)=binsum(ibin,2)/real(ptsum(ibin)) ! Bin-averaged c1**2
           binsum(ibin,3)=binsum(ibin,3)/real(ptsum(ibin)) ! Bin-averaged c2**2
           binsum(ibin,4)=binsum(ibin,4)/real(ptsum(ibin)) ! Bin-averaged r_xy
        end if
     end do
     min=MINVAL(binsum(:,1))
     do ibin=1,nbins
        if (ptsum(ibin).ge.nmin) then
29         format(i4,1x,3f9.2,2x,i5,3f8.3,4f8.2)
           write(*,29) ibin,gapmin+real(ibin)*xint-xint/2.,binsum(ibin,1),  &
                binsum(ibin,1)-min,int(ptsum(ibin)),binsum(ibin,2),binsum(ibin,3),binsum(ibin,4)
        end if
     end do !ibin
  end if !nfiles >1


     write(*,*) 
     write(*,*) 
     write(*,30)
     write(*,31)
30   format('# Part 4: Termodynamic integration:')
31   format('# lambda(1)      dGti    sum(dGti) ')
     do ifile=1,nfiles
        write (*,23) &
             FEP(ifile)%lambda(1),dgti(ifile-1),dgtisum(ifile)
     end do


     write(*,*) 
     write(*,*) 
     write(*,32)
     write(*,33)
32   format('# Part 5: Overlap sampling Lu et al:')
33   format('# lambda(1)      dG    sum(dG) ')
     do ifile=1,nfiles
        write (*,23) &
             FEP(ifile)%lambda(1),dglu(ifile-1),dglusum(ifile)
     end do


     write(*,*) 
     write(*,*) 
     write(*,33)
     write(*,34)
34   format('# Part 6: BAR Bennet:')
35   format('# lambda(1)      dG    sum(dG) ')
     do ifile=1,nfiles
        write (*,23) &
             FEP(ifile)%lambda(1),dgbar(ifile-1),dgbarsum(ifile)
     end do



  !clean up
  do ifile=1,nfiles
     deallocate(FEP(ifile)%v, FEP(ifile)%vg, FEP(ifile)%gap, &
          FEP(ifile)%c1, FEP(ifile)%c2)
     if(nnoffd > 0) deallocate(FEP(ifile)%r)
  end do
  deallocate(FEP)

  deallocate(Hij,d,e,STAT=ERR)




contains

  !----------------------------------------------------------------------------!
  !!  subroutine: startup  
  !!  Startup  
  !----------------------------------------------------------------------------!
subroutine startup
  integer :: i

  print '(a)',  '--------------------------------------------------------------------------------'
  print '(4a)', 'Welcome to ', program_name, ' version: ', program_version
  print '(a)',  ' '
  print '(2a)', 'This version was compiled using: ', compiler_version()
  print '(a)',  ' '
  print '(a)',  'And using the following compiler options: '
!  write (output_unit, *, delim='quote') options
  write ( *, '( A /)' ) trim ( compiler_options() )
!  print '(a)',  trim(compiler_options())
  print '(a)',  ' '
  print '(a)',  'For command line options type qfep --help  or qcalc -h at the terminal.'
  print '(a)',  ' '
  print '(a)',  'If you are using the interactive mode and want to quit type "quit" or Ctrl-C.' 
  print '(a)',  '--------------------------------------------------------------------------------'

end subroutine startup


  !----------------------------------------------------------------------------!
  !!  subroutine: prompt
  !----------------------------------------------------------------------------!
  subroutine prompt (outtxt)
    character(*) outtxt
#if defined (__osf__)
    !prompt to STDERR using unit 5=STDIN on OSF/1=DEC UNIX
    integer, parameter                      ::      f=5
    !write (f,'($,a)') outtxt
#elseif defined (_WIN32)
    !open the ERR file on Win32
    integer, save :: f
    if(f==0) then
       f=17
       open(unit=f,file='ERR')
    end if
    !write (f1,'($,a)') outtxt
#else
    !otherwise prompt to STDOUT
    integer, parameter                  ::      f=6
    !write (f2,'($,a)') outtxt
#endif
    write (f,'($,a)') outtxt
  end subroutine prompt

  !----------------------------------------------------------------------------!
  !!  subroutine: commandlineoptions  
  !----------------------------------------------------------------------------!
  subroutine commandlineoptions
  do k = 1, command_argument_count()
    call get_command_argument(k, arg)
    select case (arg)
    case ('-v', '--version')
      print '(3a)', program_name, ' version ', program_version
      stop
    case ('-h', '--help')
      call print_help()
      stop
    case default
      print '(a,a,/)', 'Unrecognized command-line option: ', arg
      call print_help()
      stop
    end select
  end do
  end subroutine commandlineoptions

  !----------------------------------------------------------------------------!
  !!  subroutine: print_help  
  !----------------------------------------------------------------------------!
  subroutine print_help()
    print '(a)', 'usage:'
    print '(a)', 'qfep [OPTION]'
    print '(a)', '  or'
    print '(a)', 'qfep < inputfile.inp > outputfile.out'
    print '(a)', ''
    print '(a)', 'Without options, qfep goes into interactive mode.'
    print '(a)', ''
    print '(a)', 'qfep [OPTION]:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
  end subroutine print_help

  
  subroutine tred2(A,N,NP,D,E)
    !------------------------------------------------------------
    ! This subroutine reduces a symmetric matrix to tridiagonal
    ! form. The tridiagonal matrix can further be diagonalized by
    ! the subroutine tqli.
    ! These subroutines were copied by Karin Kolmodin 20 Nov. 1997
    ! from http://rsc.anu.au/HWS/COURSES/MATHMETH/node70.html
    ! and rewritten in f90.
    !------------------------------------------------------------ 
    real(8),dimension(:)   :: D,E
    real(8),dimension(:,:) :: A
    integer                :: NP,I,J,K,L,N,ERR
    real(8)                :: SCALE,F,G,H,HH

    IF(N.GT.1)THEN
       DO I=N,2,-1  
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
             DO K=1,L
                SCALE=SCALE+ABS(A(I,K))
             end do
             IF(SCALE.EQ.0.)THEN
                E(I)=A(I,L)
             ELSE
                DO K=1,L
                   A(I,K)=A(I,K)/SCALE
                   H=H+A(I,K)**2
                end do
                F=A(I,L)
                G=-SIGN(SQRT(H),F)
                E(I)=SCALE*G
                H=H-F*G
                A(I,L)=F-G
                F=0.
                DO J=1,L
                   A(J,I)=A(I,J)/H
                   G=0.
                   DO K=1,J
                      G=G+A(J,K)*A(I,K)
                   end do
                   IF(L.GT.J)THEN
                      DO K=J+1,L
                         G=G+A(K,J)*A(I,K)
                      end do
                   END IF
                   E(J)=G/H
                   F=F+E(J)*A(I,J)
                end do
                HH=F/(H+H)
                DO J=1,L
                   F=A(I,J)
                   G=E(J)-HH*F
                   E(J)=G
                   DO K=1,J
                      A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                   end do
                end do
             END IF
          ELSE
             E(I)=A(I,L)
          END IF
          D(I)=H
       end do
    END IF
    D(1)=0.
    E(1)=0.
    DO I=1,N
       L=I-1
       IF(D(I).NE.0.)THEN
          DO J=1,L
             G=0.
             DO K=1,L
                G=G+A(I,K)*A(K,J)
             end do
             DO K=1,L
                A(K,J)=A(K,J)-G*A(K,I)
             end do
          end do
       END IF
       D(I)=A(I,I)
       A(I,I)=1.
       IF(L.GE.1)THEN
          DO J=1,L
             A(I,J)=0.
             A(J,I)=0.
          end do
       END IF
    end do
  end subroutine tred2
  !-----------------------------------

  subroutine tqli(D,E,N,NP,Z)
    !------------------------------------------------------------ 
    ! This subroutine diagonalizes a tridiagonal matrix which has 
    ! been prepared by the subroutine tred2.
    ! These subroutines were copied by Karin Kolmodin 20 Nov. 1997
    ! from http://rsc.anu.au/HWS/COURSES/MATHMETH/node70.html
    ! and rewritten in f90
    !------------------------------------------------------------  
    implicit none
    real(8),dimension(:)   :: D,E
    real(8),dimension(:,:) :: Z
    integer                :: I,N,NP,K,L,M,ITER
    real(8)                :: DD,G,R,S,C,P,F,B

    IF (N.GT.1) THEN
       DO I=2,N
          E(I-1)=E(I)
       end do
       E(N)=0.
       DO L=1,N
          ITER=0
1         DO M=L,N-1
             DD=ABS(D(M))+ABS(D(M+1))
             IF (ABS(E(M))+DD.EQ.DD) GO TO 2
          end do
          M=N
2         IF(M.NE.L)THEN
             IF(ITER.EQ.30)STOP 'too many iterations'
             ITER=ITER+1
             G=(D(L+1)-D(L))/(2.*E(L))
             R=SQRT(G**2+1.)
             G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
             S=1.
             C=1.
             P=0.
             DO I=M-1,L,-1
                F=S*E(I)
                B=C*E(I)
                IF(ABS(F).GE.ABS(G))THEN
                   C=G/F
                   R=SQRT(C**2+1.)
                   E(I+1)=F*R
                   S=1./R
                   C=C*S
                ELSE
                   S=F/G
                   R=SQRT(S**2+1.)
                   E(I+1)=G*R
                   C=1./R  
                   S=S*C
                END IF
                G=D(I+1)-P
                R=(D(I)-G)*S+2.*C*B
                P=S*R
                D(I+1)=G+P
                G=C*R-B
                DO K=1,N
                   F=Z(K,I+1)
                   Z(K,I+1)=S*Z(K,I)+C*F
                   Z(K,I)=C*Z(K,I)-S*F
                end do
             end do
             D(L)=D(L)-P
             E(L)=G
             E(M)=0.
             GO TO 1
          END IF
       end do
    END IF
    RETURN
  end subroutine tqli
  !-----------------------------------

end program qfep
