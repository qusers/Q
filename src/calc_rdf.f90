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
!  calc_rdf.f90
!  by Martin Ander and Martin Almlof
!  (multiple atoms in first mask and PBC)
!  Calculates RDF for a center atom with respect to all atoms in
!  target_mask, using Nbins bins from r=0 to r=rdf_radius
!------------------------------------------------------------------------------!
module calc_rdf
  use calc_base
  use maskmanip
  use parse
  implicit none

  ! constants
  integer, parameter               :: max_masks = 10
  real(8)                          :: pi

  ! module variables
  integer, private                 :: ats1,ats2,Nmasks = 0 ! # target masks

  type(mask_type), private, target :: temp_masks(2) ! masks for making pair lists

  type rdf_calc_type
    integer                       :: atoms_in_mask1 ! center atom idx
    integer                       :: frames = 0     ! # frames used
    integer                       :: Nbins          ! Divide rdf_radius into Nbins bins
    real(8)                       :: rdf_radius     ! Calculate RDF within this radius
    real(8),pointer               :: bins(:)        ! RDF bins
  end type rdf_calc_type
  type(rdf_calc_type), private, save :: rdf_calcs(max_masks) ! stores data for indivdual RDF calcs

  type rdf_list_type
    integer                       :: number_of_pairs
    integer, pointer              :: atom1(:), atom2(:)
  end type rdf_list_type
  type(rdf_list_type), private, allocatable :: rdf_listan(:)


contains

!------------------------------------------------------------------------------!
!>  subroutine: rdf_initialize
!>
!------------------------------------------------------------------------------!
  subroutine rdf_initialize
    allocate(rdf_listan(max_masks))
  end subroutine rdf_initialize


!------------------------------------------------------------------------------!
!>  subroutine: rdf_finalize
!>  where the rdf is finally computed and printed to
!>  standard output and to the rdf.dat file.
!------------------------------------------------------------------------------!
  subroutine rdf_finalize(i)
    ! arguments
    integer                       :: i        ! calculation #
    ! locals
    integer                       :: j        ! bin number
    real(8)                       :: shellvol ! shell volume
    real(8)                       :: totvol   ! total sphere volume
    real(kind=dp)                 :: shellrad
    real(kind=dp)                 :: normalize
    real(kind=dp)                 :: bulknormalize
!    pi = 4. * atan(1.)                       ! pi = 3.1416...
    pi = 4.d0 * datan(1.d0)

    ! write output header 1
    write (*,3) i
3   format('Results from RDF calculation' , i3 , ', normalized with respect to shell volume and frame count')

    ! write output header 2
    write (*,4)
4   format('    Bin  r_max    RDF')

    open(10, file='rdf.dat')
    ! normalize the contents of each bin with respect to shell
    ! volume and # processed frames and print the results
    ! The volume of the shell is given by:
    !  4/3 * pi * (r + dr)^3 - 4/3 * pi * r^3 or
    !  sometimes approximated by 4*pi*r^2*dr

    do j = 1, rdf_calcs(i)%Nbins
      shellrad = j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins
      shellvol = 4 * pi / 3 * ((j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3 - &
      ((j-1) * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3)
      normalize = (shellvol / rdf_calcs(i)%frames /  rdf_calcs(i)%atoms_in_mask1)
      !volume = 4 * pi * ((j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**2
      totvol = 4/3 * pi * (j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3
      bulknormalize = (rdf_calcs(i)%bins(rdf_calcs(i)%Nbins) / &
      (4 * pi / 3 * ((rdf_calcs(i)%Nbins * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3 - &
      ((rdf_calcs(i)%Nbins-1) * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins)**3))&
      / rdf_calcs(i)%frames / rdf_calcs(i)%atoms_in_mask1)
      rdf_calcs(i)%bins(j) = (rdf_calcs(i)%bins(j) / shellvol / &
      rdf_calcs(i)%frames / rdf_calcs(i)%atoms_in_mask1)  * (1/bulknormalize)
      !rdf_calcs(i)%bins(j) = rdf_calcs(i)%bins(j) / normalize
      !write (*,5) j , j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins , rdf_calcs(i)%bins(j)
      write(*,100) j , j * rdf_calcs(i)%rdf_radius / rdf_calcs(i)%Nbins , rdf_calcs(i)%bins(j)
      write(10, 100, advance='yes') j , shellrad , rdf_calcs(i)%bins(j), shellvol, totvol, bulknormalize
    end do
!5   format(i7 , f10.6 , f10.6)
100 format(i7 , f10.6 , ' ',  f10.6, ' ',  f14.6, ' ', f14.6, ' ', f14.6 )

    ! finalize and deallocate
    deallocate(rdf_calcs(i)%bins)                   ! deallocate bins
    deallocate(rdf_listan(i)%atom1)
    deallocate(rdf_listan(i)%atom2)
  end subroutine rdf_finalize


!------------------------------------------------------------------------------!
!>  function: rdf_add
!>  qcalc calls this function to set up one RDF calculation,
!>  which is then added to the list of calculations in qcalc
!------------------------------------------------------------------------------!
  integer function rdf_add(desc)
    ! arguments
    character(*)                   :: desc

    ! locals
    integer                        :: ats

    ! check the # of calculations
    if(Nmasks == max_masks) then
      write(*,10) max_masks
      return
    end if
10  format('Sorry, the maximum number of RDF calculations is ',i2)

    ! start adding new calculation
    Nmasks = Nmasks + 1               ! set current calculation
    write(*, 1011)
1011 format('Enter mask for first atom set in RDF calculation')
    call mask_initialize(temp_masks(1))
    ats1 =  maskmanip_make(temp_masks(1))
    if(ats1 == 0 ) then
      call mask_finalize(temp_masks(1))
      rdf_add = 0
      return
    end if

    write(*, 1012)
1012 format('Enter mask for second atom set in RDF calculation')
    call mask_initialize(temp_masks(2))
    ats2 =  maskmanip_make(temp_masks(2))

    if(ats2 == 0) then
      call mask_finalize(temp_masks(1))
      call mask_finalize(temp_masks(2))
      rdf_add = 0
      return
    end if


    !create rdf_listan
    allocate(rdf_listan(Nmasks)%atom1(ats1*ats2))
    allocate(rdf_listan(Nmasks)%atom2(ats1*ats2))
    !rdf_make_list returns number of pairs  
    rdf_listan(Nmasks)%number_of_pairs = rdf_make_list(temp_masks(1),temp_masks(2),rdf_listan(Nmasks))

    call mask_finalize(temp_masks(1))
    call mask_finalize(temp_masks(2))

    ! set the parameters for this calculation
    rdf_calcs(Nmasks)%atoms_in_mask1 = ats1
    rdf_calcs(Nmasks)%rdf_radius = get_real_arg('RDF radius      :')        ! get RDF radius
    rdf_calcs(Nmasks)%Nbins      = get_int_arg('Number of bins  :')         ! get # bins
    allocate(rdf_calcs(Nmasks)%bins(rdf_calcs(Nmasks)%Nbins))               ! allocate bins
    rdf_calcs(Nmasks)%bins(:) = 0                                           ! set bins(:) = 0

    rdf_add = Nmasks
    write(desc, 25) ats1
25  format('RDF calculation for ',i6,' atoms')
  end function rdf_add


!------------------------------------------------------------------------------!
!>  subroutine: rdf_calc
!>  Calculates RDF data from the current trajectory frame
!------------------------------------------------------------------------------!
  subroutine rdf_calc(i)
    !arguments
    integer, intent(in)     :: i

    !locals
    integer                 :: j       ! atom index
    integer                 :: binidx  ! bin index
    real                    :: dist    ! distance

    !       do j = 1, nat_pro
    !               if(target_masks(i)%mask(j))     then                            ! if atom j in this target_mask
    !                       dist = r(rdf_calcs(i)%center, j)                        ! calculate distance between center atom and atom j
    !                       if(dist < rdf_calcs(i)%rdf_radius) then                 ! if atom j is within rdf_radius of center atom
    !                               binidx = dist / rdf_calcs(i)%rdf_radius * rdf_calcs(i)%Nbins + 1
    !                               rdf_calcs(i)%bins(binidx) = rdf_calcs(i)%bins(binidx) + 1
    !                       end if
    !               end if
    !       end do


    do j = 1, rdf_listan(i)%number_of_pairs
      dist = r(rdf_listan(i)%atom1(j),rdf_listan(i)%atom2(j))
      if(dist < rdf_calcs(i)%rdf_radius) then
        binidx = dist / rdf_calcs(i)%rdf_radius * rdf_calcs(i)%Nbins + 1
        rdf_calcs(i)%bins(binidx) = rdf_calcs(i)%bins(binidx) + 1
      end if
    end do


    rdf_calcs(i)%frames = rdf_calcs(i)%frames + 1 ! update # processed frames

    write(*,100, advance='no')      rdf_calcs(i)%atoms_in_mask1
100 format(i10)
  end subroutine rdf_calc


!------------------------------------------------------------------------------!
!>  subroutine: rdf_heading
!>
!------------------------------------------------------------------------------!
  subroutine rdf_heading(i)
    integer                        :: i

    write(*,'(a)', advance='no') 'RDF ats'
  end subroutine rdf_heading


!------------------------------------------------------------------------------!
!>  function: r
!>  The pair distance function
!>  r(a,b) calculates the distance between atoms a and b
!------------------------------------------------------------------------------!
  real function r(a,b)
    ! arguments
    integer                        :: a,b !atom indices from mask
    ! locals
    real,dimension(3)              :: veca
    real,dimension(3)              :: vecb
    real,dimension(3)              :: delta

!   xyz vector for atom a =  xin(3*a-2:3*a)
!   xyz vector for atom b =  xin(3*b-2:3*b)
    veca = xin(3*a-2:3*a)
    vecb = xin(3*b-2:3*b)
    delta = veca - vecb
!    delta = xin(3*a-2:3*a)-xin(3*b-2:3*b)

    ! if PBC then adjust lengths according to periodicity - MA
    if( use_PBC ) then
      delta(1) = delta(1) - boxlength(1)*nint( delta(1)/boxlength(1) )
      delta(2) = delta(2) - boxlength(2)*nint( delta(2)/boxlength(2) )
      delta(3) = delta(3) - boxlength(3)*nint( delta(3)/boxlength(3) )
    end if

    r = sqrt(dot_product(delta,delta)) ! Euclidean distance


! In case there's a bug tracking need:
!    write(*,100, advance='yes')  a, b
!100 format('atom number a = ',i10, ' atom number b = ', i10)
!    write(*,200, advance='yes')  veca
!200 format(f14.6, ' coordinates atom A = ')
!    write(*,300, advance='yes')  vecb
!300 format(f14.6, ' coordinates atom B = ')
!    write(*,400, advance='yes')   delta
!400 format(' delta =  ', f14.6 )
!    write(*,500, advance='yes')   r
!500 format(' distance = ', f14.6)

  end function r


!------------------------------------------------------------------------------!
!>  function: rdf_make_list
!>  rdf_make_list makes a list of all atom-atom pairs in
!>  the masks and returns the number of pairs
!------------------------------------------------------------------------------!
  integer function rdf_make_list(mask1, mask2, rdf_list)
    integer                        :: i,j,k,l,m,nat1,nat2
    integer, allocatable           :: group1(:), group2(:)
    type(mask_type), intent(in)    :: mask1, mask2
    type(rdf_list_type)            :: rdf_list

    l = 1
    m = 1


    !make the list of atoms in first mask
    nat1 = mask1%included 
    allocate(group1(nat1))
    do j = 1, nat_pro
      if (mask1%MASK(j)) then
        group1(l) = j
        l = l+1
      end if
    end do

    !make the list of atoms in second mask
    nat2 = mask2%included 
    allocate(group2(nat2))
    do j = 1, nat_pro
      if (mask2%MASK(j)) then
        group2(m) = j
        m = m+1
      end if
    end do

    i=1
    do j = 1, nat1
      do k = 1, nat2
        if (group1(j) .ne. group2(k)) then
          rdf_list%atom1(i) = group1(j)
          rdf_list%atom2(i) = group2(k)
          i = i+1
        end if
      end do
    end do

    rdf_make_list = i-1    !set the return value

  end function rdf_make_list

end module calc_rdf
