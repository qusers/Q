!------------------------------------------------------------------------------!
!  Q version 6.0.1                                                             !
!  code authors: johan aqvist, martin almlof, martin ander, jens carlson,      !
!  isabella feierberg, peter hanspers, anders kaplan, karin kolmodin,          !
!  petra wennerstrom, kajsa ljunjberg, john marelius, martin nervall,          !
!  johan sund, ake sandgren, alexandre barrozo, masoud kazemi, paul bauer,     !
!  miha purg, irek szeler, mauricio esguerra                                   !
!  latest update: august 29, 2017                                              !
!------------------------------------------------------------------------------!


module calc_com
!!-------------------------------------------------------------------------------
!!  copyright (c) 2017 johan aqvist, john marelius, shina caroline lynn kamerlin
!!  and paul bauer
!!  **module calc_com**
!!  by martin almlof
!!  calculates the center of mass coordinates of whatever atom mask is specified
!!-------------------------------------------------------------------------------
  use calc_base
  use maskmanip
  implicit none

  !constants
  integer, parameter               :: max_masks = 10
  !module variables
  type(mask_type), private, target :: masks(max_masks)
  integer, private                 :: nmasks = 0
  type com_coord_type
    real, pointer                  :: x(:), y(:), z(:), mass(:)
  end type com_coord_type
  type(com_coord_type), private    :: coords_mass(max_masks)

  type coord_type
    real, pointer                  :: xyz(:)
  end type coord_type

  type(coord_type), private        :: coords(max_masks)

  type mass_ave_type
    real                           :: x,y,z
  end type mass_ave_type

  type(mass_ave_type), private     :: mass_ave(max_masks)
  real,private                     :: tot_mass(max_masks)


contains

subroutine com_initialize
end subroutine com_initialize

subroutine com_finalize(i)
  integer                          :: i
  call mask_finalize(masks(i))
end subroutine com_finalize


integer function com_add(desc)
  !arguments
  character(*)                     :: desc
  character(len=80)                :: line
  integer                          :: readstat
  integer                          :: ats,j
  if(nmasks == max_masks) then
    write(*,10) max_masks
    return
  end if
10 format('sorry, the maximum number of com calculations is ',i2)

  !add a new com mask
  nmasks = nmasks + 1
  call mask_initialize(masks(nmasks))
  ats =  maskmanip_make(masks(nmasks))
  !discard if no atoms in mask
  if(ats == 0) then
    call mask_finalize(masks(nmasks))
    nmasks = nmasks - 1
    com_add = 0
    return
  end if

  allocate(coords(nmasks)%xyz(3*ats))
  allocate(coords_mass(nmasks)%x(ats), coords_mass(nmasks)%y(ats), coords_mass(nmasks)%z(ats), coords_mass(nmasks)%mass(ats))

  coords_mass(nmasks)%x(:) = 0
  coords_mass(nmasks)%y(:) = 0
  coords_mass(nmasks)%z(:) = 0
  coords_mass(nmasks)%mass(:) = 0
  coords(nmasks)%xyz(:) = 0

  call com_put_mass(nmasks)
  com_add = nmasks
  write(desc, 20) masks(nmasks)%included
20 format('center of mass position ',i6,' atoms')
end function com_add


subroutine com_calc(i)
  !arguments
  integer, intent(in)              :: i

  !locals
  if(i < 1 .or. i > nmasks) return
  call mask_get(masks(i), xin, coords(i)%xyz)

  !split coords into x, y, and z coords
  coords_mass(i)%x = coords(i)%xyz(1::3)
  coords_mass(i)%y = coords(i)%xyz(2::3)
  coords_mass(i)%z = coords(i)%xyz(3::3)

  !calculate center of mass
  mass_ave(i)%x = dot_product(coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
  mass_ave(i)%y = dot_product(coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
  mass_ave(i)%z = dot_product(coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)

  write(*,100, advance='no') mass_ave(i)%x,mass_ave(i)%y,mass_ave(i)%z
100 format(3f9.4)
end subroutine com_calc


subroutine com_put_mass(i)
  integer                          :: k,j,i,at
  real                             :: mass

  if(i < 1 .or. i > nmasks) return

  tot_mass(i) = 0
  !put in masses into coords_mass
  k=1
  do j = 1, nat_pro
    if (masks(i)%mask(j)) then
      mass = iaclib(iac(j))%mass
      coords_mass(i)%mass(k) = mass
      tot_mass(i) = tot_mass(i) + mass
      k = k+1
    end if
  end do

  write(*,168) "total mass: ",tot_mass(i)
168 format(a,f10.3)

end subroutine com_put_mass


subroutine com_heading(i)
  integer                          :: i

  write(*,'(a)', advance='no') '    x        y        z    '
end subroutine com_heading


end module calc_com
