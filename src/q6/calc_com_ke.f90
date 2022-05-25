!------------------------------------------------------------------------------!
!  Q version 6.0.1                                                             !
!  code authors: johan aqvist, martin almlof, martin ander, jens carlson,      !
!  isabella feierberg, peter hanspers, anders kaplan, karin kolmodin,          !
!  petra wennerstrom, kajsa ljunjberg, john marelius, martin nervall,          !
!  johan sund, ake sandgren, alexandre barrozo, masoud kazemi, paul bauer,     !
!  miha purg, irek szeler, mauricio esguerra                                   !
!  latest update: august 29, 2017                                              !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!  copyright (c) 2017 johan aqvist, john marelius, shina caroline lynn kamerlin
!!  and paul bauer
!  calc_kineticenergy.f90
!  by martin almlof
!  calculates the kinetic energy of the center of mass of whatever atom
!  mask is specified
!------------------------------------------------------------------------------!
module calc_com_ke
  use calc_base
  use maskmanip
! use lapackmini
  implicit none

!constants
  integer, parameter               :: max_masks = 10
!module variables
  integer, parameter, private      :: conversion_factor = 2390.0574   ! gram/mol**^2/fs^2  -->  kcal/mol
  integer, private                 :: frames(max_masks), apa
  real(8), allocatable             :: kineticenergy(:)
  type(mask_type), private, target :: masks(max_masks)
  integer, private                 :: nmasks = 0

  type com_ke_coord_type
    real, pointer                  :: x(:), y(:), z(:), mass(:)
  end type com_ke_coord_type
  type(com_ke_coord_type), private :: coords_mass(max_masks), prev_coords_mass(max_masks)

  type com_ke_velocity_type
    real, pointer                  :: x(:), y(:), z(:)
  end type com_ke_velocity_type
  type(com_ke_velocity_type), private :: velocity(max_masks), rel_coords(max_masks), prev_rel_coords(max_masks),rad_vec(max_masks,3) !rel_coords is not velocities

  type coord_type
    real, pointer                  :: xyz(:)
  end type coord_type
  type(coord_type), private        :: coords(max_masks), prev_coords(max_masks)

  type dp_type
    real, pointer                  :: dp(:)
  end type dp_type
  type(dp_type), private           :: dp_vect(max_masks)

  type mass_ave_type
    real                          :: x,y,z
  end type mass_ave_type
  type(mass_ave_type), private    :: mass_ave(max_masks), prev_mass_ave(max_masks) , ang_momentum(max_masks,3)

  type eigen_stuff_type
    real                          :: evalue(3),evector(3,3)
  end type eigen_stuff_type
  type(eigen_stuff_type), private :: eigen_stuff(max_masks)


  real,private                    :: tot_mass(max_masks), ke_rot(max_masks,3)
  logical,private                 :: first_frame(max_masks) = .true.
!       real,private                    :: previous_mass_center(3,max_masks)
  real, private                   :: frame_length = 0

contains

subroutine com_ke_initialize
end subroutine com_ke_initialize


subroutine com_ke_finalize(i)
  integer                         :: i

  call mask_finalize(masks(i))
end subroutine com_ke_finalize


integer function com_ke_add(desc)
  !arguments
  character(*)                            ::      desc
  character(len=80)                       ::      line
  integer                                         ::      readstat
  integer                                         ::      ats,j
  if(nmasks == max_masks) then
    write(*,10) max_masks
    return
  end if
10 format('sorry, the maximum number of com_ke calculations is ',i2)


  !get frame time length
  if (frame_length == 0) then
    write(*,'(a)', advance='no') 'enter time span between each trajectory frame (fs): '
    read(*,'(a)', iostat=readstat) line
    read(line, *, iostat=readstat) frame_length
    if(readstat /= 0 .or. frame_length <= 0) then
      write(*, 900)
900   format('>>>>> error: invalid time span.')
      com_ke_add = 0
      return
    end if
  end if

  !add a new com_ke mask
  nmasks = nmasks + 1
  call mask_initialize(masks(nmasks))
  ats =  maskmanip_make(masks(nmasks))
  !discard if no atoms in mask
  if(ats == 0) then
    call mask_finalize(masks(nmasks))
    nmasks = nmasks - 1
    com_ke_add = 0
    return
  end if

  allocate(coords(nmasks)%xyz(3*ats), prev_coords(nmasks)%xyz(3*ats))
  allocate(coords_mass(nmasks)%x(ats), coords_mass(nmasks)%y(ats), coords_mass(nmasks)%z(ats), coords_mass(nmasks)%mass(ats))
  allocate(prev_coords_mass(nmasks)%x(ats), prev_coords_mass(nmasks)%y(ats))
  allocate(prev_coords_mass(nmasks)%z(ats), prev_coords_mass(nmasks)%mass(ats))
  allocate(velocity(nmasks)%x(ats), velocity(nmasks)%y(ats), velocity(nmasks)%z(ats))
  allocate(rel_coords(nmasks)%x(ats), rel_coords(nmasks)%y(ats), rel_coords(nmasks)%z(ats))
  allocate(dp_vect(nmasks)%dp(ats) ,prev_rel_coords(nmasks)%x(ats))
  allocate(prev_rel_coords(nmasks)%y(ats), prev_rel_coords(nmasks)%z(ats))

  do j=1,3
    allocate(rad_vec(nmasks,j)%x(ats),rad_vec(nmasks,j)%y(ats),rad_vec(nmasks,j)%z(ats))
  end do



  coords_mass(nmasks)%x(:) = 0
  coords_mass(nmasks)%y(:) = 0
  coords_mass(nmasks)%z(:) = 0
  coords_mass(nmasks)%mass(:) = 0
  coords(nmasks)%xyz(:) = 0

  frames(nmasks) = 0
  call com_ke_put_mass(nmasks)
  com_ke_add = nmasks
  write(desc, 20) masks(nmasks)%included
20 format('center of mass kinetic energy for ',i6,' atoms')
end function com_ke_add


subroutine com_ke_calc(i)
  !arguments
  integer, intent(in)     :: i
  integer                 :: info, ipiv(3,3),j
  double precision        :: a(3,3), b(3), w(30), k(6), c(6)

  !locals
  real(8)                 :: ke, ixx, ixy, ixz, iyy, iyz, izz, tot_ke_rot

  if(i < 1 .or. i > nmasks) return

  frames(i)=frames(i) + 1

  if (first_frame(i)) then
    call mask_get(masks(i), xin, coords(i)%xyz)
    prev_coords(i)%xyz = coords(i)%xyz
    first_frame(i) = .false.
  else
    prev_coords(i)%xyz = coords(i)%xyz
    call mask_get(masks(i), xin, coords(i)%xyz)
  end if

  !split coords into x, y, and z coords
  coords_mass(i)%x = coords(i)%xyz(1::3)
  coords_mass(i)%y = coords(i)%xyz(2::3)
  coords_mass(i)%z = coords(i)%xyz(3::3)
  prev_coords_mass(i)%x = prev_coords(i)%xyz(1::3)
  prev_coords_mass(i)%y = prev_coords(i)%xyz(2::3)
  prev_coords_mass(i)%z = prev_coords(i)%xyz(3::3)


  !calculate center of mass
  mass_ave(i)%x = dot_product(coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
  mass_ave(i)%y = dot_product(coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
  mass_ave(i)%z = dot_product(coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)
  prev_mass_ave(i)%x = dot_product(prev_coords_mass(i)%x(:),coords_mass(i)%mass)/tot_mass(i)
  prev_mass_ave(i)%y = dot_product(prev_coords_mass(i)%y(:),coords_mass(i)%mass)/tot_mass(i)
  prev_mass_ave(i)%z = dot_product(prev_coords_mass(i)%z(:),coords_mass(i)%mass)/tot_mass(i)


  !create coordinate set relative mass center
  rel_coords(i)%x = coords_mass(i)%x - mass_ave(i)%x
  rel_coords(i)%y = coords_mass(i)%y - mass_ave(i)%y
  rel_coords(i)%z = coords_mass(i)%z - mass_ave(i)%z
  prev_rel_coords(i)%x = prev_coords_mass(i)%x - prev_mass_ave(i)%x
  prev_rel_coords(i)%y = prev_coords_mass(i)%y - prev_mass_ave(i)%y
  prev_rel_coords(i)%z = prev_coords_mass(i)%z - prev_mass_ave(i)%z


  !calculate moment of inertia tensor
  ixx = dot_product( (rel_coords(i)%y)**2 + (rel_coords(i)%z)**2, coords_mass(i)%mass)
  iyy = dot_product( (rel_coords(i)%x)**2 + (rel_coords(i)%z)**2, coords_mass(i)%mass)
  izz = dot_product( (rel_coords(i)%y)**2 + (rel_coords(i)%x)**2, coords_mass(i)%mass)
  ixy = -1._8 * sum( (rel_coords(i)%y) * (rel_coords(i)%x) * coords_mass(i)%mass )
  ixz = -1._8 * sum( (rel_coords(i)%x) * (rel_coords(i)%z) * coords_mass(i)%mass )
  iyz = -1._8 * sum( (rel_coords(i)%y) * (rel_coords(i)%z) * coords_mass(i)%mass )


  !calculate individual atom (relative?) velocities
  velocity(i)%x = (rel_coords(i)%x - prev_rel_coords(i)%x) / frame_length
  velocity(i)%y = (rel_coords(i)%y - prev_rel_coords(i)%y) / frame_length
  velocity(i)%z = (rel_coords(i)%z - prev_rel_coords(i)%z) / frame_length

  !       write (*,*) ' '
  k(:) = (/ixx,ixy,iyy,ixz,iyz,izz/)

  !       write(*,'(6f18.3)'), k(:)
  call eigen(k,a,3,0)

  do j = 1, 3
    !               write (*,*) a(1,j),a(2,j),a(3,j),k(j*(j+1)/2)
    !               write (*,*) info
    eigen_stuff(i)%evalue(j) = k(j*(j+1)/2)
    eigen_stuff(i)%evector(:,j) = a(:,j)
  !               write (*,'(4f23.8)') eigen_stuff(i)%evalue(j), eigen_stuff(i)%evector(:,j)

  end do

  !eigen_stuff(i)%evector are the normalized principal axes of rotation
  !vector of dot_product(principal axis,atom coordinate)

  do j=1,3
    dp_vect(i)%dp = eigen_stuff(i)%evector(1,j) * &
      rel_coords(i)%x + eigen_stuff(i)%evector(2,j) * &
      rel_coords(i)%y + eigen_stuff(i)%evector(3,j) * &
      rel_coords(i)%z

    rad_vec(i,j)%x = rel_coords(i)%x - dp_vect(i)%dp*eigen_stuff(i)%evector(1,j)
    rad_vec(i,j)%y = rel_coords(i)%y - dp_vect(i)%dp*eigen_stuff(i)%evector(2,j)
    rad_vec(i,j)%z = rel_coords(i)%z - dp_vect(i)%dp*eigen_stuff(i)%evector(3,j)

  end do
  !now rad_vec(i,j)%x,y,z contains all the vectors of each atom to the principal axis j

  !time to take the cross-product of each atoms rad_vec with it's linear momentum vector

  !       do j=1,3
  !               write (*,'(3f20.12)') sum(rad_vec(i,j)%x),sum(rad_vec(i,j)%y),sum(rad_vec(i,j)%z)
  !       end do

  do j= 1,3
    ang_momentum(i,j)%x = sum(rad_vec(i,j)%y*velocity(i)%z*coords_mass(i)%mass - &
      rad_vec(i,j)%z*velocity(i)%y*coords_mass(i)%mass)
    ang_momentum(i,j)%y = sum(rad_vec(i,j)%z*velocity(i)%x*coords_mass(i)%mass - &
      rad_vec(i,j)%x*velocity(i)%z*coords_mass(i)%mass)
    ang_momentum(i,j)%z = sum(rad_vec(i,j)%x*velocity(i)%y*coords_mass(i)%mass - &
      rad_vec(i,j)%y*velocity(i)%x*coords_mass(i)%mass)
  !               write (*,'(3f23.15)') ang_momentum(i,j)%x,ang_momentum(i,j)%y,ang_momentum(i,j)%z
  end do

  tot_ke_rot = 0
  !rotational kinetic energy = lj^2/(2ij)
  do j=1,3
    ke_rot(i,j) = ((ang_momentum(i,j)%x)**2+(ang_momentum(i,j)%y)**2+ &
      (ang_momentum(i,j)%z)**2)/(2*eigen_stuff(i)%evalue(j))*conversion_factor
    tot_ke_rot = tot_ke_rot + ke_rot(i,j)
  end do

  !       write (*,'(3f23.15)') ke_rot(i,1),ke_rot(i,2),ke_rot(i,3)
  !       calculate the kinetic energy
  ke = 0.5 * tot_mass(i) * ((prev_mass_ave(i)%x-mass_ave(i)%x)**2+ &
    (prev_mass_ave(i)%y-mass_ave(i)%y)**2+(prev_mass_ave(i)%z- &
    mass_ave(i)%z)**2) / frame_length**2 * conversion_factor

  write(*,100, advance='no') ke, tot_ke_rot
100 format(2f12.3)
end subroutine com_ke_calc


subroutine com_ke_put_mass(i)
  integer                                         ::      k,j,i,at
  real                                            ::      mass

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

end subroutine com_ke_put_mass

subroutine com_ke_heading(i)
  integer                                         ::      i

  write(*,'(a)', advance='no') 'trans  rot (kcal/mol)'
end subroutine com_ke_heading

end module calc_com_ke
