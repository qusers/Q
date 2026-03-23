module settle
!!-------------------------------------------------------------------------------
!!  SETTLE constraint algorithm for rigid 3-site water molecules
!!
!!  Miyamoto S, Kollman PA (1992). "SETTLE: An analytical version of the SHAKE
!!  and RATTLE algorithm for rigid water models." J Comput Chem 13(8):952-962.
!!-------------------------------------------------------------------------------
  implicit none

  ! SETTLE geometry parameters (precomputed in init_settle)
  real(8) :: settle_dOH       ! target O-H distance
  real(8) :: settle_dHH       ! target H-H distance
  real(8) :: settle_mO        ! oxygen mass
  real(8) :: settle_mH        ! hydrogen mass
  real(8) :: settle_ra        ! height of triangle from O to HH midpoint
  real(8) :: settle_rb        ! half of H-H distance
  real(8) :: settle_rc        ! distance from COM to O along ra axis
  real(8) :: settle_wo        ! mass weight for oxygen: mO / (mO + 2*mH)
  real(8) :: settle_wh        ! mass weight for hydrogen: mH / (mO + 2*mH)
  logical :: use_settle = .false.
  integer :: settle_nwat = 0  ! number of water molecules

contains

subroutine setup_settle(dOH, angle_HOH, massO, massH)
!!-------------------------------------------------------------------------------
!!  Precompute SETTLE geometry constants from water parameters.
!!
!!  dOH       : equilibrium O-H bond length
!!  angle_HOH : equilibrium H-O-H angle in radians
!!  massO     : oxygen mass
!!  massH     : hydrogen mass
!!-------------------------------------------------------------------------------
  real(8), intent(in) :: dOH, angle_HOH, massO, massH
  real(8) :: total_mass

  settle_dOH = dOH
  settle_mO  = massO
  settle_mH  = massH

  ! H-H distance from law of cosines
  settle_dHH = 2.0d0 * dOH * sin(angle_HOH / 2.0d0)

  ! ra = height of isoceles triangle: distance from O to midpoint of HH
  settle_ra = dOH * cos(angle_HOH / 2.0d0)

  ! rb = half the H-H distance
  settle_rb = settle_dHH / 2.0d0

  ! mass weights
  total_mass = massO + 2.0d0 * massH
  settle_wo  = massO / total_mass
  settle_wh  = massH / total_mass

  ! rc = distance from COM to the HH baseline along the principal axis
  ! In the canonical frame: O at (0, ra, 0), H at (±rb, 0, 0)
  ! COM_y = wo*ra, so rc = wo*ra
  ! Then: O in COM frame at y = ra - rc = 2*wh*ra
  !        H in COM frame at y = -rc = -wo*ra
  settle_rc = settle_ra * settle_wo

end subroutine setup_settle


subroutine settle_positions(nwat, nat_solute, xx, x, winv)
!!-------------------------------------------------------------------------------
!!  Apply SETTLE constraints to all water molecules.
!!
!!  nwat       : number of water molecules
!!  nat_solute : number of solute atoms (waters start at nat_solute+1)
!!  xx         : old (pre-step) positions, dimension(3*natom)
!!  x          : unconstrained positions to be corrected, dimension(3*natom)
!!  winv       : inverse masses, dimension(natom) — only used for mass ratio
!!
!!  Water atoms are assumed to be sequential triplets: O, H1, H2
!!  starting at atom index nat_solute+1.
!!-------------------------------------------------------------------------------
  integer, intent(in) :: nwat, nat_solute
  real(8), intent(in) :: xx(:)
  real(8), intent(inout) :: x(:)
  real, intent(in) :: winv(:)  ! single precision, matching md.f90

  integer :: w, iO, iH1, iH2, i3O, i3H1, i3H2
  real(8) :: xO0(3), xH10(3), xH20(3)  ! old positions
  real(8) :: xO1(3), xH11(3), xH21(3)  ! unconstrained positions
  real(8) :: com0(3), com1(3)           ! centers of mass
  real(8) :: a0(3), b0(3), c0(3)       ! old frame vectors
  real(8) :: a1(3), b1(3), c1(3)       ! unconstrained frame vectors
  real(8) :: n0(3), n1(3), n2(3)       ! principal axes of old triangle
  real(8) :: a1p(3), b1p(3), c1p(3)    ! unconstrained in old frame
  real(8) :: sinphi, cosphi, sinpsi, cospsi
  real(8) :: xa1p, ya1p, za1p, xb1p, yb1p, zb1p, xc1p, yc1p, zc1p
  real(8) :: alpha, beta, gamma, a2b2, sinalpha, cosalpha
  real(8) :: tmp1, tmp2, tmp3
  real(8) :: ya2, xb2, yb2, xc2, yc2
  real(8) :: invlen, total_mass_inv
  real(8) :: ra, rb, rc
  real(8) :: xH1_new(3), xH2_new(3)
  real(8) :: xO_new(3)

  ra = settle_ra
  rb = settle_rb
  rc = settle_rc

  total_mass_inv = 1.0d0 / (settle_mO + 2.0d0 * settle_mH)

  do w = 1, nwat
    ! atom indices (1-based)
    iO  = nat_solute + 3*(w-1) + 1
    iH1 = iO + 1
    iH2 = iO + 2

    ! coordinate array indices (1-based, interleaved xyz)
    i3O  = iO*3  - 3
    i3H1 = iH1*3 - 3
    i3H2 = iH2*3 - 3

    ! extract old positions
    xO0  = xx(i3O+1  : i3O+3)
    xH10 = xx(i3H1+1 : i3H1+3)
    xH20 = xx(i3H2+1 : i3H2+3)

    ! extract unconstrained positions
    xO1  = x(i3O+1  : i3O+3)
    xH11 = x(i3H1+1 : i3H1+3)
    xH21 = x(i3H2+1 : i3H2+3)

    ! center of mass of old positions
    com0 = settle_wo * xO0 + settle_wh * (xH10 + xH20)

    ! center of mass of unconstrained positions (preserved by SETTLE)
    com1 = settle_wo * xO1 + settle_wh * (xH11 + xH21)

    ! vectors from COM to atoms in old configuration
    a0 = xO0  - com0
    b0 = xH10 - com0
    c0 = xH20 - com0

    ! vectors from COM to atoms in unconstrained configuration
    a1 = xO1  - com1
    b1 = xH11 - com1
    c1 = xH21 - com1

    ! build orthonormal basis from old triangle
    ! n0 = unit vector along (b0+c0)/2 - this is the axis from COM toward O (y-like)
    ! Actually, we use the Miyamoto-Kollman convention:
    ! n1 = unit vector along b0-c0 (x-axis, H1-H2 direction)
    ! n0 = unit vector perpendicular to n1 in the plane of the triangle (y-axis)
    ! n2 = n0 x n1 (z-axis, out of plane)

    n1 = b0 - c0
    invlen = 1.0d0 / sqrt(dot_product(n1, n1))
    n1 = n1 * invlen

    ! n0 is in the plane, perpendicular to n1
    ! a0 is approximately along y; subtract its n1 component
    n0 = a0 - dot_product(a0, n1) * n1
    invlen = 1.0d0 / sqrt(dot_product(n0, n0))
    n0 = n0 * invlen

    ! n2 = n0 x n1
    n2(1) = n0(2)*n1(3) - n0(3)*n1(2)
    n2(2) = n0(3)*n1(1) - n0(1)*n1(3)
    n2(3) = n0(1)*n1(2) - n0(2)*n1(1)

    ! project unconstrained displacements onto old frame
    ! a1p, b1p, c1p are the unconstrained COM-relative vectors in (n1, n0, n2) frame
    ! Convention: x = n1 component, y = n0 component, z = n2 component
    xa1p = dot_product(a1, n1)
    ya1p = dot_product(a1, n0)
    za1p = dot_product(a1, n2)
    xb1p = dot_product(b1, n1)
    yb1p = dot_product(b1, n0)
    zb1p = dot_product(b1, n2)
    xc1p = dot_product(c1, n1)
    yc1p = dot_product(c1, n0)
    zc1p = dot_product(c1, n2)

    ! solve for the constrained positions in the old frame
    ! Step 1: solve for sinphi, cosphi (rotation around x-axis to put a1 in xy plane)
    tmp1 = ya1p**2 + za1p**2
    if (tmp1 > 1.0d-30) then
      tmp1 = 1.0d0 / sqrt(tmp1)
      sinphi = za1p * tmp1
      cosphi = ya1p * tmp1
    else
      sinphi = 0.0d0
      cosphi = 1.0d0
    end if

    ! rotate b1p and c1p by -phi around x-axis
    yb1p =  cosphi * dot_product(b1, n0) + sinphi * dot_product(b1, n2)
    zb1p = -sinphi * dot_product(b1, n0) + cosphi * dot_product(b1, n2)
    yc1p =  cosphi * dot_product(c1, n0) + sinphi * dot_product(c1, n2)
    zc1p = -sinphi * dot_product(c1, n0) + cosphi * dot_product(c1, n2)
    ya1p = sqrt(ya1p**2 + za1p**2)  ! za1p is now 0

    ! Step 2: solve for sinpsi, cospsi (rotation around y-axis)
    ! This brings b1 and c1 into the xy-plane
    sinpsi = (zb1p - zc1p) / (2.0d0 * rb)
    ! Guard against numerical issues
    if (abs(sinpsi) > 1.0d0) sinpsi = sign(1.0d0, sinpsi)
    cospsi = sqrt(1.0d0 - sinpsi**2)

    ! Step 3: solve for alpha (rotation angle around z-axis)
    ! The constrained triangle is parameterized by alpha. We find alpha that
    ! minimizes mass-weighted squared displacement from unconstrained positions:
    !   d/dalpha [ sum_i m_i |r_constrained_i - r_unconstrained_i|^2 ] = 0
    ! gives: A*sin(alpha) + B*cos(alpha) = 0
    ! Coefficients simplified via COM constraints (yb1p+yc1p = -wO/wH*ya1p):
    !   A = mO*ra*ya1p + mH*xb2*(xb1p - xc1p)
    !   B = -mO*ra*xa1p + mH*xb2*(yb1p - yc1p)

    ya2 = ra - rc
    xb2 =  rb * cospsi
    xc2 = -rb * cospsi

    alpha = settle_mO * ra * ya1p + settle_mH * xb2 * (xb1p - xc1p)
    beta  = -settle_mO * ra * xa1p + settle_mH * xb2 * (yb1p - yc1p)
    a2b2  = alpha*alpha + beta*beta

    if (a2b2 > 0.0d0) then
      tmp1 = 1.0d0 / sqrt(a2b2)
      cosalpha =  alpha * tmp1
      sinalpha = -beta  * tmp1
    else
      cosalpha = 1.0d0
      sinalpha = 0.0d0
    end if

    ! Now construct the constrained positions in the rotated frame:
    ! a2'' = (0, ya2, 0)
    ! b2'' = (xb2, -rc, -zb1p + ...)
    ! Actually, the z-coordinates are recovered from the psi rotation

    ! The full constrained positions in the old principal frame are:
    ! Undo psi rotation (around y-axis) and phi rotation (around x-axis):
    ! a2' = (xa2, ya2*cosphi, -ya2*sinphi) but xa2 = 0 for oxygen
    ! Wait, I need to be more careful here.

    ! Let me restart with the clean SETTLE algorithm from M&K eq. 1-17.
    ! The constrained positions in the canonical molecular frame:
    ya2 = ra - rc
    yb2 = -rc
    xb2 =  rb * cospsi
    xc2 = -rb * cospsi

    ! These are in the frame after phi, psi, alpha rotations.
    ! We now need to undo: alpha (around z), phi (around x), then transform back to lab.

    ! Undo alpha rotation (around z-axis, so only affects x,y):
    ! For oxygen: x=0, y=ya2 → x'= 0*cos(a)+ya2*sin(a), y'= -0*sin(a)+ya2*cos(a)
    ! Wait, the forward rotation was: x_rot = x*cos(a) - y*sin(a), y_rot = x*sin(a) + y*cos(a)
    ! Inverse: x = x_rot*cos(a) + y_rot*sin(a), y = -x_rot*sin(a) + y_rot*cos(a)

    ! Undo alpha rotation:
    a1p(1) =  0.0d0 * cosalpha + ya2 * sinalpha   ! = ya2*sinalpha
    a1p(2) = -0.0d0 * sinalpha + ya2 * cosalpha    ! = ya2*cosalpha = ya1p
    a1p(3) = 0.0d0  ! oxygen stays in xy-plane

    b1p(1) =  xb2 * cosalpha + yb2 * sinalpha
    b1p(2) = -xb2 * sinalpha + yb2 * cosalpha
    b1p(3) = zb1p  ! z-coordinate preserved through alpha rotation

    c1p(1) =  xc2 * cosalpha + yb2 * sinalpha    ! yc2 = yb2 = -rc
    c1p(2) = -xc2 * sinalpha + yb2 * cosalpha
    c1p(3) = zc1p  ! z-coordinate preserved through alpha rotation

    ! Undo phi rotation (around x-axis, affects y,z):
    ! Forward: y_rot = y*cos(phi) + z*sin(phi), z_rot = -y*sin(phi) + z*cos(phi)
    ! Inverse: y = y_rot*cos(phi) - z_rot*sin(phi), z = y_rot*sin(phi) + z_rot*cos(phi)

    ! For oxygen:
    a1(1) = a1p(1)
    a1(2) = a1p(2) * cosphi - a1p(3) * sinphi
    a1(3) = a1p(2) * sinphi + a1p(3) * cosphi

    ! For H1:
    b1(1) = b1p(1)
    b1(2) = b1p(2) * cosphi - b1p(3) * sinphi
    b1(3) = b1p(2) * sinphi + b1p(3) * cosphi

    ! For H2:
    c1(1) = c1p(1)
    c1(2) = c1p(2) * cosphi - c1p(3) * sinphi
    c1(3) = c1p(2) * sinphi + c1p(3) * cosphi

    ! Transform back to lab frame: constrained = com1 + a1(1)*n1 + a1(2)*n0 + a1(3)*n2
    xO_new  = com1 + a1(1)*n1 + a1(2)*n0 + a1(3)*n2
    xH1_new = com1 + b1(1)*n1 + b1(2)*n0 + b1(3)*n2
    xH2_new = com1 + c1(1)*n1 + c1(2)*n0 + c1(3)*n2

    x(i3O+1  : i3O+3)  = xO_new
    x(i3H1+1 : i3H1+3) = xH1_new
    x(i3H2+1 : i3H2+3) = xH2_new
  end do

end subroutine settle_positions

end module settle
