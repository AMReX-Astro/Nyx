! Module providing mesh interaction routines for AGN feedback.
!
! Paul Ricker (2/2013)

module agn_mesh

!==============================================================================

use agn_geometry
use agn_params

implicit none

public

! Derived type to hold Cartesian mesh coordinate information

type agn_mesh_data
  integer                       :: Nx, Ny, Nz          ! # of zones in each coordinate
  double precision, allocatable :: x(:), y(:), z(:)    ! zone center arrays
  double precision, allocatable :: dx(:), dy(:), dz(:) ! zone spacing arrays
  double precision              :: xmin, xmax, ymin, ymax, zmin, zmax
                                           ! min, max values for each coord.
end type agn_mesh_data

contains

!==============================================================================

! average_over_sphere: Average a mesh quantity over a spherical volume of a
!                      given radius about a given center. The radius is
!                      specified in units of the mean zone spacing.
!
! arguments: A                     Mesh quantity to be averaged
!            m                     Mesh data structure
!            x, y, z               Coordinates of sphere center
!            R                     Radius of sphere

double precision function average_over_sphere(A, A_l1, A_l2, A_l3, A_h1, A_h2, A_h3, &
                                              x, y, z, R, dx)

integer         , intent(in) :: A_l1, A_l2, A_l3, A_h1, A_h2, A_h3
double precision, intent(in) :: A(A_l1:,A_l2:,A_l3:)
double precision, intent(in) :: x, y, z, R
double precision, intent(in) :: dx(:)

integer             :: i, j, k, Nr, Nmu, Nphi, is, js, ks
double precision    :: d_zone, rs, mus, phis, xs, ys, zs, Asum, dr, dmu, dphi, &
                       sint, sinp, cosp

Asum = 0.

Nr   = max(2, int(2*R))
Nmu  = max(5, int(1+2*R))
Nphi = max(10, int(5*R))

d_zone = (dx(1) + dx(2) + dx(3)) / 3.d0
dr   = R*d_zone / (Nr-1)
dmu  = 2. / (Nmu-1)
dphi = 2.*pi / (Nphi-1)

do k = 1, Nphi
  phis = (k-1)*dphi
  sinp = dsin(phis)
  cosp = dcos(phis)
  do j = 1, Nmu
    mus = (j-1)*dmu - 1.d0
    sint = sqrt(1.d0- mus**2)
    do i = 1, Nr
      rs = (i-1)*dr
      xs = x + rs*sint*cosp
      ys = y + rs*sint*sinp
      zs = z + rs*mus

      is = floor(xs/dx(1))
      js = floor(ys/dx(2))
      ks = floor(zs/dx(3))

      if (is .lt. A_l1 .or. is .gt. A_h1 .or. &
          js .lt. A_l2 .or. js .gt. A_h2 .or. &
          ks .lt. A_l3 .or. ks .gt. A_h3) then
          print *,'ACK -- sphere going outside of array bounds '
          print *,'DX         ',dx(:)
          print *,'(x,y,z) is ',xs,ys,zs
          print *,'(i,j,k) is ',is,js,ks
          print *,'   A_lo is ',A_l1,A_l2,A_l3
          print *,'   A_hi is ',A_h1,A_h2,A_h3
          stop
!     else 
!        print *,'ADDING TO ',is,js,ks
      end if

      Asum = Asum + A(is,js,ks)
    enddo
  enddo
enddo

average_over_sphere = Asum / (Nr*Nmu*Nphi)

print*,'SUM ',Asum, Nr, Nmu, Nphi, average_over_sphere  

return
end function average_over_sphere

!------------------------------------------------------------------------------

! map_sphere_onto_mesh: Map a function onto a (possibly) non-Cartesian region
!                       of a Cartesian mesh.  The routine abstracts out the size
!                       and shape of the geometrical region.
!
! arguments: x, y, z            Characteristic coordinates of region
!            R                  Radius of sphere
!            norm               Weighting function normalization (eg. mass)
!            A                  Mesh quantity to map (eg. density)

subroutine map_sphere_onto_mesh(x, y, z, R, norm, &
                                A, A_l1, A_l2, A_l3, A_h1, A_h2, A_h3, &
                                dx)

double precision, intent(in   ) :: x, y, z, R, norm
double precision, intent(in   ) :: dx(:)

integer         , intent(in   ) :: A_l1, A_l2, A_l3, A_h1, A_h2, A_h3
double precision, intent(inout) :: A(A_l1:,A_l2:,A_l3:) 

double precision                :: vol_sphere

integer, parameter  :: Nsub = 3
integer             :: i, j, k, ii, jj, kk, region
double precision    :: Aavg, f, xx, yy, zz, dxx, dyy, dzz
double precision    :: xgl, ygl, zgl, xgr, ygr, zgr

vol_sphere = 4.d0/3.d0*pi*R**3
Aavg = norm / vol_sphere

dxx = dx(1) / Nsub
dyy = dx(2) / Nsub
dzz = dx(3) / Nsub

do k = A_l1, A_h1
  zgl = (dble(k)-0.5d0)*dx(3) - z
  zgr = zgl + dx(3)

  do j = A_l2, A_h2
    ygl = (dble(j)-0.5d0)*dx(2) - y
    ygr = ygl + dx(2)

    do i = A_l3, A_h3
      xgl = (dble(i)-0.5d0)*dx(1) - x
      xgr = xgl + dx(1)

      region = loc_sphere(xgl, xgr, ygl, ygr, zgl, zgr, R)

      if (region /= OUTSIDE) then
        f = 0.
        do kk = 1, Nsub
          zz = zgl + (kk-0.5)*dzz
          do jj = 1, Nsub
            yy = ygl + (jj-0.5)*dyy
            do ii = 1, Nsub
              xx = xgl + (i-0.5)*dxx
              if (loc_sphere(xx, xx, yy, yy, zz, zz, R) /= OUTSIDE) &
                f = f + Aavg
            enddo
          enddo
        enddo
        A(i,j,k) = A(i,j,k) + f/Nsub**3
      endif

    enddo
  enddo
enddo

return
end subroutine map_sphere_onto_mesh

!------------------------------------------------------------------------------

! map_cylinder_onto_mesh: Map a cylinder onto a (possibly) non-Cartesian region
!                         of a Cartesian mesh.  The routine abstracts out the size
!                         and shape of the geometrical region.
!
! arguments: x, y, z            Characteristic coordinates of region
!            weight             Weighting function
!            wparams            List of weight function parameters (passed)
!            norm               Weighting function normalization (eg. mass)
!            location           Function that determines where a point is with
!                               respect to mapping region
!            volume             Function that returns volume of mapping region
!            A                  Mesh quantity to map (eg. density)
!            m                  Mesh data object

subroutine map_cylinder_onto_mesh(x, y, z, R, height, nxc, nyc, nzc, norm, &
                                  A, A_l1, A_l2, A_l3, A_h1, A_h2, A_h3, &
                                  weight, wparams, dx)

double precision   , intent(in   ) :: x, y, z, norm
double precision   , intent(in   ) :: R, height, nxc, nyc, nzc
double precision   , intent(in   ) :: wparams(:)
integer            , intent(in   ) :: A_l1, A_l2, A_l3, A_h1, A_h2, A_h3
double precision   , intent(inout) :: A(A_l1:,A_l2:,A_l3:) 

double precision                   :: weight, vol_cylinder
double precision                   :: dx(:)

integer, parameter  :: Nsub = 3
integer             :: i, j, k, ii, jj, kk, region
double precision    :: Aavg, f, xx, yy, zz, dxx, dyy, dzz
double precision    :: xgl, ygl, zgl, xgr, ygr, zgr

vol_cylinder = pi * R**2 * height
Aavg = norm / vol_cylinder

dxx = dx(1) / Nsub
dyy = dx(2) / Nsub
dzz = dx(3) / Nsub

do k = A_l3, A_h3
  zgl = (dble(k)-0.5d0)*dx(3) - z
  zgr = zgl + dx(3)
  do j = A_l2, A_h2
    ygl = (dble(j)-0.5d0)*dx(2) - y
    ygr = ygl + dx(2)
    do i = A_l1, A_h1
      xgl = (dble(i)-0.5d0)*dx(1) - x
      xgr = xgl + dx(1)

      region = loc_cylinder(xgl, xgr, ygl, ygr, zgl, zgr, R, height, nxc, nyc, nzc)

      if (region /= OUTSIDE) then
        f = 0.
        do kk = 1, Nsub
          zz = zgl + (kk-0.5)*dzz
          do jj = 1, Nsub
            yy = ygl + (jj-0.5)*dyy
            do ii = 1, Nsub
              xx = xgl + (i-0.5)*dxx
              if (loc_cylinder(xx, xx, yy, yy, zz, zz, R, height, nxc, nyc, nzc) /= OUTSIDE) &
                f = f + Aavg*weight(xx, yy, zz, wparams)
            enddo
          enddo
        enddo
        A(i,j,k) = A(i,j,k) + f/Nsub**3
      endif

    enddo
  enddo
enddo

return
end subroutine map_cylinder_onto_mesh

end module agn_mesh
