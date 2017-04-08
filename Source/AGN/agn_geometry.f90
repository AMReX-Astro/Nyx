! Module providing geometry-related routines for AGN feedback.
!
! Paul Ricker (2/2013)

module agn_geometry

use bl_constants_module, only: M_PI

!==============================================================================

implicit none

public

! Constants for geometric shapes

integer, parameter :: INSIDE = 0, STRADDLES = 1, OUTSIDE = 2

contains

!------------------------------------------------------------------------------

! loc_sphere: For a Cartesian zone with given left and right coordinate values,
!             determine whether it lies entirely inside, entirely outside, or
!             partially inside and outside a given sphere.
!
! arguments: xl, yl, zl             Left-face coordinates of zone relative to
!                                   sphere's center
!            xr, yr, zr             Right-face coordinates of zone relative to
!                                   sphere's center
!            params(1)              Parameter array containing radius of sphere

integer function loc_sphere(xl, xr, yl, yr, zl, zr, R)

double precision, intent(in) :: xl, xr, yl, yr, zl, zr, R

double precision             :: rl2, rr2, R2

R2  = R**2
rl2 = xl**2 + yl**2 + zl**2
rr2 = xr**2 + yr**2 + zr**2

if (rr2 <= R2) then
  loc_sphere = INSIDE
else if (rl2 > R2) then
  loc_sphere = OUTSIDE
else
  loc_sphere = STRADDLES
endif

return
end function loc_sphere

!------------------------------------------------------------------------------

! loc_cylinder: For a Cartesian zone with given left and right coordinate
!               values, determine whether it lies entirely inside, entirely
!               outside, or partially inside and outside a given cylinder.
!
! arguments: xl, yl, zl             Left-face coordinates of zone relative to
!                                   cylinder's origin point
!            xr, yr, zr             Right-face coordinates of zone relative to
!                                   cylinder's origin point
!            params(5)              Parameter array containing radius, height,
!                                   and orientation direction cosines of
!                                   cylinder

integer function loc_cylinder(xl, xr, yl, yr, zl, zr, R, h, nxc, nyc, nzc)

double precision, intent(in) :: xl, xr, yl, yr, zl, zr
double precision, intent(in) :: R, h, nxc, nyc, nzc

double precision             :: R2, xldotn, xrdotn, rl2, rr2
double precision             :: xlproj, xrproj, ylproj, yrproj, zlproj, zrproj

R2  = R**2

xldotn = xl*nxc + yl*nyc + zl*nzc
xrdotn = xr*nxc + yr*nyc + zr*nzc

xlproj = xl - xldotn*nxc
xrproj = xr - xrdotn*nxc
ylproj = yl - xldotn*nyc
yrproj = yr - xrdotn*nyc
zlproj = zl - xldotn*nzc
zrproj = zr - xrdotn*nzc

rl2 = xlproj**2 + ylproj**2 + zlproj**2
rr2 = xrproj**2 + yrproj**2 + zrproj**2

if ((xldotn >= 0.) .and. (xrdotn <= h) .and. (rr2 <= R2)) then
  loc_cylinder = INSIDE
else if ((xrdotn < 0.) .or. (xldotn > h) .or. (rl2 > R2)) then
  loc_cylinder = OUTSIDE
else
  loc_cylinder = STRADDLES
endif

return
end function loc_cylinder

!------------------------------------------------------------------------------

! weight_uniform: A uniform weighting function.
!
! arguments: x, y, z            Coordinates of sample point (ignored)
!            params(:)          Weighting function parameters (ignored)

double precision function weight_uniform(x, y, z, params)

double precision, intent(in) :: x, y, z, params(:)

weight_uniform = 1.

return
end function weight_uniform

!------------------------------------------------------------------------------

! weight_cattaneo: Cattaneo & Teyssier weighting function.
!
! arguments: x, y, z            Coordinates of sample point
!            params(:)          Weighting function parameters: R, h, nx, ny, nz, rscale

double precision function weight_cattaneo(x, y, z, R, height, nx, ny, nz, rscale)

double precision, intent(in) :: x, y, z, R, height, nx, ny, nz, rscale

double precision             :: xdotn, xproj, yproj, zproj, r2
double precision             :: vol_cylinder

xdotn = x*nx + y*ny + z*nz
xproj = x - xdotn*nx
yproj = y - xdotn*ny
zproj = z - xdotn*nz
r2    = xproj**2 + yproj**2 + zproj**2

vol_cylinder = M_PI * R**2 * height

if ((xdotn >= 0.) .and. (xdotn <= height) .and. (r2 <= R**2)) then
  weight_cattaneo = dexp(-r2/(2*rscale**2)) * &
                    xdotn/height**2 / (2*M_PI*rscale**2) /(1.d0-dexp(-0.5*(R/rscale**2))) * &
                    vol_cylinder
else
  weight_cattaneo = 0.
endif

return
end function weight_cattaneo

!==============================================================================

end module agn_geometry
