! Module providing routines implementing AGN feedback subgrid models.
!
! Paul Ricker (2/2013)

module agn_params

!==============================================================================

use fundamental_constants_module

implicit none

public

! Derived type to hold AGN particle properties

type agn_particle
  double precision :: m, mlast      ! current mass, mass at last outburst
  double precision :: x, y, z       ! Cartesian coordinates
  double precision :: vx, vy, vz    ! Cartesian velocity components
  double precision :: nx, ny, nz    ! angle cosines (for jet orientation)
end type agn_particle

! Switch constants used to choose accretion and feedback models

integer, parameter :: AGN_ACC_ALPHABONDI = 1, AGN_ACC_BETABONDI = 2, &
                      AGN_ACC_POPE = 3
integer, parameter :: AGN_FB_BUBBLE = 1, AGN_FB_JET = 2, AGN_FB_FIXBUB = 3

! Physical/mathematical constants

! These are now defined in fundamental_constants_module
! real, parameter :: G = 6.67428e-8, c = 2.99792458e10
! real, parameter :: m_p = 1.672621637e-24

double precision, parameter :: sigmaT = 6.6524e-25 ! for pure H

double precision, parameter :: kb = 1.38d-16
double precision, parameter :: Msun = 2.d33
double precision, parameter :: Mpc = 3.0857d24, kms = 1.d5

! Small quantities

double precision, parameter :: agn_small_rho = 1.E-35
double precision, parameter :: agn_small_e   = 1.E7

! Parameters for accretion models

double precision :: agn_ac_accradius         = 2.

!   Alpha-Bondi (Sijacki et al. 2008)

double precision :: agn_ac_ab_alpha          = 1.

!   Beta model (Booth & Schaye 2009)

double precision :: agn_ac_beta_nhstar       = 0.1

!   Pope model (Pope et al. 2007)

! Parameters for feedback models

double precision :: agn_fb_twomode_threshold = 0.01
double precision :: agn_fb_gasdepfac         = 0.1
double precision :: agn_fb_depradius         = 2.
double precision :: agn_fb_fbeff             = 0.1
double precision :: agn_fb_qsoeff            = 0.05

!   Bubbles (Sijacki)

double precision :: agn_fb_bub_R0            = 30.
double precision :: agn_fb_bub_E0            = 1.e55
double precision :: agn_fb_bub_rho0          = 1.e4
double precision :: agn_fb_bub_mecheff       = 0.2
double precision :: agn_fb_bub_delta         = 0.0001

!   Jets (Cattaneo & Teyssier)

double precision :: agn_fb_jet_massload      = 100.
double precision :: agn_fb_jet_cylradius     = 0.
double precision :: agn_fb_jet_cylheight     = 0.

!   Fixed bubbles

double precision :: agn_fb_fixbub_radius     = 30.

contains

!------------------------------------------------------------------------------

! agn_eddington_rate: Return Eddington accretion rate for a given mass black
! hole.

double precision function agn_eddington_rate(Mbh)

use fundamental_constants_module, only: m_p,c_light,Gconst

double precision, intent(in) :: Mbh

agn_eddington_rate =  4.d0*pi*Gconst*m_p/(sigmaT*c_light*agn_fb_fbeff) * Mbh

return
end function agn_eddington_rate

!==============================================================================

end module agn_params
