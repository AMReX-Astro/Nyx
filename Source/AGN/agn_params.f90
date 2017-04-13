! Module providing routines implementing AGN feedback subgrid models.
!
! Paul Ricker (2/2013)

module agn_params

!==============================================================================

    use fundamental_constants_module

    use amrex_fort_module, only : rt => amrex_real

    implicit none

public

! Derived type to hold AGN particle properties

type agn_particle
  real(rt) :: m, mlast      ! current mass, mass at last outburst
  real(rt) :: x, y, z       ! Cartesian coordinates
  real(rt) :: vx, vy, vz    ! Cartesian velocity components
  real(rt) :: nx, ny, nz    ! angle cosines (for jet orientation)
end type agn_particle

! Switch constants used to choose accretion and feedback models

integer, parameter :: AGN_ACC_ALPHABONDI = 1, AGN_ACC_BETABONDI = 2, &
                      AGN_ACC_POPE = 3
integer, parameter :: AGN_FB_BUBBLE = 1, AGN_FB_JET = 2, AGN_FB_FIXBUB = 3

! Physical/mathematical constants

! These are now defined in fundamental_constants_module
real(rt), parameter :: G = 6.67428e-8, c = 2.99792458e10
real(rt), parameter :: m_p = 1.672621637e-24

real(rt), parameter :: sigmaT = 6.6524e-25 ! for pure H

real(rt), parameter :: kb = 1.38d-16
real(rt), parameter :: Msun = 2.d33
real(rt), parameter :: Mpc = 3.0857d24, kms = 1.d5

! Small quantities

real(rt), parameter :: agn_small_rho = 1.E-35
real(rt), parameter :: agn_small_e   = 1.E7

! Parameters for accretion models

real(rt) :: agn_ac_accradius         = 2.

!   Alpha-Bondi (Sijacki et al. 2008)

real(rt) :: agn_ac_ab_alpha          = 1.

!   Beta model (Booth & Schaye 2009)

real(rt) :: agn_ac_beta_nhstar       = 0.1

!   Pope model (Pope et al. 2007)

! Parameters for feedback models

real(rt) :: agn_fb_twomode_threshold = 0.01
real(rt) :: agn_fb_gasdepfac         = 0.1
real(rt) :: agn_fb_depradius         = 2.
real(rt) :: agn_fb_fbeff             = 0.1
real(rt) :: agn_fb_qsoeff            = 0.05

!   Bubbles (Sijacki)

real(rt) :: agn_fb_bub_R0            = 30.
real(rt) :: agn_fb_bub_E0            = 1.e55
real(rt) :: agn_fb_bub_rho0          = 1.e4
real(rt) :: agn_fb_bub_mecheff       = 0.2
real(rt) :: agn_fb_bub_delta         = 0.0001

!   Jets (Cattaneo & Teyssier)

real(rt) :: agn_fb_jet_massload      = 100.
real(rt) :: agn_fb_jet_cylradius     = 0.
real(rt) :: agn_fb_jet_cylheight     = 0.

!   Fixed bubbles

real(rt) :: agn_fb_fixbub_radius     = 30.

contains

!------------------------------------------------------------------------------

! agn_eddington_rate: Return Eddington accretion rate for a given mass black
! hole.

real(rt) function agn_eddington_rate(Mbh)

use fundamental_constants_module, only: m_p,c_light,Gconst

real(rt), intent(in) :: Mbh

agn_eddington_rate =  4.d0*pi*Gconst*m_p/(sigmaT*c_light*agn_fb_fbeff) * Mbh

return
end function agn_eddington_rate

!==============================================================================

end module agn_params
