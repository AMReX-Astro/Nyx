!
! Nyx code units are defined as: 
!      Length [Mpc]
!      Velocity [km/s]
!      Mass [M_sun]
!      (+ Kelvin for temperature)
!
! Fundamental constants are taken from Phys. Rev. D 86 (2012), see 
!          pdg.lbl.gov/2012/reviews/rpp2012-rev-phys-constants.pdf
!          pdg.lbl.gov/2012/reviews/rpp2012-rev-astrophysical-constants.pdf
!

module fundamental_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! Pi
  real(rt), parameter :: pi = 3.141592653589793238d0

  ! Relation of our code units & CGS: 
  real(rt), parameter :: M_unit = 1.9884d33     ! M_sun
  real(rt), parameter :: L_unit = 3.0856776d24  ! Mpc
  real(rt), parameter :: V_unit = 1.d5          ! km/s
  real(rt), parameter :: T_unit = L_unit/V_unit ! time unit

  !
  ! Fundamental constants
  !
  real(rt), parameter :: Gconst = 6.6743e-8_rt * &            ! Newton [g^-1*s^-2*cm^3]
                                  M_unit*T_unit**2/L_unit**3

  real(rt), parameter :: k_B    = 1.3806504e-16_rt * &        ! Boltzmann [g*cm^2/s^2*K]
                                  T_unit**2/(M_unit*L_unit**2)

  real(rt), parameter :: hbar   = 1.054571628e-27_rt * &      ! Planck/2pi [g*cm^2/s]
                                  T_unit/(M_unit*L_unit**2)

  real(rt), parameter :: n_A    = 6.02214179e23_rt * M_unit   ! Avogadro's number [mol^-1]

  real(rt), parameter :: m_proton = 1.6726231e-24_rt / M_unit ! Proton mass [g]

  real(rt), parameter :: c_light = 2.99792458e10_rt / V_unit  ! Speed of light [cm/s] 

  real(rt), parameter :: Hubble_const = 100._rt              ! Hubble constant / h

  !
  ! Useful quantities and conversions
  !
  real(rt), parameter :: mp_over_kb = m_proton/k_B

  real(rt), parameter :: density_to_cgs = M_unit / L_unit**3

  ! For internal energy
  real(rt), parameter :: e_to_cgs = V_unit**2 

  ! For source terms we convert [erg/(s*cm^3) = g/(s^3*cm)] into code units
  real(rt), parameter :: heat_from_cgs = L_unit*(T_unit**3 / M_unit)


  ! Keep these for now for JF H/C files (not used elsewhere)
  ! energy per volume from code units to cgs
  real(rt), parameter :: epv_to_cgs = e_to_cgs * density_to_cgs
  real(rt), parameter :: Mpc_to_km = 3.08568025d19
  ! energy per volume per time from cosmo to cgs
  real(rt), parameter :: epvpt_to_cgs = epv_to_cgs / T_unit


end module fundamental_constants_module
