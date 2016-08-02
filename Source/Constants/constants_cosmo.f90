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

  implicit none

  ! Relation of our code units & CGS: 
  double precision, parameter :: M_unit = 1.9884d33     ! M_sun
  double precision, parameter :: L_unit = 3.0856776d24  ! Mpc
  double precision, parameter :: V_unit = 1.d5          ! km/s
  double precision, parameter :: T_unit = L_unit/V_unit ! time unit

  !
  ! Fundamental constants
  !
  double precision, parameter :: Gconst = 6.6743d-8 * &       ! Newton [g^-1*s^-2*cm^3]
                                           M_unit*T_unit**2/L_unit**3

  double precision, parameter :: k_B    = 1.3806504d-16 * &   ! Boltzmann [g*cm^2/s^2*K]
                                           T_unit**2/(M_unit*L_unit**2)

  double precision, parameter :: hbar   = 1.054571628d-27 * & ! Planck/2pi [g*cm^2/s]
                                           T_unit/(M_unit*L_unit**2)

  double precision, parameter :: n_A    = 6.02214179d23 * &   ! Avogadro's number [mol^-1]
                                           M_unit

  double precision, parameter :: m_proton = 1.6726231d-24 / & ! Proton mass [g]
                                           M_unit

  double precision, parameter :: c_light = 2.99792458d10 / &  ! Speed of light [cm/s] 
                                           V_unit

  double precision, parameter :: Hubble_const = 1.0d2 ! Hubble constant / h

  !
  ! Useful quantities and conversions
  !
  double precision, parameter :: mp_over_kb = m_proton/k_B

  double precision, parameter :: density_to_cgs = M_unit / L_unit**3

  ! For internal energy
  double precision, parameter :: e_to_cgs = V_unit**2 

  ! For source terms we convert [erg/(s*cm^3) = g/(s^3*cm)] into code units
  double precision, parameter :: heat_from_cgs = L_unit*(T_unit**3 / M_unit)


  ! Keep these for now for JF H/C files (not used elsewhere)
  ! energy per volume from code units to cgs
  double precision, parameter :: epv_to_cgs = e_to_cgs * density_to_cgs
  double precision, parameter :: Mpc_to_km = 3.08568025d19
  ! energy per volume per time from cosmo to cgs
  double precision, parameter :: epvpt_to_cgs = epv_to_cgs / T_unit


end module fundamental_constants_module
