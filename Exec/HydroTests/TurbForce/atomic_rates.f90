! *************************************************************************************
! Dummy module for tabulates cooling and UV heating rates.
! *************************************************************************************

module atomic_rates_module

  implicit none

  ! Routine which acts like a class constructor
  public  :: tabulate_rates, interp_to_this_z

 ! Photo- rates (from file)
  integer         , parameter, private :: NCOOLFILE=1
  double precision, dimension(NCOOLFILE), public :: lzr
  double precision, dimension(NCOOLFILE), public :: rggh0, rgghe0, rgghep
  double precision, dimension(NCOOLFILE), public :: reh0, rehe0, rehep

  ! Other rates (from equations)
  integer, parameter, public :: NCOOLTAB=1
  double precision, dimension(NCOOLTAB+1), public :: AlphaHp, AlphaHep, AlphaHepp, Alphad
  double precision, dimension(NCOOLTAB+1), public :: GammaeH0, GammaeHe0, GammaeHep
  double precision, dimension(NCOOLTAB+1), public :: BetaH0, BetaHe0, BetaHep, Betaff1, Betaff4
  double precision, dimension(NCOOLTAB+1), public :: RecHp, RecHep, RecHepp

  double precision, public, save :: this_z, ggh0, gghe0, gghep, eh0, ehe0, ehep
 
  double precision, parameter, public :: MPROTON = 1.6726231d-24, BOLTZMANN = 1.3806e-16

  double precision, parameter, public :: TCOOLMIN = 0.0d0, TCOOLMAX = 9.0d0  ! in log10
  double precision, parameter, public :: deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB

  ! Note that XHYDROGEN can be set by a call to set_xhydrogen which now
  ! lives in set_method_params.
  double precision, public :: XHYDROGEN = 0.76d0
  double precision, public :: YHELIUM   = 7.8947368421d-2  ! (1.0d0-XHYDROGEN)/(4.0d0*XHYDROGEN)

  contains

      subroutine tabulate_rates()

      lzr(:) = 0.0d0
      
      rggh0(:) = 0.0d0 
      rgghe0(:) = 0.0d0 
      rgghep(:) = 0.0d0                            
      reh0(:) = 0.0d0  
      rehe0(:) = 0.0d0  
      rehep(:) = 0.0d0

      AlphaHp(:) = 0.0d0
      AlphaHep(:) = 0.0d0
      AlphaHepp(:) = 0.0d0
      Alphad(:) = 0.0d0

      GammaeH0(:) = 0.0d0
      GammaeHe0(:) = 0.0d0
      GammaeHep(:) = 0.0d0
      
      BetaH0(:) = 0.0d0
      BetaHe0(:) = 0.0d0
      BetaHep(:) = 0.0d0
      Betaff1(:) = 0.0d0
      Betaff4(:) = 0.0d0

      RecHp(:) = 0.0d0
      RecHep(:) = 0.0d0
      RecHepp(:) = 0.0d0

      end subroutine tabulate_rates

      ! ****************************************************************************

      subroutine interp_to_this_z(z)

      double precision, intent(in) :: z

      this_z = z

      ggh0  = 0.0d0
      gghe0 = 0.0d0
      gghep = 0.0d0
      eh0   = 0.0d0
      ehe0  = 0.0d0
      ehep  = 0.0d0
 
      end subroutine interp_to_this_z


end module atomic_rates_module

! *************************************************************************************
! This must live outside of atomic_rates module so it can be called by the C++
! *************************************************************************************

subroutine init_this_z(comoving_a)

use atomic_rates_module

double precision, intent(in   ) :: comoving_a
double precision                :: z

    z = 1.d0/comoving_a - 1.d0
    call interp_to_this_z(z)

end subroutine init_this_z
