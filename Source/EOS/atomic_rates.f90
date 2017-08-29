! *************************************************************************************
! Tabulates cooling and UV heating rates. 
!
!     Units are CGS, temperature is in K
!
!     Two sets of rates, as in:
!      1) Katz, Weinberg & Hernquist, 1996: Astrophysical Journal Supplement v. 105, p.19
!      2) Lukic et al. 2015: Monthly Notices of the Royal Astronomical Society, v. 446, p.3697
! 
!     NOTE: This is executed only once per run, and rates are ugly, thus 
!           execution efficiency is not important, but readability of 
!           the code is. -- Zarija
!
! *************************************************************************************

module atomic_rates_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! Routine which acts like a class constructor
  public  :: tabulate_rates, interp_to_this_z

  ! Photo- rates (from file)
  integer   , parameter          , private :: NCOOLFILE=301
  real(rt), dimension(NCOOLFILE), public :: lzr
  real(rt), dimension(NCOOLFILE), public :: rggh0, rgghe0, rgghep
  real(rt), dimension(NCOOLFILE), public :: reh0, rehe0, rehep

  ! Other rates (from equations)
  integer, parameter, public :: NCOOLTAB=2000
  real(rt), dimension(NCOOLTAB+1), public :: AlphaHp, AlphaHep, AlphaHepp, Alphad
  real(rt), dimension(NCOOLTAB+1), public :: GammaeH0, GammaeHe0, GammaeHep
  real(rt), dimension(NCOOLTAB+1), public :: BetaH0, BetaHe0, BetaHep, Betaff1, Betaff4
  real(rt), dimension(NCOOLTAB+1), public :: RecHp, RecHep, RecHepp

  real(rt), public, save :: this_z, ggh0, gghe0, gghep, eh0, ehe0, ehep
 
  real(rt), parameter, public :: TCOOLMIN = 0.0d0, TCOOLMAX = 9.0d0  ! in log10
  real(rt), parameter, public :: TCOOLMIN_R = 10.0d0**TCOOLMIN, TCOOLMAX_R = 10.0d0**TCOOLMAX
  real(rt), parameter, public :: deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB

  real(rt), parameter, public :: MPROTON = 1.6726231d-24, BOLTZMANN = 1.3806e-16

  ! Note that XHYDROGEN can be set by a call to set_xhydrogen which now
  ! lives in set_method_params.
  real(rt), public :: XHYDROGEN = 0.76d0
  real(rt), public :: YHELIUM   = 7.8947368421d-2  ! (1.0d0-XHYDROGEN)/(4.0d0*XHYDROGEN)

  contains

      subroutine tabulate_rates()
      integer :: i
      logical, parameter :: Katz96=.false.
      real(rt), parameter :: t3=1.0d3, t5=1.0d5, t6=1.0d6
      real(rt) :: t, U, E, y, sqrt_t, corr_term
      logical, save :: first=.true.

      !$OMP CRITICAL(TREECOOL_READ)
      if (first) then

         first = .false.

         ! Read in photoionization rates and heating from a file
         open(unit=11,file='TREECOOL_middle',status='old')
         do i = 1, NCOOLFILE
            read(11,*) lzr(i), rggh0(i), rgghe0(i), rgghep(i), &
                                reh0(i),  rehe0(i),  rehep(i)
         end do
         close(11)

         ! Initialize cooling tables
         t = 10.0d0**TCOOLMIN
         if (Katz96) then
            do i = 1, NCOOLTAB+1
               ! Rates are as in Katz et al. 1996
               sqrt_t = dsqrt(t)      
               corr_term    = 1.d0 / (1.0d0 + sqrt_t/dsqrt(t5))

               ! Recombination rates
               ! Alphad: dielectronic recombination rate of singly ioniozed helium
               Alphad(i)    = 1.90d-03/(t*sqrt_t) * dexp(-4.7d5/t) * (1.0d0+0.3d0*dexp(-9.4d4/t))
               AlphaHp(i)   = 8.40d-11/sqrt_t * (t/t3)**(-0.2d0) / (1.0d0 + (t/t6)**0.7d0)
               AlphaHep(i)  = 1.50d-10 * t**(-0.6353d0)
               AlphaHepp(i) = 3.36d-10/sqrt_t * (t/t3)**(-0.2d0) / (1.0d0 + (t/t6)**0.7d0)

               ! Collisional ionization rates
               GammaeH0(i)  = 5.85d-11*sqrt_t * dexp(-157809.1d0/t) * corr_term
               GammaeHe0(i) = 2.38d-11*sqrt_t * dexp(-285335.4d0/t) * corr_term
               GammaeHep(i) = 5.68d-12*sqrt_t * dexp(-631515.0d0/t) * corr_term

               ! Collisional ionization & excitation cooling rates
               BetaH0(i)  = 7.5d-19 * dexp(-118348.0d0/t) * corr_term + 2.171d-11*GammaeH0(i)
               BetaHe0(i) = 3.941d-11 * GammaeHe0(i)
               BetaHep(i) = 5.54d-17 * t**(-0.397d0) * dexp(-473638.0d0/t) * corr_term + &
                            8.715d-11 * GammaeHep(i)

               ! Recombination cooling rates
               RecHp(i)   = 1.036d-16 * t * AlphaHp(i)
               RecHep(i)  = 1.036d-16 * t * AlphaHep(i) + 6.526d-11 * Alphad(i)
               RecHepp(i) = 1.036d-16 * t * AlphaHepp(i)

               ! Free-free cooling rate
               Betaff1(i) = 1.42d-27 * sqrt_t * (1.1d0 + 0.34d0*dexp(-(5.5d0 - dlog10(t))**2 / 3.0d0))
               Betaff4(i) = Betaff1(i)

               t = t*10.0d0**deltaT
            enddo
         else
            do i = 1, NCOOLTAB+1
               ! Rates are as in Lukic et al.
               sqrt_t = dsqrt(t)

               ! Recombination rates
               ! Alphad: dielectronic recombination rate of singly ioniozed helium
               Alphad(i)    = 1.90d-03/(t*sqrt_t) * dexp(-4.7d5/t) * (1.0d0+0.3d0*dexp(-9.4d4/t))
               AlphaHp(i)   = 7.982d-11 / (dsqrt(t/3.148d0)* (1.0d0+dsqrt(t/3.148d0))**0.252 * &
                                           (1.0d0+dsqrt(t/7.036d5))**1.748)
               if (t .le. 1.0d6) then
                  AlphaHep(i)  = 3.294d-11 / (dsqrt(t/15.54d0)* (1.0d0+dsqrt(t/15.54d0))**0.309 * &
                                              (1.0d0+dsqrt(t/3.676d7))**1.691)
               else
                  AlphaHep(i)  = 9.356d-10 / (dsqrt(t/4.266d-2)* (1.0d0+dsqrt(t/4.266d-2))**0.2108 * &
                                              (1.0d0+dsqrt(t/4.677d6))**1.7892)
               endif
               AlphaHepp(i) = 1.891d-10 / (dsqrt(t/9.37d0)* (1.0d0+dsqrt(t/9.37d0))**0.2476 * &
                                           (1.0d0+dsqrt(t/2.774d6))**1.7524)

               ! Collisional ionization rates
               E = 13.6d0
               U = 1.16045d4*E/t
               GammaeH0(i)  = 2.91d-8*U**0.39*dexp(-U) / (0.232d0+U)
               E = 24.6d0
               U = 1.16045d4*E/t
               GammaeHe0(i) = 1.75d-8*U**0.35*dexp(-U) / (0.18d0+U)
               E = 54.4d0
               U = 1.16045d4*E/t
               GammaeHep(i) = 2.05d-9*(1.0d0+dsqrt(U))*U**0.25*dexp(-U) / (0.265d0+U)

               ! Collisional ionization & excitation cooling rates
               corr_term = 1.d0 / (1.0d0 + sqrt_t/dsqrt(5.0d7))
               y = dlog(t)
               if (t .le. 1.0d5) then
                  BetaH0(i)  = 1.0d-20 * dexp( 2.137913d2 - 1.139492d2*y + 2.506062d1*y**2 - &
                                               2.762755d0*y**3 + 1.515352d-1*y**4 - &
                                               3.290382d-3*y**5 - 1.18415d5/t )
               else
                  BetaH0(i)  = 1.0d-20 * dexp( 2.7125446d2 - 9.8019455d1*y + 1.400728d1*y**2 - &
                                               9.780842d-1*y**3 + 3.356289d-2*y**4 - &
                                               4.553323d-4*y**5 - 1.18415d5/t )
               endif
               BetaHe0(i) = 9.38d-22 * sqrt_t * dexp(-285335.4d0/t) * corr_term
               BetaHep(i) = (5.54d-17 * t**(-0.397d0) * dexp(-473638.0d0/t) + & 
                             4.85d-22 * sqrt_t * dexp(-631515.0d0/t) )*corr_term

               ! Recombination cooling rates
               RecHp(i)   = 2.851d-27 * sqrt_t * (5.914d0-0.5d0*dlog(t)+1.184d-2*t**(1.0d0/3.0d0))
               RecHep(i)  = 1.55d-26 * t**0.3647 + 1.24d-13/(t*sqrt_t) * dexp(-4.7d5/t) * & 
                                                                         (1.0d0+0.3d0*dexp(-9.4d4/t))
               RecHepp(i) = 1.14d-26 * sqrt_t * (6.607d0-0.5d0*dlog(t)+7.459d-3*t**(1.0d0/3.0d0))

               ! Free-free cooling rate
               if (t .le. 3.2d5) then
                  Betaff1(i) = 1.426d-27 * sqrt_t * (0.79464d0 + 0.1243d0*dlog10(t))
               else
                  Betaff1(i) = 1.426d-27 * sqrt_t * (2.13164d0 - 0.1240d0*dlog10(t))
               endif

               if (t/4.0d0 .le. 3.2d5) then
                  Betaff4(i) = 1.426d-27 * sqrt_t * 4.0d0*(0.79464d0 + 0.1243d0*dlog10(t))
               else
                  Betaff4(i) = 1.426d-27 * sqrt_t * 4.0d0*(2.13164d0 - 0.1240d0*dlog10(t))
               endif
               
               t = t*10.0d0**deltaT
            enddo
         endif  ! Katz rates

      end if  ! first_call
      !$OMP END CRITICAL(TREECOOL_READ)

      end subroutine tabulate_rates

      ! ****************************************************************************

      subroutine interp_to_this_z(z)

      real(rt), intent(in) :: z
      real(rt) :: lopz, fact
      integer :: i, j

      this_z = z
      lopz   = dlog10(1.0d0 + z)

      if (lopz .ge. lzr(NCOOLFILE)) then
         ggh0  = 0.0d0
         gghe0 = 0.0d0
         gghep = 0.0d0
         eh0   = 0.0d0
         ehe0  = 0.0d0
         ehep  = 0.0d0
         return
      endif

      if (lopz .le. lzr(1)) then
         j = 1
      else
         do i = 2, NCOOLFILE
            if (lopz .lt. lzr(i)) then
               j = i-1
               exit
            endif
         enddo
      endif

      fact  = (lopz-lzr(j))/(lzr(j+1)-lzr(j))

      ggh0  = rggh0(j)  + (rggh0(j+1)-rggh0(j))*fact
      gghe0 = rgghe0(j) + (rgghe0(j+1)-rgghe0(j))*fact
      gghep = rgghep(j) + (rgghep(j+1)-rgghep(j))*fact
      eh0   = reh0(j)   + (reh0(j+1)-reh0(j))*fact
      ehe0  = rehe0(j)  + (rehe0(j+1)-rehe0(j))*fact
      ehep  = rehep(j)  + (rehep(j+1)-rehep(j))*fact

      end subroutine interp_to_this_z

end module atomic_rates_module

! *************************************************************************************
! This must live outside of atomic_rates module so it can be called by the C++
! *************************************************************************************

subroutine fort_init_this_z(comoving_a) &
    bind(C, name="fort_init_this_z")

    use amrex_fort_module, only : rt => amrex_real
    use atomic_rates_module

    implicit none

    real(rt), intent(in   ) :: comoving_a
    real(rt)                :: z

    z = 1.d0/comoving_a - 1.d0
    call interp_to_this_z(z)

end subroutine fort_init_this_z
