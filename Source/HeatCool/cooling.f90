! Calculates cooling (H & He) + UV heating rates. 
!
!     Working units are CGS here, temperature is in K 
!

module heating_cooling_module

  implicit none

  public :: hc_rates

  contains

      subroutine hc_rates(z, R_in, e_in, t, ne, energy, prnt_d)

      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, & 
                                              heat_from_cgs
      use eos_module, only: iterate_ne
      use atomic_rates_module, ONLY: TCOOLMIN, TCOOLMAX, NCOOLTAB, deltaT, &
                                     MPROTON, XHYDROGEN, &
                                     AlphaHp, AlphaHep, AlphaHepp, Alphad, &
                                     GammaeH0, GammaeHe0, GammaeHep, &
                                     BetaH0, BetaHe0, BetaHep, Betaff1, Betaff4, &
                                     RecHp, RecHep, RecHepp, &
                                     eh0, ehe0, ehep

      double precision, intent(in   ) :: z, R_in, e_in
      double precision, intent(inout) :: t, ne
      double precision, intent(  out) :: energy
      logical, intent(in)             :: prnt_d ! for diagnostics print

      double precision, parameter :: compt_c = 1.01765467d-37, T_cmb = 2.725d0

      double precision :: logT, tmp, fhi, flo
      double precision :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      double precision :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      double precision :: lambda_c, lambda_ff, lambda, heat
      double precision :: rho, U
      double precision :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer :: j


     ! Converts from code units to CGS
      rho = R_in * density_to_cgs * (1.0d0+z)**3
        U = e_in * e_to_cgs
      nh  = rho*XHYDROGEN/MPROTON

      ! Get gas temperature and individual ionization species
      call iterate_ne(z, U, t, nh, ne, nh0, nhp, nhe0, nhep, nhepp)

      ! Convert species to CGS units: 
      ne    = nh * ne
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logT = dlog10(t)
      if (logT .ge. TCOOLMAX) then ! Only free-free and Compton cooling are relevant
         lambda_ff = 1.42d-27 * dsqrt(t) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logT)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne
         lambda_c  = compt_c*T_cmb**4*ne*(t - T_cmb*(1.0d0+z))*(1.0d0 + z)**4

         energy = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z)**4
         ne     = ne / nh
         return
      endif

      ! Temperature floor
      if (logT .le. TCOOLMIN) logT = TCOOLMIN + 0.5d0*deltaT

      ! Interpolate rates
      tmp = (logT-TCOOLMIN)/deltaT
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1

      ahp   = flo*AlphaHp  (j) + fhi*AlphaHp  (j+1)
      ahep  = flo*AlphaHep (j) + fhi*AlphaHep (j+1)
      ahepp = flo*AlphaHepp(j) + fhi*AlphaHepp(j+1)
      ad    = flo*Alphad   (j) + fhi*Alphad   (j+1)
      geh0  = flo*GammaeH0 (j) + fhi*GammaeH0 (j+1)
      gehe0 = flo*GammaeHe0(j) + fhi*GammaeHe0(j+1)
      gehep = flo*GammaeHep(j) + fhi*GammaeHep(j+1)
      bh0   = flo*BetaH0   (j) + fhi*BetaH0   (j+1)
      bhe0  = flo*BetaHe0  (j) + fhi*BetaHe0  (j+1)
      bhep  = flo*BetaHep  (j) + fhi*BetaHep  (j+1)
      bff1  = flo*Betaff1  (j) + fhi*Betaff1  (j+1)
      bff4  = flo*Betaff4  (j) + fhi*Betaff4  (j+1)
      rhp   = flo*RecHp    (j) + fhi*RecHp    (j+1)
      rhep  = flo*RecHep   (j) + fhi*RecHep   (j+1)
      rhepp = flo*RecHepp  (j) + fhi*RecHepp  (j+1)

      ! Cooling: 
      lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + &
                 rhp*nhp + rhep*nhep + rhepp*nhepp + &
                 bff1*(nhp+nhep) + bff4*nhepp ) * ne

      lambda_c = compt_c*T_cmb**4*ne*(t - T_cmb*(1.0d0+z))*(1.0d0 + z)**4   ! Compton cooling
      lambda = lambda + lambda_c

      ! Heating terms
      heat = nh0*eh0 + nhe0*ehe0 + nhep*ehep

      ! Convert back to code units
      ne     = ne / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z)**4

      end subroutine hc_rates

end module heating_cooling_module
