
subroutine f_rhs(num_eq, time, e_in, energy, rpar, ipar)

      use amrex_fort_module, only : rt => amrex_real
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

      use vode_aux_module       , only: z_vode, rho_vode, T_vode, ne_vode, i_vode, j_vode, k_vode

      integer, intent(in)             :: num_eq, ipar
      real(rt), intent(inout) :: e_in(num_eq)
      real(rt), intent(in   ) :: time
      real(rt), intent(in   ) :: rpar
      real(rt), intent(  out) :: energy

      real(rt), parameter :: compt_c = 1.01765467d-37, T_cmb = 2.725d0

      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
      integer :: j

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
        U = e_in(1) * e_to_cgs
      nh  = rho*XHYDROGEN/MPROTON

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         print *,'AT              ',i_vode,j_vode,k_vode
         call bl_pd_abort("TOO BIG TIME IN F_RHS")
      end if

      ! Get gas temperature and individual ionization species
      call iterate_ne(z_vode, U, T_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)

      ! Convert species to CGS units: 
      ne_vode = nh * ne_vode
      nh0   = nh * nh0
      nhp   = nh * nhp
      nhe0  = nh * nhe0
      nhep  = nh * nhep
      nhepp = nh * nhepp

      logT = dlog10(T_vode)
      if (logT .ge. TCOOLMAX) then ! Only free-free and Compton cooling are relevant
         lambda_ff = 1.42d-27 * dsqrt(T_vode) * (1.1d0 + 0.34d0*dexp(-(5.5d0 - logT)**2 / 3.0d0)) &
                              * (nhp + 4.0d0*nhepp)*ne_vode
         lambda_c  = compt_c*T_cmb**4 * ne_vode * (T_vode - T_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4

         energy  = (-lambda_ff -lambda_c) * heat_from_cgs/(1.0d0+z_vode)**4

         ! Convert to the actual term to be used in e_out = e_in + dt*energy
         energy  = energy / rho_vode * (1.0d0+z_vode)
         ne_vode = ne_vode / nh
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
                 bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode

      lambda_c = compt_c*T_cmb**4*ne_vode*(T_vode - T_cmb*(1.0d0+z_vode))*(1.0d0 + z_vode)**4   ! Compton cooling
      lambda = lambda + lambda_c

      ! Heating terms
      heat = nh0*eh0 + nhe0*ehe0 + nhep*ehep

      ! Convert back to code units
      ne_vode     = ne_vode / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z_vode)**4

      ! Convert to the actual term to be used in e_out = e_in + dt*energy
      a = 1.d0 / (1.d0 + z_vode)
      energy = energy / rho_vode / a

end subroutine f_rhs

subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer         , intent(in   ) :: neq, ml, mu, nrpd, ipar
  real(rt), intent(in   ) :: y(neq), rpar, t
  real(rt), intent(  out) :: pd(neq,neq)

  ! Should never get here, we are using a numerical Jacobian
  print *,'IN JAC ROUTINE'
  stop

end subroutine jac
