#ifndef __f_rhs_struct__
#define __f_rhs_struct__


#include <atomic_rates_data.H>
#include <constants_cosmo.H>
#include <eos_hc.H>

struct RhsData {
	// Set once in init
	amrex::Real	uvb_density_A;
	amrex::Real	uvb_density_B;
	int         flash_h;
	int         flash_he;
	amrex::Real zhi_flash;
	amrex::Real zheii_flash;
	amrex::Real T_zhii;
	amrex::Real T_zheii;
	amrex::Real H_reion_z;
	amrex::Real He_reion_z;
	int         inhomogeneous_on;

	// Filled from Nyx.cpp
	amrex::Real gamma_minus_1;
	amrex::Real h_species;

	// Change at most twice per Nyx coarse time-step
	amrex::Real a;
    int         JH_vode;
    int         JHe_vode;
	amrex::Real* T_vode;
	amrex::Real* ne_vode;
	amrex::Real* rho_vode;
	amrex::Real* rho_init_vode;
	amrex::Real* rho_src_vode;
	amrex::Real* rhoe_src_vode;
	amrex::Real* e_src_vode;
};

extern RhsData* f_rhs_data;

void
ode_eos_setup(RhsData* f_rhs_data, amrex::Real gamma_minus_1, amrex::Real h_species)
{

	f_rhs_data->gamma_minus_1 = gamma_minus_1;
	f_rhs_data->h_species = h_species;
	ParmParse pp_nyx("nyx");
	pp_nyx.query("inhomo_reion", inhomogeneous_on);

    f_rhs_data->uvb_density_A=1.0;
	f_rhs_data->uvb_density_B=0.0;
	f_rhs_data->zhi_flash=-1.0;
	f_rhs_data->zheii_flash=-1.0;
	f_rhs_data->T_zhi=0.0;
	f_rhs_data->T_zheii=0.0;
	f_rhs_data->flash_h=0;
	f_rhs_data->flash_he=0;
	f_rhs_data->inhomogeneous_on=0;
	pp_nyx.query("uvb_density_A"            , f_rhs_data->uvb_density_A);
	pp_nyx.query("uvb_density_B"            , f_rhs_data->uvb_density_B);
	pp_nyx.query("reionization_zHI_flash"   , f_rhs_data->zhi_flash);
	pp_nyx.query("reionization_zHeII_flash" , f_rhs_data->zheii_flash);
	pp_nyx.query("reionization_T_zHI"       , f_rhs_data->T_zhi);
	pp_nyx.query("reionization_T_zHeII"     , f_rhs_data->T_zheii);

	// Set options in reion_aux_module
    //   Hydrogen reionization
    if (f_rhs_data->zhi_flash > 0.0)
	{
        if (f_rhs_data->inhomogeneous_on > 0)
		{
			f_rhs_data->flash_h = 0;
		}
		else
		{
			f_rhs_data->flash_h = 1;
		}
	}
    else
	{
        f_rhs_data->flash_h = 0;
	}

	//   Helium reionization
	if (f_rhs_data->zheii_flash > 0.0)
		f_rhs_data->flash_he = 1;
	else
		f_rhs_data->flash_he = 0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_initialize(/*Add ijk, idx*/,/*Add Array4 for state4, diag_eos4*/)
{

    // Flash reionization?
    if ((flash_h == 1) && (z > zhi_flash))
	   JH_vode = 0;
    else
		JH_vode = 1;

    if ((flash_he == 1) && (z > zheii_flash))
		JHe_vode = 0;
    else
		JHe_vode = 1;

    if (flash_h ) H_reion_z  = zhi_flash;
    if (flash_he) He_reion_z = zheii_flash;
	/*
    if (inhomogeneous_on)
	{
	   H_reion_z = diag_eos(i,j,k,ZHI_COMP);
       if (z > H_reion_z)
		   JH_vode = 0;
       else
		   JH_vode = 1;
	}
	*/
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_finalize(amrex::Real& e_out, RhsData* f_rhs_data/*Add ijk, idx*/,/*Add Array4 for state4, diag_eos4*/)
{
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;
    amrex::Real rho, U, a, mu;
    amrex::Real uvb_density_A, uvb_density_B, mean_rhob, JH_vode, JHe_vode, gamma_minus_1;
    const amrex::Real h_species_in=.76;
    uvb_density_A=1.0;
	uvb_density_B=0.0;
	gamma_minus_1=2.0/3.0;
	/*
    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    if (inhomogeneous_on) then
       STOP "Do not currently support inhomogenous_on with box"
       !H_reion_z = diag_eos(i,j,k,ZHI_COMP)
       if (z .gt. H_reion_z) then
          JH_vode = 0
       else
          JH_vode = 1
       endif
    endif*/
    JH_vode = 1.0;
	JHe_vode = 1.0;

    T_vode=rpar[0];
    ne_vode=rpar[1];
    rho_vode=rpar[2];
    z_vode=rpar[3];
	a=1/(z_vode+1.0);
	if(e_out < 0.e0)
	{
		T_vode  = 10.0;
		ne_vode = 0.0;
		mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+ne_vode);
		e_out  = T_vode / (gamma_minus_1 * mp_over_kb * mu);
	}

    nyx_eos_T_given_Re_device(gamma_minus_1, h_species_in, JH_vode, JHe_vode, &T_vode, &ne_vode,
                                               rho_vode, e_out, a);
	/*
      ! Instanteneous heating from reionization:
      T_H = 0.0d0
      if (inhomogeneous_on .or. flash_h) then
         if ((H_reion_z  .lt. z) .and. (H_reion_z  .ge. z_end)) T_H  = (1.0d0 - species(2))*max((T_zhi-T_vode), 0.0d0)
      endif

      T_He = 0.0d0
      if (flash_he) then
         if ((He_reion_z .lt. z) .and. (He_reion_z .ge. z_end)) T_He = (1.0d0 - species(5))*max((T_zheii-T_vode), 0.0d0)
      endif

      if ((T_H .gt. 0.0d0) .or. (T_He .gt. 0.0d0)) then
         T_vode = T_vode + T_H + T_He                            ! For simplicity, we assume
         ne_vode = 1.0d0 + YHELIUM                              !    completely ionized medium at
         if (T_He .gt. 0.0d0) ne_vode = ne_vode + YHELIUM        !    this point.  It's a very minor
         mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_vode)   !    detail compared to the overall approximation.
         e_out  = T_vode / (gamma_minus_1 * mp_over_kB * mu)
         call nyx_eos_T_given_Re_device(JH_vode, JHe_vode, T_vode, ne_vode, rho_vode, e_out(1), a, species)
      endif
    */
    rpar[0]=T_vode;
    rpar[1]=ne_vode;
    rpar[2]=rho_vode;
}

subroutine f_rhs_split(num_eq, time, y_in, yp_out, rpar, ipar)

      use amrex_error_module, only : amrex_abort
      use amrex_fort_module, only : rt => amrex_real
      use fundamental_constants_module, only: e_to_cgs, density_to_cgs, & 
                                              heat_from_cgs
      use eos_module, only: iterate_ne
      use atomic_rates_module, ONLY: TCOOLMIN, TCOOLMAX, NCOOLTAB, deltaT, &
                                     MPROTON, XHYDROGEN, &
                                     uvb_density_A, uvb_density_B, mean_rhob, &
                                     BetaH0, BetaHe0, BetaHep, Betaff1, Betaff4, &
                                     RecHp, RecHep, RecHepp, &
                                     eh0, ehe0, ehep

      use vode_aux_module       , only: z_vode, rho_vode, T_vode, ne_vode, &
                                        JH_vode, JHe_vode, i_vode, j_vode, k_vode, fn_vode, NR_vode, &
                                        rho_init_vode, e_src_vode, rho_src_vode

      integer, intent(in)             :: num_eq, ipar
      real(rt), intent(inout) :: y_in(num_eq)
      real(rt), intent(in   ) :: time
      real(rt), intent(in   ) :: rpar
      real(rt), intent(  out) :: yp_out(num_eq)

      real(rt), parameter :: compt_c = 1.01765467d-37, T_cmb = 2.725d0


      real(rt) :: e_in(1)
      real(rt) :: energy
      real(rt) :: rho_in
      real(rt) :: logT, tmp, fhi, flo
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp
      real(rt) :: lambda_c, lambda_ff, lambda, heat
      real(rt) :: rho, U, a, rho_heat
      real(rt) :: nh, nh0, nhp, nhe0, nhep, nhepp
#ifdef AMREX_USE_CUDA_FORTRAN
     attributes(managed) :: nh, nh0, nhp, nhe0, nhep, nhepp
#endif
      integer :: j

      fn_vode=fn_vode+1;

      e_in = y_in(1)
      rho_vode = y_in(2)
!      rho_vode = rho_init_vode + time * rho_src_vode

      if (e_in(1) .lt. 0.d0) &
         e_in(1) = tiny(e_in(1))

     ! Converts from code units to CGS
      rho = rho_vode * density_to_cgs * (1.0d0+z_vode)**3
        U = e_in(1) * e_to_cgs
      nh  = rho*XHYDROGEN/MPROTON

      if (time .gt. 1) then
         print *,'TIME INTO F_RHS ',time
         print *,'AT              ',i_vode,j_vode,k_vode
         call amrex_abort("TOO BIG TIME IN F_RHS")
      end if

      ! Get gas temperature and individual ionization species
      ! testing different memory structures
!      NR_vode=0
      call iterate_ne(JH_vode, JHe_vode, z_vode, U, T_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp)

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
         yp_out(1) = energy + e_src_vode
         yp_out(2) = rho_src_vode
         return
      end if

      ! Temperature floor
      if (logT .le. TCOOLMIN)  logT = TCOOLMIN + 0.5d0*deltaT

      ! Interpolate rates
      tmp = (logT-TCOOLMIN)/deltaT
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1

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
      heat = JH_vode*nh0*eh0 + JH_vode*nhe0*ehe0 + JHe_vode*nhep*ehep
      rho_heat = uvb_density_A * (rho_vode/mean_rhob)**uvb_density_B
      heat = rho_heat*heat

      ! Convert back to code units
      ne_vode     = ne_vode / nh
      energy = (heat - lambda)*heat_from_cgs/(1.0d0+z_vode)**4

      ! Convert to the actual term to be used in e_out = e_in + dt*energy
      a = 1.d0 / (1.d0 + z_vode)
      energy = (energy) / rho_vode / a

      yp_out(1) = energy + e_src_vode
      yp_out(2) = rho_src_vode


end subroutine f_rhs_split
#endif
