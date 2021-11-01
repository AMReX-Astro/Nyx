#ifndef __f_rhs__
#define __f_rhs__

#include <atomic_rates_data.H>
#include <constants_cosmo.H>
#include <eos_hc.H>
#include <AMReX_MultiFab.H>

static amrex::Vector<void*> ptr_lst;
//static amrex::Arena* Managed_Arena;

static void* sunalloc(size_t mem_size)
{
  amrex::MultiFab::updateMemUsage ("Sunalloc", mem_size, nullptr);
  amrex::MultiFab::updateMemUsage ("All", mem_size, nullptr);
  void * ptr = (void*) The_Arena()->alloc(mem_size);
  ptr_lst.push_back(ptr);
  return ptr;
}

static void sunfree(void* ptr)
{
  size_t mem_size = dynamic_cast<CArena*>(The_Arena())->sizeOf(ptr);
  ptr_lst.erase(std::remove_if(ptr_lst.begin(), ptr_lst.end(), [ptr](void* x) { return x == ptr; }));
  The_Arena()->free(ptr);
  amrex::MultiFab::updateMemUsage ("Sunalloc", -mem_size, nullptr);
  amrex::MultiFab::updateMemUsage ("All", -mem_size, nullptr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ode_eos_finalize(amrex::Real& e_out, amrex::Real* rpar, int neq, AtomicRates* atomic_rates, amrex::Real h_species_in)
{
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;
    amrex::Real rho, U, a, mu;
    amrex::Real uvb_density_A, uvb_density_B, mean_rhob, JH_vode, JHe_vode, gamma_minus_1;
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
                Real YHELIUM = (1.0-h_species_in)/(4.0*h_species_in);
                T_vode  = 10.0;
                ne_vode = 0.0;
                mu     = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+ne_vode);
                e_out  = T_vode / (gamma_minus_1 * mp_over_kb * mu);
        }

    nyx_eos_T_given_Re_device(atomic_rates, gamma_minus_1, h_species_in, JH_vode, JHe_vode, &T_vode, &ne_vode,
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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int
f_rhs_rpar(amrex::Real const time, amrex::Real& e_in, amrex::Real& energy, amrex::Real* rpar, AtomicRates* atomic_rates, amrex::Real h_species_in)
{

    const amrex::Real compt_c = 1.01765467e-37;
    const amrex::Real T_cmb = 2.725e0;
    amrex::Real deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB;

    amrex::Real logT, tmp, fhi, flo;
    amrex::Real ggh0, gghe0, gghep, eh0, ehe0, ehep;
    amrex::Real ahp, ahep, ahepp, ad, geh0, gehe0, gehep;
    amrex::Real bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp;
    amrex::Real lambda_c, lambda_ff, lambda, heat;
    amrex::Real rho, U, a, rho_heat;
    amrex::Real nh, nh0, nhp, nhe0, nhep, nhepp;
    int j;
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;
    amrex::Real uvb_density_A, uvb_density_B, JH_vode, JHe_vode, gamma_minus_1;
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

    if(e_in<=0 || amrex::isnan(e_in)) e_in=std::numeric_limits<amrex::Real>::min();

    // Converts from code units to CGS
    rho = rho_vode * density_to_cgs * (1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode));
    U = e_in * e_to_cgs;
    nh  = rho*h_species_in/MPROTON;

    // Get gas temperature and individual ionization species
    // testing different memory structures
    iterate_ne_device(atomic_rates, JH_vode, JHe_vode, z_vode, U, T_vode, nh, ne_vode, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1, h_species_in);
    // Convert species to CGS units: 
    ne_vode = nh * ne_vode;
    nh0   = nh * nh0;
    nhp   = nh * nhp;
    nhe0  = nh * nhe0;
    nhep  = nh * nhep;
    nhepp = nh * nhepp;

    logT = std::log10(T_vode);
    if (logT >= TCOOLMAX) // Only free-free and Compton cooling are relevant
    {
        lambda_ff = 1.42e-27 * std::sqrt(T_vode) * (1.1e0 + 0.34e0*std::exp(-(5.5e0 - logT)*(5.5e0 - logT) / 3.0e0))
            * (nhp + 4.0e0*nhepp)*ne_vode;
        lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));

        energy  = (-lambda_ff -lambda_c) * heat_from_cgs/((1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode)));

        // Convert to the actual term to be used in e_out = e_in + dt*energy
        energy  = energy / rho_vode * (1.0e0+amrex::Math::abs(z_vode));
        ne_vode = ne_vode / nh;
        rpar[0]=T_vode;
        rpar[1]=ne_vode;
        rpar[2]=rho_vode;
        return 0;
    }

    // Temperature floor
    if (logT <= TCOOLMIN)  logT = TCOOLMIN + 0.5e0*deltaT;

    // Interpolate rates
    // Interpolate rates
    tmp = (logT-TCOOLMIN)/deltaT;
    j = amrex::Math::floor(tmp); // replacing int(tmp)
    fhi = tmp - j;
    flo = 1.0e0 - fhi;
    j = j + 1; // F90 arrays start with 1, using Array1D<Real,1,size>

    bh0   = flo*atomic_rates->BetaH0   (j) + fhi*atomic_rates->BetaH0   (j+1);
    bhe0  = flo*atomic_rates->BetaHe0  (j) + fhi*atomic_rates->BetaHe0  (j+1);
    bhep  = flo*atomic_rates->BetaHep  (j) + fhi*atomic_rates->BetaHep  (j+1);
    bff1  = flo*atomic_rates->Betaff1  (j) + fhi*atomic_rates->Betaff1  (j+1);
    bff4  = flo*atomic_rates->Betaff4  (j) + fhi*atomic_rates->Betaff4  (j+1);
    rhp   = flo*atomic_rates->RecHp    (j) + fhi*atomic_rates->RecHp    (j+1);
    rhep  = flo*atomic_rates->RecHep   (j) + fhi*atomic_rates->RecHep   (j+1);
    rhepp = flo*atomic_rates->RecHepp  (j) + fhi*atomic_rates->RecHepp  (j+1);

    // Cooling: 
    lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + 
               rhp*nhp + rhep*nhep + rhepp*nhepp + 
               bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode;

    lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));
    lambda = lambda + lambda_c;

    interp_to_this_z(atomic_rates, z_vode, ggh0, gghe0, gghep, eh0, ehe0, ehep);
    // Heating terms
    heat = JH_vode*nh0*eh0 + JH_vode*nhe0*ehe0 + JHe_vode*nhep*ehep;
    
    rho_heat = uvb_density_A * std::pow((rho_vode/atomic_rates->mean_rhob),uvb_density_B);
    heat = rho_heat*heat;

    // Convert back to code units
    ne_vode     = ne_vode / nh;
    energy  = (heat - lambda)*heat_from_cgs / ((1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode)));

    // Convert to the actual term to be used in e_out = e_in + dt*energy
    a = 1.e0 / (1.e0 + amrex::Math::abs(z_vode));
    energy = energy / rho_vode / a;
    rpar[0]=T_vode;
    rpar[1]=ne_vode;
    rpar[2]=rho_vode;
    return 0;
}
#endif
