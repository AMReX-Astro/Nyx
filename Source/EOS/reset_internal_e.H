#include "Nyx.H"
#include "IndexDefines.H"
#include <constants_cosmo.H>
#ifdef HEATCOOL
#include <eos_hc.H>
#include <atomic_rates_data.H>
#else
#include <eos_adiabatic.H>
#endif

using namespace amrex;

#ifndef NO_HYDRO
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void reset_internal_e (const int i,
                       const int j,
                       const int k,
                       amrex::Array4<amrex::Real> const& u,
                       amrex::Array4<amrex::Real> const& d,
                       amrex::Array4<amrex::Real> const& r,
                       AtomicRates* atomic_rates,
                       const amrex::Real comoving_a,
                       const amrex::Real gamma_minus_1,
                       const amrex::Real h_species,
                       const amrex::Real small_temp,
                       const int interp)
{
    Real rhoInv = 1.0 / u(i,j,k,Density_comp);
    Real Up     = u(i,j,k,Xmom_comp) * rhoInv;
    Real Vp     = u(i,j,k,Ymom_comp) * rhoInv;
    Real Wp     = u(i,j,k,Zmom_comp) * rhoInv;
    Real ke     = 0.5 * u(i,j,k,Density_comp) * (Up*Up + Vp*Vp + Wp*Wp);

    Real rho_eint = u(i,j,k,Eden_comp) - ke;
    Real dummy_pres = 1e200;
    Real eint_new = -1e200;

    // Reset (e from e) if it's greater than 0.01% of big E.
    if (rho_eint > 0.0 && rho_eint / u(i,j,k,Eden_comp) > 1.0e-6 && interp == 0)
    {
        // Create reset source so u(i,j,k,Eint_comp) = u(i,j,k,Eint_comp) + r(i,j,k) = rho_eint
        r(i,j,k) = rho_eint - u(i,j,k,Eint_comp);
        u(i,j,k,Eint_comp) = rho_eint;
    }
    // If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
    else if (u(i,j,k,Eint_comp) > 0.0)
    {
        // e is not updated, so reset source is zero
        r(i,j,k) += 0.0;
        u(i,j,k,Eden_comp) = u(i,j,k,Eint_comp) + ke;
    }
    // If not resetting and little e is negative ...
    else if (u(i,j,k,Eint_comp) <= 0.0)
    {
        nyx_eos_given_RT(atomic_rates, gamma_minus_1, h_species, &eint_new, &dummy_pres, u(i,j,k,Density_comp), small_temp, 
                         d(i,j,k,Ne_comp),comoving_a);

        // Create reset source so u(i,j,k,Eint_comp) = u(i,j,k,Eint_comp) + r(i,j,k) = u(i,j,k,Density_comp) * eint_new
        r(i,j,k) = u(i,j,k,Density_comp) *  eint_new - u(i,j,k,Eint_comp);
        u(i,j,k,Eint_comp) = u(i,j,k,Density_comp) *  eint_new;

        u(i,j,k,Eden_comp) = u(i,j,k,Eint_comp) + ke;
    }

    // Scale reset source by 1/rho so src is in terms of e
    //           r(i,j,k) = r(i,j,k) / u(i,j,k,Density_comp)
}

#endif
