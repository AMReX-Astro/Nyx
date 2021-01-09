#include <Nyx.H>
#include <Hydro.H>
#include <constants_cosmo.H>

using namespace amrex;

using std::string;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
floor_density(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& state,
  AtomicRates* atomic_rates,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  amrex::Real const a_new,
  const amrex::Real gamma_minus_1,
  const amrex::Real small_dens,
  const amrex::Real small_temp,
  const amrex::Real h_species_in)
{
    //
    //  Reset negative density to small_dens, set (rho e) and (rho E) from small_temp  and zero out momenta
    //
    if (state(i,j,k,Density_comp) < small_dens)
    {
#ifndef CONST_SPECIES
        for (int n = 0; n < NumSpec; ++n)
            state(i,j,k,FirstSpec_comp+n) *= small_dens / state(i,j,k,Density_comp);
#endif

        state(i,j,k,Density_comp ) = small_dens;

        state(i,j,k,Xmom_comp) = 0.0;
        state(i,j,k,Ymom_comp) = 0.0;
        state(i,j,k,Zmom_comp) = 0.0;

        amrex::Real dummy_pres = 1e200;
        amrex::Real eint_new = -1e200;
        amrex::Real Ne = 0.0;

        // Re-create "e" from {small_dens, small_temp}
        nyx_eos_given_RT(atomic_rates, gamma_minus_1, h_species_in, &eint_new, &dummy_pres, state(i,j,k,Density_comp),
                         small_temp, Ne,a_new);

        // Define (rho e) with small_dens and the new "e" 
        state(i,j,k,Eint_comp) = state(i,j,k,Density_comp) *  eint_new;

        // Here we effectively zero the KE so set (rho E) = (rho e)
        state(i,j,k,Eden_comp) = state(i,j,k,Eint_comp);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_mu_for_enforce_min(
  const int i,
  const int j,
  const int k,
  const int icomp,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& mu_x,
  amrex::Array4<amrex::Real> const& mu_y,
  amrex::Array4<amrex::Real> const& mu_z,
  const amrex::Real small_value)
{
    // NOTE: This is called on a tile box grown by one, and we have already done a FillBoundary/FillPatch.
    //       We include ghost cells in our tests for state(i,j,k,icomp) > small_value so that we can create 
    //       fluxes here regardless of which side of the grid boundary has the cell needing fixing.

    // NOTE: We trigger a correction when state(i,j,k,icomp) < small_value, but if this is true, we adjust so that 
    //       state(i,j,k,icomp) == 1.01 * small_value ... just to give ourselves some breathing room 
    Real target = 1.01 * small_value;

    if (state(i,j,k,icomp) < small_value)
    {
        // std::cout << "UPDATING " << IntVect(i,j,k) << " WITH INITIAL VALUE " << state(i,j,k,icomp) << std::endl;
        // std::cout << " and x neighbors " <<  state(i+1,j,k,icomp) << " " <<  state(i-1,j,k,icomp) << std::endl;
        // std::cout << " and y neighbors " <<  state(i,j+1,k,icomp) << " " <<  state(i,j-1,k,icomp) << std::endl;
        // std::cout << " and z neighbors " <<  state(i,j,k+1,icomp) << " " <<  state(i,j,k-1,icomp) << std::endl;

        Real total_need  = target - state(i,j,k,icomp);

        // We only let each cell contribute 1/6 of its difference from target;
        //    this guarantees that the redistribution can't make a cell with (state > target)
        //    go under target
        Real avail_from_ihi = amrex::max((state(i+1,j,k,icomp) - target) / 6.0, 0.0);
        Real avail_from_ilo = amrex::max((state(i-1,j,k,icomp) - target) / 6.0, 0.0);
        Real avail_from_jhi = amrex::max((state(i,j+1,k,icomp) - target) / 6.0, 0.0);
        Real avail_from_jlo = amrex::max((state(i,j-1,k,icomp) - target) / 6.0, 0.0);
        Real avail_from_khi = amrex::max((state(i,j,k+1,icomp) - target) / 6.0, 0.0);
        Real avail_from_klo = amrex::max((state(i,j,k-1,icomp) - target) / 6.0, 0.0);

        // std::cout << " avail x " << avail_from_ihi << " " << avail_from_ilo << std::endl;
        // std::cout << " avail y " << avail_from_jhi << " " << avail_from_jlo << std::endl;
        // std::cout << " avail z " << avail_from_khi << " " << avail_from_klo << std::endl;

        Real total_avail = avail_from_ihi + avail_from_ilo +
                           avail_from_jhi + avail_from_jlo +
                           avail_from_khi + avail_from_klo; 

        // std::cout << " total avail / total need " << total_avail << " " << total_need << std::endl;

        Real fac;
        if (total_need < total_avail)
            fac = total_need / total_avail;
        else
            fac = 1.0;

        Real from_ihi = fac * avail_from_ihi;
        Real from_ilo = fac * avail_from_ilo;
        Real from_jhi = fac * avail_from_jhi;
        Real from_jlo = fac * avail_from_jlo;
        Real from_khi = fac * avail_from_khi;
        Real from_klo = fac * avail_from_klo;

        // mu should always be positive - we view this as a "diffusive" operator

        if (from_ihi > 0)
        {
            mu_x(i+1,j,k) =  from_ihi / (state(i+1,j,k,icomp) - state(i  ,j,k,icomp));
            if (mu_x(i+1,j,k) < 0.) amrex::Abort("mu_x(i+1,j,k) < 0");
        }
        if (from_ilo > 0)
        {
            mu_x(i  ,j,k) = -from_ilo / (state(i  ,j,k,icomp) - state(i-1,j,k,icomp));
            if (mu_x(i  ,j,k) < 0.) amrex::Abort("mu_x(i  ,j,k) < 0");
        }
        if (from_jhi > 0)
        {
            mu_y(i,j+1,k) =  from_jhi  / (state(i,j+1,k,icomp) - state(i,j  ,k,icomp));
            if (mu_y(i,j+1,k) < 0.) amrex::Abort("mu_y(i,j+1,k) < 0");
        }
        if (from_jlo > 0)
        {
            mu_y(i,j  ,k) = -from_jlo  / (state(i,j  ,k,icomp) - state(i,j-1,k,icomp));
            if (mu_y(i,j  ,k) < 0.) amrex::Abort("mu_y(i,j  ,k) < 0");
        }
        if (from_khi > 0)
        {
            mu_z(i,j,k+1) =  from_khi  / (state(i,j,k+1,icomp) - state(i,j,k  ,icomp));
            if (mu_z(i,j,k+1) < 0.) amrex::Abort("mu_z(i,j,k+1) < 0");
        }
        if (from_klo > 0)
        {
            mu_z(i,j,k  ) = -from_klo  / (state(i,j,k  ,icomp) - state(i,j,k-1,icomp));
            if (mu_z(i,j,k  ) < 0.) amrex::Abort("mu_z(i,j,k  ) < 0");
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
create_update_for_minimum(
  const int i,
  const int j,
  const int k,
#ifndef CONST_SPECIES
  const int icomp,
#endif
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& mu_x,
  amrex::Array4<amrex::Real> const& mu_y,
  amrex::Array4<amrex::Real> const& mu_z,
  amrex::Array4<amrex::Real> const& update 
#ifndef CONST_SPECIES
 ,const int FirstSpec_comp, const int NumSpec
#endif
                                             )
{
    for (int n = 0; n < state.nComp(); n++)
    {
        if (n == Density_comp || n == Xmom_comp || n == Ymom_comp || n == Zmom_comp || n == Eden_comp || n == Eint_comp)
        {
            update(i,j,k,n) = mu_x(i+1,j,k) * (state(i+1,j,k,n) - state(i  ,j,k,n))
                             -mu_x(i  ,j,k) * (state(i  ,j,k,n) - state(i-1,j,k,n))
                             +mu_y(i,j+1,k) * (state(i,j+1,k,n) - state(i,j  ,k,n))
                             -mu_y(i,j  ,k) * (state(i,j  ,k,n) - state(i,j-1,k,n))
                             +mu_z(i,j,k+1) * (state(i,j,k+1,n) - state(i,j,k  ,n))
                             -mu_z(i,j,k  ) * (state(i,j,k  ,n) - state(i,j,k-1,n));
        } else {
            update(i,j,k,n) = 0.0;
        }
    }

#ifndef CONST_SPECIES
    // Enforce that X is unchanged with change in rho
    if (icomp == Density_comp)
        for (int n = 0; n < NumSpec; ++n)
        {
            state(i,j,k,FirstSpec_comp+n) *= (state(i,j,k,Density_comp) + update(i,j,k,Density_comp)) / state(i,j,k,Density_comp);
        }
#endif
}
