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
  const int NumSpec,
  amrex::Real const a_new,
  const amrex::Real gamma_minus_1,
  const amrex::Real small_dens,
  const amrex::Real small_temp)
{
    //
    //  Reset negative density to small_dens, set (rho e) and (rho E) from small_temp  and zero out momenta
    //
    if (state(i,j,k,URHO) < small_dens)
    {
        for (int n = 0; n < NumSpec; ++n)
        {
            state(i,j,k,FirstSpec+n) *= small_dens / state(i,j,k,URHO);
        }

        state(i,j,k,URHO ) = small_dens;

        state(i,j,k,UMX) = 0.0;
        state(i,j,k,UMY) = 0.0;
        state(i,j,k,UMZ) = 0.0;

        amrex::Real dummy_pres = 1e200;
        amrex::Real eint_new = -1e200;
        amrex::Real Ne = 0.0;
        amrex::Real h_species = 0.76;

        // Re-create "e" from {small_dens, small_temp}
        nyx_eos_given_RT(atomic_rates, gamma_minus_1, h_species, &eint_new, &dummy_pres, state(i,j,k,URHO),
			 small_temp, Ne,a_new);

        // Define (rho e) with small_dens and the new "e" 
        state(i,j,k,UEINT) = state(i,j,k,URHO) *  eint_new;

        // Here we effectively zero the KE so set (rho E) = (rho e)
        state(i,j,k,UEDEN) = state(i,j,k,UEINT);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_mu_for_enforce_min_density(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& mu_x,
  amrex::Array4<amrex::Real> const& mu_y,
  amrex::Array4<amrex::Real> const& mu_z,
  const amrex::Real small_dens)
{
    // NOTE: This is called on a tile box grown by one, and we have already done a FillBoundary/FillPatch.
    //       We include ghost cells in our tests for rho > small_dens so that we can create fluxes here
    //       regardless of which side of the grid boundary has the cell needing fixing.

    // NOTE: We trigger a correction when rho < small_dens, but if this is true, we adjust so that 
    //       rho == 1.1 * small_dens ... just to give ourselves some breathing room 
    Real target_density = 1.01 * small_dens;

    if (state(i,j,k,URHO) < small_dens)
    {
        // std::cout << "UPDATING " << IntVect(i,j,k) << " WITH INITIAL RHO " << state(i,j,k,URHO) << std::endl;
        // std::cout << " and x neighbors " <<  state(i+1,j,k,URHO) << " " <<  state(i-1,j,k,URHO) << std::endl;
        // std::cout << " and y neighbors " <<  state(i,j+1,k,URHO) << " " <<  state(i,j-1,k,URHO) << std::endl;
        // std::cout << " and z neighbors " <<  state(i,j,k+1,URHO) << " " <<  state(i,j,k-1,URHO) << std::endl;

        Real total_need  = target_density - state(i,j,k,URHO);

        // We only let each cell contribute 1/6 of its difference from target_density;
        //    this guarantees that the redistribution can't make a cell with (rho > target_density)
        //    go under target_density
        Real avail_from_ihi = amrex::max((state(i+1,j,k,URHO) - target_density) / 6.0, 0.0);
        Real avail_from_ilo = amrex::max((state(i-1,j,k,URHO) - target_density) / 6.0, 0.0);
        Real avail_from_jhi = amrex::max((state(i,j+1,k,URHO) - target_density) / 6.0, 0.0);
        Real avail_from_jlo = amrex::max((state(i,j-1,k,URHO) - target_density) / 6.0, 0.0);
        Real avail_from_khi = amrex::max((state(i,j,k+1,URHO) - target_density) / 6.0, 0.0);
        Real avail_from_klo = amrex::max((state(i,j,k-1,URHO) - target_density) / 6.0, 0.0);

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
            mu_x(i+1,j,k) =  from_ihi / (state(i+1,j,k,URHO) - state(i  ,j,k,URHO));
            if (mu_x(i+1,j,k) < 0.) amrex::Abort("mu_x(i+1,j,k) < 0");
        }
        if (from_ilo > 0)
        {
            mu_x(i  ,j,k) = -from_ilo / (state(i  ,j,k,URHO) - state(i-1,j,k,URHO));
            if (mu_x(i  ,j,k) < 0.) amrex::Abort("mu_x(i  ,j,k) < 0");
        }
        if (from_jhi > 0)
        {
            mu_y(i,j+1,k) =  from_jhi  / (state(i,j+1,k,URHO) - state(i,j  ,k,URHO));
            if (mu_y(i,j+1,k) < 0.) amrex::Abort("mu_y(i,j+1,k) < 0");
        }
        if (from_jlo > 0)
        {
            mu_y(i,j  ,k) = -from_jlo  / (state(i,j  ,k,URHO) - state(i,j-1,k,URHO));
            if (mu_y(i,j  ,k) < 0.) amrex::Abort("mu_y(i,j  ,k) < 0");
        }
        if (from_khi > 0)
        {
            mu_z(i,j,k+1) =  from_khi  / (state(i,j,k+1,URHO) - state(i,j,k  ,URHO));
            if (mu_z(i,j,k+1) < 0.) amrex::Abort("mu_z(i,j,k+1) < 0");
        }
        if (from_klo > 0)
        {
            mu_z(i,j,k  ) = -from_klo  / (state(i,j,k  ,URHO) - state(i,j,k-1,URHO));
            if (mu_z(i,j,k  ) < 0.) amrex::Abort("mu_z(i,j,k  ) < 0");
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
create_update_for_minimum_density(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& mu_x,
  amrex::Array4<amrex::Real> const& mu_y,
  amrex::Array4<amrex::Real> const& mu_z,
  amrex::Array4<amrex::Real> const& update,
  const int FirstSpec, const int NumSpec)
{
    for (int n = 0; n < state.nComp(); n++)
    {
        if (n == URHO or n == UMX or n == UMY or n == UMZ or n == UEDEN or n == UEINT) 
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

    // Enforce that X is unchanged with change in rho
    for (int n = 0; n < NumSpec; ++n)
    {
        state(i,j,k,FirstSpec+n) *= (state(i,j,k,URHO) + update(i,j,k,URHO)) / state(i,j,k,URHO);
    }
}
