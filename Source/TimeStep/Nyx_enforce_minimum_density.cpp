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

void
Nyx::enforce_minimum_density( MultiFab& S_old, MultiFab& S_new,
                              amrex::Real dt, amrex::Real a_old, amrex::Real a_new)
{
    BL_PROFILE("Nyx::enforce_minimum_density()");

    if (verbose)
      amrex::Print() << "Enforce minimum density... " << std::endl;

    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));

    int lnum_spec    = NumSpec;
    Real lsmall_dens = small_dens;
    Real lgamma_minus_1 = gamma - 1.0;
    Real lsmall_temp = small_temp;
    auto atomic_rates = atomic_rates_glob;

//  if (S_new.min(Density) < small_dens)
    {
        for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            // Only update on valid cells
            const amrex::Box& bx = mfi.tilebox();
            auto const& uout = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                floor_density(i, j, k, uout, atomic_rates, lnum_spec, a_new, lgamma_minus_1, lsmall_dens, lsmall_temp);
            });
        }
    }
}
