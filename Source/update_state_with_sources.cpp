#include <Nyx.H>
#include <Nyx_F.H>
#include <Hydro.H>
#include <constants_cosmo.H>

using namespace amrex;

using std::string;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
do_enforce_minimum_density(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& state,
  const int NumSpec,
  amrex::Real const a_new,
  const amrex::Real gamma_minus_1,
  const amrex::Real small_dens,
  const amrex::Real small_temp)
{
    //Reset negative density to small_dens and zero out momenta:
    if(state(i,j,k,URHO) < small_dens)
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

        // HACK HACK -- we can't yet do this call
        // Re-create "e" from {small_dens, small_temp}
        nyx_eos_given_RT(gamma_minus_1, h_species, &eint_new, &dummy_pres, state(i,j,k,URHO),
						 small_temp, Ne,a_new);

        // Define (rho e) with small_dens and the new "e" 
        state(i,j,k,UEINT) = state(i,j,k,URHO) *  eint_new;

        // Here we effectively zero the KE so set (rho E) = (rho e)
        state(i,j,k,UEDEN) = state(i,j,k,UEINT);
    }
}

void
Nyx::update_state_with_sources( MultiFab& S_old, MultiFab& S_new,
                                MultiFab& ext_src_old, MultiFab& hydro_source,
                                MultiFab& grav_vector, MultiFab& divu_cc,
                                amrex::Real dt, amrex::Real a_old, amrex::Real a_new)
{
    BL_PROFILE("Nyx::update_state_with_sources()");

    if (verbose)
      amrex::Print() << "Updating state with the hydro sources ... " << std::endl;

    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));

    const amrex::Real a_half = 0.5 * (a_old + a_new);
    const amrex::Real a_half_inv = 1 / a_half;
    const amrex::Real a_oldsq = a_old * a_old;
    const amrex::Real a_newsq = a_new * a_new;
    const amrex::Real a_new_inv = 1.0 / a_new;
    const amrex::Real a_newsq_inv = 1.0 / a_newsq;
    const amrex::Real dt_a_new    = dt / a_new;

    int lnum_spec    = NumSpec;
    Real lsmall_dens = small_dens;
	Real lgamma_minus_1 = gamma - 1.0;
	Real lsmall_temp = small_temp;

        ////This set of dt should be used for Saxpy dt like setup
    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {

        const amrex::Box& bx = mfi.tilebox();
        auto const& uin = S_old.array(mfi);
        auto const& uout = S_new.array(mfi);
        auto const& hydro_src = hydro_source.array(mfi);
        auto const& src = ext_src_old.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            for (int n = 0; n < uout.nComp(); ++n) 
            {
                        if(n==URHO)
                        {
                                uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) 
                                        + dt *  src(i,j,k,n) * a_half_inv;
                        }
                        else if(n>=UMX&&n<=UMZ)
                        {
                                uout(i,j,k,n) = a_old * uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n);
                                uout(i,j,k,n) = uout(i,j,k,n) * a_new_inv;
                        }
                        else if(n==UEDEN)
                        {
                                uout(i,j,k,n) =  a_oldsq * uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n);
                                uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
                        }
                        else if(n==UEINT)
                        {
                                uout(i,j,k,n) =  a_oldsq*uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n) ;
                                uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
                        }
                        else
                        {
                                uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n) * a_half_inv;
                        }
            }
        });

        //Unclear whether this should be part of previous ParallelFor
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            do_enforce_minimum_density(i, j, k, uout, lnum_spec, a_new, lgamma_minus_1, lsmall_dens, lsmall_temp);
        });

    }

#ifndef CONST_SPECIES
    enforce_nonnegative_species(S_new);
#endif

    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const amrex::Box& bx = mfi.tilebox();
        auto const& uin = S_old.array(mfi);
        auto const& uout = S_new.array(mfi);
        auto const& grav = grav_vector.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {


                const amrex::Real rhoInv = 1.0 / uout(i,j,k,URHO);
                const amrex::Real vx = uout(i, j, k, UMX);
                const amrex::Real vy = uout(i, j, k, UMY);
                const amrex::Real vz = uout(i, j, k, UMZ);

                // **** Start Diagnostics ****
                //Multiplies by rhoInv
                const amrex::Real old_ke = 0.5 * rhoInv * (vx * vx + vy * vy + vz * vz);
                const amrex::Real old_rhoeint = uout(i,j,k,UEDEN) - old_ke;
                // ****   End Diagnostics ****

                const amrex::Real rho = uin(i, j, k, URHO);

                const amrex::Real SrU = rho * grav(i,j,k,0);
                const amrex::Real SrV = rho * grav(i,j,k,1);
                const amrex::Real SrW = rho * grav(i,j,k,2);

                // We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
                uout(i,j,k,UMX)   += SrU * dt_a_new;
                uout(i,j,k,UMY)   += SrV * dt_a_new;
                uout(i,j,k,UMZ)   += SrW * dt_a_new;

                if (grav_source_type == 1)
                {
                    // This does work (in 1-d)
                    // Src = rho u dot g, evaluated with all quantities at t^n
                    const amrex::Real SrE = uin(i,j,k,UMX) * grav(i,j,k,0) +
                        uin(i,j,k,UMY) * grav(i,j,k,1) +
                        uin(i,j,k,UMZ) * grav(i,j,k,2);
                    uout(i,j,k,UEDEN) = (a_newsq*uout(i,j,k,UEDEN) + SrE * (dt*a_half)) * a_newsq_inv;
                }
                else if (grav_source_type == 3)
                {
                    //Multiplies by rhoInv
                    const amrex::Real new_ke = 0.5 * rhoInv * (vx * vx + vy * vy + vz * vz);
                    uout(i,j,k,UEDEN) = old_rhoeint + new_ke;
                                }
                else
                    amrex::Abort("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type");


        });
    }
}
