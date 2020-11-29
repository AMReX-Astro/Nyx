#include <Nyx.H>
#include <Hydro.H>
#include <constants_cosmo.H>

using namespace amrex;

using std::string;

void
Nyx::update_state_with_sources( MultiFab& S_old, MultiFab& S_new,
                                MultiFab& ext_src_old, MultiFab& hydro_source,
                                MultiFab& grav_vector, 
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
    }

    // Enforce minimum density over the whole MultiFab
    enforce_minimum_density(S_old, S_new, dt, a_old, a_new);

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

                // Src = rho u dot g, evaluated with all quantities at t^n
                const amrex::Real SrE = uin(i,j,k,UMX) * grav(i,j,k,0) +
                                        uin(i,j,k,UMY) * grav(i,j,k,1) +
                                        uin(i,j,k,UMZ) * grav(i,j,k,2);
                uout(i,j,k,UEDEN) = (a_newsq*uout(i,j,k,UEDEN) + SrE * (dt*a_half)) * a_newsq_inv;
        });
    }
}
