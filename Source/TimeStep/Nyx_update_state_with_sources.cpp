#include <Nyx.H>
#include <Hydro.H>
#include <constants_cosmo.H>

using namespace amrex;

using std::string;

void
Nyx::update_state_with_sources( MultiFab& S_old, MultiFab& S_new,
                                MultiFab& ext_src_old, MultiFab& hydro_source,
                                MultiFab& grav_vector, 
#ifdef SDC
                                MultiFab& reset_e_src,
#endif
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
#ifdef SDC
        auto const& reset_src = reset_e_src.array(mfi);
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            for (int n = 0; n < uout.nComp(); ++n) 
            {
                if(n==Density_comp)
                {
                        uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) 
                                + dt *  src(i,j,k,n) * a_half_inv;
                }
                else if(n>=Xmom_comp&&n<=Zmom_comp)
                {
                        uout(i,j,k,n) = a_old * uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n);
                        uout(i,j,k,n) = uout(i,j,k,n) * a_new_inv;
                }
                else if(n==Eden_comp)
                {
                        uout(i,j,k,n) =  a_oldsq * uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n);
                        uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
                }
                else if(n==Eint_comp)
                {
                        uout(i,j,k,n) =  a_oldsq*uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n) ;
                        uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
//                        reset_src(i,j,k) = uout(i,j,k,n);
                }
                else
                {
                        uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n) * a_half_inv;
                }
            }
        });
    }

    // Enforce minimum density over the whole MultiFab
    enforce_minimum_density(S_old, S_new, 
#ifdef SDC
                            hydro_source,
                            reset_e_src, 
#endif
                            a_new);

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
            const amrex::Real rho = uin(i, j, k, Density_comp);

            const amrex::Real SrU = rho * grav(i,j,k,0);
            const amrex::Real SrV = rho * grav(i,j,k,1);
            const amrex::Real SrW = rho * grav(i,j,k,2);

            // We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
            uout(i,j,k,Xmom_comp)   += SrU * dt_a_new;
            uout(i,j,k,Ymom_comp)   += SrV * dt_a_new;
            uout(i,j,k,Zmom_comp)   += SrW * dt_a_new;

            // Src = rho u dot g, evaluated with all quantities at t^n
            const amrex::Real SrE = uin(i,j,k,Xmom_comp) * grav(i,j,k,0) +
                                    uin(i,j,k,Ymom_comp) * grav(i,j,k,1) +
                                    uin(i,j,k,Zmom_comp) * grav(i,j,k,2);
            uout(i,j,k,Eden_comp) = (a_newsq*uout(i,j,k,Eden_comp) + SrE * (dt*a_half)) * a_newsq_inv;
        });
    }
}
