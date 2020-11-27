#include <Nyx.H>
#include <Nyx_enforce_minimum_density.H>
#include <constants_cosmo.H>

using namespace amrex;

void
Nyx::enforce_minimum_density( MultiFab& S_old, MultiFab& S_new,
                              Real dt, Real a_old, Real a_new)
{
    BL_PROFILE("Nyx::enforce_minimum_density()");

    if (verbose)
      amrex::Print() << "Enforce minimum density... " << std::endl;

    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));

    if (S_new.min(Density) < small_dens)
    {
        if (enforce_min_density_type == "floor")
        {
            enforce_minimum_density_floor(S_old, S_new, dt, a_old, a_new);

        } else if (enforce_min_density_type == "conservative") {
            enforce_minimum_density_cons(S_old, S_new, dt, a_old, a_new);

        } else {
            amrex::Abort("Don't know this enforce_min_density_type");
        }
    }
}

void
Nyx::enforce_minimum_density_floor( MultiFab& S_old, MultiFab& S_new,
                                    Real dt, Real a_old, Real a_new )
{
    int lnum_spec    = NumSpec;
    Real lsmall_dens = small_dens;
    Real lgamma_minus_1 = gamma - 1.0;
    Real lsmall_temp = small_temp;
    auto atomic_rates = atomic_rates_glob;

    //
    //  Reset negative density to small_dens, set (rho e) and (rho E) from small_temp  and zero out momenta
    //
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

void
Nyx::enforce_minimum_density_cons ( MultiFab& S_old, MultiFab& S_new,
                                    Real dt, Real a_old, Real a_new )
{
    bool debug = false;

    int lnum_spec    = NumSpec;
    int lfirst_spec  = FirstSpec;
    Real lsmall_dens = small_dens;
  
    Real cur_time = state[State_Type].curTime();

    // We need to define this temporary because S_new only has one ghost cell and we need two.
    MultiFab Sborder;
    Sborder.define(grids, S_new.DistributionMap(), S_new.nComp(), 2);

    // Define face-based coefficients to be defined when enforcing minimum density 
    //     then used to enjoy the updates of all the other variables
    // The ghost face space is only needed as temp space; we only use "valid" faces...
    MultiFab mu_x(amrex::convert(grids,IntVect(1,0,0)), dmap, 1, 1);
    MultiFab mu_y(amrex::convert(grids,IntVect(0,1,0)), dmap, 1, 1);
    MultiFab mu_z(amrex::convert(grids,IntVect(0,0,1)), dmap, 1, 1);

    Real rho_old_min_before = S_old.min(URHO);
    Real  ru_old_min_before = S_old.min(UMX);
    Real  rv_old_min_before = S_old.min(UMY);
    Real  rw_old_min_before = S_old.min(UMZ);
    Real  re_old_min_before = S_old.min(UEINT);
    Real  rE_old_min_before = S_old.min(UEDEN);

    Real rho_new_min_before = S_new.min(URHO);
    Real  ru_new_min_before = S_new.min(UMX);
    Real  rv_new_min_before = S_new.min(UMY);
    Real  rw_new_min_before = S_new.min(UMZ);
    Real  re_new_min_before = S_new.min(UEINT);
    Real  rE_new_min_before = S_new.min(UEDEN);

    Real rho_old_sum_before = S_old.sum(0);
    Real rho_new_sum_before = S_new.sum(0);

    Real rho_new_min_after;
    Real  re_new_min_after;
    Real  ru_new_min_after;
    Real  rv_new_min_after;
    Real  rw_new_min_after;
    Real  rE_new_min_after;

    Real rho_new_sum_after;

    Real rho_new_min = rho_new_min_before;

    bool too_low = (rho_new_min < small_dens);

    int iter = 0;

    if (S_new.contains_nan())
       amrex::Abort("NaN in enforce_minimum_density before we start iterations");

    // 10 is an arbitrary limit here -- just to make sure we don't get stuck here somehow 
    while (too_low and iter < 10)
    {
        // First make sure that all ghost cells are updated because we use them in defining fluxes
        // Note that below we update S_new, not Sborder, so we must FillPatch each time.
        FillPatch(*this, Sborder, 2, cur_time, State_Type, Density, Sborder.nComp());

        // Initialize to zero; these will only be non-zero at a face across which density is passed...
        mu_x.setVal(0.);
        mu_y.setVal(0.);
        mu_z.setVal(0.);

        // This will hold the update to each cell due to enforcing minimum density in a conservative way
        MultiFab update(grids , dmap, Sborder.nComp(), 0);
        update.setVal(0.);

        for (amrex::MFIter mfi(Sborder, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const amrex::Box& gbx = mfi.growntilebox(1);
            auto const& sbord = Sborder.array(mfi);
            auto const& mu_x_arr = mu_x.array(mfi);
            auto const& mu_y_arr = mu_y.array(mfi);
            auto const& mu_z_arr = mu_z.array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
              compute_mu_for_enforce_min_density(i, j, k, sbord, mu_x_arr, mu_y_arr, mu_z_arr, lsmall_dens);
            });
        }

        for (amrex::MFIter mfi(Sborder, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            // Only update on valid cells
            const amrex::Box& bx = mfi.tilebox();
            auto const& sbord_arr = Sborder.array(mfi);
            auto const& mu_x_arr  = mu_x.array(mfi);
            auto const& mu_y_arr  = mu_y.array(mfi);
            auto const& mu_z_arr  = mu_z.array(mfi);
            auto const& upd_arr   = update.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
              create_update_for_minimum_density(i, j, k, sbord_arr, 
                                                mu_x_arr, mu_y_arr, mu_z_arr, upd_arr, 
                                                lfirst_spec, lnum_spec);
            });
        }

        S_new.plus(update,0,S_new.nComp(),0);

        if (S_new.contains_nan())
        {
           amrex::Print() << "Doing iteration iter " << std::endl; 
           amrex::Abort("   and finding NaN in enforce_minimum_density");
        }

        // This is used to decide whether to continue the iteration
        rho_new_min_after = S_new.min(URHO);

        if (debug) 
        {
             ru_new_min_after = S_new.min(UMX);
             rv_new_min_after = S_new.min(UMY);
             rw_new_min_after = S_new.min(UMZ);
             re_new_min_after = S_new.min(UEINT);
             rE_new_min_after = S_new.min(UEDEN);
            amrex::Print() << "After " << iter+1 << " iterations " << std::endl;
            amrex::Print() << "  MIN OF rho: old / new / new new " << 
                rho_old_min_before << " " << rho_new_min_before << " " << rho_new_min_after << std::endl;
            amrex::Print() << "  MIN OF  ru: old / new / new new " << 
                ru_old_min_before << " " <<  ru_new_min_before << " " << ru_new_min_after << std::endl;
            amrex::Print() << "  MIN OF  rv: old / new / new new " << 
                rv_old_min_before << " " << rv_new_min_before << " " << rv_new_min_after << std::endl;
            amrex::Print() << "  MIN OF  rw: old / new / new new " << 
                rw_old_min_before << " " << rw_new_min_before << " " << rw_new_min_after << std::endl;
            amrex::Print() << "  MIN OF  re: old / new / new new " << 
                re_old_min_before << " " << re_new_min_before << " " << re_new_min_after << std::endl;
            amrex::Print() << "  MIN OF  rE: old / new / new new " << 
                rE_old_min_before << " " << rE_new_min_before << " " << rE_new_min_after << std::endl;
        }

        too_low = (rho_new_min_after < small_dens);
        iter++;

    } // iter

    rho_new_sum_after = S_new.sum(0);

    amrex::Print() << "After " << iter << " iterations " << std::endl;
    amrex::Print() << "  SUM OF rho_old / rho_new / new rho_new " << 
            rho_old_sum_before << " " << rho_new_sum_before << " " << rho_new_sum_after << std::endl;

    if (rho_new_min_after < small_dens)
       amrex::Abort("Not able to enforce small_dens this way after all");
}
