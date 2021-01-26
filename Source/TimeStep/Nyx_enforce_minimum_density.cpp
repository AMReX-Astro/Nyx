#include <Nyx.H>
#include <Nyx_enforce_minimum_density.H>
#include <constants_cosmo.H>

using namespace amrex;

void
Nyx::enforce_minimum_density( MultiFab& S_old, MultiFab& S_new,
#ifdef SDC
                              MultiFab& hydro_source,
                              MultiFab& reset_e_src,
#endif
                              Real a_new)
{
    BL_PROFILE("Nyx::enforce_minimum_density()");

    if (verbose)
      amrex::Print() << "Enforce minimum density... " << std::endl;

    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));

    if (S_new.min(Density_comp) < small_dens)
    {
        if (enforce_min_density_type == "floor")
        {
            enforce_minimum_density_floor(S_new, a_new);

        } else if (enforce_min_density_type == "conservative") {
#ifdef SDC
            enforce_minimum_density_cons(S_old, S_new, reset_e_src);
#else
            enforce_minimum_density_cons(S_old, S_new);
#endif
        } else {
            amrex::Abort("Don't know this enforce_min_density_type");
        }
#ifdef SDC
    //
    //  Reset hydro_src = A_rho for SDC solve (This assumes ext_src_old(rho)=0)
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(hydro_source,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Only update on valid cells
            const amrex::Box& bx = mfi.tilebox();
            auto const& hydro_src = hydro_source.array(mfi);
            auto const& uin  = S_old.array(mfi);
            auto const& uout = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                hydro_src(i,j,k,Density_comp) = uout(i,j,k,Density_comp) - uin(i,j,k,Density_comp);
            });

        }
#endif
    }
}

void
Nyx::enforce_minimum_density_floor( MultiFab& S_new, Real a_new_in)
{
    Real lsmall_dens = small_dens;
    Real lgamma_minus_1 = gamma - 1.0;
    Real lsmall_temp = small_temp;
    auto atomic_rates = atomic_rates_glob;

#ifndef CONST_SPECIES
    int lnum_spec    = NumSpec;
#endif

    Real l_h_species = h_species;

    //
    //  Reset negative density to small_dens, set (rho e) and (rho E) from small_temp  and zero out momenta
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Only update on valid cells
        const amrex::Box& bx = mfi.tilebox();
        auto const& uout = S_new.array(mfi);

         amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
         {
             floor_density(i, j, k, uout, atomic_rates,
#ifndef CONST_SPECIES
                           lnum_spec,
#endif
                           a_new_in, lgamma_minus_1, lsmall_dens, lsmall_temp, l_h_species);
         });
    }
}

void
#ifdef SDC
Nyx::enforce_minimum_density_cons ( MultiFab& S_old, MultiFab& S_new, MultiFab& reset_e_src)
#else
Nyx::enforce_minimum_density_cons ( MultiFab& S_old, MultiFab& S_new)
#endif
{
    bool debug = false;

#ifndef CONST_SPECIES
    int lnum_spec    = NumSpec;
    int lfirst_spec  = FirstSpec_comp;
#endif
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

    Real rho_old_min_before = S_old.min(Density_comp);
    Real  ru_old_min_before = S_old.min(Xmom_comp);
    Real  rv_old_min_before = S_old.min(Ymom_comp);
    Real  rw_old_min_before = S_old.min(Zmom_comp);
    Real  re_old_min_before = S_old.min(Eint_comp);
    Real  rE_old_min_before = S_old.min(Eden_comp);

    Real rho_new_min_before = S_new.min(Density_comp);
    Real  ru_new_min_before = S_new.min(Xmom_comp);
    Real  rv_new_min_before = S_new.min(Ymom_comp);
    Real  rw_new_min_before = S_new.min(Zmom_comp);
    Real  re_new_min_before = S_new.min(Eint_comp);
    Real  rE_new_min_before = S_new.min(Eden_comp);

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
    {
        for (MFIter mfi(S_new,MFItInfo().UseDefaultStream()); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.growntilebox();
            for (int i = 0; i < S_new[mfi].nComp(); i++)
            {
                IntVect p_nan(D_DECL(-10, -10, -10));
//                if (ParallelDescriptor::IOProcessor())
//                    std::cout << "enforce_minimum_density: testing component " << i << " for NaNs" << std::endl;
                bool has_nan=S_new[mfi].contains_nan<RunOn::Device>(bx,Density_comp+i,1,p_nan);
                if (has_nan)
                {
                    std::cout<<"nans in comp "<<i<<" at "<<p_nan<<std::flush<<std::endl;
                }
            }
        }
        if (S_new.contains_nan(Density_comp, S_new.nComp(), 0))
            amrex::Abort("NaN in enforce_minimum_density before we start iterations");
    }

    // 10 is an arbitrary limit here -- just to make sure we don't get stuck here somehow
    while (too_low && iter < 10)
    {
        // First make sure that all ghost cells are updated because we use them in defining fluxes
        // Note that below we update S_new, not Sborder, so we must FillPatch each time.
        FillPatch(*this, Sborder, 2, cur_time, State_Type, Density_comp, Sborder.nComp());

        // Initialize to zero; these will only be non-zero at a face across which density is passed...
        mu_x.setVal(0.);
        mu_y.setVal(0.);
        mu_z.setVal(0.);

        // This will hold the update to each cell due to enforcing minimum density in a conservative way
        MultiFab update(grids , dmap, Sborder.nComp(), 0);
        update.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& gbx = mfi.growntilebox(1);
            auto const& sbord = Sborder.array(mfi);
            auto const& mu_x_arr = mu_x.array(mfi);
            auto const& mu_y_arr = mu_y.array(mfi);
            auto const& mu_z_arr = mu_z.array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
              compute_mu_for_enforce_min(i, j, k, Density_comp, sbord, mu_x_arr, mu_y_arr, mu_z_arr, lsmall_dens);
            });
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
              create_update_for_minimum(i, j, k,
#ifndef CONST_SPECIES
                                        Density_comp,
#endif
                                        sbord_arr, mu_x_arr, mu_y_arr, mu_z_arr, upd_arr
#ifndef CONST_SPECIES
                                       ,lfirst_spec, lnum_spec
#endif
                                                              );
            });
        }

        S_new.plus(update,0,S_new.nComp(),0);

#ifdef SDC
        MultiFab::Copy(reset_e_src,update,Eint_comp,0,1,0);
#endif

        if (S_new.contains_nan())
        {
           amrex::Print() << "Doing iteration iter " << std::endl;
           amrex::Abort("   and finding NaN in enforce_minimum_density");
        }

        // This is used to decide whether to continue the iteration
        rho_new_min_after = S_new.min(Density_comp);

        // This is used to decide whether to call enforce_min_energy_cons at the end of this routine
        re_new_min_after = S_new.min(Eint_comp);

        if (debug)
        {
             ru_new_min_after = S_new.min(Xmom_comp);
             rv_new_min_after = S_new.min(Ymom_comp);
             rw_new_min_after = S_new.min(Zmom_comp);
             rE_new_min_after = S_new.min(Eden_comp);
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
