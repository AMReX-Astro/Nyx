
#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>

#ifdef GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

using std::string;

#ifndef NO_HYDRO
void
Nyx::just_the_hydro (Real time,
                     Real dt,
                     Real a_old,
                     Real a_new)
{
    BL_PROFILE("Nyx::just_the_hydro()");

    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();
    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  D_old        = get_old_data(DiagEOS_Type);
    MultiFab&  S_new        = get_new_data(State_Type);
    MultiFab&  D_new        = get_new_data(DiagEOS_Type);

    if (std::abs(time-prev_time) > (1.e-10*cur_time) )
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "just_the_hydro:  prev_time = " << prev_time << std::endl;
            std::cout << "just_the_hydro:       time = " <<      time << std::endl;
        }
        amrex::Abort("time should equal prev_time in just_the_hydro!");
    }

#ifndef NDEBUG
    if (S_old.contains_nan(Density, S_old.nComp(), 0))
    {
        for (int i = 0; i < S_old.nComp(); i++)
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "just_the_hydro: testing component " << i << " for NaNs" << std::endl;
            if (S_old.contains_nan(Density+i,1,0))
                amrex::Abort("S_old has NaNs in this component");
        }
    }
#endif

    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

    if (do_reflux && level < finest_level)
    {
        //
        // Set reflux registers to zero.
        //
        get_flux_reg(level+1).setVal(0);
    }
    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister* fine    = 0;
    FluxRegister* current = 0;

    if (do_reflux && level < finest_level)
        fine = &get_flux_reg(level+1);
    if (do_reflux && level > 0)
        current = &get_flux_reg(level);

    const Real* dx     = geom.CellSize();
    Real        courno = -1.0e+200;

    MultiFab ext_src_old(grids, dmap, NUM_STATE, 3);
    ext_src_old.setVal(0);

    if (add_ext_src && ParallelDescriptor::IOProcessor())
    {
       if (strang_split)
       {
          std::cout << "Source terms are handled with strang splitting" << std::endl; 
       } else {
          std::cout << "Source terms are handled with predictor/corrector" << std::endl; 
       }
    }

    if (add_ext_src && !strang_split)
    {
#ifndef NO_OLD_SRC
        get_old_source(prev_time, dt, ext_src_old);
#endif //#ifndef NO_OLD_SRC
        ext_src_old.FillBoundary();
    }

    // Define the gravity vector so we can pass this to ca_umdrv.
    MultiFab grav_vector(grids, dmap, BL_SPACEDIM, 3);
    grav_vector.setVal(0.);

#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary(geom.periodicity());
#endif

    MultiFab fluxes[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
        fluxes[j].setVal(0.0);
    }

    BL_ASSERT(NUM_GROW == 4);

    Real  e_added = 0;
    Real ke_added = 0;

    // Create FAB for extended grid values (including boundaries) and fill.
    MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    FillPatch(*this, S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);
    MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), 2, NUM_GROW);
    FillPatch(*this, D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, 2);

    if (add_ext_src && strang_split) 
        strang_first_step(time,dt,S_old_tmp,D_old_tmp);

#ifdef _OPENMP
#pragma omp parallel
#endif
       {
       FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM];
       Real cflLoc = -1.e+200;

       for (MFIter mfi(S_old_tmp,true); mfi.isValid(); ++mfi)
       {

        const Box& bx        = mfi.tilebox();

        FArrayBox& state     = S_old_tmp[mfi];
        FArrayBox& dstate    = D_old_tmp[mfi];
        FArrayBox& stateout  = S_new[mfi];

#ifdef SHEAR_IMPROVED
        FArrayBox& am_tmp = AveMom_tmp[mfi];
#endif

        Real se  = 0;
        Real ske = 0;

        // Allocate fabs for fluxes.
        for (int i = 0; i < BL_SPACEDIM ; i++) {
            const Box &bxtmp = amrex::surroundingNodes(bx, i);
            flux[i].resize(bxtmp, NUM_STATE);
            u_gdnv[i].resize(amrex::grow(bxtmp, 1), 1);
            u_gdnv[i].setVal(1.e200);
        }

        fort_advance_gas
            (&time, bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(state),
             BL_TO_FORTRAN(stateout),
             BL_TO_FORTRAN(u_gdnv[0]),
             BL_TO_FORTRAN(u_gdnv[1]),
             BL_TO_FORTRAN(u_gdnv[2]),
             BL_TO_FORTRAN(ext_src_old[mfi]),
             BL_TO_FORTRAN(grav_vector[mfi]),
             dx, &dt,
             BL_TO_FORTRAN(flux[0]),
             BL_TO_FORTRAN(flux[1]),
             BL_TO_FORTRAN(flux[2]),
             &cflLoc, &a_old, &a_new, &se, &ske, &print_fortran_warnings, &do_grav);

        for (int i = 0; i < BL_SPACEDIM; ++i) {
          fluxes[i][mfi].copy(flux[i], mfi.nodaltilebox(i));
        }

         e_added += se;
        ke_added += ske;
       } // end of MFIter loop

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
       {
        courno = std::max(courno, cflLoc);
       }

       } // end of omp parallel region

       // We copy old Temp and Ne to new Temp and Ne so that they can be used
       //    as guesses when we next need them.
       MultiFab::Copy(D_new,D_old,0,0,2,0);

       if (do_reflux) {
         if (current) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
             current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1);
           }
         }
         if (fine) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
	         fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
           }
         }
       }

#ifdef GRAVITY
    if (verbose > 1)
    {
        Real added[2] = {e_added,ke_added};

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(added, 2, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            const Real vol = D_TERM(dx[0],*dx[1],*dx[2]);

            e_added  = vol*added[0];
            ke_added = vol*added[1];

            const Real sum_added = std::abs(e_added) + std::abs(ke_added);

            if (sum_added > 0)
            {
                std::cout << "Gravitational work at level "
                          << level
                          << " is "
                          << std::abs( e_added)/sum_added*100
                          << " % into (rho e) and "
                          << std::abs(ke_added)/sum_added*100
                          << " % into (KE) " << '\n';
            }
        }
    }
#endif /*GRAVITY*/

    grav_vector.clear();

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level
                      << " IS " << courno << '\n';

        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
    {
        for (int i = 0; i < S_new.nComp(); i++)
        {
            if (S_new.contains_nan(Density + i, 1, 0))
            {
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "BAD -- S_new component " << Density + i << 
                    " has NaNs after the hydro update" << std::endl;
                amrex::Abort();
            }
        }
    }
#endif

    if (add_ext_src && !strang_split)
    {
        get_old_source(prev_time, dt, ext_src_old);
        // Must compute new temperature in case it is needed in the source term
        // evaluation
        compute_new_temp();

        // Compute source at new time (no ghost cells needed)
        MultiFab ext_src_new(grids, dmap, NUM_STATE, 0);
        ext_src_new.setVal(0);

        get_new_source(prev_time, cur_time, dt, ext_src_new);

        time_center_source_terms(S_new, ext_src_old, ext_src_new, dt);

        compute_new_temp();
    } // end if (add_ext_src && !strang_split)

#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
        amrex::Abort("S_new has NaNs before the second strang call");
#endif

    // This returns updated (rho e), (rho E), and Temperature
    if (add_ext_src && strang_split)
        strang_second_step(cur_time,dt,S_new,D_new);

#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
        amrex::Abort("S_new has NaNs after the second strang call");
#endif
}
#endif
