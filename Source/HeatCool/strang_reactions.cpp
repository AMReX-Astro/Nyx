
#include <Nyx.H>

using namespace amrex;
using std::string;

#ifdef HEATCOOL    
void
Nyx::strang_first_step (Real time, Real dt, MultiFab& S_old, MultiFab& D_old)
{
    BL_PROFILE("Nyx::strang_first_step()");
    const Real strt_time = ParallelDescriptor::second();
    Real half_dt = 0.5*dt;

    const Real a = get_comoving_a(time);

    const Real z = 1.0/a - 1.0;
    if(heat_cool_type != 11)
        amrex::Abort("Invalid heating cooling type");

    if(strang_grown_box != 1)
      {
            int ierr=integrate_state_vec(S_old,       D_old,       a, half_dt);
            S_old.FillBoundary(geom.periodicity());
            D_old.FillBoundary(geom.periodicity());
            if(ierr)
              amrex::Abort("error out of integrate_state_vec");
      }
    else
      {
            if(use_typical_steps)
                amrex::ParallelDescriptor::ReduceLongMax(new_max_sundials_steps);
            int ierr=integrate_state_grownvec(S_old,       D_old,       a, half_dt);
            // Not sure how to fill patches
            //    FillPatch(*this, S_old, NUM_GROW, time, State_Type, 0, NUM_STATE);
            if(ierr)
              amrex::Abort("error out of integrate_state_vec");
      }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
          std::cout << "Nyx::strang_first_step() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}
#else
void
Nyx::strang_first_step (Real /*time*/, Real /*dt*/, MultiFab& /*S_old*/, MultiFab& /*D_old*/)
{}
#endif

#ifdef HEATCOOL    
void
Nyx::strang_second_step (Real time, Real dt, MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::strang_second_step()");
    const Real strt_time = ParallelDescriptor::second();
    Real half_dt = 0.5*dt;
    int  min_iter = 100000;
    int  max_iter =      0;

    int min_iter_grid;
    int max_iter_grid;

    // Set a at the half of the time step in the second strang
    const Real a = get_comoving_a(time-half_dt);

    MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
    reset_e_src.setVal(0.0);
    reset_internal_energy(S_new,D_new,reset_e_src);
    compute_new_temp     (S_new,D_new);

    if(heat_cool_type== 11)
      {
          if(use_typical_steps)
              amrex::ParallelDescriptor::ReduceLongMax(old_max_sundials_steps);
          int ierr=integrate_state_vec(S_new,       D_new,       a, half_dt);
          if(ierr)
              amrex::Abort("error out of integrate_state_box");
      }
    else
            amrex::Abort("Invalid heating cooling type");

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
          std::cout << "Nyx::strang_second_step() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}
#else
void
Nyx::strang_second_step (Real /*time*/, Real /*dt*/, MultiFab& /*S_new*/, MultiFab& /*D_new*/)
{}
#endif
