
#include <Nyx.H>

using namespace amrex;
using std::string;

void
Nyx::sdc_reactions (MultiFab& S_old, MultiFab& S_new, MultiFab& D_new, 
                    MultiFab& hydro_src, MultiFab& IR,
                    Real delta_time, Real a_old, Real a_new, int sdc_iter)
{
    BL_PROFILE("Nyx::sdc_reactions()");
#ifdef HEATCOOL
    const Real* dx = geom.CellSize();

    // First reset internal energy before call to compute_temp
    MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
    reset_e_src.setVal(0.0);

    reset_internal_energy(S_new,D_new,reset_e_src);
    compute_new_temp     (S_new,D_new);
    
    Gpu::LaunchSafeGuard lsg(true);
    const Real strt_time = ParallelDescriptor::second();

    const Real z = 1.0/a_old - 1.0;

    if(heat_cool_type != 11)
        amrex::Abort("Invalid heating cooling type");

    if(use_typical_steps)
        amrex::ParallelDescriptor::ReduceLongMax(new_max_sundials_steps);
    int ierr=integrate_state_struct(S_old,S_new, D_new, hydro_src, IR, reset_e_src,  a_old, a_new, delta_time, sdc_iter);
    if(ierr)
        amrex::Abort("error out of integrate_state_vec");

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
#endif
}
