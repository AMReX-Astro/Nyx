#include <Nyx.H>
#include <Nyx_F.H>

using namespace amrex;

#ifndef NO_HYDRO

#ifndef AGN
void
Nyx::get_old_source (Real      old_time,
                     Real      dt,
                     MultiFab& ext_src)
{
    const Real strt     = ParallelDescriptor::second();
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real a        = get_comoving_a(old_time);
    const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);

    // We need to define these temporary multifabs because S_old and D_old only have one ghost cell.
    MultiFab Sborder, Dborder;

    Sborder.define(grids, S_old.nComp(), 4, Fab_allocate);
    Dborder.define(grids, D_old.nComp(), 4, Fab_allocate);

    FillPatch(*this, Sborder, 4, old_time, State_Type, Density, Sborder.nComp());
    FillPatch(*this, Dborder, 4, old_time, DiagEOS_Type, 0, 2);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(CA_EXT_SRC, ca_ext_src)
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(Sborder[mfi]), BL_TO_FORTRAN(Sborder[mfi]),
             BL_TO_FORTRAN(Dborder[mfi]), BL_TO_FORTRAN(Dborder[mfi]),
             BL_TO_FORTRAN(ext_src[mfi]),
             prob_lo, dx, &old_time, &z, &dt);

        // The formulae in subroutine ctoprim assume that the source term for density is zero
        // Here we abort if it is non-zero.
        if (ext_src[mfi].norm(0,Density,1) != 0)
        {
            std::cout << "The source terms for density are non-zero" << std::endl;
            amrex::Error();
        }
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Nyx::get_old_source() time = " << end << '\n';
    }

}

void
Nyx::get_new_source (Real      old_time,
                     Real      new_time,
                     Real      dt,
                     MultiFab& ext_src)
{
    const Real strt     = ParallelDescriptor::second();
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real a        = get_comoving_a(new_time);
    const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);
    MultiFab& S_new = get_old_data(State_Type);
    MultiFab& D_new = get_old_data(DiagEOS_Type);

    // We need to define these temporary multifabs because S_old and D_old only have one ghost cell.
    MultiFab Sborder_old, Dborder_old;
    MultiFab Sborder_new, Dborder_new;

    Sborder_old.define(grids, S_old.nComp(), 4, Fab_allocate);
    Dborder_old.define(grids, D_old.nComp(), 4, Fab_allocate);

    Sborder_new.define(grids, S_new.nComp(), 4, Fab_allocate);
    Dborder_new.define(grids, D_new.nComp(), 4, Fab_allocate);

    FillPatch(*this, Sborder_old, 4, old_time, State_Type  , Density, Sborder_old.nComp());
    FillPatch(*this, Sborder_new, 4, new_time, State_Type  , Density, Sborder_new.nComp());
    FillPatch(*this, Dborder_old, 4, old_time, DiagEOS_Type, 0      , 2);
    FillPatch(*this, Dborder_new, 4, new_time, DiagEOS_Type, 0      , 2);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(CA_EXT_SRC, ca_ext_src)
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(Sborder_old[mfi]), BL_TO_FORTRAN(Sborder_new[mfi]),
             BL_TO_FORTRAN(Dborder_old[mfi]), BL_TO_FORTRAN(Dborder_new[mfi]),
             BL_TO_FORTRAN(ext_src[mfi]),
             prob_lo, dx, &new_time, &z, &dt);
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Nyx::get_new_source() time = " << end << '\n';
    }

}

void
Nyx::time_center_source_terms (MultiFab& S_new,
                               MultiFab& ext_src_old,
                               MultiFab& ext_src_new,
                               Real      dt)
{
    const Real strt     = ParallelDescriptor::second();

    // Subtract off half of the old source term, and add half of the new.
    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    Real a_old = get_comoving_a(prev_time);
    Real a_new = get_comoving_a(cur_time);

    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(TIME_CENTER_SOURCES, time_center_sources)
            (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(ext_src_old[mfi]), BL_TO_FORTRAN(ext_src_new[mfi]),
             &a_old, &a_new, &dt, &print_fortran_warnings);
    }

#if 0
    if (heat_cool_type == 1)  
    {
        MultiFab& S_old = get_old_data(State_Type);
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            BL_FORT_PROC_CALL(ADJUST_HEAT_COOL, adjust_heat_cool)
                (bx.loVect(), bx.hiVect(), 
                 BL_TO_FORTRAN(S_old[mfi]), BL_TO_FORTRAN(S_new[mfi]),
                 BL_TO_FORTRAN(ext_src_old[mfi]), BL_TO_FORTRAN(ext_src_new[mfi]),
                 &a_old, &a_new, &dt);
        }
    }
#endif

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Nyx::time_center_sources() time = " << end << '\n';
    }
}
#endif
#endif
