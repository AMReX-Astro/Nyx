#include <Nyx.H>

using namespace amrex;

#ifndef NO_HYDRO

#ifndef AGN

void Nyx::integrate_state_force(
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& state_arr,
    amrex::Array4<amrex::Real> const& diag_eos_arr,
    const amrex::Real* /*dx*/,
    const amrex::Real /*time*/,
    const amrex::Real a,
    const amrex::Real half_dt)
{

}

void Nyx::ext_src_force(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& /*old_state*/,
    amrex::Array4<const amrex::Real> const& /*new_state*/,
    amrex::Array4<const amrex::Real> const& /*old_diag*/,
    amrex::Array4<      amrex::Real> const& /*new_diag*/,
    amrex::Array4<      amrex::Real> const& /*src*/,
    const amrex::Real* problo,
    const amrex::Real* dx,
    const amrex::Real time,
    const amrex::Real z,
    const amrex::Real dt)
{
}

void
Nyx::get_old_source (Real      old_time,
                     Real      /*dt*/,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_old_source()");
    //const Real* dx      = geom.CellSize();
    //const Real* prob_lo = geom.ProbLo();
    //const Real a        = get_comoving_a(old_time);
    //const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);

    // We need to define these temporary multifabs because S_old and D_old only have one ghost cell.
    MultiFab Sborder, Dborder;

    Sborder.define(grids, S_old.DistributionMap(), S_old.nComp(), 4);
    Dborder.define(grids, D_old.DistributionMap(), D_old.nComp(), 4);

    FillPatch(*this, Sborder, 4, old_time, State_Type, Density_comp, Sborder.nComp());
    FillPatch(*this, Dborder, 4, old_time, DiagEOS_Type, 0, D_old.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Abort("C++ interface for external source not written");

        // The formulae in subroutine ctoprim assume that the source term for density is zero
        // Here we abort if it is non-zero.
        Real norm_density = ext_src[mfi].norm<RunOn::Device>(0,Density_comp,1);
        amrex::Gpu::streamSynchronize();
        if (norm_density != 0)
        {
            std::cout << "The source terms for density are non-zero" << std::endl;
            amrex::Error();
        }
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());
}

void
Nyx::get_new_source (Real      old_time,
                     Real      new_time,
                     Real      /*dt*/,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_new_source()");
    //const Real* dx      = geom.CellSize();
    //const Real* prob_lo = geom.ProbLo();
    //const Real a        = get_comoving_a(new_time);
    //const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);
    MultiFab& S_new = get_old_data(State_Type);
    MultiFab& D_new = get_old_data(DiagEOS_Type);

    // We need to define these temporary multifabs because S_old and D_old only have one ghost cell.
    MultiFab Sborder_old, Dborder_old;
    MultiFab Sborder_new, Dborder_new;

    Sborder_old.define(grids, S_old.DistributionMap(), S_old.nComp(), 4);
    Dborder_old.define(grids, D_old.DistributionMap(), D_old.nComp(), 4);

    Sborder_new.define(grids, S_new.DistributionMap(), S_new.nComp(), 4);
    Dborder_new.define(grids, D_new.DistributionMap(), D_new.nComp(), 4);

    FillPatch(*this, Sborder_old, 4, old_time, State_Type, Density_comp, Sborder_old.nComp());
    FillPatch(*this, Sborder_new, 4, new_time, State_Type, Density_comp, Sborder_new.nComp());
    FillPatch(*this, Dborder_old, 4, old_time, DiagEOS_Type, 0      , Dborder_old.nComp());
    FillPatch(*this, Dborder_new, 4, new_time, DiagEOS_Type, 0      , Dborder_new.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Abort("C++ interface for external source not written");
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());
}
#endif
#endif
