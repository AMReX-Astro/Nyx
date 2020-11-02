#include <Nyx.H>
#include <Nyx_F.H>

using namespace amrex;

#ifndef NO_HYDRO

#ifndef AGN

void integrate_state_force(
  amrex::Box const& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& diag_eos,
  const amrex::Real* dx,
  const amrex::Real time,
  const amrex::Real a,
  const amrex::Real half_dt)
{

}

void ext_src_force(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& old_state,
  amrex::Array4<const amrex::Real> const& new_state,
  amrex::Array4<const amrex::Real> const& old_diag,
  amrex::Array4<amrex::Real> const& new_diag,
  amrex::Array4<amrex::Real> const& src,
  const amrex::Real* problo,
  const amrex::Real* dx,
  const amrex::Real time,
  const amrex::Real z,
  const amrex::Real dt)
{
	// Make a copy of the state so we can evolve it then throw it away
    amrex::FArrayBox tmp_state_fab(bx, QVAR);
    amrex::Elixir tmp_state_eli = tmp_state_fab.elixir();
	auto const& tmp_state = tmp_state_fab.array();

	//	tmp_state_fab.copy<RunOn::Device>(new_state);
    amrex::ParallelFor(bx, QVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
			tmp_state(i,j,k,n) = new_state(i,j,k,n);
		});
	amrex::Real a = 1.0 / (1.0+z);

    // Note that when we call this routine to compute the "old" source,
    //      both "old_state" and "new_state" are acutally the "old" state.
    // When we call this routine to compute the "new" source,
    //      both "old_state" is in fact the "old" state and
    //           "new_state" is in fact the "new" state

	amrex::Real half_dt = 0.5 * dt;
    integrate_state_force(bx, tmp_state, new_diag,
						  dx, time, a, half_dt);

	// Recall that this routine is called from a tiled MFIter 
	//  !   For old source: lo(:), hi(:) are the bounds of the growntilebox(src.nGrow9))
    //   For new source: lo(:), hi(:) are the bounds of the      tilebox, e.g. valid region only
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
			src(i,j,k,UMX)   = (tmp_state(i,j,k,UMX)   - new_state(i,j,k,UMX)) * a / half_dt;
			src(i,j,k,UMY)   = (tmp_state(i,j,k,UMY)   - new_state(i,j,k,UMY)) * a / half_dt;
			src(i,j,k,UMZ)   = (tmp_state(i,j,k,UMZ)   - new_state(i,j,k,UMZ)) * a / half_dt;
			src(i,j,k,UEINT) = (tmp_state(i,j,k,UEINT) - new_state(i,j,k,UEINT)) * a / half_dt;
			src(i,j,k,UEDEN) = (tmp_state(i,j,k,UEDEN) - new_state(i,j,k,UEDEN)) * a / half_dt;
	 });
}

void
Nyx::get_old_source (Real      old_time,
                     Real      dt,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_old_source()");
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real a        = get_comoving_a(old_time);
    const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);

    // We need to define these temporary multifabs because S_old and D_old only have one ghost cell.
    MultiFab Sborder, Dborder;

    Sborder.define(grids, S_old.DistributionMap(), S_old.nComp(), 4);
    Dborder.define(grids, D_old.DistributionMap(), D_old.nComp(), 4);

    FillPatch(*this, Sborder, 4, old_time, State_Type, Density, Sborder.nComp());
    FillPatch(*this, Dborder, 4, old_time, DiagEOS_Type, 0, D_old.nComp());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old, MFItInfo().SetDynamic(true).EnableTiling()); mfi.isValid(); ++mfi)
    {
		// We explicitly want to fill the ghost regions of the ext_src array
        const Box& bx = mfi.growntilebox(ext_src.nGrow());
        ext_src_force
            (bx,
             Sborder.array(mfi), Sborder.array(mfi),
             Dborder.array(mfi), Dborder.array(mfi),
             ext_src.array(mfi),
             prob_lo, dx, old_time, z, dt);
		ext_src[mfi].setVal(0.0);

        // The formulae in subroutine ctoprim assume that the source term for density is zero
        // Here we abort if it is non-zero.
        Real norm_density = ext_src[mfi].norm<RunOn::Device>(0,Density,1);
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
                     Real      dt,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_new_source()");
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

    Sborder_old.define(grids, S_old.DistributionMap(), S_old.nComp(), 4);
    Dborder_old.define(grids, D_old.DistributionMap(), D_old.nComp(), 4);

    Sborder_new.define(grids, S_new.DistributionMap(), S_new.nComp(), 4);
    Dborder_new.define(grids, D_new.DistributionMap(), D_new.nComp(), 4);

    FillPatch(*this, Sborder_old, 4, old_time, State_Type  , Density, Sborder_old.nComp());
    FillPatch(*this, Sborder_new, 4, new_time, State_Type  , Density, Sborder_new.nComp());
    FillPatch(*this, Dborder_old, 4, old_time, DiagEOS_Type, 0      , Dborder_old.nComp());
    FillPatch(*this, Dborder_new, 4, new_time, DiagEOS_Type, 0      , Dborder_new.nComp());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old, MFItInfo().SetDynamic(true).EnableTiling()); mfi.isValid(); ++mfi)
    {
        // We explicitly only want to fill the valid region
        const Box& bx = mfi.tilebox();
        ext_src_force
            (bx,
             Sborder_old.array(mfi), Sborder_new.array(mfi),
             Dborder_old.array(mfi), Dborder_new.array(mfi),
             ext_src.array(mfi),
             prob_lo, dx, new_time, z, dt);
		ext_src[mfi].setVal(0.0);
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());
}
#endif
#endif
