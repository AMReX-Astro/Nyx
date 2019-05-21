#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;

using std::string;

void
Nyx::update_state_with_sources( MultiFab& S_old, MultiFab& S_new, 
                                MultiFab& ext_src_old, MultiFab& hydro_src, 
                                MultiFab& grav, MultiFab& divu_cc,
                                amrex::Real dt, amrex::Real a_old, amrex::Real a_new)
{
    BL_PROFILE("Nyx::update_state_with_sources()");
    amrex::Print() << "Updating state with the hydro sources ... " << std::endl;

    if(hydro_convert)
    {
      int print_fortran_warnings_tmp=print_fortran_warnings;
      int do_grav_tmp=do_grav;

    FArrayBox sum_state;
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

      FArrayBox* fab_S_old = S_old.fabPtr(mfi);
      FArrayBox* fab_S_new = S_new.fabPtr(mfi);
      FArrayBox* fab_ext_src_old = ext_src_old.fabPtr(mfi);
      FArrayBox* fab_hydro_src = hydro_src.fabPtr(mfi);
      FArrayBox* fab_divu_cc = divu_cc.fabPtr(mfi);
      FArrayBox* fab_grav = grav.fabPtr(mfi);
      
      const Box& bx = mfi.tilebox();
      const Box& obx = amrex::grow(bx, 1);

      sum_state.resize(obx, AMREX_SPACEDIM);
      Elixir elix_s = sum_state.elixir();
      const auto fab_sum_state = sum_state.array();

	AMREX_LAUNCH_DEVICE_LAMBDA(bx,tbx,
	{
	  ca_fort_update_state (
		  tbx.loVect(), tbx.hiVect(),
		  BL_TO_FORTRAN(*fab_S_old),
		  BL_TO_FORTRAN(*fab_S_new),
		  BL_TO_FORTRAN(*fab_ext_src_old),
		  BL_TO_FORTRAN(*fab_hydro_src),
		  BL_TO_FORTRAN(*fab_divu_cc),
		  fab_sum_state.p,&((fab_sum_state).begin.x),amrex::GpuArray<int,3>{(fab_sum_state).end.x-1,(fab_sum_state).end.y-1,(fab_sum_state).end.z-1}.data(),
		  &dt, &a_old, &a_new, &print_fortran_warnings_tmp);

	  // Note this increments S_new, it doesn't add source to S_old
	  // However we create the source term using rho_old
	  if (do_grav_tmp)
	    ca_fort_add_grav_source (
		    tbx.loVect(), tbx.hiVect(),
		    BL_TO_FORTRAN(*fab_S_old),
		    BL_TO_FORTRAN(*fab_S_new),
		    BL_TO_FORTRAN(*fab_grav),
		    &dt, &a_old, &a_new);
	});
    }

    }
    else
      {
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        fort_update_state (
             bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(S_old[mfi]),
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(ext_src_old[mfi]),
             BL_TO_FORTRAN(hydro_src[mfi]),
             BL_TO_FORTRAN(divu_cc[mfi]),
             &dt, &a_old, &a_new, &print_fortran_warnings);

        // Note this increments S_new, it doesn't add source to S_old
        // However we create the source term using rho_old
        if (do_grav)
           fort_add_grav_source (
                bx.loVect(), bx.hiVect(), 
                BL_TO_FORTRAN(S_old[mfi]),
                BL_TO_FORTRAN(S_new[mfi]),
                BL_TO_FORTRAN(grav[mfi]),
                &dt, &a_old, &a_new);
   }
      }
}
