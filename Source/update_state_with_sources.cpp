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
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Updating state with the hydro sources ... " << std::endl;
    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));
    if(verbose>1) {
        std::cout<<"hydro_src norm2(0)"<<hydro_src.norm2(0)<<std::endl;
	std::cout<<"hydro_src norm2(Eint)"<<hydro_src.norm2(Eint)<<std::endl;
	std::cout<<"hydro_src norm2(Eint)"<<hydro_src.norm2(Eden)<<std::endl;
}


    if(hydro_convert)
    {
      Gpu::LaunchSafeGuard lsg(true);
      int print_fortran_warnings_tmp=print_fortran_warnings;
      int do_grav_tmp=do_grav;

      FArrayBox sum_state, divu_cc_small;
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

      const auto fab_ext_src_old = ext_src_old.array(mfi);
      const auto fab_hydro_src = hydro_src.array(mfi);
      const auto fab_grav = grav.array(mfi);
      const auto fab_S_old = S_old.array(mfi);
      const auto fab_S_new = S_new.array(mfi);
      
      const Box& bx = mfi.tilebox();
      const Box& obx = amrex::grow(bx, 1);

      sum_state.resize(obx, AMREX_SPACEDIM);
      Elixir elix_s = sum_state.elixir();
      const auto fab_sum_state = sum_state.array();

      divu_cc_small.resize(bx, 1);
      Elixir elix_divu_cc_small = divu_cc_small.elixir();
      const auto fab_divu_cc = divu_cc_small.array();

      S_old[mfi].prefetchToDevice();
      S_new[mfi].prefetchToDevice();
      ext_src_old[mfi].prefetchToDevice();
      hydro_src[mfi].prefetchToDevice();
      divu_cc_small.prefetchToDevice();
      grav[mfi].prefetchToDevice();

      sum_state.prefetchToDevice();

	AMREX_LAUNCH_DEVICE_LAMBDA(bx,tbx,
	{
	  ca_fort_update_state (
		  tbx.loVect(), tbx.hiVect(),
		  BL_ARR4_TO_FORTRAN(fab_S_old),
		  BL_ARR4_TO_FORTRAN(fab_S_new),
		  BL_ARR4_TO_FORTRAN(fab_ext_src_old),
		  BL_ARR4_TO_FORTRAN(fab_hydro_src),
		  BL_ARR4_TO_FORTRAN(fab_divu_cc),
		  BL_ARR4_TO_FORTRAN_3D(fab_sum_state),
		  &dt, &a_old, &a_new, &print_fortran_warnings_tmp);

	  // Note this increments S_new, it doesn't add source to S_old
	  // However we create the source term using rho_old
	  if (do_grav_tmp)
	    ca_fort_add_grav_source (
		    tbx.loVect(), tbx.hiVect(),
		    BL_ARR4_TO_FORTRAN(fab_S_old),
		    BL_ARR4_TO_FORTRAN(fab_S_new),
		    BL_ARR4_TO_FORTRAN(fab_grav),
		    &dt, &a_old, &a_new);
	});
    }

    }
    else
      {
	Gpu::synchronize();
	Gpu::LaunchSafeGuard lsg(false);
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
   if(verbose>1) {
	std::cout<<"S_new norm2(0)"<<S_new.norm2(0)<<std::endl;
	std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eint)<<std::endl;
	std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eden)<<std::endl;
}

}
