
#include "Nyx.H"
#include "Nyx_F.H"


#ifdef GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

using std::string;

#ifndef NO_HYDRO
void
Nyx::strang_hydro_fuse (Real time,
                   Real dt,
                   Real a_old,
                   Real a_new)
{
    BL_PROFILE("Nyx::strang_hydro_fuse()");

    BL_ASSERT(NUM_GROW == 4);

    Gpu::LaunchSafeGuard lsg(true);

    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();
    bool use_grav_zero = false;
    bool use_evolving_a = true;
    
    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  S_new        = get_new_data(State_Type);

    MultiFab&  D_old        = get_old_data(DiagEOS_Type);
    MultiFab&  D_new        = get_new_data(DiagEOS_Type);

#ifndef NDEBUG
    if (std::abs(time-prev_time) > (1.e-10*cur_time) )
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "strang_hydro:  prev_time = " << prev_time << std::endl;
            std::cout << "strang_hydro:       time = " <<      time << std::endl;
        }
        amrex::Abort("time should equal prev_time in strang_hydro!");
    }

    {
      amrex::Gpu::Device::streamSynchronize();
      if (S_old.contains_nan(Density, S_old.nComp(), 0))
	{
	  for (int i = 0; i < S_old.nComp(); i++)
	    {
	      if (ParallelDescriptor::IOProcessor())
                std::cout << "strang_hydro: testing component " << i << " for NaNs" << std::endl;
	      if (S_old.contains_nan(Density+i,1,0))
                amrex::Abort("S_old has NaNs in this component");
	    }
	}
    }
#endif
    
    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

    MultiFab ext_src_old(grids, dmap, NUM_STATE, NUM_GROW);
	//    std::unique_ptr<MultiFab> ext_src_old;

    //assume user-provided source is not CUDA
    if (add_ext_src)
      {
	//	ext_src_old.reset(new MultiFab(grids, dmap, NUM_STATE, NUM_GROW));
#ifndef GPU_COMPATIBLE_PROBLEM
	amrex::Gpu::Device::streamSynchronize();
	Gpu::LaunchSafeGuard lsg(false);
#endif
	get_old_source(prev_time, dt, ext_src_old);
      }

    // Define the gravity vector 
    MultiFab grav_vector(grids, dmap, BL_SPACEDIM, NUM_GROW);
    //in case do_grav==0
    grav_vector.setVal(0);

#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary(geom.periodicity());
#endif
    amrex::Gpu::Device::streamSynchronize();

    BL_PROFILE_VAR("Nyx::strang_hydro()::old_tmp_patch",old_tmp);
    // Create FAB for extended grid values (including boundaries) and fill.
    MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);

    FillPatch(*this, S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);
    FillPatch(*this, D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, D_old.nComp());

    BL_PROFILE_VAR_STOP(old_tmp);

    std::unique_ptr<MultiFab> divu_cc;
    std::unique_ptr<MultiFab> hydro_src;

#ifndef NDEBUG
    {
	amrex::Gpu::Device::streamSynchronize();
	if (S_old_tmp.contains_nan(Density, S_old_tmp.nComp(), 0)) {
          {
	    for (int i = 0; i < S_old_tmp.nComp(); i++)
	      {
	        if (ParallelDescriptor::IOProcessor())
		  std::cout << "strang_hydro: testing component " << i << " for NaNs" << std::endl;
		if (S_old_tmp.contains_nan(Density+i,1,0))
		  amrex::Abort("S_old_tmp has NaNs in this component before first strang");
              }
	    amrex::Abort("S_new has NaNs before the second strang call");
          }
      }
#endif

    bool   init_flux_register = true;
    bool add_to_flux_register = true;

    //    MultiFab hydro_src(grids, dmap, NUM_STATE, 0);

    if(hydro_convert)
      ctu_hydro_fuse(time,dt,a_old,a_new,S_old_tmp,D_old_tmp,
				 ext_src_old,*hydro_src,grav_vector,
				 init_flux_register, add_to_flux_register);
    else
      {
	amrex::Abort("Behavior not defined for fuse and hydro_convert=0");
      }
	
    D_old_tmp.clear();
    S_old_tmp.clear();

#ifndef NDEBUG
    {
    amrex::Gpu::Device::streamSynchronize();
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
      {
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();
        for (int i = 0; i < S_new[mfi].nComp(); i++)
        {
	  IntVect p_nan(D_DECL(-10, -10, -10));
            if (ParallelDescriptor::IOProcessor())
                std::cout << "strang_hydro: testing component " << i << " for NaNs" << std::endl;
            if (S_new[mfi].contains_nan(bx,Density+i,1,p_nan))
	      {
		std::cout<<"nans"<<p_nan<<std::flush<<std::endl;
		amrex::Abort("S_new has NaNs in this component after hydro");
	      }
        }
        amrex::Abort("S_new has NaNs before the second strang call");
    }
      }
#endif

    grav_vector.clear();

    if (add_ext_src)
    {
#ifndef GPU_COMPATIBLE_PROBLEM
      {
        amrex::Gpu::Device::streamSynchronize();
	Gpu::LaunchSafeGuard lsg(false);
        get_old_source(prev_time, dt, ext_src_old);
      }
#else
        get_old_source(prev_time, dt, ext_src_old);
#endif

        // Must compute new temperature in case it is needed in the source term evaluation
        compute_new_temp(S_new,D_new);

        // Compute source at new time (no ghost cells needed)
        MultiFab ext_src_new(grids, dmap, NUM_STATE, 0);

#ifndef GPU_COMPATIBLE_PROBLEM
      {
        amrex::Gpu::Device::streamSynchronize();
	Gpu::LaunchSafeGuard lsg(false);
        get_new_source(prev_time, cur_time, dt, ext_src_new);
      }
#else
        get_new_source(prev_time, cur_time, dt, ext_src_new);
#endif

        time_center_source_terms(S_new, ext_src_old, ext_src_new, dt);
	ext_src_old.clear();
        compute_new_temp(S_new,D_new);
    } // end if (add_ext_src)


#ifndef NDEBUG
    amrex::Gpu::Device::streamSynchronize();
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
      {
        for (int i = 0; i < S_new.nComp(); i++)
        {
	  IntVect p_nan(D_DECL(-10, -10, -10));
            if (ParallelDescriptor::IOProcessor())
                std::cout << "strang_hydro: testing component " << i << " for NaNs" << std::endl;
            if (S_new.contains_nan(Density+i,1,0))
	      {
		amrex::Print()<<p_nan<<std::endl;
                amrex::Abort("S_old has NaNs in this component");
	      }
        }
        amrex::Abort("S_new has NaNs before the second strang call");
      }
#endif

#ifdef HEATCOOL
    if(verbose) {
      amrex::Print()<<"Before second strang:"<<std::endl;
      amrex::Arena::PrintUsage();
    }
    // This returns updated (rho e), (rho E), and Temperature
    strang_second_step(cur_time,dt,S_new,D_new);
#endif

#ifndef NDEBUG
    amrex::Gpu::Device::streamSynchronize();
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
        amrex::Abort("S_new has NaNs after the second strang call");
#endif

    amrex::Gpu::Device::streamSynchronize();

}

#endif
