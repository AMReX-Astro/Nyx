
#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;
using std::string;

void
Nyx::strang_first_step (Real time, Real dt, MultiFab& S_old, MultiFab& D_old)
{
    BL_PROFILE("Nyx::strang_first_step()");
    
    Gpu::LaunchSafeGuard lsg(true);
    const Real strt_time = ParallelDescriptor::second();
    Real half_dt = 0.5*dt;

    const Real a = get_comoving_a(time);

#ifndef FORCING
    {
      const Real z = 1.0/a - 1.0;
      fort_interp_to_this_z(&z);
      /*      int neq=1;
      AMREX_LAUNCH_DEVICE_LAMBDA(neq,i,
				 {
				   ca_interp_to_this_z(&z);
				 });*/
    }
#endif

    /////////////////////Consider adding ifdefs for whether CVODE is compiled in for these statements
    if(heat_cool_type == 3 || heat_cool_type==5 || heat_cool_type==7 || heat_cool_type==9)
      {
	Gpu::streamSynchronize();
	Gpu::LaunchSafeGuard lsg(false);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
	  {
	    // Note that this "bx" includes the grow cells 
	    const Box& bx = mfi.growntilebox(S_old.nGrow());

	    int  min_iter = 100000;
	    int  max_iter =      0;

	    integrate_state
                (bx.loVect(), bx.hiVect(), 
                 BL_TO_FORTRAN(S_old[mfi]),
                 BL_TO_FORTRAN(D_old[mfi]),
                 &a, &half_dt, &min_iter, &max_iter);

#ifdef AMREX_DEBUG
	    /*        if (S_old[mfi].contains_nan())
		      amrex::Abort("state has NaNs after the first strang call");*/
#endif

    }

      }
    else if(heat_cool_type == 4)
      {
	int ierr=integrate_state_grownexact(S_old, D_old, a, half_dt);
	if(ierr)
	  amrex::Abort("error out of integrate_state_exact");

      }
    else if(strang_grown_box != 1)
      {
	if(heat_cool_type== 10)
	  {
	    //#ifdef CVODE_LIBS
	    int ierr=integrate_state_box(S_old,       D_old,       a, half_dt);
	    S_old.FillBoundary(geom.periodicity());
	    D_old.FillBoundary(geom.periodicity());
	    if(ierr)
	      amrex::Abort("error out of integrate_state_box");
	  }
	else if(heat_cool_type== 11)
	  {
	    //#ifdef CVODE_LIBS
	    int ierr=integrate_state_vec(S_old,       D_old,       a, half_dt);
	    S_old.FillBoundary(geom.periodicity());
	    D_old.FillBoundary(geom.periodicity());
	    // Not sure how to fill patches
	    //    FillPatch(*this, S_old, NUM_GROW, time, State_Type, 0, NUM_STATE);
	    if(ierr)
	      amrex::Abort("error out of integrate_state_vec");
	  }
	else if(heat_cool_type== 12)
	  {
	    //#ifdef CVODE_LIBS
	    int ierr=integrate_state_cell(S_old,       D_old,       a, half_dt);
	    S_old.FillBoundary(geom.periodicity());
	    D_old.FillBoundary(geom.periodicity());
	    if(ierr)
	      amrex::Abort("error out of integrate_state_cell");
	  }
	else
	  amrex::Abort("Invalid heating cooling type");
      }
    else
      {
	if(heat_cool_type== 10)
	  {
	    //#ifdef CVODE_LIBS
	    int ierr=integrate_state_grownbox(S_old,       D_old,       a, half_dt);
	    if(ierr)
	      amrex::Abort("error out of integrate_state_box");
	  }
	else if(heat_cool_type== 11)
	  {
	    //#ifdef CVODE_LIBS
	    if(use_typical_steps)
	        amrex::ParallelDescriptor::ReduceLongMax(new_max_sundials_steps);
	    int ierr=integrate_state_grownvec(S_old,       D_old,       a, half_dt);
	    // Not sure how to fill patches
	    //    FillPatch(*this, S_old, NUM_GROW, time, State_Type, 0, NUM_STATE);
	    if(ierr)
	      amrex::Abort("error out of integrate_state_vec");
	  }
	else if(heat_cool_type== 12)
	  {
	    //#ifdef CVODE_LIBS
	    int ierr=integrate_state_growncell(S_old,       D_old,       a, half_dt);
	    if(ierr)
	      amrex::Abort("error out of integrate_state_cell");
	  }
	else
	  amrex::Abort("Invalid heating cooling type");
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

void
Nyx::strang_second_step (Real time, Real dt, MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::strang_second_step()");

    Gpu::LaunchSafeGuard lsg(true);
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

#ifndef FORCING
    {
      const Real z = 1.0/a - 1.0;
      fort_interp_to_this_z(&z);
      /*      int neq=1;
      AMREX_LAUNCH_DEVICE_LAMBDA(neq,i,
      {
      ca_interp_to_this_z(&z);
      });*/
    }
#endif

    /////////////////////Consider adding ifdefs for whether CVODE is compiled in for these statements
    if(heat_cool_type == 3 || heat_cool_type==5 || heat_cool_type==7 || heat_cool_type==9)
      {
	Gpu::streamSynchronize();
	Gpu::LaunchSafeGuard lsg(false);
#ifdef _OPENMP
#pragma omp parallel private(min_iter_grid,max_iter_grid) reduction(min:min_iter) reduction(max:max_iter)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        // Here bx is just the valid region
        const Box& bx = mfi.tilebox();

        min_iter_grid = 100000;
        max_iter_grid =      0;

        integrate_state
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(D_new[mfi]),
             &a, &half_dt, &min_iter_grid, &max_iter_grid);

	/*        if (S_new[mfi].contains_nan<RunOn::Device>(bx,0,S_new.nComp()))
        {
            amrex::Gpu::streamSynchronize();
            std::cout << "NANS IN THIS GRID " << bx << std::endl;
	    }*/

        min_iter = std::min(min_iter,min_iter_grid);
        max_iter = std::max(max_iter,max_iter_grid);
    }

      }
    else if(heat_cool_type == 4)
      {
	int ierr=integrate_state_exact(S_new, D_new, a, half_dt);
	if(ierr)
	  amrex::Abort("error out of integrate_state_exact");
      }
    else if(heat_cool_type== 10)
      {
	//#ifdef CVODE_LIBS
    int ierr=integrate_state_box(S_new,       D_new,       a, half_dt);
    if(ierr)
      amrex::Abort("error out of integrate_state_box");
      }
    else if(heat_cool_type== 11)
      {
	//#ifdef CVODE_LIBS
	  if(use_typical_steps)
	      amrex::ParallelDescriptor::ReduceLongMax(old_max_sundials_steps);
    int ierr=integrate_state_vec(S_new,       D_new,       a, half_dt);
    if(ierr)
      amrex::Abort("error out of integrate_state_box");
      }
    else if(heat_cool_type== 12)
      {
	//#ifdef CVODE_LIBS
    int ierr=integrate_state_cell(S_new,       D_new,       a, half_dt);
    if(ierr)
      amrex::Abort("error out of integrate_state_cell");
      }
    else
            amrex::Abort("Invalid heating cooling type");

    if(heat_cool_type == 3 || heat_cool_type==5 || heat_cool_type==7 || heat_cool_type==9)
    {
        ParallelDescriptor::ReduceIntMax(max_iter);
	ParallelDescriptor::ReduceIntMin(min_iter);
    }

    if (heat_cool_type == 1)
        if (ParallelDescriptor::IOProcessor())
	  {
	    ParallelDescriptor::ReduceIntMax(max_iter);
	    ParallelDescriptor::ReduceIntMin(min_iter);
            std::cout << "Min/Max Number of Iterations in Second Strang: " << min_iter << " " << max_iter << std::endl;
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
	  std::cout << "Nyx::strang_second_step() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }
}
