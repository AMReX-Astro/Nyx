// @todo: deprecate windows includes

#include <winstd.H>

#include <iostream>
#include <iomanip>
#include <sstream>

#ifndef WIN32
#include <unistd.h>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

#ifdef IN_SITU
#include <boxlib_in_situ_analysis.H>
#elif IN_TRANSIT
#include <InTransitAnalysis.H>
#endif

const int NyxHaloFinderSignal = 42;

// This anonymous namespace defines the workflow of the sidecars when running
// in in-transit mode.
namespace
{
#ifdef IN_TRANSIT
  static int STATIC_SIGNAL_HANDLER (int in_signal) {

    BL_ASSERT(NyxHaloFinderSignal != ParallelDescriptor::SidecarQuitSignal);

    int out_signal = in_signal;

    if (in_signal == NyxHaloFinderSignal)
    {
      if (ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars got the halo finder signal!" << std::endl;

      MultiFab mf;
      Geometry geom;
      int time_step;
      // Receive the necessary data for doing analysis.
      MultiFab::SendMultiFabToSidecars(&mf);
      Geometry::SendGeometryToSidecars(&geom);
      ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());

      // Here Reeber constructs the local-global merge trees and computes the
      // halo locations.
      runInSituAnalysis(mf, geom, time_step);

      if (ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars completed halo finding analysis." << std::endl;
    }
    else if (in_signal != ParallelDescriptor::SidecarQuitSignal)
    {
      std::ostringstream ss_error_msg;
      ss_error_msg << "Unknown signal sent to sidecars: -----> " << in_signal << " <-----" << std::endl;
      BoxLib::Error(const_cast<const char*>(ss_error_msg.str().c_str()));
    }

    return out_signal;
  }

  static void STATIC_INIT () {
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "Initializing Reeber on sidecars ... ";
    initInSituAnalysis();
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "done." << std::endl;
    ParallelDescriptor::AddSignalHandler(STATIC_SIGNAL_HANDLER);
  }

  static void STATIC_CLEAN () {
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars all done!" << std::endl;
  }

  static void RunAtStatic () {
    BoxLib::ExecOnInitialize(STATIC_INIT);
    BoxLib::ExecOnFinalize(STATIC_CLEAN);
  };
#endif
}


int
main (int argc, char* argv[])
{
#ifdef IN_TRANSIT
    RunAtStatic();
    // Olders version of Reeber give wrong results if the # of sidecars is not
    // a power of 2. This bug has been fixed in newer versions.
    const int nSidecarProcs(256);
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
#endif

    BoxLib::Initialize(argc, argv);
#ifdef IN_TRANSIT
    if (ParallelDescriptor::InSidecarGroup()) return 0;
#endif


    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }
    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

// If we use sidecars (i.e., IN_TRANSIT), then this initialization gets called
// via BoxLib::ExecOnInitialize() through the anonymous namespace defined in
// Nyx.cpp
#ifdef IN_SITU
      initInSituAnalysis();
#endif

    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int max_step;
    Real strt_time;
    Real stop_time;
    ParmParse pp;

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);

    if (strt_time < 0.0)
    {
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0)
    {
        BoxLib::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;

    amrptr->init(strt_time, stop_time);
    //
    // If we set the regrid_on_restart flag and if we are *not* going to take
    // a time step then we want to go ahead and regrid here.
    //
    if (amrptr->RegridOnRestart())
    {
        if (    (amrptr->levelSteps(0) >= max_step ) ||
                ( (stop_time >= 0.0) &&
                  (amrptr->cumTime() >= stop_time)  )    )
        {
            //
            // Regrid only!
            //
            amrptr->RegridOnly(amrptr->cumTime());
        }
    }

    const Real time_before_main_loop = ParallelDescriptor::second();
    while (amrptr->okToContinue()
           && (amrptr->levelSteps(0) < max_step || max_step < 0)
           && (amrptr->cumTime() < stop_time || stop_time < 0.0))

    {
        //
        // Do a timestep.
        //
        amrptr->coarseTimeStep(stop_time);
    }
    const Real time_without_init = ParallelDescriptor::second() - time_before_main_loop;

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    delete amrptr;
    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, IOProc);

#ifdef IN_TRANSIT
    int signal = ParallelDescriptor::SidecarQuitSignal;
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
    ParallelDescriptor::Bcast(&signal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
#endif

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << dRunTime2 << std::endl;
        std::cout << "Run time w/o init = " << time_without_init << " sec" << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(dRunTime2);

    BoxLib::Finalize();

    return 0;
}
