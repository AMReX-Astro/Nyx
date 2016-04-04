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
#include <Geometry.H>
#include <MultiFab.H>

#ifdef IN_SITU
#include <boxlib_in_situ_analysis.H>
#elif IN_TRANSIT
//#include <InTransitAnalysis.H>
#endif

#ifdef FAKE_REEBER
#include <unistd.h>
namespace
{
  void runInSituAnalysis(const MultiFab& simulation_data, const Geometry &geometry, int time_step)
  {
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "<<||||||||||>> _in runInSituAnalysis:  faking it." << std::endl;
    }
    BoxLib::USleep(0.42);  // ---- seconds
  }

  void initInSituAnalysis()
  {
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "|||||||||||||| _in initInSituAnalysis:  faking it." << std::endl;
      std::cout << "|||||||||||||| nProcs (np all comp sidecar) ="
		<< ParallelDescriptor::NProcs() << "  "
		<< ParallelDescriptor::NProcsAll() << "  "
		<< ParallelDescriptor::NProcsComp() << "  "
		<< ParallelDescriptor::NProcsSidecar() << "  "
                << std::endl;
    }
  }

}
#endif

const int NyxHaloFinderSignal(42);
const int resizeSignal(43);

// This anonymous namespace defines the workflow of the sidecars when running
// in in-transit mode.
namespace
{
#ifdef IN_TRANSIT
  static void ResizeSidecars(int newSize, Amr *amrPtr) {
    // ---- everyone meets here
    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << ParallelDescriptor::MyProcAll() << ":  _in ResizeSidecars::newSize = "
                << newSize << std::endl;
    }
    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    ParallelDescriptor::SetNProcsSidecar(newSize);
  }


  static int SidecarEventLoop() {

#ifdef BL_USE_MPI
    BL_ASSERT(NyxHaloFinderSignal != ParallelDescriptor::SidecarQuitSignal);

    bool finished(false);
    int sidecarSignal(-1);
    int myProcAll(ParallelDescriptor::MyProcAll());
    const int quitSignal(ParallelDescriptor::SidecarQuitSignal);

    while ( ! finished) {
        // Receive the sidecarSignal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter());

        switch(sidecarSignal) {
          case NyxHaloFinderSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars got the halo finder sidecarSignal!" << std::endl;
	    }
            MultiFab mf;
            Geometry geom;
            int time_step;
	    ParallelDescriptor::Barrier();
            // Receive the necessary data for doing analysis.
            MultiFab::SendMultiFabToSidecars(&mf);
            Geometry::SendGeometryToSidecars(&geom);
            ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());

            // Here Reeber constructs the local-global merge trees and computes the
            // halo locations.
            runInSituAnalysis(mf, geom, time_step);

            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars completed halo finding analysis." << std::endl;
	    }
	  }
          break;

          case resizeSignal:
	  {
	    int newSize(-1);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "_in sidecars:  Sidecars received the resize sidecarSignal." << std::endl;
            }
	    finished = true;
	  }
          break;

          case quitSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the quit sidecarSignal." << std::endl;
            }
            finished = true;
	  }
          break;

          default:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecars received bad sidecarSignal = " << sidecarSignal << std::endl;
            }
	  }
          break;
        }

    }
    if(ParallelDescriptor::IOProcessor()) {
      if(sidecarSignal == resizeSignal) {
        std::cout << "===== Sidecars exiting for resize. =====" << std::endl;
      }
      if(sidecarSignal == quitSignal) {
        std::cout << "===== Sidecars quitting. =====" << std::endl;
      }
    }
    return sidecarSignal;
#endif
  }


  static void SidecarInit() {
    if(ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor()) {
      std::cout << "Initializing Reeber on sidecars ... " << std::endl;
    }
    initInSituAnalysis();
  }
#endif
}




int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc, argv);


    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }
    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

#ifdef IN_SITU
      initInSituAnalysis();
#endif

    //
    // Don't start timing until all CPUs are ready to go.
    //
    ParallelDescriptor::Barrier("Starting main.");

    BL_COMM_PROFILE_NAMETAG("main TOP");

    int MPI_IntraGroup_Broadcast_Rank;
    int myProcAll(ParallelDescriptor::MyProcAll());
    int nSidecarProcs(0), nSidecarProcsFromParmParse(-3);
    int prevSidecarProcs(-2), maxSidecarProcs(0);
    int sidecarSignal(NyxHaloFinderSignal);
    bool resizeSidecars(false);


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

    int how(-1);
    pp.query("how",how);

#ifdef IN_TRANSIT
    pp.query("nSidecars", nSidecarProcsFromParmParse);
    pp.query("maxSidecarProcs", maxSidecarProcs);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecarProcs from parmparse = " << nSidecarProcsFromParmParse << std::endl;
    }
    if(prevSidecarProcs != nSidecarProcs) {
      resizeSidecars = true;
    }
    prevSidecarProcs = nSidecarProcs;
    if(nSidecarProcsFromParmParse >= 0) {
      if(nSidecarProcsFromParmParse >= ParallelDescriptor::NProcsAll()) {
        BoxLib::Abort("**** Error:  nSidecarProcsFromParmParse >= nProcs");
      }
      nSidecarProcs = nSidecarProcsFromParmParse;
    }

#endif

    if (strt_time < 0.0)
    {
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0)
    {
        BoxLib::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

#ifdef IN_TRANSIT
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    //if(nSidecarProcs > 0) {
      SidecarInit();
    //}
#endif

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 0" << std::endl;

DistributionMapping::InitProximityMap();
DistributionMapping::Initialize();
    //Amr *amrptr = new Amr;
    Amr *amrptr;
    if(ParallelDescriptor::InCompGroup()) {
ParallelDescriptor::Barrier();
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 1" << std::endl;

      amrptr = new Amr;
ParallelDescriptor::Barrier();
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 2" << std::endl;

      amrptr->init(strt_time,stop_time);
ParallelDescriptor::Barrier();
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 3" << std::endl;

    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "************** sizeof(Amr)      = " << sizeof(Amr) << std::endl;
      //std::cout << "************** sizeof(*amrptr)  = " << sizeof(*amrptr) << std::endl;
      std::cout << "************** sizeof(AmrLevel) = " << sizeof(AmrLevel) << std::endl;
      //std::cout << "************** sizeof(Nyx)      = " << sizeof(Nyx) << std::endl;
      //amrptr->PrintData(std::cout);
    }


ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 4" << std::endl;

    bool finished(false);

    while ( ! finished) {


    if(ParallelDescriptor::InSidecarGroup()) {

      int returnCode = SidecarEventLoop();
      if(returnCode == ParallelDescriptor::SidecarQuitSignal) {
        finished = true;
      }

    } else {

ParallelDescriptor::Barrier();
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 5" << std::endl;

    int sendstep(0), rtag(0);

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
ParallelDescriptor::Barrier();
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 10" << std::endl;

    if (amrptr->okToContinue()
           && (amrptr->levelSteps(0) < max_step || max_step < 0)
           && (amrptr->cumTime() < stop_time || stop_time < 0.0))

    {
        // ---- Do a timestep.
        amrptr->coarseTimeStep(stop_time);
    } else {
      finished = true;
    }
    const Real time_without_init = ParallelDescriptor::second() - time_before_main_loop;

    }

#ifdef IN_TRANSIT
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "************** nSidecarProcs prevSidecarProcs = " << nSidecarProcs << "  " << prevSidecarProcs << std::endl;
    }
    if(prevSidecarProcs != nSidecarProcs) {
      resizeSidecars = true;
    } else {
      resizeSidecars = false;
    }

    if(amrptr->levelSteps(0) > 1) {
      resizeSidecars = true;
    }
    if(resizeSidecars && ! finished) {

    if(ParallelDescriptor::InCompGroup()) {
      // ---- stop the sidecars
      if(nSidecarProcs > 0) {
        int sidecarSignal(resizeSignal);
        ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                  ParallelDescriptor::CommunicatorInter());
      }
    }

      // ---- test resizing the sidecars
      prevSidecarProcs = nSidecarProcs;
      nSidecarProcs = BoxLib::Random_int(ParallelDescriptor::NProcsAll()/2);
      nSidecarProcs = std::min(nSidecarProcs, maxSidecarProcs);
      ParallelDescriptor::Bcast(&nSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "NNNNNNNN new nSidecarProcs = " << nSidecarProcs << std::endl;
        std::cout << "NNNNNNNN     maxSidecarProcs = " << maxSidecarProcs << std::endl;
      }

      MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    if(nSidecarProcs < prevSidecarProcs) {
      ResizeSidecars(nSidecarProcs, amrptr);
    }

    if(ParallelDescriptor::InCompGroup()) {
      if(nSidecarProcs > prevSidecarProcs) {
        amrptr->AddProcsToSidecar(nSidecarProcs, prevSidecarProcs);
      } else {
      }
    }

    if(nSidecarProcs < prevSidecarProcs) {
      amrptr->AddProcsToComp(nSidecarProcs, prevSidecarProcs);

      amrptr->RedistributeGrids(how);
    }

    if(nSidecarProcs > prevSidecarProcs) {
      ResizeSidecars(nSidecarProcs, amrptr);
    }
    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "@@@@@@@@ after resize sidecars:  restarting event loop." << std::endl;
    }
    if(prevSidecarProcs != nSidecarProcs) {
      resizeSidecars = true;
    } else {
      resizeSidecars = false;
    }

    }

    if(finished) {
      if(ParallelDescriptor::InCompGroup()) {
        // ---- stop the sidecars
        if(nSidecarProcs > 0) {
          int sidecarSignal(ParallelDescriptor::SidecarQuitSignal);
          ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                    ParallelDescriptor::CommunicatorInter());
        }
      }
    }

#endif

    }  // ---- end while( ! finished)



    if(ParallelDescriptor::InCompGroup()) {
      // Write final checkpoint and plotfile
      if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
      }

      if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
      }
    }

#ifdef IN_TRANSIT
    if(nSidecarProcs > 0) {
      // ---- stop the sidecars
      sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
      MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
      ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                ParallelDescriptor::CommunicatorInter());
    }
#endif

    ParallelDescriptor::SetNProcsSidecar(0);

    delete amrptr;


    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, IOProc);

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
