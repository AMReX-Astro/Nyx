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
#include <MemInfo.H>

#ifdef REEBER
#ifdef IN_SITU
#include <boxlib_in_situ_analysis.H>
#elif defined IN_TRANSIT
#include <InTransitAnalysis.H>
#endif
#endif

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
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
const int GimletSignal(55);
const int quitSignal(-44);


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

    while ( ! finished) {
        // Receive the sidecarSignal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter());

        switch(sidecarSignal) {
          case NyxHaloFinderSignal:
	  {
#ifdef REEBER
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars got the halo finder sidecarSignal!" << std::endl;
	    }
	    BoxArray bac;
            Geometry geom;
            int time_step, nComp(0), nGhost(0);

            // Receive the necessary data for doing analysis.
            ParallelDescriptor::Bcast(&nComp, 1, 0, ParallelDescriptor::CommunicatorInter());
	    BoxArray::RecvBoxArray(bac);
            MultiFab mf(bac, nComp, nGhost);

            MultiFab *mfSource = 0;
            MultiFab *mfDest = &mf;
            int srcComp(0), destComp(0);
            int srcNGhost(0), destNGhost(0);
            MPI_Comm commInter(ParallelDescriptor::CommunicatorInter());
            MPI_Comm commSrc(ParallelDescriptor::CommunicatorComp());
            MPI_Comm commDest(ParallelDescriptor::CommunicatorSidecar());
            bool isSrc(false);

            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost,
                                commSrc, commDest, commInter,
                                isSrc);

            Geometry::SendGeometryToSidecars(&geom);
            ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());

            // Here Reeber constructs the local-global merge trees and computes the
            // halo locations.
            runInSituAnalysis(mf, geom, time_step);

            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars completed halo finding analysis." << std::endl;
	    }
#else
      BoxLib::Abort("Nyx received halo finder signal but not compiled with Reeber");
#endif
	  }
          break;

          case GimletSignal:
	  {
#ifdef GIMLET
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars got the halo finder GimletSignal!" << std::endl;
	    }
	    BoxArray bac;
            Geometry geom;
            int time_step;
            Real new_a, omega_m, omega_b, omega_l, comoving_h;

	    BoxArray::RecvBoxArray(bac);

            MultiFab *mfSource = 0;
            MultiFab *mfDest = 0;
            int srcComp(0), destComp(0), nComp(1), nGhost(0);
            int srcNGhost(0), destNGhost(0);
            MPI_Comm commInter(ParallelDescriptor::CommunicatorInter());
            MPI_Comm commSrc(ParallelDescriptor::CommunicatorComp());
            MPI_Comm commDest(ParallelDescriptor::CommunicatorSidecar());
            bool isSrc(false);

BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  bac = " << bac << std::endl;
	    // ---- we should probably combine all of these into one MultiFab
BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  before density:   np npa npc nps = "
          << ParallelDescriptor::NProcs() << "  "
          << ParallelDescriptor::NProcsAll() << "  "
          << ParallelDescriptor::NProcsComp() << "  "
          << ParallelDescriptor::NProcsSidecar() << "  "
          << std::endl;
Array<int> dmA(bac.size() + 1, 0);
for(int i(0); i < dmA.size() - 1; ++i) {
  dmA[i] = i % ParallelDescriptor::NProcs();
}
dmA[dmA.size() - 1] = ParallelDescriptor::MyProc();
DistributionMapping dm(dmA);
            MultiFab density(bac, nComp, nGhost, dm);
BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  after density" << std::endl;
	    MultiFab temperature(bac, nComp, nGhost, dm);
	    MultiFab e_int(bac, nComp, nGhost, dm);
	    MultiFab dm_density(bac, nComp, nGhost, dm);
	    MultiFab xmom(bac, nComp, nGhost, dm);
	    MultiFab ymom(bac, nComp, nGhost, dm);
	    MultiFab zmom(bac, nComp, nGhost, dm);

BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  after mf " << std::endl;

            //MultiFab::SendMultiFabToSidecars(&density);
	    mfDest = &density;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);
BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  after copyInter density " << std::endl;

            //MultiFab::SendMultiFabToSidecars(&temperature);
	    mfDest = &temperature;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

            //MultiFab::SendMultiFabToSidecars(&e_int);
	    mfDest = &e_int;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

            //MultiFab::SendMultiFabToSidecars(&dm_density);
	    mfDest = &dm_density;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

            //MultiFab::SendMultiFabToSidecars(&xmom);
	    mfDest = &xmom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

            //MultiFab::SendMultiFabToSidecars(&ymom);
	    mfDest = &ymom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

            //MultiFab::SendMultiFabToSidecars(&zmom);
	    mfDest = &zmom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, isSrc);

BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  after copyInter " << std::endl;

            Geometry::SendGeometryToSidecars(&geom);
BoxLib::USleep(ParallelDescriptor::MyProcAll());
std::cout << ParallelDescriptor::MyProcAll() << "::_in GimletSignal:  after copyGeom " << std::endl;

            ParallelDescriptor::Bcast(&new_a, 1, 0, commInter);
            ParallelDescriptor::Bcast(&omega_m, 1, 0, commInter);
            omega_l = 1.0 - omega_m;
            ParallelDescriptor::Bcast(&omega_b, 1, 0, commInter);
            ParallelDescriptor::Bcast(&comoving_h, 1, 0, commInter);
            ParallelDescriptor::Bcast(&time_step, 1, 0, commInter);

            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "===== Sidecars got everything ..." << std::endl;
            }

      Real time1 = ParallelDescriptor::second();
      do_analysis(omega_b, omega_m, omega_l, comoving_h, new_a, density, temperature,
                  e_int, dm_density, xmom, ymom, zmom, geom, time_step);
      Real dtime = ParallelDescriptor::second() - time1;
      ParallelDescriptor::ReduceRealMax(dtime, ParallelDescriptor::IOProcessorNumber());
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl << "===== Time for Gimlet in-transit to post-process (sec): "
	          << dtime << " sec" << std::endl << std::flush;
      }
      ParallelDescriptor::Barrier();
#else
      BoxLib::Abort("Nyx received gimlet finder signal but not compiled with gimlet");
#endif
	  }
          break;

          case resizeSignal:
	  {
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
#ifdef IN_SITU
    if(ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor()) {
      std::cout << "Initializing Reeber on sidecars ... " << std::endl;
    }
    initInSituAnalysis();
#endif
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
    int nSidecarProcs(0), nSidecarProcsFromParmParse(-3);
    int prevSidecarProcs(0), minSidecarProcs(0), maxSidecarProcs(0);
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
    pp.query("minSidecarProcs", minSidecarProcs);
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
    nSidecarProcs = std::min(nSidecarProcs, maxSidecarProcs);
    nSidecarProcs = std::max(nSidecarProcs, minSidecarProcs);

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
    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    //if(nSidecarProcs > 0) {
      SidecarInit();
    //}
#endif

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
BoxLib::USleep(myProcAll/10.0);
std::cout << myProcAll << ":: _here 0" << std::endl;

    Amr *amrptr = new Amr;
    amrptr->init(strt_time,stop_time);

#if BL_USE_MPI
        // ---- initialize nyx memory monitoring
        MemInfo *mInfo = MemInfo::GetInstance();
        mInfo->LogSummary("MemInit  ");
#endif

    // ---- set initial sidecar size
    ParallelDescriptor::Bcast(&nSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "NNNNNNNN new nSidecarProcs = " << nSidecarProcs << std::endl;
      std::cout << "NNNNNNNN     maxSidecarProcs = " << maxSidecarProcs << std::endl;
    }

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
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "@@@@@@@@ after resize sidecars:  restarting event loop." << std::endl;
    }

#ifdef BL_USE_MPI
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
#endif
    //if(ParallelDescriptor::InCompGroup()) {
      //amrptr->init(strt_time,stop_time);
    //}

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "************** sizeof(Amr)      = " << sizeof(Amr) << std::endl;
      std::cout << "************** sizeof(AmrLevel) = " << sizeof(AmrLevel) << std::endl;
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

      if(ParallelDescriptor::InCompGroup()) {    // ---- stop the sidecars
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
      nSidecarProcs = std::max(nSidecarProcs, minSidecarProcs);
      ParallelDescriptor::Bcast(&nSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "NNNNNNNN new nSidecarProcs = " << nSidecarProcs << std::endl;
        std::cout << "NNNNNNNN     minSidecarProcs = " << minSidecarProcs << std::endl;
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
      if(ParallelDescriptor::InCompGroup()) {    // ---- stop the sidecars
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
