// @todo: deprecate windows includes


#include <iostream>
#include <iomanip>
#include <sstream>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#ifdef BL_USE_MPI
#include <MemInfo.H>
#endif
#include <Nyx.H>

#ifdef REEBER
#include <ReeberAnalysis.H> // This actually works both in situ and in-transit.
#endif

#include "Nyx_output.H"

std::string inputs_name = "";

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#include <fftw3-mpi.h>
#include <MakeFFTWBoxes.H>
#endif

using namespace amrex;

const int NyxHaloFinderSignal(42);
const int resizeSignal(43);
const int GimletSignal(55);
const int quitSignal(-44);


// This anonymous namespace defines the workflow of the sidecars when running
// in in-transit mode.
namespace
{
  static void ResizeSidecars(int newSize) {
#ifdef BL_USE_MPI
    // ---- everyone meets here
    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << ParallelDescriptor::MyProcAll() << ":  _in ResizeSidecars::newSize = "
                << newSize << std::endl;
    }
    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    ParallelDescriptor::SetNProcsSidecars(newSize);
#endif /* BL_USE_MPI */
  }


  static int SidecarEventLoop() {

#ifdef BL_USE_MPI
    BL_ASSERT(NyxHaloFinderSignal != quitSignal);

    bool finished(false);
    int sidecarSignal(-1);
    int whichSidecar(0);  // ---- this is sidecar zero
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "SSSSSSSS:  Starting SidecarEventLoop." << std::endl;
    }

    while ( ! finished) {
        // ---- Receive the sidecarSignal from the compute group.
        if(ParallelDescriptor::IOProcessor()) {
          std::cout << "SSSSSSSS:  waiting for signal from comp..." << std::endl;
        }
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter(whichSidecar));

        switch(sidecarSignal) {
          case NyxHaloFinderSignal:
          {
#ifdef REEBER
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars got the halo finder sidecarSignal!" << std::endl;
            }

            Geometry geom;
            Geometry::SendGeometryToSidecar(&geom, whichSidecar);

            int time_step, nComp(0), nGhost(0);
            int do_analysis;

            // Receive the necessary data for doing analysis.
            ParallelDescriptor::Bcast(&nComp, 1, 0, ParallelDescriptor::CommunicatorInter(whichSidecar));

            // Get desired box array and distribution mapping from Reeber
            BoxArray ba;
            DistributionMapping dm;
            getAnalysisDecomposition(geom, ParallelDescriptor::NProcsSidecar(whichSidecar), ba, dm);

            MultiFab mf(ba, dm, nComp + 1, nGhost);

            MultiFab *mfSource = 0;
            MultiFab *mfDest = &mf;
            int srcComp(0), destComp(1);
            int srcNGhost(0), destNGhost(0);
            const MPI_Comm &commSrc = ParallelDescriptor::CommunicatorComp();
            const MPI_Comm &commDest = ParallelDescriptor::CommunicatorSidecar();
            const MPI_Comm &commInter = ParallelDescriptor::CommunicatorInter(whichSidecar);
            const MPI_Comm &commBoth = ParallelDescriptor::CommunicatorBoth(whichSidecar);
            bool isSrc(false);

            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost,
                                commSrc, commDest, commInter, commBoth,
                                isSrc);

            mf.setVal(0, 0, 1, 0);
            for (int comp = 1; comp < mf.nComp(); ++comp)
                MultiFab::Add(mf, mf, comp, 0, 1, 0);


            ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter(whichSidecar));
            ParallelDescriptor::Bcast(&do_analysis, 1, 0, ParallelDescriptor::CommunicatorInter(whichSidecar));

            // Here Reeber constructs the local-global merge trees and computes the
            // halo locations.
            runReeberAnalysis(mf, geom, time_step, bool(do_analysis));

            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars completed halo finding analysis." << std::endl;
            }
#else
            amrex::Abort("Nyx received halo finder signal but not compiled with Reeber");
#endif /* REEBER */
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

            BoxArray::RecvBoxArray(bac, whichSidecar);

            MultiFab *mfSource = 0;
            MultiFab *mfDest = 0;
            int srcComp(0), destComp(0), nComp(1), nGhost(0);
            int srcNGhost(0), destNGhost(0);
            const MPI_Comm &commSrc = ParallelDescriptor::CommunicatorComp();
            const MPI_Comm &commDest = ParallelDescriptor::CommunicatorSidecar();
            const MPI_Comm &commInter = ParallelDescriptor::CommunicatorInter(whichSidecar);
            const MPI_Comm &commBoth = ParallelDescriptor::CommunicatorBoth(whichSidecar);
            bool isSrc(false);

            // ---- we should probably combine all of these into one MultiFab
            MultiFab density(bac, nComp, nGhost);
            MultiFab temperature(bac, nComp, nGhost);
            MultiFab e_int(bac, nComp, nGhost);
            MultiFab dm_density(bac, nComp, nGhost);
            MultiFab xmom(bac, nComp, nGhost);
            MultiFab ymom(bac, nComp, nGhost);
            MultiFab zmom(bac, nComp, nGhost);

            mfDest = &density;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &temperature;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &e_int;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &dm_density;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &xmom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &ymom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            mfDest = &zmom;
            MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                                srcNGhost, destNGhost, commSrc, commDest, commInter, commBoth, isSrc);

            Geometry::SendGeometryToSidecar(&geom, whichSidecar);

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
            amrex::Abort("Nyx received Gimlet signal but not compiled with Gimlet");
#endif /* GIMLET */
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
#else
    return 0;
#endif /* BL_USE_MPI */
    }

// The following function does not seem to be used
//    static void SidecarInit() {
//#ifdef REEBER
//      if(ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor()) {
//        std::cout << "Initializing Reeber on sidecars ... " << std::endl;
//      } else if (ParallelDescriptor::IOProcessor()) {
//        std::cout << "Initializing Reeber in situ ... " << std::endl;
//      }
//      // Reeber reads its ParmParse stuff here so we don't need to do any of it
//      // in Nyx proper.
//      reeber_int = initReeberAnalysis();
//      std::cout << "SidecarInit(): reeber_int = " << reeber_int << std::endl;
//#endif
//    }
}




int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);


    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }
    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    //
    // Don't start timing until all CPUs are ready to go.
    //
    ParallelDescriptor::Barrier("Starting main.");

    BL_COMM_PROFILE_NAMETAG("main TOP");

#ifdef BL_USE_MPI
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
    int nSidecarProcsFromParmParse(-3);
    Nyx::nSidecarProcs = 0;
    int prevSidecarProcs(0);
    int sidecarSignal(NyxHaloFinderSignal);
    int resizeSidecars(false);  // ---- instead of bool for bcast
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

    int how(-1);
    pp.query("how",how);

#ifdef BL_USE_MPI
    // Set up sidecars if user wants them.
    if (pp.query("nSidecars", nSidecarProcsFromParmParse))
    {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "nSidecarProcs from parmparse = " << nSidecarProcsFromParmParse << std::endl;
      }
      resizeSidecars = !(prevSidecarProcs == Nyx::nSidecarProcs);
      prevSidecarProcs = Nyx::nSidecarProcs;
      if(nSidecarProcsFromParmParse >= 0) {
        if(nSidecarProcsFromParmParse >= ParallelDescriptor::NProcsAll()) {
          amrex::Abort("**** Error:  nSidecarProcsFromParmParse >= nProcs");
        }
        Nyx::nSidecarProcs = nSidecarProcsFromParmParse;
      }
    }
#endif

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0)
    {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    // Reeber has to do some initialization.
#ifdef REEBER
    reeber_int = initReeberAnalysis();
#endif

    if (Nyx::nSidecarProcs > 0)
      Nyx::forceParticleRedist = true;

    Amr *amrptr = new Amr;
    amrptr->init(strt_time,stop_time);

#if BL_USE_MPI
    // ---- initialize nyx memory monitoring
    MemInfo *mInfo = MemInfo::GetInstance();
    mInfo->LogSummary("MemInit  ");
#endif

    // ---- set initial sidecar size
    ParallelDescriptor::Bcast(&Nyx::nSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "IIIIIIII new nSidecarProcs = " << Nyx::nSidecarProcs << std::endl;
    }

#ifdef BL_USE_MPI
    if(Nyx::nSidecarProcs < prevSidecarProcs) {
      ResizeSidecars(Nyx::nSidecarProcs);
      amrptr->AddProcsToComp(Nyx::nSidecarProcs, prevSidecarProcs);
      amrptr->RedistributeGrids(how);
    } else if (Nyx::nSidecarProcs > prevSidecarProcs) {
      if(ParallelDescriptor::InCompGroup()) {
        amrptr->AddProcsToSidecar(Nyx::nSidecarProcs, prevSidecarProcs);
      }
      ResizeSidecars(Nyx::nSidecarProcs);
    }
#endif
    const Real time_before_main_loop = ParallelDescriptor::second();

#ifdef USE_CVODE
    Nyx::alloc_simd_vec();
#endif

    bool finished(false);

    while ( ! finished) {

      Nyx::forceParticleRedist = true;

      if(ParallelDescriptor::InSidecarGroup()) {  // ------------------- start sidecars
#ifdef BL_USE_MPI

        int returnCode = SidecarEventLoop();
        if(returnCode == quitSignal) {
          finished = true;
        }
        resizeSidecars = (returnCode == resizeSignal);

#endif
      } else {  // ----------------------------------------------------- start comp

        // If we set the regrid_on_restart flag and if we are *not* going to take
        // a time step then we want to go ahead and regrid here.
        //
        if (amrptr->RegridOnRestart()) {
          if (    (amrptr->levelSteps(0) >= max_step ) ||
                  ( (stop_time >= 0.0) &&
                    (amrptr->cumTime() >= stop_time)  )    )
          {
              // Regrid only!
              amrptr->RegridOnly(amrptr->cumTime());
          }
        }

        if (amrptr->okToContinue()
             && (amrptr->levelSteps(0) < max_step || max_step < 0)
             && (amrptr->cumTime() < stop_time || stop_time < 0.0))

        {
          amrptr->coarseTimeStep(stop_time);          // ---- Do a timestep.
        } else {
          finished = true;
#ifdef BL_USE_MPI
          resizeSidecars = false;
#endif
        }

#ifdef BL_USE_MPI
        if((finished || resizeSidecars) && Nyx::nSidecarProcs > 0 && prevSidecarProcs > 0) {
          // ---- stop the sidecars
          int sidecarSignal(-1);
          if(finished) {
            sidecarSignal = quitSignal;
          } else if(resizeSidecars) {
            sidecarSignal = resizeSignal;
          }
          int whichSidecar(0);
          ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                    ParallelDescriptor::CommunicatorInter(whichSidecar));
        }
#endif

      }  // ---------------- end start comp



#ifdef BL_USE_MPI
      if(resizeSidecars) {    // ---- both comp and sidecars are here
        ParallelDescriptor::Bcast(&prevSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
        ParallelDescriptor::Bcast(&Nyx::nSidecarProcs, 1, 0, ParallelDescriptor::CommunicatorAll());
        if(ParallelDescriptor::InCompGroup()) {
          if(ParallelDescriptor::IOProcessor()) {
            std::cout << "NNNNNNNN new nSidecarProcs    = " << Nyx::nSidecarProcs    << std::endl;
            std::cout << "NNNNNNNN     prevSidecarProcs = " << prevSidecarProcs << std::endl;
          }
        }
        Nyx::forceParticleRedist = true;

        if(Nyx::nSidecarProcs < prevSidecarProcs) {
          ResizeSidecars(Nyx::nSidecarProcs);
        }

        if(Nyx::nSidecarProcs > prevSidecarProcs) {
          if(ParallelDescriptor::InCompGroup()) {
            amrptr->AddProcsToSidecar(Nyx::nSidecarProcs, prevSidecarProcs);
          }
        }

        if(Nyx::nSidecarProcs < prevSidecarProcs) {

          amrptr->AddProcsToComp(Nyx::nSidecarProcs, prevSidecarProcs);
          amrptr->RedistributeGrids(how);
        }

        if(Nyx::nSidecarProcs > prevSidecarProcs) {
          ResizeSidecars(Nyx::nSidecarProcs);
        }
        if(ParallelDescriptor::IOProcessor()) {
          std::cout << "@@@@@@@@ after resize sidecars:  restarting event loop." << std::endl;
        }
      }
#endif
    }  // ---- end while( ! finished)

#ifdef USE_CVODE
    Nyx::dealloc_simd_vec();
#endif

    const Real time_without_init = ParallelDescriptor::second() - time_before_main_loop;
    if (ParallelDescriptor::IOProcessor()) std::cout << "Time w/o init: " << time_without_init << std::endl;


    if(ParallelDescriptor::InCompGroup()) {
      // Write final checkpoint and plotfile
      if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
      }
      if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
      }
    }

#ifdef BL_USE_MPI
    if(Nyx::nSidecarProcs > 0) {    // ---- stop the sidecars
      sidecarSignal = quitSignal;
      int whichSidecar(0);
      ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                ParallelDescriptor::CommunicatorInter(whichSidecar));
    }

    ParallelDescriptor::SetNProcsSidecars(0);
#endif

#ifdef AMREX_USE_CATALYST
    amrptr->computeInSitu();
#endif

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
    }

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(dRunTime2);

    amrex::Finalize();

    return 0;
}
