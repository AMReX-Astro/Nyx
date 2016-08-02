// @todo: deprecate windows includes

#include <winstd.H>

#include <iostream>
#include <iomanip>

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
#endif

#include "Nyx_output.H"

std::string inputs_name = "";

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
