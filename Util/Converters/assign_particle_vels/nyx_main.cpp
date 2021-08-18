
#include <iostream>
#include <iomanip>
#include <sstream>

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
#ifdef REEBER_HIST
#include <ReeberAnalysis.H> // This actually works both in situ and in-transit.
#endif
#endif

#include <AMReX_PlotFileUtil.H>
#include <Nyx_output.H>
#include <Derive.H>

std::string inputs_name = "";

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#include <fftw3-mpi.h>
#include <MakeFFTWBoxes.H>
#endif

#ifdef HENSON
#include <henson/context.h>
#include <henson/data.h>
#endif

using namespace amrex;

const int NyxHaloFinderSignal(42);
const int resizeSignal(43);
const int GimletSignal(55);
const int quitSignal(-44);

amrex::LevelBld* getLevelBld ();

    const std::string dm_chk_particle_file("DM");
    const std::string agn_chk_particle_file("AGN");
    const std::string npc_chk_particle_file("NPC");

void
nyx_main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {

    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }
    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);
    
    ParallelDescriptor::Barrier("Starting main.");

    std::cout << std::setprecision(10);

    int max_step;
    Real stop_time;
    ParmParse pp;

    max_step  = -1;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("stop_time", stop_time);

    if (max_step < 0 && stop_time < 0.0)
    {
        amrex::Abort("**** Error: either max_step or stop_time has to be positive!");
    }

    // Reeber has to do some initialization.
#ifdef REEBER
#ifdef REEBER_HIST
    reeber_int = initReeberAnalysis();
#endif
#endif

    // We hard-wire the initial time to 0
    Real strt_time =  0.0;

    Amr *amrptr = new Amr(getLevelBld());
    amrptr->init(strt_time,stop_time);

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");

    }
    amrex::Finalize();
}
/*
int
main (int argc, char* argv[])
{
    nyx_main(argc, argv);
    return 0;
}
*/
