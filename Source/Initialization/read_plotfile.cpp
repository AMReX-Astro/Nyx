//
// This is the version that reads input from PlotFiles.
//
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>

#include <Nyx.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>

using namespace amrex;

void
Nyx::ReadPlotFile (bool               first,
                   const std::string& file, bool& rhoe_infile)
{
    amrex::Print() << "Reading data from plotfile: " << file << std::endl;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(file, fileType);

    if ( ! dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData&                 amrData    = dataServices.AmrDataRef();
    const Vector<std::string> plotnames  = amrData.PlotVarNames();

    // Sanity checks
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (amrData.ProbLo()[i] != geom.ProbLo()[i])
        {
            std::cout << "AmrData.prob_lo(i) = " << amrData.ProbLo()[i] << std::endl;
            std::cout << "   geom.prob_lo(i) = " <<    geom.ProbLo()[i] << std::endl;
            amrex::Error("prob_lo from plotfile doesn't match prob_lo from inputs file");
        }
        if (amrData.ProbHi()[i] != geom.ProbHi()[i])
        {
            std::cout << "AmrData.prob_hi(i) = " << amrData.ProbHi()[i] << std::endl;
            std::cout << "   geom.prob_hi(i) = " <<    geom.ProbHi()[i] << std::endl;
            amrex::Error("prob_hi from plotfile doesn't match prob_hi from inputs file");
        }
    }

    if (amrData.ProbDomain()[level] != parent->getLevel(level).Domain())
    {
        std::cout << "AmrData.domain = " << amrData.ProbDomain()[level] << std::endl;
        std::cout << "   geom.domain = " << parent->getLevel(level).Domain() << std::endl;
        amrex::Error("prob_domain from plotfile doesn't match prob_domain from inputs file");
    }

    if (amrData.boxArray(level) != grids)
    {
        std::cout << "AmrData.boxArray = " << amrData.boxArray(level) << std::endl;
        std::cout << "   grids         = " << grids << std::endl;
        amrex::Error("boxArray from plotfile doesn't match grids ");
    }

    if (amrData.FinestLevel() != parent->finestLevel()) {
        amrex::Error("finest_level from plotfile doesn't match finest_level from inputs file");
    }

    const int Nlev = parent->finestLevel() + 1;

#ifndef NO_HYDRO
    int iDensity_comp = -1, iXmom_comp = -1, iYmom_comp = -1, iZmom_comp = -1, iTEMP = -1, iNE = -1, iDensity_compE = -1;

    for (int i = 0; i < plotnames.size(); ++i)
    {
        if (plotnames[i] == "density")                  iDensity_comp  = i;
        if (plotnames[i] == "xmom")                     iXmom_comp   = i;
        if (plotnames[i] == "ymom")                     iYmom_comp   = i;
        if (plotnames[i] == "zmom")                     iZmom_comp   = i;
        if (plotnames[i] == "rho_e")                    iDensity_compE = i;
        if (plotnames[i] == "Temp")                     iTEMP  = i;
        if (plotnames[i] == "Ne")                       iNE    = i;
    }

    if ( iDensity_comp < 0 ) amrex::Abort("\"density\" is not in the plotfile");
    if ( iXmom_comp  < 0 ) amrex::Abort("\"xmom\" is not in the plotfile");
    if ( iYmom_comp  < 0 ) amrex::Abort("\"ymom\" is not in the plotfile");
    if ( iZmom_comp  < 0 ) amrex::Abort("\"zmom\" is not in the plotfile");

    if ( iDensity_compE < 0 )
    {
        if ( iTEMP < 0 ) amrex::Abort("\"rho_e\" nor \"Temp\" is in the plotfile");
        if ( iNE   < 0 ) amrex::Abort("\"rho_e\" nor \"Ne\" is in the plotfile");

        if (iNE != iTEMP + 1) amrex::Abort("We assume Ne = Temp +1");
    }
    else {rhoe_infile = true;}

    Real time = amrData.Time();
    parent->setCumTime(time);

    // Note that we only store this for level 0.
    nsteps_from_plotfile = amrData.LevelSteps()[0];

    //
    // Read density and momentum
    //
    for (int lev = 0; lev < Nlev; ++lev)
    {

        MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);
        const Box bx = grids.minimalBox();

        S_new.copy(amrData.GetGrids(lev,iDensity_comp,bx),0,Density,1);
        amrData.FlushGrids(iDensity_comp);

        S_new.copy(amrData.GetGrids(lev,iXmom_comp,bx),0,Xmom,1);
        amrData.FlushGrids(iXmom_comp);

        S_new.copy(amrData.GetGrids(lev,iYmom_comp,bx),0,Ymom,1);
        amrData.FlushGrids(iYmom_comp);

        S_new.copy(amrData.GetGrids(lev,iZmom_comp,bx),0,Zmom,1);
        amrData.FlushGrids(iZmom_comp);

        if (rhoe_infile)
        {
            S_new.copy(amrData.GetGrids(lev,iDensity_compE,bx),0,Eint_comp,1);
            amrData.FlushGrids(iDensity_compE);
        }
    }

    amrex::Print() << "Successfully read state data" << std::endl;

    //
    // Read temperature and Ne if there is no rho_e in the file
    //
    if ( ! rhoe_infile)
    {
        for (int lev = 0; lev < Nlev; ++lev)
        {
            MultiFab& D_new = parent->getLevel(lev).get_new_data(DiagEOS_Type);
            const Box bx = grids.minimalBox();

            D_new.copy(amrData.GetGrids(lev,iTEMP,bx),0,Temp_comp,1);
            amrData.FlushGrids(iTEMP);

            D_new.copy(amrData.GetGrids(lev,iNE,bx),0,Ne_comp,1);
            amrData.FlushGrids(iNE);

            amrex::Print() << "D_new.max " << D_new.norm0() << std::endl;;
        }

        amrex::Print() << "Successfully read temperature and Ne" << std::endl;
    }
#endif
}
