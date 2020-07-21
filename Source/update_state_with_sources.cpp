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

    int ng = 0;

    amrex::MultiFab::Copy(S_new, S_old, 0, 0, S_new.nComp(), ng);
    amrex::MultiFab::Saxpy(S_new, 1.0, hydro_src, 0, 0, S_new.nComp(), ng);

    if(verbose>1) {
        std::cout<<"S_new norm2(0)"<<S_new.norm2(0)<<std::endl;
        std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eint)<<std::endl;
        std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eden)<<std::endl;
}

}
