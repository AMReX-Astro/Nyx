#include <Nyx.H>
#include <Nyx_F.H>

using namespace amrex;

using std::string;

void
Nyx::compute_hydro_sources(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new,
                           MultiFab& S_border, MultiFab& D_border, 
                           MultiFab& ext_src_old, MultiFab& hydro_src, 
                           MultiFab& grav_vector, MultiFab& divu_cc,
                           bool init_flux_register, bool add_to_flux_register) 
{
    amrex::Print() << "Computing the hydro sources ... " << std::endl;

    Nyx::construct_hydro_source(S_border, ext_src_old, hydro_src, grav_vector, 
                                a_old, a_new, dt,
                                init_flux_register,
                                add_to_flux_register);
}
