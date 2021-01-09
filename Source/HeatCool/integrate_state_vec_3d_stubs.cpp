#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>
#include <Nyx.H>

using namespace amrex;

/* Private function to check function return values */
static int check_retval(void* /*flagvalue*/, const char* /*funcname*/, int /*opt*/);

int Nyx::integrate_state_vec
  (amrex::MultiFab& /*S_old*/,
   amrex::MultiFab& /*D_old*/,
   const Real& /*a*/, const Real& /*delta_time*/)
{
  amrex::Abort("Using stubs file for heat_cool_type=11");
  return 0;
}

int Nyx::integrate_state_grownvec
  (amrex::MultiFab& /*S_old*/,
   amrex::MultiFab& /*D_old*/,
   const Real& /*a*/, const Real& /*delta_time*/)
{
  amrex::Abort("Using stubs file for heat_cool_type=11");
  return 0;
}
