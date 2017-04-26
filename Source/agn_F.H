#include <AMReX_BLFort.H>

extern "C"
{
    void nyx_compute_overlap(const amrex::Real*, int ns, int np, int* my_id, 
                             const amrex::Real* ghosts, int ng, const amrex::Real* dx); }
//                           const amrex::Real* ghosts, int ng, const amrex::Real* dx,
//                           const amrex::Real* density, const int* dlo, const int* dhi); }
