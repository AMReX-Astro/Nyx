#include <AMReX_BLFort.H>

extern "C"
{
    void nyx_compute_overlap(const amrex::Real*, int ns, int np, int* my_id, 
                             const amrex::Real* ghosts, int ng, const amrex::Real* dx);

    void agn_particle_velocity(const amrex::Real*, int ns, int np, 
                               const amrex::Real* state_old, const int* sold_lo, const int* sold_hi,
                               const amrex::Real* state_new, const int* snew_lo, const int* snew_hi,
                               int& start_comp, const amrex::Real* dx);
}
