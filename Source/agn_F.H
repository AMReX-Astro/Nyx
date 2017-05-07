#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void nyx_compute_overlap(const int* np, void* particles, 
                             const int* ng, void*    ghosts, const amrex::Real* dx);

    void agn_merge_particles(const int* np, void* particles, 
                             const int* ng, void*    ghosts, const amrex::Real* dx);

    void agn_particle_velocity(const int* np, void* particles,
                               const amrex::Real* state_old, const int* sold_lo, const int* sold_hi,
                               const amrex::Real* state_new, const int* snew_lo, const int* snew_hi,
                               const amrex::Real* dx, const int* add_energy);

    void agn_accrete_mass(const int* np, void* particles,
                          const amrex::Real* state, const int* s_lo, const int* s_hi,
                          const amrex::Real* eps_rad, const amrex::Real* dt, const amrex::Real* dx);
#ifdef __cplusplus
}
#endif
