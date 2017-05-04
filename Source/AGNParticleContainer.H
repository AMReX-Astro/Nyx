
#ifndef _AGNParticleContainer_H_
#define _AGNParticleContainer_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

#include "NyxParticleContainer.H"

namespace {

struct GhostCommTag {

    GhostCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const GhostCommTag& l, const GhostCommTag& r) {
    return (l.proc_id < r.proc_id ||
           (l.proc_id == r.proc_id && l.grid_id < r.grid_id) ||
           (l.proc_id == r.proc_id && l.grid_id == r.grid_id && l.tile_id < r.tile_id ));
}

}

class AGNParticleContainer
    : public NyxParticleContainer<5>
{
public:

    using MyParIter = amrex::ParIter<5>;
    using PairIndex = std::pair<int, int>;
    using GhostCommMap = std::map<GhostCommTag, amrex::Array<char> >;

    AGNParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<5>(amr)
    { }

    virtual ~AGNParticleContainer () {}

    const int NumberOfParticles(MyParIter& pti) { return pti.GetArrayOfStructs().size(); }

    void AddOneParticle (int lev, int grid, int tile, amrex::Real mass, 
                         amrex::Real x, amrex::Real y, amrex::Real z)
    {
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];
        AddOneParticle(particle_tile, mass, x, y, z);
    }

    void AddOneParticle (ParticleTileType& particle_tile, amrex::Real mass, 
                         amrex::Real x, amrex::Real y, amrex::Real z)
    {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        p.pos(0) = x;
        p.pos(1) = y;
        p.pos(2) = z;

        // Set mass 
        p.rdata(0) = mass;

        // Zero initial velocity
        p.rdata(1) = 0.;
        p.rdata(2) = 0.;
        p.rdata(3) = 0.;

        // Zero initial energy
        p.rdata(4) = 0.;

        particle_tile.push_back(p);
    }
    ///
    /// Invalidate particles in cells that are already occupied
    ///
    void ComputeOverlap(int lev);

    ///
    /// Invalidate particles that have been merged with other particles
    ///
    void Merge(int lev);

    ///
    /// Compute the momentum that has been removed from the gas in order to define the particle velocity
    ///
    void ComputeParticleVelocity(int lev, amrex::MultiFab& state_old, amrex::MultiFab& state_new, int add_energy);

    ///
    /// Accrete mass from the grid onto the existing AGN particles
    ///
    void AccreteMass(int lev, amrex::MultiFab& state, amrex::Real eps_rad, amrex::Real dt);

    ///
    /// This fills the ghost buffers for each tile with the proper data
    ///
    void fillGhosts(int lev);

    ///
    /// Each tile clears its ghosts, freeing the memory
    ///
    void clearGhosts(int lev);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev);

private:

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packGhostParticle(int lev,
                           const amrex::IntVect& neighbor_cell,
                           const amrex::BaseFab<int>& mask,
                           const ParticleType& p,
                           GhostCommMap& ghosts_to_comm);

    ///
    /// Perform the MPI communication neccesary to fill ghost buffers
    ///
    void fillGhostsMPI(GhostCommMap& ghosts_to_comm);

    const size_t pdata_size = 8*sizeof(RealType); // we communicate the position for position(3), mass, velocity(3) and energy
    amrex::FabArray<amrex::BaseFab<int> > mask;
    std::map<PairIndex, amrex::Array<char> > ghosts;
};

#endif /* _AGNParticleContainer_H_ */
