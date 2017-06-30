
#ifndef _AGNParticleContainer_H_
#define _AGNParticleContainer_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

#include "NyxParticleContainer.H"

namespace {

struct NeighborCommTag {

    NeighborCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const NeighborCommTag& l, const NeighborCommTag& r) {
    return (l.proc_id < r.proc_id ||
           (l.proc_id == r.proc_id && l.grid_id < r.grid_id) ||
           (l.proc_id == r.proc_id && l.grid_id == r.grid_id && l.tile_id < r.tile_id ));
}

}

class AGNParticleContainer
    : public NyxParticleContainer<3+BL_SPACEDIM>
{
public:

    using MyParIter = amrex::ParIter<3+BL_SPACEDIM>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, amrex::Array<char> >;

    AGNParticleContainer (amrex::Amr* amr, int nghost)
      : NyxParticleContainer<3+BL_SPACEDIM>(amr),
      mask_defined(false),
      ng(nghost)
    {
      real_comp_names.clear();
      real_comp_names.push_back("mass");
      real_comp_names.push_back("xvel");
      real_comp_names.push_back("yvel");
      real_comp_names.push_back("zvel");
      real_comp_names.push_back("energy");
      real_comp_names.push_back("mdot");
    }

    virtual ~AGNParticleContainer () {}

    const int NumberOfParticles(MyParIter& pti) { return pti.GetArrayOfStructs().size(); }

    virtual void moveKickDrift (amrex::MultiFab& acceleration,
                                int level,
                                amrex::Real timestep,
                                amrex::Real a_old = 1.0,
                                amrex::Real a_half = 1.0,
                                int where_width = 0);

    virtual void moveKick      (amrex::MultiFab& acceleration,
                                int level,
                                amrex::Real timestep,
                                amrex::Real a_new = 1.0,
                                amrex::Real a_half = 1.0);

  void AddOneParticle (int lev,
                       int grid,
                       int tile,
                       amrex::Real mass, 
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z)
    {
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];
        AddOneParticle(particle_tile, mass, x, y, z);
    }

    void AddOneParticle (ParticleTileType& particle_tile,
                         amrex::Real mass, 
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z)
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

        // Zero initial mdot
        p.rdata(5) = 0.;

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
    void ComputeParticleVelocity(int lev,
                                 amrex::MultiFab& state_old,
                                 amrex::MultiFab& state_new,
                                 int add_energy);

    ///
    /// Accrete mass from the grid onto the existing AGN particles
    ///
    void AccreteMass(int lev,
                     amrex::MultiFab& state,
                     amrex::MultiFab& density_lost,
                     amrex::Real dt);

    ///
    ///  Release energy if it exceeds thermal feedback threshold.
    ///
  void ReleaseEnergy(int lev,
                     amrex::MultiFab& state,
                     amrex::MultiFab& D_new,
                     amrex::Real a);

    ///
    /// This fills the ghost buffers for each tile with the proper data
    ///
    void fillNeighbors(int lev);

    ///
    /// Each tile clears its ghosts, freeing the memory
    ///
    void clearNeighbors(int lev);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev);

private:

    ///
    /// Apply periodic shift to particle position, so naive distance calculation
    /// between neighbors and regular particles will be correct.
    ///
    void applyPeriodicShift(int lev, ParticleType& p, const amrex::IntVect& neighbor_cell);

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packNeighborParticle(int lev,
                           const amrex::IntVect& neighbor_cell,
                           const amrex::BaseFab<int>& mask,
                           const ParticleType& p,
                           NeighborCommMap& ghosts_to_comm);

    ///
    /// Perform the MPI communication neccesary to fill ghost buffers
    ///
    void fillNeighborsMPI(NeighborCommMap& ghosts_to_comm);

    void defineMask();

    const size_t pdata_size = sizeof(ParticleType);
    int ng;
    amrex::FabArray<amrex::BaseFab<int> > mask;
    bool mask_defined;
    std::map<PairIndex, amrex::Array<char> > ghosts;
};

#endif /* _AGNParticleContainer_H_ */
