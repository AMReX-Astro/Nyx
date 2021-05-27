#include <AGNParticleContainer.H>
#include <AMReX_RealVect.H>
#include <agn_F.H>

using namespace amrex;

using std::cout;
using std::endl;

void
AGNParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
                                     int                    lev,
                                     amrex::Real            dt,
                                     amrex::Real            a_old,
                                     amrex::Real            a_half,
                                     int                    where_width)
{
    BL_PROFILE("AGNParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const auto dxi = Geom(lev).InvCellSizeArray();
    const GpuArray<Real,AMREX_SPACEDIM> plo = Geom(lev).ProbLoArray();
    const Periodicity& periodic = Geom(lev).periodicity();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        const IntVect& ng = acceleration.nGrowVect();
        ac_ptr = new amrex::MultiFab(this->ParticleBoxArray(lev),
                                     this->ParticleDistributionMap(lev),
                                     acceleration.nComp(),acceleration.nGrow());
        ac_ptr->setVal(0.);
        if(acceleration.boxArray() == ac_ptr->boxArray())//this->finestLevel() == 0)
        {
            ac_ptr->Redistribute(acceleration,0,0,acceleration.nComp(),ng);
            ac_ptr->FillBoundary();
        }
        else
        {
            ac_ptr->ParallelCopy(acceleration,0,0,acceleration.nComp(),ng,ng);
            ac_ptr->FillBoundary();
        }
    }

    int do_move = 1;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        ParticleType* pstruct = particles().data();
        int Np = particles.size();

        const FArrayBox& accel_fab= ((*ac_ptr)[pti.index()]);
        Array4<amrex::Real const> accel= accel_fab.array();

        if (Np > 0)
        {
            amrex::ParallelFor(Np,
                               [=] AMREX_GPU_HOST_DEVICE ( long i)
                               {
                                   update_agn_particle_single(pstruct[i],
                                                              accel,
                                                              plo,dxi,dt,a_old, a_half,do_move);
                               });
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;

    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        amrex::ParticleLocData pld;
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetArrayOfStructs();
            const int   n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell.
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->Where(p, pld, lev, lev, where_width))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.id() == amrex::GhostParticleID)
                    {
                        p.id() = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.id() << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
AGNParticleContainer::moveKick (MultiFab&       acceleration,
                                int             lev,
                                Real            dt,
                                Real            a_new,
                                Real            a_half)
{
    BL_PROFILE("AGNParticleContainer::moveKick()");

    const auto dxi = Geom(lev).InvCellSizeArray();
    const Periodicity& periodic = Geom(lev).periodicity();
    const GpuArray<Real,AMREX_SPACEDIM> plo = Geom(lev).ProbLoArray();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new MultiFab(ParticleBoxArray(lev),
                                  ParticleDistributionMap(lev),
                                  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary(periodic);
    }

    int do_move = 0;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        ParticleType* pstruct = particles().data();
        int Np = particles.size();

        const FArrayBox& accel_fab= ((*ac_ptr)[pti.index()]);
        Array4<amrex::Real const> accel= accel_fab.array();

        if (Np > 0)
        {
            amrex::ParallelFor(Np,
                               [=] AMREX_GPU_HOST_DEVICE ( long i)
                               {
                                   update_agn_particle_single(pstruct[i],
                                                              accel,
                                                              plo,dxi,dt,a_half,a_new,do_move);
                               });
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;
}

void AGNParticleContainer::ComputeOverlap(int lev)
{
    BL_PROFILE("AGNParticleContainer::ComputeOverlap()");
    Vector<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = neighbors[lev][index].size() / pdata_size;

        nyx_compute_overlap(&Np, particles.data(), 
                            &Ng, neighbors[lev][index].dataPtr(), dx);

    }
}

void AGNParticleContainer::Merge(int lev)
{
    BL_PROFILE("AGNParticleContainer::Merge()");
    Vector<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = neighbors[lev][index].size() / pdata_size;

        agn_merge_particles(&Np, particles.data(), 
                            &Ng, neighbors[lev][index].dataPtr(), dx);
    }
}

void AGNParticleContainer::ComputeParticleVelocity(int lev,
                                                   amrex::MultiFab& state_old, 
                                                   amrex::MultiFab& state_new,
                                                   int add_energy)
{
    BL_PROFILE("AGNParticleContainer::ComputeParticleVelocity()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state_old.FillBoundary(periodic);
    state_new.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& soldbox = state_old[pti].box();
        const Box& snewbox = state_new[pti].box();

        agn_particle_velocity(&Np, particles.data(), 
                              state_old[pti].dataPtr(), 
                              soldbox.loVect(), soldbox.hiVect(),
                              state_new[pti].dataPtr(), 
                              snewbox.loVect(), snewbox.hiVect(),
                              dx, &add_energy);
    }
}

void AGNParticleContainer::AccreteMass(int lev,
                                       amrex::MultiFab& state,
                                       amrex::MultiFab& density_lost,
                                       amrex::Real dt)
{
    BL_PROFILE("AGNParticleContainer::AccreteMass()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& sbox = state[pti].box();

        agn_accrete_mass(&Np, particles.data(),
                         state[pti].dataPtr(),
                         density_lost[pti].dataPtr(),
                         sbox.loVect(), sbox.hiVect(),
                         &dt, dx);
    }
}

void AGNParticleContainer::ReleaseEnergy(int lev, amrex::MultiFab& state, amrex::MultiFab& D_new, amrex::Real a)
{
    BL_PROFILE("AGNParticleContainer::ReleaseEnergy()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state.FillBoundary(periodic);
    D_new.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& sbox = state[pti].box();
        const Box& Dbox = D_new[pti].box();
        agn_release_energy(&Np, particles.data(), 
                           state[pti].dataPtr(), 
                           sbox.loVect(), sbox.hiVect(),
                           D_new[pti].dataPtr(),
                           Dbox.loVect(), Dbox.hiVect(),
                           &a, dx); 
    }
}

void AGNParticleContainer::writeAllAtLevel(int lev)
{
  BL_PROFILE("AGNParticleContainer::writeAllAtLevel()");
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
      Print() << "There are " << Np  << " AGN particles at level " << lev
              << " in grid " << pti.index() << std::endl;
      for (unsigned i = 0; i < Np; ++i)
        {
          const ParticleType& p = particles[i];
          const IntVect& iv = Index(p, lev);

          int id = p.id();
          int cpu = p.cpu();
          RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
          Real mass = p.rdata(0);
          RealVect uvw(p.rdata(1), p.rdata(2), p.rdata(3));
          Real energy = p.rdata(4);
          Real mdot = p.rdata(5);

          Print() << "[" << i << "]: id " << id << " cpu " << cpu
                  << " mass " << mass
                  << " index " << iv
                  << " position " << xyz
                  << " velocity " << uvw
                  << " energy " << energy
                  << " mdot " << mdot
                  << endl;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void AGNParticleContainer::
update_agn_particle_single (ParticleType&  p,
                            amrex::Array4<amrex::Real const> const& acc,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                            const amrex::Real& dt, const amrex::Real& a_prev,
                            const amrex::Real& a_cur, const int& do_move)
{
    if (p.id() <= 0) return;

    amrex::Real half_dt       = 0.5 * dt;
    amrex::Real a_cur_inv    = 1.0 / a_cur;
    amrex::Real dt_a_cur_inv = dt * a_cur_inv;

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

    int i = amrex::Math::floor(lx);
    int j = amrex::Math::floor(ly);
    int k = amrex::Math::floor(lz);

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
        amrex::Real val = 0.0;
        for (int kk = 0; kk<=1; ++kk)
        {
            for (int jj = 0; jj <= 1; ++jj)
            {
                for (int ii = 0; ii <= 1; ++ii)
                {
                    val += sx[amrex::Math::abs(ii-1)]*
                           sy[amrex::Math::abs(jj-1)]*
                           sz[amrex::Math::abs(kk-1)]*acc(i-ii,j-jj,k-kk,d);
                }
            }
        }

        p.rdata(d+1)=a_prev*p.rdata(d+1)+half_dt * val;
        p.rdata(d+1)*=a_cur_inv;
    }

    if (do_move == 1)
    {
        for (int comp=0; comp < AMREX_SPACEDIM; ++comp) {
            p.pos(comp) = p.pos(comp) + dt_a_cur_inv * p.rdata(comp+1);
        }
    }
}
