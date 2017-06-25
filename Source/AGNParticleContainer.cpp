#include "AGNParticleContainer.H"
#include "AMReX_RealVect.H"
#include "agn_F.H"

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

    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
					 this->m_gdb->ParticleDistributionMap(lev),
					 acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary(periodic);
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 1;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_agn_particles(&Np, particles.data(),
                                (*ac_ptr)[pti].dataPtr(),
                                ac_box.loVect(), ac_box.hiVect(),
                                plo,dx,dt,a_old,a_half,&do_move);
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

    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

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

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 0;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_agn_particles(&Np, particles.data(),
                                (*ac_ptr)[pti].dataPtr(),
                                ac_box.loVect(), ac_box.hiVect(),
                                plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

void AGNParticleContainer::ComputeOverlap(int lev)
{
    BL_PROFILE("AGNParticleContainer::ComputeOverlap()");
    Array<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = ghosts[index].size() / pdata_size;

        nyx_compute_overlap(&Np, particles.data(), 
                            &Ng, ghosts[index].dataPtr(), dx);

    }
}

void AGNParticleContainer::Merge(int lev)
{
    BL_PROFILE("AGNParticleContainer::Merge()");
    Array<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = ghosts[index].size() / pdata_size;

        agn_merge_particles(&Np, particles.data(), 
                            &Ng, ghosts[index].dataPtr(), dx);
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

void AGNParticleContainer::defineMask() 
{
    BL_PROFILE("AGNParticleContainer::defineMask()");
    const int lev = 0;
    const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
    const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);
    const Geometry& gm = m_gdb->Geom(lev);

    mask.define(ba, dm, 2, ng);
    mask.setVal(-1, ng);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        mask.setVal(grid_id, box, 0, 1);
        mask.setVal(tile_id, box, 1, 1);
    }

    mask.FillBoundary(gm.periodicity());
    mask_defined = true;
}

void AGNParticleContainer::fillNeighbors(int lev) 
{
    BL_PROFILE("AGNParticleContainer::fillNeighbors()");
    BL_ASSERT(lev == 0);
    if (!mask_defined) defineMask();

    NeighborCommMap ghosts_to_comm;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& tile_box = pti.tilebox();
        const IntVect& lo = tile_box.smallEnd();
        const IntVect& hi = tile_box.bigEnd();
        
        Box shrink_box = pti.tilebox();
        shrink_box.grow(-ng);
        
        auto& particles = pti.GetArrayOfStructs();
        for (unsigned i = 0; i < pti.numParticles(); ++i) {
            const ParticleType& p = particles[i];
            const IntVect& iv = Index(p, lev);
            
            // if the particle is more than one cell away from 
            // the tile boundary, it's not anybody's neighbor
            if (shrink_box.contains(iv)) continue;
            
            // shift stores whether we are near the tile boundary in each direction.
            // -1 means lo, 1 means hi, 0 means not near the boundary
            IntVect shift = IntVect::TheZeroVector();
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                if (iv[idim] == lo[idim])
                    shift[idim] = -ng;
                else if (iv[idim] == hi[idim])
                    shift[idim] = ng;
            }
            
            // Based on the value of shift, we add the particle to a map to be sent
            // to the neighbors. A particle can be sent to up to 3 neighbors in 2D
            // and up to 7 in 3D, depending on whether is near the tile box corners,
            // edges, or just the faces. First, add the particle for the "face" neighbors
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                if (shift[idim] == 0) continue;
                IntVect neighbor_cell = iv;
                neighbor_cell.shift(idim, shift[idim]);
                BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                packNeighborParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
            }
            
            // Now add the particle to the "edge" neighbors
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                for (int jdim = 0; jdim < idim; ++jdim) {
                    if (shift[idim] != 0 and shift[jdim] != 0) {
                        IntVect neighbor_cell = iv;
                        neighbor_cell.shift(idim, shift[idim]);
                        neighbor_cell.shift(jdim, shift[jdim]);
                        BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                        packNeighborParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
                    }
                }
            }
            
            // Finally, add the particle for the "vertex" neighbors (only relevant in 3D)
            if (shift[0] != 0 and shift[1] != 0 and shift[2] != 0) {
                IntVect neighbor_cell = iv;
                neighbor_cell.shift(shift);
                BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                packNeighborParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
            }
        }
    }
    
    fillNeighborsMPI(ghosts_to_comm);
}

void AGNParticleContainer::clearNeighbors(int lev) 
{
    BL_PROFILE("AGNParticleContainer::clearNeighbors()");
    ghosts.clear();
}

void AGNParticleContainer::applyPeriodicShift(int lev, ParticleType& p,
                                              const IntVect& neighbor_cell) 
{
    BL_PROFILE("AGNParticleContainer::applyPeriodicShift()");
    const Periodicity& periodicity = Geom(lev).periodicity();
    if (not periodicity.isAnyPeriodic()) return;

    const Box& domain = Geom(lev).Domain();
    const IntVect& lo = domain.smallEnd();
    const IntVect& hi = domain.bigEnd();
    const RealBox& prob_domain = Geom(lev).ProbDomain();

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        if (not periodicity.isPeriodic(dir)) continue;
        if (neighbor_cell[dir] < lo[dir]) {
            p.pos(dir) += prob_domain.length(dir);
        }
        else if (neighbor_cell[dir] > hi[dir]) {
            p.pos(dir) -= prob_domain.length(dir);
        }
    }
}

void AGNParticleContainer::packNeighborParticle(int lev,
                                             const IntVect& neighbor_cell,
                                             const BaseFab<int>& mask,
                                             const ParticleType& p,
                                             NeighborCommMap& ghosts_to_comm) 
{
    BL_PROFILE("AGNParticleContainer::packNeighborParticle()");
    const int neighbor_grid = mask(neighbor_cell, 0);
    if (neighbor_grid >= 0) {
        const int who = ParticleDistributionMap(lev)[neighbor_grid];
        const int MyProc = ParallelDescriptor::MyProc();
        const int neighbor_tile = mask(neighbor_cell, 1);
        PairIndex dst_index(neighbor_grid, neighbor_tile);
        ParticleType particle = p;
        applyPeriodicShift(lev, particle, neighbor_cell);
        if (who == MyProc) {
            size_t old_size = ghosts[dst_index].size();
            size_t new_size = ghosts[dst_index].size() + pdata_size;
            ghosts[dst_index].resize(new_size);
            std::memcpy(&ghosts[dst_index][old_size], &particle, pdata_size);
        } else {
            NeighborCommTag tag(who, neighbor_grid, neighbor_tile);
            Array<char>& buffer = ghosts_to_comm[tag];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + pdata_size;
            buffer.resize(new_size);
            std::memcpy(&buffer[old_size], &particle, pdata_size);
        }
    }
}

void AGNParticleContainer::fillNeighborsMPI(NeighborCommMap& ghosts_to_comm) 
{
    BL_PROFILE("AGNParticleContainer::fillNeighborsMPI()");
#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    
    // count the number of tiles to be sent to each proc
    std::map<int, int> tile_counts;
    for (const auto& kv: ghosts_to_comm) {
        tile_counts[kv.first.proc_id] += 1;
    }
    
    // flatten all the data for each proc into a single buffer
    // once this is done, each dst proc will have an Array<char>
    // the buffer will be packed like:
    // ntiles, gid1, tid1, size1, data1....  gid2, tid2, size2, data2... etc. 
    std::map<int, Array<char> > send_data;
    for (const auto& kv: ghosts_to_comm) {
        Array<char>& buffer = send_data[kv.first.proc_id];
        buffer.resize(sizeof(int));
        std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
    }
    
    for (auto& kv : ghosts_to_comm) {
        int data_size = kv.second.size();
        Array<char>& buffer = send_data[kv.first.proc_id];
        size_t old_size = buffer.size();
        size_t new_size = buffer.size() + 2*sizeof(int) + sizeof(int) + data_size;
        buffer.resize(new_size);
        char* dst = &buffer[old_size];
        std::memcpy(dst, &(kv.first.grid_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.tile_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &data_size,          sizeof(int)); dst += sizeof(int);
        if (data_size == 0) continue;
        std::memcpy(dst, &kv.second[0], data_size);
        Array<char>().swap(kv.second);
    }
    
    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Array<long> snds(NProcs, 0), rcvs(NProcs, 0);
    long num_snds = 0;
    for (const auto& kv : send_data) {
        num_snds      += kv.second.size();
        snds[kv.first] = kv.second.size();
    }
    ParallelDescriptor::ReduceLongMax(num_snds);
    if (num_snds == 0) return;
    
    // communicate that information
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(rcvs[MyProc] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());
    
    Array<int> RcvProc;
    Array<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Array<MPI_Status>  stats(nrcvs);
    Array<MPI_Request> rreqs(nrcvs);
    
    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.
    Array<char> recvdata(TotRcvBytes);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }
    
    // unpack the received data and put them into the proper ghost buffers
    if (nrcvs > 0) {
        BL_MPI_REQUIRE( MPI_Waitall(nrcvs, rreqs.data(), stats.data()) );
        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, gid, tid, size;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);
                
                if (size == 0) continue;
                
                PairIndex dst_index(gid, tid);
                size_t old_size = ghosts[dst_index].size();
                size_t new_size = ghosts[dst_index].size() + size;
                ghosts[dst_index].resize(new_size);
                std::memcpy(&ghosts[dst_index][old_size], buffer, size); buffer += size;
            }
        }
    }
#endif
}

void AGNParticleContainer::writeAllAtLevel(int lev)
{
  BL_PROFILE("AGNParticleContainer::writeAllAtLevel()");
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
      pout() << "There are " << Np  << " AGN particles at level " << lev
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

          pout() << "[" << i << "]: id " << id << " cpu " << cpu
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
