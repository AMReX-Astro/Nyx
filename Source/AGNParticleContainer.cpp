#include "AGNParticleContainer.H"
#include "AMReX_RealVect.H"
#include "agn_F.H"

using namespace amrex;

using std::cout;
using std::endl;

void AGNParticleContainer::ComputeOverlap(int lev)
{
    Array<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        size_t Np = particles.size();

        my_id.resize(Np);
        for (int i = 0; i < Np; ++i ) {
          my_id[i] = particles[i].id();
        }

        int nstride = particles.dataShape().first;
        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = ghosts[index].size() / pdata_size;

//      const Box& dbox = density_to_subtract[pti].box();

        nyx_compute_overlap(particles.data(), nstride, Np, my_id.dataPtr(),
                           (RealType*) ghosts[index].dataPtr(), Ng, dx);
//                         density_to_subtract[pti].dataPtr(), 
//                         dbox.loVect(), dbox.hiVect());

        for (int i = 0; i < Np; ++i ) {
          particles[i].id() = my_id[i];
        }
    }
}

void AGNParticleContainer::fillGhosts(int lev) {
    GhostCommMap ghosts_to_comm;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& tile_box = pti.tilebox();
        const IntVect& lo = tile_box.smallEnd();
        const IntVect& hi = tile_box.bigEnd();
        
        Box shrink_box = pti.tilebox();
        shrink_box.grow(-1);
        
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
                    shift[idim] = -1;
                else if (iv[idim] == hi[idim])
                    shift[idim] = 1;
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
                packGhostParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
            }
            
            // Now add the particle to the "edge" neighbors
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                for (int jdim = 0; jdim < idim; ++jdim) {
                    if (shift[idim] != 0 and shift[jdim] != 0) {
                        IntVect neighbor_cell = iv;
                        neighbor_cell.shift(idim, shift[idim]);
                        neighbor_cell.shift(jdim, shift[jdim]);
                        BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                        packGhostParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
                    }
                }
            }
            
#if (BL_SPACEDIM == 3)
            // Finally, add the particle for the "vertex" neighbors (only relevant in 3D)
            if (shift[0] != 0 and shift[1] != 0 and shift[2] != 0) {
                IntVect neighbor_cell = iv;
                neighbor_cell.shift(shift);
                BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                packGhostParticle(lev, neighbor_cell, mask[pti], p, ghosts_to_comm);
            }
#endif
        }
    }
    
    fillGhostsMPI(ghosts_to_comm);
}

void AGNParticleContainer::clearGhosts(int lev) {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        auto& ghost_particles = ghosts[std::make_pair(grid_id, tile_id)];
        Array<char>().swap(ghost_particles);
    }
}

void AGNParticleContainer::packGhostParticle(int lev,
                                             const IntVect& neighbor_cell,
                                             const BaseFab<int>& mask,
                                             const ParticleType& p,
                                             GhostCommMap& ghosts_to_comm) {
    const int neighbor_grid = mask(neighbor_cell, 0);
    if (neighbor_grid >= 0) {
        const int who = ParticleDistributionMap(lev)[neighbor_grid];
        const int MyProc = ParallelDescriptor::MyProc();
        const int neighbor_tile = mask(neighbor_cell, 1);
        PairIndex dst_index(neighbor_grid, neighbor_tile);
        if (who == MyProc) {
            size_t old_size = ghosts[dst_index].size();
            size_t new_size = ghosts[dst_index].size() + pdata_size;
            ghosts[dst_index].resize(new_size);
            std::memcpy(&ghosts[dst_index][old_size], &p, pdata_size);
        } else {
            GhostCommTag tag(who, neighbor_grid, neighbor_tile);
            Array<char>& buffer = ghosts_to_comm[tag];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + pdata_size;
            buffer.resize(new_size);
            std::memcpy(&buffer[old_size], &p, pdata_size);
        }
    }
}

void AGNParticleContainer::fillGhostsMPI(GhostCommMap& ghosts_to_comm) {

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
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
      cout << "AGN particles: " << Np << " << at level " << lev << endl;
      for (unsigned i = 0; i < Np; ++i)
        {
          const ParticleType& p = particles[i];
          const IntVect& iv = Index(p, lev);

          RealVect xyz(p.pos(0), p.pos(1), p.pos(2));

          cout << "[" << i << "]: id " << p.id()
               << " mass " << p.rdata(0)
               << " index " << iv
               << " position " << xyz << endl;
        }
    }
}
