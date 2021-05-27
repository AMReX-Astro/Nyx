#include <stdint.h>

#include <DarkMatterParticleContainer.H>

using namespace amrex;

/// These are helper functions used when initializing from a morton-ordered
/// binary particle file.
namespace {

  inline uint64_t split(unsigned int a) {
    uint64_t x = a & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8)  & 0x100f00f00f00f00f;
    x = (x | x << 4)  & 0x10c30c30c30c30c3;
    x = (x | x << 2)  & 0x1249249249249249;
    return x;
  }
  
  inline uint64_t get_morton_index(unsigned int x,
                                   unsigned int y,
                                   unsigned int z) {
    uint64_t morton_index = 0;
    morton_index |= split(x) | ( split(y) << 1) | (split(z) << 2);
    return morton_index;
  }  

  struct BoxMortonKey {
    uint64_t morton_id;
    int box_id;
  };

  struct by_morton_id { 
    bool operator()(const BoxMortonKey &a, const BoxMortonKey &b) { 
      return a.morton_id < b.morton_id;
    }
  };

  std::string get_file_name(const std::string& base, int file_num) {
    std::stringstream ss;
    ss << base << file_num;
    return ss.str();
  }

  struct ParticleMortonFileHeader {
    long NP;
    int  DM;
    int  NX;
    int  SZ;
    int  NF;
  };
  
  void ReadHeader(const std::string& dir,
                  const std::string& file,
                  ParticleMortonFileHeader& hdr) {
    std::string header_filename = dir;
    header_filename += "/";
    header_filename += file;
    
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(header_filename, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);

    HdrFile >> hdr.NP;
    HdrFile >> hdr.DM;
    HdrFile >> hdr.NX;
    HdrFile >> hdr.SZ;
    HdrFile >> hdr.NF;    
  }

}

void
DarkMatterParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
                                            int                    lev,
                                            amrex::Real            dt,
                                            amrex::Real            a_old,
                                            amrex::Real            a_half,
                                            int                    where_width)
{
    BL_PROFILE("DarkMatterParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;
    const auto dxi              = Geom(lev).InvCellSizeArray();

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

    const GpuArray<Real,AMREX_SPACEDIM> plo = Geom(lev).ProbLoArray();

    int do_move = 1;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        ParticleType* pstruct = particles().data();
        const long np = pti.numParticles();
        int grid    = pti.index();

        const FArrayBox& accel_fab= ((*ac_ptr)[grid]);
        Array4<amrex::Real const> accel= accel_fab.array();

        int nc=AMREX_SPACEDIM;
        amrex::ParallelFor(np,
                           [=] AMREX_GPU_HOST_DEVICE ( long i)
                           {
                             update_dm_particle_single(pstruct[i],nc,
                                                       accel,
                                                       plo,dxi,dt,a_old, a_half,do_move);
                           });
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
#pragma omp parallel for private(pld) if (Gpu::notInLaunchRegion())
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
                        int grid = kv.first.first;
                        
                        
                        std::cout << "Oops -- removing particle " << p << " " << this->Index(p, lev) << " " << lev << " " << (this->m_gdb->ParticleBoxArray(lev))[grid] << " " << where_width << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
DarkMatterParticleContainer::moveKick (MultiFab&       acceleration,
                                       int             lev,
                                       Real            dt,
                                       Real            a_new,
                                       Real            a_half) 
{
    BL_PROFILE("DarkMatterParticleContainer::moveKick()");

    const auto dxi              = Geom(lev).InvCellSizeArray();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
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

    const GpuArray<Real,AMREX_SPACEDIM> plo = Geom(lev).ProbLoArray();

    int do_move = 0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        ParticleType* pstruct = particles().data();
        const long np = pti.numParticles();
        int grid    = pti.index();

        const FArrayBox& accel_fab= ((*ac_ptr)[grid]);
        Array4<amrex::Real const> accel= accel_fab.array();

        int nc=AMREX_SPACEDIM;
        amrex::ParallelFor(np,
                           [=] AMREX_GPU_HOST_DEVICE ( long i)
                           {
                             update_dm_particle_single(pstruct[i],nc,
                                                       accel,
                                                       plo,dxi,dt,a_half,a_new,do_move);
                           });
    }

    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void update_dm_particle_single (amrex::ParticleContainer<4, 0>::SuperParticleType&  p,
                                const int nc,
                                amrex::Array4<amrex::Real const> const& acc,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                const amrex::Real& dt, const amrex::Real& a_prev,
                                const amrex::Real& a_cur, const int& do_move)
{
    amrex::Real half_dt       = 0.5 * dt;
    amrex::Real a_cur_inv    = 1.0 / a_cur;
    amrex::Real dt_a_cur_inv = dt * a_cur_inv;

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;
    
    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));
    
    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;
    
    amrex::Real sx[] = {amrex::Real(1.)-xint, xint};
    amrex::Real sy[] = {amrex::Real(1.)-yint, yint};
    amrex::Real sz[] = {amrex::Real(1.)-zint, zint};

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
           for (int comp=0; comp < nc; ++comp) {
             p.pos(comp) = p.pos(comp) + dt_a_cur_inv * p.rdata(comp+1);
           }
         }

}

void
DarkMatterParticleContainer::InitCosmo1ppcMultiLevel(
                        MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                        const Real particleMass, int disp_idx, int vel_idx, 
                        BoxArray &baWhereNot, int lev, int nlevs)
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Vector<int> calls;

    calls.resize(nlevs);

    calls[lev]++;

    if (calls[lev] > 1) return;

    Vector<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(nlevs);

    ParticleType p;
    Real         disp[AMREX_SPACEDIM];
    Real         vel[AMREX_SPACEDIM];
    
    Real        mean_disp[AMREX_SPACEDIM]={D_DECL(0,0,0)};


    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();
        ParticleLocData pld;
        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
                    IntVect indices(D_DECL(ix, jx, kx));
                    totalcount++;
                    if (baWhereNot.contains(indices)) 
                    {
                       continue;
                    }

                    for (int n = 0; n < AMREX_SPACEDIM; n++)
                    {
                        disp[n] = myFab(indices,disp_idx+n);
                        //
                        // Start with homogeneous distribution (for 1 p per cell in the center of the cell),
                        //
                        p.pos(n) = geom.ProbLo(n) + 
                            (indices[n]+Real(0.5))*dx[n];
                        if(disp[n]*disp_fac[n]>dx[n]/2.0)
                          outcount[n]++;
                        if(disp[n]*disp_fac[n]<-dx[n]/2.0)
                          outcountminus[n]++;
                        mean_disp[n]+=fabs(disp[n]);
                        //
                        // then add the displacement (input values weighted by domain length).
                        //
                        p.pos(n) += disp[n] * disp_fac[n];

                        //
                        // Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
                        p.rdata(n+1) = vel[n] * vel_fac[n];
                    }
                    //
                    // Set the mass of the particle from the input value.
                    //
                    p.rdata(0)  = particleMass;
                    p.id()      = ParticleType::NextID();
                    p.cpu()     = MyProc;
        
                    if (!this->Where(p, pld))
                    {
                        this->PeriodicShift(p);

                        if (!this->Where(p, pld))
                            amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                    }

                    BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
                    //handle particles that ran out of this level into a finer one. 
                    if (baWhereNot.contains(pld.m_cell))
                    {
                      outside_counter++;
                      ParticleType newp[8];
                      ParticleLocData new_pld;
                      for (int i=0;i<8;i++)
                      {
                          newp[i].rdata(0)   = particleMass/8.0;
                          newp[i].id()       = ParticleType::NextID();
                          newp[i].cpu()      = MyProc;
                          for (int dim=0;dim<AMREX_SPACEDIM;dim++)
                          {
                              newp[i].pos(dim)=p.pos(dim)+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
                              newp[i].rdata(dim+1)=p.rdata(dim+1);
                          }
                          
                          if (!this->Where(newp[i], new_pld))
                          {
                              this->PeriodicShift(newp[i]);
                              
                              if (!this->Where(newp[i], new_pld))
                                  amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                          }
                          particles[new_pld.m_lev][std::make_pair(new_pld.m_grid, 
                                                                  new_pld.m_tile)].push_back(newp[i]);
                      }
                      
                    }
                    
                    //
                    // Add it to the appropriate PBox at the appropriate level.
                    //
                    else
                        particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

/*
  Particle deposition
*/

void
DarkMatterParticleContainer::AssignDensityAndVels (Vector<std::unique_ptr<MultiFab> >& mf, int lev_min) const
{
     AssignDensity(mf, lev_min, AMREX_SPACEDIM+1);
}

void 
DarkMatterParticleContainer::InitFromBinaryMortonFile(const std::string& particle_directory,
                                                      int /*nextra*/, int skip_factor) {
  BL_PROFILE("DarkMatterParticleContainer::InitFromBinaryMortonFile");
  
  ParticleMortonFileHeader hdr;
  ReadHeader(particle_directory, "Header", hdr);    
  
  uint64_t num_parts = hdr.NP;
  int DM             = hdr.DM;
  int NX             = hdr.NX;
  int float_size     = hdr.SZ;
  int num_files      = hdr.NF;
  size_t psize       = (DM + NX) * float_size;
  
  std::string particle_file_base = particle_directory + "/particles.";
  std::vector<std::string> file_names;
  for (int i = 0; i < num_files; ++i)
    file_names.push_back(get_file_name(particle_file_base, i));
  
  const int lev = 0;
  const BoxArray& ba = ParticleBoxArray(lev);
  int num_boxes = ba.size();
  uint64_t num_parts_per_box  = num_parts / num_boxes;
  uint64_t num_parts_per_file = num_parts / num_files;
  uint64_t num_bytes_per_file = num_parts_per_file * psize;
  
  std::vector<BoxMortonKey> box_morton_keys(num_boxes);
  for (int i = 0; i < num_boxes; ++i) {
    const Box& box = ba[i];
    unsigned int x = box.smallEnd(0);
    unsigned int y = box.smallEnd(1);
    unsigned int z = box.smallEnd(2);
    box_morton_keys[i].morton_id = get_morton_index(x, y, z);
    box_morton_keys[i].box_id = i;
  }
  
  std::sort(box_morton_keys.begin(), box_morton_keys.end(), by_morton_id());
  
  std::vector<int> file_indices(num_boxes);
  for (int i = 0; i < num_boxes; ++i)
    file_indices[box_morton_keys[i].box_id] = i;
  
  ParticleType p;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {  // no tiling
    const int grid = mfi.index();
    const int tile = mfi.LocalTileIndex();      
    auto& particles = GetParticles(lev);
    
    uint64_t start    = file_indices[grid]*num_parts_per_box;
    uint64_t stop     = start + num_parts_per_box;

    int file_num      = start / num_parts_per_file;
    uint64_t seek_pos = (start * psize ) % num_bytes_per_file;
    std::string file_name = file_names[file_num];
    
    std::ifstream ifs;
    ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
    if (!ifs ) {
      amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
      amrex::Abort();
    } 

    ifs.seekg(seek_pos, std::ios::beg);
    
    for (uint64_t i = start; i < stop; ++i) {
      int next_file = i / num_parts_per_file;
      if (next_file != file_num) {
        file_num = next_file;
        file_name = file_names[file_num];
        ifs.close();
        ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
        if (!ifs ) {
          amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
          amrex::Abort();
        }
      }

      Vector<float> fpos(DM);
      Vector<float> fextra(NX);
      ifs.read((char*)&fpos[0],   DM*sizeof(float));
      ifs.read((char*)&fextra[0], NX*sizeof(float));
      
      if ( (i - start) % skip_factor == 0 ) {
        AMREX_D_TERM(p.pos(0) = fpos[0];,
                     p.pos(1) = fpos[1];,
                     p.pos(2) = fpos[2];);
        
        for (int comp = 0; comp < NX; comp++)
          p.rdata(AMREX_SPACEDIM+comp) = fextra[comp];
        
        p.rdata(AMREX_SPACEDIM) *= skip_factor;
        
        p.id()  = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        particles[std::make_pair(grid, tile)].push_back(p);
      }
    }    
  }
  
  Redistribute();
}

