#include <stdint.h>

#include <DarkMatterParticleContainer.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_ParticleHDF5.H>
#endif

using namespace amrex;
#include <constants_cosmo.H>
#include <LightConeParticle.H>

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

struct ShellFilter
{
    GpuArray<Real, AMREX_SPACEDIM> m_plo, m_phi, m_center;
    Real m_radius_inner, m_radius_outer, m_z, m_t, m_dt;
    Box m_domain;

    ShellFilter (const GpuArray<Real, AMREX_SPACEDIM>& plo,
                 const GpuArray<Real, AMREX_SPACEDIM>& phi,
                 const GpuArray<Real, AMREX_SPACEDIM>& center,
                 Real radius_inner,
                 Real radius_outer,
                 Real z,
                 Real t,
                 Real dt,
                 const Box& domain)
        : m_plo(plo), m_phi(phi), m_center(center), m_radius_inner(radius_inner), m_radius_outer(radius_outer), m_z(z), m_t(t), m_dt(dt), m_domain(domain)
    {}
#if 1
    //c_light is the constant for speed of light [cm/s]
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        int result=0;
        if(m_radius_inner<=0 || m_radius_outer<=0)
            return 0;

        if(src.m_aos[i].id()>0) {

            Real xlen, ylen, zlen;
            
            Real lenx = m_phi[0]-m_plo[0];
            Real leny = m_phi[1]-m_plo[1];
            Real lenz = m_phi[2]-m_plo[2];
            int maxind[3];
            maxind[0] = floor((m_radius_outer+lenx*0.5)/lenx);
            maxind[1] = floor((m_radius_outer+leny*0.5)/leny);
            maxind[2] = floor((m_radius_outer+lenz*0.5)/lenz);

            //printf("Value is %d\n", maxind);
            for(int idir=-maxind[0];idir<=maxind[0];idir++)
              for(int jdir=-maxind[1];jdir<=maxind[1];jdir++)
                for(int kdir=-maxind[2];kdir<=maxind[2];kdir++)
                    {
                        xlen = src.m_aos[i].rdata(0+1+3)+(idir)*(m_phi[0]-m_plo[0]) - m_center[0];
                        ylen = src.m_aos[i].rdata(1+1+3)+(jdir)*(m_phi[1]-m_plo[1]) - m_center[1];
                        zlen = src.m_aos[i].rdata(2+1+3)+(kdir)*(m_phi[2]-m_plo[2]) - m_center[2];
                        Real mag = sqrt(xlen*xlen+ylen*ylen+zlen*zlen);
                        result+=int(mag>m_radius_inner && mag<m_radius_outer);
                    }
        }
        return result;
    }
};

#else
    //c_light is the constant for speed of light [cm/s]
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    bool operator() (const SrcData& src, int i) const noexcept
    {
        bool result=false;
        if(m_radius_inner<=0 || m_radius_outer<=0)
            return false;
	if(src.m_aos[i].id()>0) {

			Real xlen, ylen, zlen;
            
			Real lenx = m_phi[0]-m_plo[0];
            Real leny = m_phi[1]-m_plo[1];
            Real lenz = m_phi[2]-m_plo[2];
            int maxind[3];
            maxind[0] = floor((m_radius_outer+lenx*0.5)/lenx);
            maxind[1] = floor((m_radius_outer+leny*0.5)/leny);
            maxind[2] = floor((m_radius_outer+lenz*0.5)/lenz);

            //printf("Value is %d\n", maxind);

        for(int idir=-maxind[0];idir<=maxind[0];idir++)
            for(int jdir=-maxind[1];jdir<=maxind[1];jdir++)
                for(int kdir=-maxind[2];kdir<=maxind[2];kdir++)
                    {
                        xlen = src.m_aos[i].rdata(0+1+3)+(idir)*(m_phi[0]-m_plo[0]) - m_center[0];
                        ylen = src.m_aos[i].rdata(1+1+3)+(jdir)*(m_phi[1]-m_plo[1]) - m_center[1];
                        zlen = src.m_aos[i].rdata(2+1+3)+(kdir)*(m_phi[2]-m_plo[2]) - m_center[2];
                        Real mag = sqrt(xlen*xlen+ylen*ylen+zlen*zlen);
                        result=result? true : (mag>m_radius_inner && mag<m_radius_outer);
			//     	                Print()<<xlen<<"\t"<<ylen<<"\t"<<zlen<<"\t"<<mag<<"\t"<<m_radius_inner<<"\t"<<m_radius_outer<<"\t"<<result<<std::endl;
                    }
	}
        return (result);
    }
};
#endif
struct ShellStoreFilter
{
    GpuArray<Real, AMREX_SPACEDIM> m_plo, m_phi, m_center;
    Real m_radius_inner, m_radius_outer, m_z, m_t, m_dt, m_dt_a_cur_inv;
    Box m_domain;

    ShellStoreFilter (const GpuArray<Real, AMREX_SPACEDIM>& plo,
                 const GpuArray<Real, AMREX_SPACEDIM>& phi,
                 const GpuArray<Real, AMREX_SPACEDIM>& center,
                 Real radius_inner,
                 Real radius_outer,
                 Real z,
                 Real t,
                 Real dt,
                 const Box& domain,
                 const Real dt_a_cur_inv)
        : m_plo(plo), m_phi(phi), m_center(center), m_radius_inner(radius_inner), m_radius_outer(radius_outer), m_z(z), m_t(t), m_dt(dt), m_domain(domain), m_dt_a_cur_inv(dt_a_cur_inv)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst, const SrcData& src,
                     int src_i, int dst_i, int index) const noexcept
    {
      int local_index=-1;
      int nc=AMREX_SPACEDIM;
      bool result=false;
      Real xlen = m_phi[0]-m_plo[0];
      Real ylen = m_phi[1]-m_plo[1];
      Real zlen = m_phi[2]-m_plo[2];
      int maxind[3];
      maxind[0] = floor((m_radius_outer+xlen*0.5)/xlen);
      maxind[1] = floor((m_radius_outer+ylen*0.5)/ylen);
      maxind[2] = floor((m_radius_outer+zlen*0.5)/zlen);

      for(int idir=-maxind[0];idir<=maxind[0];idir++)
        for(int jdir=-maxind[1];jdir<=maxind[1];jdir++)
          for(int kdir=-maxind[2];kdir<=maxind[2];kdir++)
            {
              xlen = src.m_aos[src_i].rdata(0+1+3)+(idir)*(m_phi[0]-m_plo[0]) - m_center[0];
              ylen = src.m_aos[src_i].rdata(1+1+3)+(jdir)*(m_phi[1]-m_plo[1]) - m_center[1];
              zlen = src.m_aos[src_i].rdata(2+1+3)+(kdir)*(m_phi[2]-m_plo[2]) - m_center[2];
              Real mag = sqrt(xlen*xlen+ylen*ylen+zlen*zlen);

              if(int(mag>m_radius_inner && mag<m_radius_outer))
                local_index++;

              if(local_index==index) {
                dst.m_aos[dst_i].rdata(0)=src.m_aos[src_i].rdata(0);
                int comp=0;
                dst.m_aos[dst_i].pos(comp) = src.m_aos[src_i].pos(comp)+(idir)*(m_phi[comp]-m_plo[comp]);
                comp=1;
                dst.m_aos[dst_i].pos(comp) = src.m_aos[src_i].pos(comp)+(jdir)*(m_phi[comp]-m_plo[comp]);
                comp=2;
                dst.m_aos[dst_i].pos(comp) = src.m_aos[src_i].pos(comp)+(kdir)*(m_phi[comp]-m_plo[comp]);
                for (int comp=0; comp < nc; ++comp) {
                  dst.m_aos[dst_i].rdata(comp+1+3)=src.m_aos[src_i].pos(comp);
                  dst.m_aos[dst_i].rdata(comp+1+3+3) = src.m_aos[src_i].pos(comp) + m_dt_a_cur_inv * src.m_aos[src_i].rdata(comp+1);
                  dst.m_aos[dst_i].rdata(comp+1)=src.m_aos[src_i].rdata(comp+1);
                              //              p2.pos(comp)=p.pos(comp);
                  dst.m_aos[dst_i].id()=src.m_aos[src_i].id();
                  dst.m_aos[dst_i].cpu()=src.m_aos[src_i].cpu();
                }
                return;
              }
              //                      Print()<<xlen<<"\t"<<ylen<<"\t"<<zlen<<"\t"<<mag<<"\t"<<m_radius_inner<<"\t"<<m_radius_outer<<"\t"<<result<<std::endl;
            }
    }
};

  //  the steps would be count number of output particles
  //  prefix sum the counts
  //  resize the destination tile
  //  then use the prefix sum to a big threaded copy
/**
 * \brief Reorders the ParticleTile into two partitions
 * left [0, num_left-1] and right [num_left, ptile.numParticles()-1]
 * and returns the number of particles in the left partition.
 *
 * The functor is_left [(ParticleTileData ptd, int index) -> bool] maps each particle to
 * either the left [return true] or the right [return false] partition.
 * It must return the same result if evaluated multiple times for the same particle.
 *
 * \param ptile the ParticleTile to partition
 * \param is_left functor to map particles to a partition
 */
template <typename DstTile, typename SrcTile, typename Index, typename N,
    std::enable_if_t<std::is_integral_v<Index>, int> foo = 0, typename Pred, typename F>
Index filterAndTransformAndCopyParticles (DstTile& dst, const SrcTile& src,
                                          Index src_start, Index dst_start, N n, Pred&& p, F&& f) noexcept
{
    long np = n;
    const auto & src_data = src.getConstParticleTileData();
 
    if (np == 0) { return 0; }
    Gpu::DeviceVector<Index> mask_vec(np);
    auto * mask=mask_vec.dataPtr();
    const auto & ptd = src;
    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
      mask[i] = int(p(src_data,i));
    });

    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());
    const int num_output = offsets[np-1];
    dst.resize(num_output);

    auto dst_data = dst.getParticleTileData();

    Index last_mask=0, last_offset=0;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    auto const* p_offsets = offsets.dataPtr();

    // f here should do the copy and transform parts (or maybe instead just do copy here, although we want the images specifically)
    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
      if(mask[i]>0)
        for(int j=0;j<mask[i];j++)
            f(dst_data,src_data, src_start+i,
              dst_start+p_offsets[src_start+i]+j,j);
    });
    return num_output;
}

template <typename PC, typename F>
void filterParticles (PC& pc, F&& f)
{
    BL_PROFILE("filterParticles");

    using ParIter = typename PC::ParIterType;
    using ParticleTileType = typename PC::ParticleTileType;

    for (int lev = 0; lev <= pc.finestLevel(); ++lev)
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& ptile = pc.ParticlesAt(lev, pti);

            ParticleTileType ptile_tmp;
            ptile_tmp.resize(ptile.size());

            auto num_output = amrex::filterParticles(ptile_tmp, ptile, f);

            ptile.swap(ptile_tmp);
            ptile.resize(num_output);
        }
    }
}


}

void
DarkMatterParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
                                            int                    lev,
                                            amrex::Real            t,
                                            amrex::Real            dt,
                                            amrex::Real            a_old,
                                            amrex::Real            a_half,
                                            int                    where_width,
                                            amrex::Real            radius_inner,
                                            amrex::Real            radius_outer)
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

    auto pc= this;
    auto ShellPC = new ParticleContainer<10,0>(pc->Geom(lev), pc->ParticleDistributionMap(lev), pc->ParticleBoxArray(lev));
    ShellPC->resizeData();
    ParticleContainer<10,0>::ParticleInitData pdata = {{}, {}, {}, {}};
    ShellPC->InitNRandomPerCell(1,pdata);
    //    ShellPC->resize(pc.TotalNumberOfParticles());
    auto geom_test=pc->Geom(lev);
    const GpuArray<Real,AMREX_SPACEDIM> phi=geom_test.ProbHiArray();
    const GpuArray<Real,AMREX_SPACEDIM> center({AMREX_D_DECL((phi[0]-plo[0])*0.5,(phi[1]-plo[1])*0.5,(phi[2]-plo[2])*0.5)});
    //const Real a_half = 0.5 * (a_old + a_new);

    auto domain=geom_test.Domain();
    auto z=1/a_old-1;
    //From write_info.cpp
    amrex::Real a_cur        = 1.0 / a_half;
    amrex::Real a_cur_inv    = 1.0 / a_cur;
    amrex::Real dt_a_cur_inv = dt * a_cur_inv;
    ShellFilter shell_filter_test(plo, phi, center, radius_inner, radius_outer, z, t, dt, domain);
    ShellStoreFilter shell_store_filter_test(plo, phi, center, radius_inner, radius_outer, z, t, dt, domain, dt_a_cur_inv);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        ParticleType* pstruct = particles().data();
        const long np = pti.numParticles();
        int grid    = pti.index();

        auto& ptile = ShellPC->DefineAndReturnParticleTile(lev, pti);
        int old_np = ptile.size();
        int num_to_add = np;
        int new_np = old_np + num_to_add;
        ptile.resize(new_np);

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
    /*
    for (MFIter mfi(*ShellPC->m_dummy_mf[0], false); mfi.isValid(); ++mfi) {
        Box grid = ParticleBoxArray(0)[mfi.index()];
        auto ind = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        RealBox grid_box (grid,dx,geom.ProbLo());
        ParticleTile<ParticleType, NArrayReal, NArrayInt, amrex::PinnedArenaAllocator> ptile_tmp;

        Gpu::Device::streamSynchronize();
    }*/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        auto& particles = (this->ParticlesAt(lev,pti)).GetArrayOfStructs();

        auto* pstruct = particles().data();
        auto& ptile = ShellPC->ParticlesAt(lev,pti);

        auto& particles2 = (ShellPC->ParticlesAt(lev,pti)).GetArrayOfStructs();
        auto* pstruct2 = particles2().data();
#ifdef AMREX_USE_GPU
        const long np = pti.numParticles();
        int grid    = pti.index();
        ParticleContainer<10,0>::ParticleTileType ptile_tmp;

        auto num_output = filterAndTransformAndCopyParticles(ptile_tmp, ptile, 0, 0, np, shell_filter_test, shell_store_filter_test);
        ptile.swap(ptile_tmp);
        ptile.resize(num_output);
    }
#else
        auto ptile_tmp = ptile;
        ParticleContainer<10,0>::ParticleTileType ptile_tmp2;
	//        ptile_tmp.resize((ShellPC->ParticlesAt(lev,pti)).size());
        ptile_tmp.resize((this->ParticlesAt(lev,pti)).size());

        const long np = pti.numParticles();
        int grid    = pti.index();

        const FArrayBox& accel_fab= ((*ac_ptr)[0]);
        Array4<amrex::Real const> accel= accel_fab.array();

        int nc=AMREX_SPACEDIM;
	    int num_output=0;
        for(int i=0;i<np;i++) {
	    	int index=0;
        	store_dm_particle_single(pstruct[i],pstruct2[i],nc,
                                     accel,plo,phi,dxi,dt,a_old,
                                     a_half,do_move, radius_inner, radius_outer,index);
       const amrex::ParticleContainer<1+AMREX_SPACEDIM+6, 0>::SuperParticleType p = pstruct2[i];
       int mask=shell_filter_test(ptile.getConstParticleTileData(),i);
       if(mask>0) {
	   for(int j=0;j<mask;j++) {
           store_dm_particle_single(pstruct[i],pstruct2[i],nc,
                                    accel,plo,phi,dxi,dt,a_old,
                                    a_half,do_move, radius_inner, radius_outer,j,true);
           ptile_tmp2.push_back(pstruct2[i]);
	   num_output++;
	   }
       }
    }
    //        auto num_output = amrex::filterParticles(ptile_tmp, ptile, shell_filter_test);
        ptile.swap(ptile_tmp2);
        ptile.resize(num_output);
}
#endif

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
                             update_dm_particle_move_single(pstruct[i],nc,
                                                       accel,
                                                       plo,dxi,dt,a_old, a_half,do_move);
                           });

    }

    auto dir=Concatenate("plt_light", int(100*(1/a_old-1)), 7);
    auto name="ShellPC";
    amrex::Vector<std::string> real_comp_names_shell;
    real_comp_names_shell.clear();
    real_comp_names_shell.push_back("mass");
    real_comp_names_shell.push_back("xvel");
    real_comp_names_shell.push_back("yvel");
    real_comp_names_shell.push_back("zvel");
    real_comp_names_shell.push_back("xposold");
    real_comp_names_shell.push_back("yposold");
    real_comp_names_shell.push_back("zposold");
    real_comp_names_shell.push_back("xposvalid");
    real_comp_names_shell.push_back("yposvalid");
    real_comp_names_shell.push_back("zposvalid");
    std::string compression = "None@0";
    if(radius_inner>0&&radius_outer>radius_inner) {
        int write_hdf5=0;
#ifdef AMREX_USE_HDF5
        write_hdf5=1;
#endif
        if(write_hdf5!=1){
            //ShellPC->WritePlotFile(dir, name, real_comp_names_shell);
		}
#ifdef AMREX_USE_HDF5
        else{
            ShellPC->WritePlotFileHDF5(dir, name, real_comp_names_shell, compression);
		}
#endif
    }
    Print()<<"After write\t"<<ShellPC->TotalNumberOfParticles()<<"\t"<<a_old<<"\t"<<do_move<<"\t"<<lev<<"\t"<<t<<"\t"<<dt<<"\t"<<a_half<<"\t"<<where_width<<"\t"<<radius_inner<<std::endl;
    //    ShellPC->amrex::ParticleContainer<7,0>::WritePlotFile(dir, name, real_comp_names_shell);
    if (ac_ptr != &acceleration) delete ac_ptr;
    
    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        if (! m_particle_locator.isValid(GetParGDB())) m_particle_locator.build(GetParGDB());
        m_particle_locator.setGeometry(GetParGDB());
        AmrAssignGrid<DenseBinIteratorFactory<Box>> assign_grid = m_particle_locator.getGridAssignor();

        amrex::ParticleLocData pld; 
        for (auto& kv : pmap) {
            AoS&  particles = kv.second.GetArrayOfStructs();
            ParticleType* pstruct = particles().data();
            const long np = particles.size();
            amrex::ParallelFor(np,
                           [=] AMREX_GPU_HOST_DEVICE ( long i)
                           {
                               //                              amrex::ParticleContainer<4, 0>::SuperParticleType&  p=pstruct[i];
                               auto&  p=pstruct[i];
                               if(p.id()>0) {
                               const auto tup = assign_grid(p, lev, lev, where_width);
                               auto p_boxes = amrex::get<0>(tup);
                               auto p_levs  = amrex::get<1>(tup);
                               if(p_boxes<0||p_levs<0) {
                                   if (p.id() == amrex::GhostParticleID)
                                   {
                                       p.id() = -1;
                                   }
                                   else
                                   {
                                       amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                                   }
                               }
                               }
                           });
            Gpu::streamSynchronize();
        }
    }
	delete ShellPC;
}

void
DarkMatterParticleContainer::moveKick (MultiFab&       acceleration,
                                       int             lev,
                                       Real            t,
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
                             update_dm_particle_move_single(pstruct[i],nc,
                                                       accel,
                                                       plo,dxi,dt,a_half,a_new,do_move);
                           });
    }

    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void update_dm_particle_single (amrex::ParticleContainer<1+AMREX_SPACEDIM, 0>::SuperParticleType&  p,
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

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void store_dm_particle_single (amrex::ParticleContainer<1+AMREX_SPACEDIM, 0>::SuperParticleType&  p,
                               amrex::ParticleContainer<1+AMREX_SPACEDIM+6, 0>::SuperParticleType&  p2,
                               const int nc,
                               amrex::Array4<amrex::Real const> const& acc,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& phi,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Real& dt, const amrex::Real& a_prev,
                               const amrex::Real& a_cur, const int& do_move, const Real& radius_inner, 
							   const Real& radius_outer, int index, const bool is_file_write)
{
    amrex::Real half_dt       = 0.5 * dt;
    amrex::Real a_cur_inv    = 1.0 / a_cur;
    amrex::Real dt_a_cur_inv = dt * a_cur_inv;
    const GpuArray<Real,AMREX_SPACEDIM> center({AMREX_D_DECL((phi[0]-plo[0])*0.5,(phi[1]-plo[1])*0.5,(phi[2]-plo[2])*0.5)});

    if (do_move == 1) {
		p2.rdata(0)=p.rdata(0);
		bool result=false;
		Real lenx = phi[0]-plo[0];
		Real leny = phi[1]-plo[1];
		Real lenz = phi[2]-plo[2];
		int maxind[3]; 
		maxind[0] = floor((radius_outer+lenx*0.5)/lenx);
		maxind[1] = floor((radius_outer+leny*0.5)/leny);
		maxind[2] = floor((radius_outer+lenz*0.5)/lenz);

		Real xlen, ylen, zlen;
		int local_index=-1;
		for(int idir=-maxind[0];idir<=maxind[0];idir++) {
			for(int jdir=-maxind[1];jdir<=maxind[1];jdir++) {
				for(int kdir=-maxind[2];kdir<=maxind[2];kdir++) {
					xlen = p.pos(0)+(idir)*(phi[0]-plo[0]) - center[0];
					ylen = p.pos(1)+(jdir)*(phi[1]-plo[1]) - center[1];
					zlen = p.pos(2)+(kdir)*(phi[2]-plo[2]) - center[2];
					Real mag = sqrt(xlen*xlen+ylen*ylen+zlen*zlen);
					result=result? true : (mag>radius_inner && mag<radius_outer);

					if(int(mag>radius_inner && mag<radius_outer)){
						local_index++;
					}

					if(local_index==index) {
						int comp=0;
						p2.pos(comp) = p.pos(comp)+(idir)*(phi[comp]-plo[comp]);
						Real x1 = p2.pos(comp);
						Real vx = p.rdata(comp+1);
						comp=1;
						p2.pos(comp) = p.pos(comp)+(jdir)*(phi[comp]-plo[comp]);
						Real y1 = p2.pos(comp);
						Real vy = p.rdata(comp+1);
						comp=2;
						p2.pos(comp) = p.pos(comp)+(kdir)*(phi[comp]-plo[comp]);
						Real z1 = p2.pos(comp);
						Real vz = p.rdata(comp+1);

						// There are two calls to this routine. Do the particle output writing only once 
						// based on the bool.
						if(is_file_write) {
							LightConeParticle tmp;
							tmp.x = x1; tmp.y = y1; tmp.z = z1;
							tmp.vx = vx; tmp.vy = vy; tmp.vz = vz;
							shell_particles.emplace_back(tmp);
						}

						for (int comp=0; comp < nc; ++comp) {
							p2.rdata(comp+1+3)=p.pos(comp);
							p2.rdata(comp+1+3+3) = p.pos(comp) + dt_a_cur_inv * p.rdata(comp+1);
							p2.rdata(comp+1)=p.rdata(comp+1);
							//              p2.pos(comp)=p.pos(comp);
							p2.id()=p.id();
							p2.cpu()=p.cpu();
						}
						return;
					}
				}
			}
		}
	}

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void update_dm_particle_move_single (amrex::ParticleContainer<1+AMREX_SPACEDIM, 0>::SuperParticleType&  p,
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
    
    Real        mean_disp[AMREX_SPACEDIM]={AMREX_D_DECL(0,0,0)};


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
                    IntVect indices(AMREX_D_DECL(ix, jx, kx));
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

