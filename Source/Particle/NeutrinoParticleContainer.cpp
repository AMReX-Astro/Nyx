#ifdef NEUTRINO_PARTICLES

#include <NeutrinoParticleContainer.H>
#include <NeutrinoParticles_K.H>

using namespace amrex;

#ifndef NEUTRINO_DARK_PARTICLES
void
NeutrinoParticleContainer::AssignRelativisticDensity (Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled,
                                                      int               lev_min,
                                                      int               ncomp,
                                                      int               finest_level,
                                                      int               ngrow) const
{    
    BL_PROFILE("NeutrinoParticleContainer::AssignRelativisticDensity()");
    
    BL_ASSERT(NStructReal >= 1);
    BL_ASSERT(NStructReal >= ncomp);
    BL_ASSERT(ncomp == AMREX_SPACEDIM+1);
    
    if (finest_level == -1) {
        finest_level = this->finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level)) {
        finest_level--;
    }

    ngrow = std::max(ngrow, 2);
    
    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; lev++)
    {        
        auto ng = lev == lev_min ? IntVect(AMREX_D_DECL(ngrow,ngrow,ngrow)) : m_gdb->refRatio(lev-1);
        mf_to_be_filled[lev].reset(new MultiFab(m_gdb->boxArray(lev),
                                                m_gdb->DistributionMap(lev), ncomp, ng));
        mf_to_be_filled[lev]->setVal(0.0);
    }
    
    // Test whether the grid structure of the boxArray is the same as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        if (!OnSameGrids(lev, *mf_to_be_filled[lev]))
        {
            all_grids_the_same = false;
            break;
        }
    }
    
    Vector<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            auto ng = lev == lev_min ? IntVect(AMREX_D_DECL(ngrow,ngrow,ngrow)) : m_gdb->refRatio(lev-1);
            mf_part[lev].reset(new MultiFab(ParticleBoxArray(lev), 
                                            ParticleDistributionMap(lev), ncomp, ng));
            mf_part[lev]->setVal(0.0);
        }
    }
    
    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;
    
    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignRelativisticDensitySingleLevel(*mf[0], 0, ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            MultiFab::ParallelCopy(*mf_to_be_filled[0],*mf[0],0,0,ncomp,0,0);
        }
        return;
    }
    
    int lo_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir}; // periodic boundaries
    int hi_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    PCInterp mapper;
    int rho_index = 0;
    
    Vector<std::unique_ptr<MultiFab> > tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const BoxArray& ba = mf[lev]->boxArray();
        const DistributionMapping& dm = mf[lev]->DistributionMap();
        tmp[lev].reset(new MultiFab(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }
    
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        
        AssignRelativisticDensitySingleLevel(*mf[lev], lev, 1, 0);

        if (lev < finest_level) {
          PhysBCFunct<BndryFuncArray> cphysbc(this->m_gdb->Geom(lev), bcs,
                                              BndryFuncArray([](Real* data,
                                                                AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                                                const int* dom_lo, const int* dom_hi,
                                                                const Real* dx, const Real* grd_lo,
                                                                const Real* time, const int* bc){}));
          PhysBCFunct<BndryFuncArray> fphysbc(this->m_gdb->Geom(lev+1), bcs,
                                              BndryFuncArray([](Real* data,
                                                                AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                                                const int* dom_lo, const int* dom_hi,
                                                                const Real* dx, const Real* grd_lo,
                                                                const Real* time, const int* bc){}));

          amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf[lev],
                                       rho_index, rho_index, ncomp,
                                       this->m_gdb->Geom(lev),
                                       this->m_gdb->Geom(lev+1),
                                       cphysbc, 0, fphysbc, 0,
                                       this->m_gdb->refRatio(lev), &mapper,
                                       bcs, 0);

        }
        
        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(*mf[lev], *mf[lev-1], 0, 1,
                                      m_gdb->refRatio(lev-1),
                                      m_gdb->Geom(lev-1),
                                      m_gdb->Geom(lev));
        }
        
        mf[lev]->plus(*tmp[lev], 0, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf[lev+1], *mf[lev], 0, ncomp, m_gdb->refRatio(lev));
    }
    
    if (!all_grids_the_same) {
        for (int lev = lev_min; lev <= finest_level; lev++) {
            MultiFab::ParallelCopy(*mf_to_be_filled[lev],*mf_part[lev],0,0,1,0,0);
        }
    }
    if (lev_min > 0) {
        int nlevels = finest_level - lev_min + 1;
        for (int i = 0; i < nlevels; i++)
            {
                mf_to_be_filled[i] = std::move(mf_to_be_filled[i+lev_min]);
            }
        mf_to_be_filled.resize(nlevels);
    }
}

//
// This is the single-level version for cell-centered density
//
void
NeutrinoParticleContainer::AssignRelativisticDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                                 int       lev,
                                                                 int       ncomp,
                                                                 int       particle_lvl_offset) const
{
    BL_PROFILE("NeutrinoParticleContainer::AssignCellDensitySingleLevel()");

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(ParticleBoxArray(lev),
                                ParticleDistributionMap(lev),
                                ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1)
       amrex::Error("Must have at least one ghost cell when in AssignRelativisticDensitySingleLevel");

#ifdef _OPENMP
    const int       ng          = mf_pointer->nGrow();
#endif
    const Real      strttime    = amrex::second();
    const Geometry& gm          = Geom(lev);
    const Real*     dx_particle = Geom(lev + particle_lvl_offset).CellSize();

    const auto dx               = Geom(lev).CellSize();
    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto plo              = Geom(lev).ProbLoArray();
    const auto pdxi             = Geom(lev + particle_lvl_offset).InvCellSizeArray();

    if (Geom(lev).isAnyPeriodic() && ! Geom(lev).isAllPeriodic()) {
        amrex::Error(
            "AssignCellDensitySingleLevel: problem must be periodic in no or all directions"
            );
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*mf_pointer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        (*mf_pointer)[mfi].setVal(0);
    }

//  using ParConstIter = ParConstIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_rho;
        for (MyConstParIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const auto pstruct = particles().data();
            const long np = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            Real* data_ptr;
            const int *lo, *hi;
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(ng);
            local_rho.resize(tile_box,ncomp);
            local_rho = 0.0;
            data_ptr = local_rho.dataPtr();
            lo = tile_box.loVect();
            hi = tile_box.hiVect();
            auto rhoarr = local_rho.array();
#else
            const Box& box = fab.box();
            data_ptr = fab.dataPtr();
            lo = box.loVect();
            hi = box.hiVect();
            auto rhoarr = fab.array();
#endif

            if (particle_lvl_offset == 0) {
                if (m_relativistic) {
                AMREX_FOR_1D( np, i,
                {
                    neutrino_deposit_relativistic_cic(pstruct[i], rhoarr, plo, dxi);
                });
                } else {
                AMREX_FOR_1D( np, i,
                {
                    amrex_deposit_cic(pstruct[i], ncomp, rhoarr, plo, dxi);
                });
                }
            } else {
                if (m_relativistic) {
                    AMREX_FOR_1D( np, i,
                    {
                        neutrino_deposit_particle_dx_relativistic_cic(pstruct[i],
                                                                      rhoarr, plo, dxi);
                    }
                } else {
                AMREX_FOR_1D( np, i,
                {
                    amrex_deposit_particle_dx_cic(pstruct[i], ncomp, rhoarr, plo, dxi, pdxi);
                });
                }
            }

#ifdef _OPENMP
            fab.atomicAdd(local_rho, tile_box, tile_box, 0, 0, ncomp);
#endif

        }
    }

    mf_pointer->SumBoundary(Geom(lev).periodicity());

    // If ncomp > 1, first divide the momenta (component n)
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    for (int n = 1; n < ncomp; n++)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*mf_pointer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        MultiFab::ParallelCopy(mf_to_be_filled,*mf_pointer,0,0,ncomp,0,0);
        delete mf_pointer;
    }

    if (m_verbose > 1) {
      Real stoptime = amrex::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "NeutrinoParticleContainer::AssignRelativisticDensitySingleLevel time: " << stoptime << '\n';
    }
}
#endif

void
NeutrinoParticleContainer::moveKick (amrex::MultiFab&       acceleration,
                                            int                    lev,
                                            amrex::Real            dt,
                                            amrex::Real            a_new,
                                            amrex::Real            a_half)

{
  if(m_relativistic)
    amrex::Error("We shouldnt actually call moveKick for neutrinos");
  else
    {
    amrex::Error("MoveKick not defined for neutrinos");
    }
}

void
NeutrinoParticleContainer::moveKickDrift (MultiFab&       acceleration,
                                       int             lev,
                                       Real            dt,
                                       Real            a_old,
                                          Real            a_half,
                                          int where_width)
{

  if(m_relativistic)
    amrex::Error("We shouldnt actually call moveKickDrift for neutrinos");
  else
    {
      amrex::Error("MoveKickDrift not defined for neutrinos");
    }

}

#endif

