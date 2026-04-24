#include <Nyx.H>
#include <Gravity.H>

using namespace amrex;

std::unique_ptr<MultiFab>
Nyx::particle_derive (const std::string& name, Real time, int ngrow)
{
#ifdef AMREX_PARTICLES
    if (Nyx::theDMPC() && name == "particle_count")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theDMPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#ifdef AGN
    else if (Nyx::theAPC() && name == "agn_particle_count")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theAPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#endif
#ifdef NEUTRINO_PARTICLES
    else if (Nyx::theNPC() && name == "neutrino_particle_count")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theNPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#endif
    else if (Nyx::theDMPC() && name == "total_particle_count")
    {
        //
        // We want the total particle count at this level or higher.
        //
        std::unique_ptr<MultiFab> derive_dat = particle_derive("particle_count", time, ngrow);
        IntVect trr(1);

        // @todo: level vs. lev
        for (int lev = level + 1; lev <= parent->finestLevel(); lev++)
        {
            auto ba = parent->boxArray(lev);
            const auto& dm = parent->DistributionMap(lev);
            MultiFab temp_dat(ba, dm, 1, 0);

            trr *= parent->refRatio(lev - 1);

            ba.coarsen(trr);
            MultiFab ctemp_dat(ba, dm, 1, 0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            Nyx::theDMPC()->Increment(temp_dat, lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab = temp_dat[mfi];
                FArrayBox& cfab = ctemp_dat[mfi];
                auto farr = temp_dat.array(mfi);
                auto carr = ctemp_dat.array(mfi);
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

                amrex::ParallelFor(fbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    auto p = IntVect(AMREX_D_DECL(i,j,k));
                    const Real val = farr(i,j,k);
                    if (val > 0)
                        carr(amrex::coarsen(p, trr)) += val;
                });
            }

            temp_dat.clear();

            MultiFab dat(grids, dmap, 1, 0);
            dat.setVal(0);
            dat.MultiFab::ParallelCopy(ctemp_dat, 0, 0, 1, 0, 0);

            MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
        }

        return derive_dat;
    }
    else if (Nyx::theDMPC() && name == "particle_mass_density")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 0, 0, 1, 0, 0);

        return derive_dat;
    }
    else if (Nyx::theDMPC() && name == "particle_x_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 1, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 1, 0, 1, 0, 0);

        return derive_dat;
    }
    else if (Nyx::theDMPC() && name == "particle_y_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 2, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 2, 0, 1, 0, 0);

        return derive_dat;
    }
    else if (Nyx::theDMPC() && name == "particle_z_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 3, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 3, 0, 1, 0, 0);

        return derive_dat;
    }
#ifdef AGN
    else if (Nyx::theAPC() && name == "agn_mass_density")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theAPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 0, 0, 1, 0, 0);

        return derive_dat;
    }
#endif
#ifdef NEUTRINO_PARTICLES
    else if (Nyx::theNPC() && name == "neutrino_mass_density")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 0, 0, 1, 0, 0);

        return derive_dat;
    }
#ifdef NEUTRINO_DARK_PARTICLES
    else if (Nyx::theNPC() && name == "neutrino_x_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 1, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 1, 0, 1, 0, 0);

        return derive_dat;
    }
    else if (Nyx::theNPC() && name == "neutrino_y_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 2, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 2, 0, 1, 0, 0);

        return derive_dat;
    }
    else if (Nyx::theNPC() && name == "neutrino_z_velocity")
    {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensityAndVels(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 3, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 3, 0, 1, 0, 0);

        return derive_dat;
    }
#else
    else if (Nyx::theNPC() && (name == "neutrino_x_velocity" || name == "neutrino_y_velocity" || name == "neutrino_z_velocity" ))
    {
        amrex::Print()<<"Returning mass density for neutrinos, since velocity not implemented for NEUTRINO_DARK_PARTICLES=FALSE"<<std::endl;
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 0, 0, 1, 0, 0);

        return derive_dat;
    }
    //////////////////////////////////////////////////////////
#endif
#endif
    else if (name == "total_density")
    {
      if (Nyx::theDMPC())
      {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensity(particle_mf);
       
        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        derive_dat->ParallelCopy(*particle_mf[level], 0, 0, 1, 0, 0);

#ifndef NO_HYDRO
        std::unique_ptr<MultiFab> gas_density = derive("density",time,0);
        MultiFab::Add(*derive_dat,*gas_density, 0, 0, 1, 0);
#endif
        return derive_dat; 
      }
      else 
      {
        return derive("density",time,0);
      }
    }
    else
#endif
    {
        return AmrLevel::derive(name, time, ngrow);
    }
}

#ifdef __cplusplus
extern "C"
{
#endif

  void dernull(const Box& /*bx*/, FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
               const FArrayBox& /*datfab*/, const Geometry& /*geomdata*/,
               Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    // This routine is used by particle_count.  Yes it does nothing.

  }

    void dermaggrav(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = std::sqrt(dat(i,j,k,0)*dat(i,j,k,0) +
                                 dat(i,j,k,1)*dat(i,j,k,1) +
                                 dat(i,j,k,2)*dat(i,j,k,2));

      });
    }

    void derdenvol(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      // Here dat contains (Density, Xmom, Ymom, Zmom_comp)
      const Real V_cell = dx[0] * dx[1] * dx[2];
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = V_cell * dat(i,j,k,0);

      });
    }

    void deroverden(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int level)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (Density, Xmom, Ymom, Zmom_comp)
      const Real over_den = Nyx::average_total_density * std::pow(Nyx::tagging_base,level+1);

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = dat(i,j,k,0) / over_den;

      });
    }

    void deroverdenzoom(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                        const FArrayBox& /*datfab*/, const Geometry& geomdata,
                        Real /*time*/, const int* /*bcrec*/, int level)
    {
      auto const der = derfab.array();

      //Assume Domain is a cube
      int idim = 0;
      int domlo = geomdata.Domain().smallEnd(idim);
      int domhi = geomdata.Domain().bigEnd(idim);

      int ref_size = domhi / (2*static_cast<int>(std::round(std::pow(2,(level+1)))));
      int center   = (domhi-domlo+1) / 2;

      auto const bx_ref = Box(IntVect(AMREX_D_DECL(amrex::max(center-ref_size+1, bx.smallEnd(0)),
                                                   amrex::max(center-ref_size+1, bx.smallEnd(1)),
                                                   amrex::max(center-ref_size+1, bx.smallEnd(2)))),
                              IntVect(AMREX_D_DECL(amrex::min(center+ref_size,   bx.bigEnd(0)),
                                                   amrex::min(center+ref_size,   bx.bigEnd(1)),
                                                   amrex::min(center+ref_size,   bx.bigEnd(2))) ));
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = 0.0;

      });
      amrex::ParallelFor(bx_ref,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = 1.0;

      });
      
    }

#ifdef __cplusplus
}
#endif

void
Nyx::compute_overdensity (const MultiFab& mf_pmd, MultiFab& mf_od)
{
    AMREX_ASSERT(mf_pmd.nComp() == 1);
    AMREX_ASSERT(mf_od.nComp() == 1);
    AMREX_ASSERT(mf_pmd.boxArray() == mf_od.boxArray());
    AMREX_ASSERT(mf_pmd.DistributionMap() == mf_od.DistributionMap());

    mf_od.define(mf_pmd.boxArray(), mf_pmd.DistributionMap(), 1, 0);

    // --------------------------------------------------
    // 1. Compute global mean density
    // --------------------------------------------------

    // FIX 1: correct signature is sum(comp, local) — no int nghost argument
    Real rho_sum = mf_pmd.sum(0, false);

    // FIX 2: BoxArray is fully replicated on every rank, so this is already
    // the global cell count. ReduceLongSum is both wrong (returns void) and
    // unnecessary (would multiply the count by nprocs).
    Long ncells = mf_pmd.boxArray().numPts();

    Real rho_mean = rho_sum / static_cast<Real>(ncells);

    // --------------------------------------------------
    // 2. Compute overdensity on GPU
    // δ = ρ / ρ̄ - 1
    // --------------------------------------------------
    for (MFIter mfi(mf_pmd, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& rho   = mf_pmd.const_array(mfi);
        auto const& delta = mf_od.array(mfi);
        Real inv_mean = 1.0_rt / rho_mean;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            delta(i,j,k) = rho(i,j,k) * inv_mean - 1.0_rt;
        });
    }
    if (mf_od.contains_nan()) {
        std::cout << "Value of rho_mean is " << rho_mean << " " << rho_sum << std::endl;
        amrex::Abort("NaNs detected in MultiFab");
    }
}

void
Nyx::compute_matter_power_spectrum(const MultiFab& mf_od)
{
    Geometry geom = parent->Geom(0);
    cMultiFab c_fft_pmd;
    {
        // Note that the complex Hermitian output array Y has (nx/2+1,ny,nz) elements.
        // Y[nx-i,j,k] = Y[i,j,k]*
        FFT::R2C<Real,FFT::Direction::forward> fft(geom.Domain());
        auto const& [ba, dm] = fft.getSpectralDataLayout();
        c_fft_pmd.define(ba,dm,1,0);
        fft.forward(mf_od, c_fft_pmd);
    }

    // For simplicity, we assume the domain is a cube, and we are not
    // going to worry about scaling.

    int nx = geom.Domain().length(0);
    int ny = geom.Domain().length(1);
    int nz = geom.Domain().length(2);
    int nk = int(std::sqrt((nx/2)*(nx/2) +
                       (ny/2)*(ny/2) +
                       (nz/2)*(nz/2))) + 1;
    Gpu::DeviceVector<Real> power_spec_mf_od_d(nk,Real(0.0));
    Gpu::DeviceVector<int> counts_d(nk, 0);
    Real* power_spec_mf_od_d_ptr = power_spec_mf_od_d.data();
    int* counts_d_ptr = counts_d.data();
    auto const& c_fft_pmd_arr = c_fft_pmd.const_arrays();
    Real Lx = geom.ProbLength(0);
    Real Ly = geom.ProbLength(1);
    Real Lz = geom.ProbLength(2);
    Real dx = geom.CellSize(0);
    Real dy = geom.CellSize(1);
    Real dz = geom.CellSize(2);
    Real kfund = 2.0 * M_PI / Lx;
    ParallelFor(c_fft_pmd, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        int ki = i;
        int kj = (j <= ny/2) ? j : ny-j;
        int kk = (k <= nz/2) ? k : nz-k;
        Real kmag = kfund * std::sqrt(ki*ki + kj*kj + kk*kk);
        int di = int(kmag / kfund);
        Real kx = 2.0 * M_PI * ki / Lx;
        Real ky = 2.0 * M_PI * kj / Ly;
        Real kz = 2.0 * M_PI * kk / Lz;

        Real wx = (kx == 0.0) ? 1.0 : sin(0.5*kx*dx)/(0.5*kx*dx);
        Real wy = (ky == 0.0) ? 1.0 : sin(0.5*ky*dy)/(0.5*ky*dy);
        Real wz = (kz == 0.0) ? 1.0 : sin(0.5*kz*dz)/(0.5*kz*dz);

        // CIC window (note the square per direction)
        Real W_CIC = (wx*wx) * (wy*wy) * (wz*wz);

        if (di < nk) {
        Real value = amrex::norm(c_fft_pmd_arr[b](i,j,k));
        // Account for Hermitian symmetry in x-direction
        // Hermitian symmetry Y[nx-i,j,k] = Y[i,j,k]*
        if ((i > 0) && (2*i != nx)) {
            // Multiply by 2 because we have +ki and -ki
            value *= Real(2.0);
        }
        value /= (W_CIC * W_CIC);
        int weight = ((i > 0) && (2*i != nx)) ? 2 : 1;
        HostDevice::Atomic::Add(power_spec_mf_od_d_ptr+di, value);
        HostDevice::Atomic::Add(counts_d_ptr + di, weight);
        }
    });

    Real* power_spec_mf_od_h_ptr = nullptr;
    int*  counts_h_ptr = nullptr;
    #ifdef AMREX_USE_GPU
        Gpu::HostVector<Real> power_spec_mf_od_h(counts_d.size());
        Gpu::HostVector<int> counts_h(counts_d.size());
        Gpu::copyAsync(Gpu::deviceToHost, power_spec_mf_od_d.begin(), power_spec_mf_od_d.end(), power_spec_mf_od_h.begin());
        Gpu::copyAsync(Gpu::deviceToHost, counts_d.begin(), counts_d.end(), counts_h.begin());
        Gpu::streamSynchronize();
        power_spec_mf_od_h_ptr = power_spec_mf_od_h.data();
        counts_h_ptr = counts_h.data();
    #else
        power_spec_mf_od_h_ptr = power_spec_mf_od_d.data();
        counts_h_ptr = counts_d.data();
    #endif

    ParallelDescriptor::ReduceRealSum(power_spec_mf_od_h_ptr, nk);
    ParallelDescriptor::ReduceIntSum(counts_h_ptr, nk);

    for (int i = 0; i < nk; ++i) {
        if (counts_h_ptr[i] > 0) {
            power_spec_mf_od_h_ptr[i] /= Real(counts_h_ptr[i]);
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        // --- Step 1: format redshift ---
        char zbuf[32];
        const Real prev_time = state[State_Type].prevTime();
        const Real a_old     = get_comoving_a(prev_time);
        Real z_old = 1/a_old-1;
        sprintf(zbuf, "%07d", int(z_old * 10000 + 0.5));  // → "0028165"

        // --- Step 2: ensure directory exists ---
        fs::path dir = "MatterPk";
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
        // --- Step 3: build full filename ---
        fs::path filepath = dir / ("spectrum_" + std::string(zbuf) + ".txt");

        // --- Step 4: open file ---
        std::ofstream ofs(filepath);

        if (!ofs) {
            throw std::runtime_error("Failed to open file: " + filepath.string());
        }
   
        Real dV = dx*dy*dz;
        Real dom_vol = Lx*Ly*Lz;
        Real dk = 2.0 * M_PI / Lx;
        for (int i = 0; i < nk; ++i) {
            Real k = dk * (i + 0.5);
            ofs << k << " " << power_spec_mf_od_h_ptr[i]* dV*dV/dom_vol << "\n";
        }
    }
}
