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
        IntVect trr(D_DECL(1, 1, 1));

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
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
                    const Real val = ffab(p);
                    if (val > 0)
                        cfab(amrex::coarsen(p, trr)) += val;
                }
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
