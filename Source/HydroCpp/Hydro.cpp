#include "Hydro.H"

/**
 *  Set up the source terms to go into the hydro.
 */
void
Nyx::construct_hydro_source(
  const amrex::MultiFab& S,
        amrex::MultiFab& sources_for_hydro,
  amrex::MultiFab& hydro_source, 
  amrex::Real time,
  amrex::Real dt,
  bool init_flux_register, bool add_to_flux_register)
{
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Computing hydro advance" << std::endl;
    }

    sources_for_hydro.FillBoundary(geom.periodicity());
    hydro_source.setVal(0);

    int nGrowF = 0;
    // Compute_hydro_sources style
    amrex::MultiFab fluxes[BL_SPACEDIM];
    int finest_level = parent->finestLevel();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    amrex::FluxRegister* fine    = 0;
    amrex::FluxRegister* current = 0;

    if(finest_level!=0)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
            fluxes[j].setVal(0.0);
        }
        if (do_reflux)
        {
            if (level < finest_level)
            {
                fine = &get_flux_reg(level+1);
                if (init_flux_register)
                    fine->setVal(0);
            }
            if (level > 0) {
                current = &get_flux_reg(level);
            }
        }
    }

    const amrex::Real* dx = geom.CellSize();

    amrex::Real dx1 = dx[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      dx1 *= dx[dir];
    }

    std::array<amrex::Real, AMREX_SPACEDIM> dxD = {AMREX_D_DECL(dx1, dx1, dx1)};
    const amrex::Real* dxDp = &(dxD[0]);

    amrex::Real courno = -1.0e+200;

    amrex::MultiFab& S_new = get_new_data(State_Type);

    // note: the radiation consup currently does not fill these
    amrex::Real E_added_flux = 0.;
    amrex::Real mass_added_flux = 0.;
    amrex::Real xmom_added_flux = 0.;
    amrex::Real ymom_added_flux = 0.;
    amrex::Real zmom_added_flux = 0.;
    amrex::Real mass_lost = 0.;
    amrex::Real xmom_lost = 0.;
    amrex::Real ymom_lost = 0.;
    amrex::Real zmom_lost = 0.;
    amrex::Real eden_lost = 0.;
    amrex::Real xang_lost = 0.;
    amrex::Real yang_lost = 0.;
    amrex::Real zang_lost = 0.;

    BL_PROFILE_VAR("Nyx::advance_hydro_pc_umdrv()", PC_UMDRV);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())               \
    reduction(+:E_added_flux,mass_added_flux)                           \
    reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux)    \
    reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost)        \
    reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)        \
    reduction(max:courno)
#endif
    {
      // amrex::IArrayBox bcMask[AMREX_SPACEDIM];
      amrex::Real cflLoc = -1.0e+200;
      int is_finest_level = (level == finest_level) ? 1 : 0;
      int flag_nscbc_isAnyPerio = (geom.isAnyPeriodic()) ? 1 : 0;
      int flag_nscbc_perio[AMREX_SPACEDIM]; // For 3D, we will know which
                                            // corners have a periodicity
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        flag_nscbc_perio[dir] =
          (amrex::DefaultGeometry().isPeriodic(dir)) ? 1 : 0;
      }

      const int* domain_lo = geom.Domain().loVect();
      const int* domain_hi = geom.Domain().hiVect();

      // Temporary Fabs needed for Hydro Computation
      for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& qbx = amrex::grow(bx, NUM_GROW + nGrowF);
        const amrex::Box& fbx = amrex::grow(bx, nGrowF);
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();

        amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
        amrex::Elixir flux_eli[AMREX_SPACEDIM];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          const amrex::Box& efbx = surroundingNodes(fbx, dir);
          flux[dir].resize(efbx, S.nComp());
          flux_eli[dir] = flux[dir].elixir();
        }

        auto const& s = S.array(mfi);
        auto const& hyd_src = hydro_source.array(mfi);

        // Resize Temporary Fabs
        amrex::FArrayBox q(qbx, QVAR), qaux(qbx, NQAUX), src_q(qbx, QVAR);
        // Use Elixir Construct to steal the Fabs metadata
        amrex::Elixir qeli = q.elixir();
        amrex::Elixir qauxeli = qaux.elixir();
        amrex::Elixir src_qeli = src_q.elixir();
        // Get Arrays to pass to the gpu.
        auto const& qarr = q.array();
        auto const& qauxar = qaux.array();
        auto const& srcqarr = src_q.array();

        BL_PROFILE_VAR("Nyx::ctoprim()", ctop);
        amrex::ParallelFor(
          qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, qarr, qauxar);
          });
        BL_PROFILE_VAR_STOP(ctop);

        // TODO GPUize NCSCBC
        // Imposing Ghost-Cells Navier-Stokes Characteristic BCs if "UserBC" are
        // used For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp.
        // 3399-3408, 2017.
        //
        // The user should provide a bcnormal routine in bc_fill_module with
        // additional optional arguments to temporary fill ghost-cells for
        // EXT_DIR and to provide target BC values. See the examples.

        // Allocate fabs for bcMask. Note that we grow in the opposite direction
        // because the Riemann solver wants a face value in a ghost-cell
#if 0 
      for (int dir = 0; dir < AMREX_SPACEDIM ; dir++)  {
        const Box& bxtmp = amrex::surroundingNodes(fbx,dir);
        Box TestBox(bxtmp);
        for(int d=0; d<AMREX_SPACEDIM; ++d) {
          if (dir!=d) TestBox.grow(d,1);
        }
        bcMask[dir].resize(TestBox,1);
        bcMask[dir].setVal(0);
      }
      
      // Becase bcMask is read in the Riemann solver in any case,
      // here we put physbc values in the appropriate faces for the non-nscbc case
      set_bc_mask(lo, hi, domain_lo, domain_hi,
                  AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                               AMREX_TO_FORTRAN(bcMask[1]),
                               AMREX_TO_FORTRAN(bcMask[2])));

      if (nscbc_adv == 1)
      {
        impose_NSCBC(lo, hi, domain_lo, domain_hi,
                     AMREX_TO_FORTRAN(*statein),
                     AMREX_TO_FORTRAN(q.fab()),
                     AMREX_TO_FORTRAN(qaux.fab()),
                     AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                            AMREX_TO_FORTRAN(bcMask[1]),
                            AMREX_TO_FORTRAN(bcMask[2])),
                     &flag_nscbc_isAnyPerio, flag_nscbc_perio, 
                     &time, dx, &dt);
      }
#endif
        BL_PROFILE_VAR("Nyx::srctoprim()", srctop);
        const auto& src_in = sources_for_hydro.array(mfi);
        amrex::ParallelFor(
          qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr);
          });
        BL_PROFILE_VAR_STOP(srctop);

        amrex::FArrayBox pradial(amrex::Box::TheUnitBox(), 1);
        if (!amrex::DefaultGeometry().IsCartesian()) {
          pradial.resize(amrex::surroundingNodes(bx, 0), 1);
        }
        amrex::Elixir pradial_eli = pradial.elixir();

#ifdef AMREX_USE_GPU
        auto device = amrex::RunOn::Gpu;
#else
        auto device = amrex::RunOn::Cpu;
#endif
        BL_PROFILE_VAR("Nyx::umdrv()", purm);
        amrex::MultiFab             volume;
        amrex::MultiFab             area[3];
        amrex::MultiFab             dLogArea[1];
        amrex::Vector< amrex::Vector<amrex::Real> > radius;

        volume.clear();
        volume.define(grids,dmap,1,NUM_GROW);
        geom.GetVolume(volume);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            area[dir].clear();
            area[dir].define(getEdgeBoxArray(dir),dmap,1,NUM_GROW);
            geom.GetFaceArea(area[dir],dir);
        }
        for (int dir = BL_SPACEDIM; dir < 3; dir++)
        {
            area[dir].clear();
            area[dir].define(grids, dmap, 1, 0);
            area[dir].setVal(0.0);
        }

        const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          flx_arr{
            AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())};
        const amrex::GpuArray<
          const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
          a{AMREX_D_DECL(
            area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))};
        pc_umdrv(
          is_finest_level, time, bx, domain_lo, domain_hi, phys_bc.lo(),
          phys_bc.hi(), s, hyd_src, qarr, qauxar, srcqarr, dx, dt, flx_arr, a,
          volume.array(mfi), cflLoc);
        BL_PROFILE_VAR_STOP(purm);

        BL_PROFILE_VAR("courno", crno);
        courno = amrex::max(courno, cflLoc);

        // Replacing YAFluxRegister or EBFluxRegister function with copy:
        //HostDevice::Atomic::Add(fluxes_fab(i,j,k,n),flux_fab(i,j,k,n));
        if(finest_level!=0)
        {
            for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
                amrex::Array4<amrex::Real> const flux_fab = (flux[idir]).array();
                amrex::Array4<amrex::Real> fluxes_fab = (fluxes[idir]).array(mfi);
                const int numcomp = NUM_STATE;
                fluxes[idir][mfi].prefetchToDevice();
                flux[idir].prefetchToDevice();

                AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
                {
                    fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
                });
            }
        }

        BL_PROFILE_VAR_STOP(crno);
      } // MFIter loop
      // These seem to check if the provided flux is a gpuptr, and use launches
      if (add_to_flux_register && finest_level!=0)
      {
        if (do_reflux) {
          if (current) {
            for (int i = 0; i < BL_SPACEDIM ; i++) {
              current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1);
            }
          }
          if (fine) {
            for (int i = 0; i < BL_SPACEDIM ; i++) {
              fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,amrex::FluxRegister::ADD);
            }
          }
        }
      }
    }   // end of OMP parallel region

    BL_PROFILE_VAR_STOP(PC_UMDRV);

#ifdef AMREX_DEBUG
//  if (print_energy_diagnostics) 
    {
      amrex::Real foo[5] = {
        E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux,
        mass_added_flux};

#ifdef AMREX_LAZY
      Lazy::QueueReduction([=]() mutable {
#endif
        amrex::ParallelDescriptor::ReduceRealSum(
          foo, 5, amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor()) {
          E_added_flux = foo[0];
          xmom_added_flux = foo[1];
          ymom_added_flux = foo[2];
          zmom_added_flux = foo[3];
          mass_added_flux = foo[4];
          amrex::Print() << "mass added from fluxes                      : "
                         << mass_added_flux << std::endl;
          amrex::Print() << "xmom added from fluxes                      : "
                         << xmom_added_flux << std::endl;
          amrex::Print() << "ymom added from fluxes                      : "
                         << ymom_added_flux << std::endl;
          amrex::Print() << "zmom added from fluxes                      : "
                         << zmom_added_flux << std::endl;
          amrex::Print() << "(rho E) added from fluxes                   : "
                         << E_added_flux << std::endl;
        }
#ifdef AMREX_LAZY
      });
#endif
    }
#endif

    if (courno > 1.0) {
      amrex::Print() << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level
                     << " IS " << courno << '\n';
    }
}

void
pc_umdrv(
  const int is_finest_level,
  const amrex::Real time,
  amrex::Box const& bx,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::Real* dx,
  const amrex::Real dt,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<amrex::Real> const& vol,
  amrex::Real cflLoc)
{
  //  Set Up for Hydro Flux Calculations
  auto const& bxg2 = grow(bx, 2);
  amrex::FArrayBox qec[AMREX_SPACEDIM];
  amrex::Elixir qec_eli[AMREX_SPACEDIM];
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const amrex::Box eboxes = amrex::surroundingNodes(bxg2, dir);
    qec[dir].resize(eboxes, NGDNV);
    qec_eli[dir] = qec[dir].elixir();
  }
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> qec_arr{
    AMREX_D_DECL(qec[0].array(), qec[1].array(), qec[2].array())};

  //  Temporary FArrayBoxes
  amrex::FArrayBox divu(bxg2, 1);
  amrex::FArrayBox pdivu(bx, 1);
  amrex::Elixir divueli = divu.elixir();
  amrex::Elixir pdiveli = pdivu.elixir();
  auto const& divarr = divu.array();
  auto const& pdivuarr = pdivu.array();

  BL_PROFILE_VAR("Nyx::umeth()", umeth);
#if AMREX_SPACEDIM == 1
  amrex::Abort("PLM isn't implemented in 1D.");
  //pc_umeth_1D(
  //  bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
  //  flx[0], qec_arr[0], a[0], pdivuarr, vol, dx, dt);
#elif AMREX_SPACEDIM == 2
  pc_umeth_2D(
    bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
    flx[0], flx[1], qec_arr[0], qec_arr[1], a[0], a[1], pdivuarr, vol, dx, dt);
#elif AMREX_SPACEDIM == 3
  pc_umeth_3D(
    bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
    flx[0], flx[1], flx[2], qec_arr[0], qec_arr[1], qec_arr[2], a[0], a[1],
    a[2], pdivuarr, vol, dx, dt);
#endif
  BL_PROFILE_VAR_STOP(umeth);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    qec_eli[dir].clear();
  }

  // divu
  AMREX_D_TERM(const amrex::Real dx0 = dx[0];, const amrex::Real dx1 = dx[1];
               , const amrex::Real dx2 = dx[2];);
  amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_divu(i, j, k, q, AMREX_D_DECL(dx0, dx1, dx2), divarr);
  });

  // consup
  amrex::Real difmag = 0.1;
  pc_consup(bx, uin, uout, flx, a, vol, divarr, pdivuarr, dx, difmag);
}

void
pc_consup(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& div,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const* del,
  amrex::Real const difmag)
{
  // Flux alterations
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::Box const& fbx = surroundingNodes(bx, dir);
    const amrex::Real dx = del[dir];
    amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_artif_visc(i, j, k, flx[dir], div, u, dx, difmag, dir);
      // Normalize Species Flux
      pc_norm_spec_flx(i, j, k, flx[dir]);
      // Make flux extensive
      pc_ext_flx(i, j, k, flx[dir], a[dir]);
    });
  }

  // Combine for Hydro Sources
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_update(i, j, k, update, flx, vol, pdivu);
  });
}
