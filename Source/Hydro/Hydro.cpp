#include <Hydro.H>
#include <Utilities.H>
#include <Nyx.H>

/**
 *  Set up the source terms to go into the hydro.
*/

void
Nyx::construct_hydro_source(
  const amrex::MultiFab& S,
        amrex::MultiFab& sources_for_hydro,
  amrex::MultiFab& hydro_source,
  amrex::MultiFab& grav_vector, 
  amrex::Real a_old,
  amrex::Real a_new,
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
    amrex::MultiFab hydro_fluxes[AMREX_SPACEDIM];
    int finest_level = parent->finestLevel();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    amrex::FluxRegister* fine    = 0;
    amrex::FluxRegister* current = 0;

    if(finest_level!=0)
    {
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            hydro_fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
            hydro_fluxes[j].setVal(0.0);
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

    amrex::Real courno = -1.0e+200;

    amrex::MultiFab& S_new = get_new_data(State_Type);

    amrex::Real E_added_flux = 0.;
    amrex::Real mass_added_flux = 0.;
    amrex::Real xmom_added_flux = 0.;
    amrex::Real ymom_added_flux = 0.;
    amrex::Real zmom_added_flux = 0.;

    BL_PROFILE_VAR("Nyx::advance_hydro_pc_umdrv()", PC_UMDRV);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())               \
    reduction(+:E_added_flux,mass_added_flux)                           \
    reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux)    \
    reduction(max:courno)
#endif
    {
      // amrex::IArrayBox bcMask[AMREX_SPACEDIM];
      amrex::Real cflLoc = -1.0e+200;

      const amrex::Real a_dot = (a_new - a_old) / dt;

      const int NumSpec_loc = QVAR - NGDNV;
      const amrex::Real gamma_minus_1_loc = gamma-1.0;

      if (S.nComp() != QVAR) amrex::Print() << "NCOMP QVAR " << S.nComp() << " " << QVAR << std::endl;
      AMREX_ALWAYS_ASSERT(S.nComp() == QVAR);

      // Temporary Fabs needed for Hydro Computation
      for ( amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, Nyx::minimize_memory!=0 ? amrex::Gpu::synchronize() : amrex::Gpu::streamSynchronize())
      {
        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& qbx = amrex::grow(bx, NUM_GROW + nGrowF);
        const amrex::Box& fbx = amrex::grow(bx, nGrowF);

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
        amrex::FArrayBox     q(qbx, QVAR);
        amrex::FArrayBox src_q(qbx, QVAR);

        // Use Elixir Construct to steal the Fabs metadata
        amrex::Elixir     qeli = q.elixir();
        amrex::Elixir src_qeli = src_q.elixir();

        // Get Arrays to pass to the gpu.
        auto const&    qarr =     q.array();
        auto const& srcqarr = src_q.array();

        BL_PROFILE_VAR("Nyx::ctoprim()", ctop);
        amrex::ParallelFor(qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            pc_ctoprim(i, j, k, s, qarr, NumSpec_loc, gamma_minus_1_loc);
        });
        BL_PROFILE_VAR_STOP(ctop);

        BL_PROFILE_VAR("Nyx::srctoprim()", srctop);
        const auto&  src_in = sources_for_hydro.array(mfi);
        const auto& grav_in = grav_vector.array(mfi);

        amrex::ParallelFor(qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            pc_srctoprim(i, j, k, qarr, grav_in, src_in, srcqarr,
                         a_dot, NumSpec_loc, gamma_minus_1_loc);
        });
        BL_PROFILE_VAR_STOP(srctop);

        BL_PROFILE_VAR("Nyx::umdrv()", purm);

        const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          flx_arr{flux[0].array(), flux[1].array(), flux[2].array()};

        pc_umdrv(
          bx, s, hyd_src, qarr, srcqarr, flx_arr, dx, dt, a_old, a_new, cfl,
          gamma, gamma_minus_1_loc, NumSpec,
          small_dens, small_pres, small, cflLoc, ppm_type);
        BL_PROFILE_VAR_STOP(purm);

        BL_PROFILE_VAR("courno", crno);
        courno = amrex::max(courno, cflLoc);

        // Replacing YAFluxRegister or EBFluxRegister function with copy:
        //HostDevice::Atomic::Add(fluxes_fab(i,j,k,n),flux_fab(i,j,k,n));
        if(finest_level!=0)
        {
            for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
                amrex::Array4<amrex::Real> const flux_fab = (flux[idir]).array();
                amrex::Array4<amrex::Real> fluxes_fab = (hydro_fluxes[idir]).array(mfi);
                const int numcomp = NUM_STATE;
                hydro_fluxes[idir][mfi].prefetchToDevice();
                flux[idir].prefetchToDevice();

                AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
                {
                    fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
                });
            }
        }

        BL_PROFILE_VAR_STOP(crno);
      } // MFIter loop
    }   // end of OMP parallel region

    if (add_to_flux_register && finest_level!=0)
    {
      if (do_reflux) {
        if (current) {
          for (int i = 0; i < AMREX_SPACEDIM ; i++) {
            current->FineAdd(hydro_fluxes[i], i, 0, 0, NUM_STATE, 1);
          }
        }
        if (fine) { // HACK
          for (int i = 0; i < AMREX_SPACEDIM ; i++) {
            fine->CrseInit(hydro_fluxes[i],i,0,0,NUM_STATE,-1.,amrex::FluxRegister::ADD);
          }
        }
      }
    }

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
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& src_q, 
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::Real* dx,
  const amrex::Real dt,
  const amrex::Real a_old,
  const amrex::Real a_new,
  const amrex::Real cfl,
  const amrex::Real gamma, const amrex::Real gamma_minus_1, 
  const int NumSpec,
  const amrex::Real small_dens, const amrex::Real small_pres, 
  const amrex::Real small, 
  amrex::Real cflLoc,
  const int ppm_type) 

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

  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> qec_arr
    {qec[0].array(), qec[1].array(), qec[2].array()};

  //  Temporary FArrayBoxes
  amrex::FArrayBox divu(bxg2, 1);
  amrex::FArrayBox pdivu(bx, 1);
  amrex::Elixir divueli = divu.elixir();
  amrex::Elixir pdiveli = pdivu.elixir();
  auto const& divarr = divu.array();
  auto const& pdivuarr = pdivu.array();

  BL_PROFILE_VAR("Nyx::umeth()", umeth);
  pc_umeth_3D(
    bx, 
    q, src_q,
    flx[0], flx[1], flx[2], 
    qec_arr[0], qec_arr[1], qec_arr[2], 
    pdivuarr, dx, dt, a_old, a_new, NumSpec, gamma, gamma_minus_1,
    small_dens, small_pres, small, ppm_type);
  BL_PROFILE_VAR_STOP(umeth);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    qec_eli[dir].clear();
  }

#if 0
  // if (limit_fluxes_on_small_dens == 1) 
  {
      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) 
      {
          const Box& nbx = amrex::surroundingNodes(bx, idir);
          const Real& dx_dir = dx[idir];
          const Real& lcfl   = cfl;
          amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              limit_hydro_fluxes_on_small_dens (i, j, k, idir, uin, q, flx[idir], small_dens, lcfl, dx_dir, dt);
            });
      }
  }
#endif

  // divu
  const amrex::Real dx0 = dx[0]; 
  const amrex::Real dx1 = dx[1];
  const amrex::Real dx2 = dx[2];
  amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
  {
    pc_divu(i, j, k, q, dx0, dx1, dx2, divarr);
  });

  // consup
  amrex::Real difmag = 0.1;
  pc_consup(bx, uin, uout, flx, divarr, pdivuarr, a_old, a_new, dx, dt, NumSpec, gamma_minus_1, difmag);
}

void
pc_consup(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  amrex::Array4<const amrex::Real> const& div,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const a_old,
  amrex::Real const a_new,
  amrex::Real const* del,
  amrex::Real const dt,
  const int NumSpec,
  amrex::Real const gamma_minus_1,
  amrex::Real const difmag)
{

  // Flux alterations
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::Box const& fbx = surroundingNodes(bx, dir);
    const amrex::Real dx = del[dir];

    const GpuArray<Real,BL_SPACEDIM>  area{ del[1]*del[2], del[0]*del[2], del[0]*del[1] };

    amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

         pc_artif_visc(i, j, k, flx[dir], div, u, dx, difmag, dir);

         // Normalize Species Flux
         normalize_species_fluxes(i, j, k, flx[dir], NumSpec);

         // Make flux extensive
         pc_ext_flx(i, j, k, flx[dir], area[dir], dt);
      });
  }

  amrex::Real vol = del[0] * del[1] * del[2];

  // Combine for Hydro Sources
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_update(i, j, k, u, update, flx, vol, pdivu, a_old, a_new, dt, gamma_minus_1);
  });

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::Box const& fbx = surroundingNodes(bx, dir);
    amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Change scaling for flux registers
        pc_ext_flx_dt(i, j, k, flx[dir], a_old, a_new);
    });
  }
}
