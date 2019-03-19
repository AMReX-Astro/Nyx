#include "Nyx.H"
#include "Nyx_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Nyx::construct_ctu_hydro_source(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new,
                           MultiFab& Sborder, MultiFab& D_border, 
                           MultiFab& ext_src_old, MultiFab& hydro_source, 
                           MultiFab& grav_vector, MultiFab& divu_cc,
                           bool init_flux_register, bool add_to_flux_register) 
{

  BL_PROFILE("Nyx::construct_ctu_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering construct_ctu_hydro_source" << std::endl << std::endl;

  // Compute_hydro_sources style
  MultiFab fluxes[BL_SPACEDIM];
  const int finest_level = parent->finestLevel();

  for (int j = 0; j < BL_SPACEDIM; j++)
    {
      fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
      fluxes[j].setVal(0.0);
    }

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister* fine    = 0;
    FluxRegister* current = 0;

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
  /*/
  fluxes.resize(3);

  for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
      fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, NUM_STATE, 0));
      (*fluxes[dir]).setVal(0);
    }

  for (int dir = BL_SPACEDIM; dir < 3; ++dir)
    fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, NUM_STATE, 0));

  mass_fluxes.resize(3);

  for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    mass_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));

  for (int dir = BL_SPACEDIM; dir < 3; ++dir)
    mass_fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, 1, 0));
  */
///
/// The primitive variable state array.
///
    amrex::MultiFab q;

///
/// The auxiliary primitive variable state array.
///
    amrex::MultiFab qaux;

///
/// The source terms in primitive form.
///
    amrex::MultiFab src_q;

    q.define(grids, dmap, QVAR, NUM_GROW);
    q.setVal(0.0);
    // for qc
    qaux.define(grids, dmap, 1, NUM_GROW);

    src_q.define(grids, dmap, NQSRC, NUM_GROW);
    //src_q.setVal(0.0);

    amrex::MultiFab csml;
    csml.define(grids, dmap, 1, NUM_GROW);
    csml.setVal(0.0);

    cons_to_prim(Sborder, q, qaux, grav_vector, ext_src_old, src_q, csml, a_old, a_new, dt);

///
/// The data.
///
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

    dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0],grids,dmap,0,NUM_GROW);
#endif
    
    hydro_source.setVal(0.0);

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Nyx::construct_ctu_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  Real mass_lost = 0.;
  Real xmom_lost = 0.;
  Real ymom_lost = 0.;
  Real zmom_lost = 0.;
  Real eden_lost = 0.;
  Real xang_lost = 0.;
  Real yang_lost = 0.;
  Real zang_lost = 0.;

  BL_PROFILE_VAR("Nyx::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp) \
                     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#else
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
#endif
  {

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).

    FArrayBox flatn;
#ifdef RADIATION
    FArrayBox flatg;
#endif
    FArrayBox dq;
    FArrayBox Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc;
    FArrayBox sm, sp;
    FArrayBox shk;
    FArrayBox qxm, qxp;
#if AMREX_SPACEDIM >= 2
    FArrayBox qym, qyp;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qzm, qzp;
#endif
    FArrayBox div;
    FArrayBox q_int;
#ifdef RADIATION
    FArrayBox lambda_int;
#endif
#if AMREX_SPACEDIM >= 2
    FArrayBox ftmp1, ftmp2;
#ifdef RADIATION
    FArrayBox rftmp1, rftmp2;
#endif
    FArrayBox qgdnvtmp1, qgdnvtmp2;
    FArrayBox ql, qr;
#endif
    FArrayBox flux[AMREX_SPACEDIM];
    FArrayBox qe[AMREX_SPACEDIM];
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qmyx, qpyx;
    FArrayBox qmzx, qpzx;
    FArrayBox qmxy, qpxy;
    FArrayBox qmzy, qpzy;
    FArrayBox qmxz, qpxz;
    FArrayBox qmyz, qpyz;
#endif
    FArrayBox pdivu;

    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      //      for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      //      q.resize(obx, 1);
      //      Elixir elix_q = q.elixir();
      
      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();

#ifdef RADIATION
      flatg.resize(obx, 1);
      Elixir elix_flatg = flatg.elixir();
#endif

      q[mfi];
      // compute the flattening coefficient
      // remove first order check
      if (use_flattening == 1) {
#ifdef RADIATION
        ca_rad_flatten(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(flatg));
#else
#pragma gpu
        ca_uflatten(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(q[mfi]),
                    BL_TO_FORTRAN_ANYD(flatn), QPRES+1);
#endif
      } else {
        flatn.setVal(1.0, obx);
      }

      const Box& xbx = amrex::surroundingNodes(bx, 0);
      const Box& gxbx = amrex::grow(xbx, 1);
#if AMREX_SPACEDIM >= 2
      const Box& ybx = amrex::surroundingNodes(bx, 1);
      const Box& gybx = amrex::grow(ybx, 1);
#endif
#if AMREX_SPACEDIM == 3
      const Box& zbx = amrex::surroundingNodes(bx, 2);
      const Box& gzbx = amrex::grow(zbx, 1);
#endif

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();

      qxm.resize(obx, QVAR);
      Elixir elix_qxm = qxm.elixir();

      qxp.resize(obx, QVAR);
      Elixir elix_qxp = qxp.elixir();

#if AMREX_SPACEDIM >= 2
      qym.resize(obx, QVAR);
      Elixir elix_qym = qym.elixir();

      qyp.resize(obx, QVAR);
      Elixir elix_qyp = qyp.elixir();
#endif

#if AMREX_SPACEDIM == 3
      qzm.resize(obx, QVAR);
      Elixir elix_qzm = qzm.elixir();

      qzp.resize(obx, QVAR);
      Elixir elix_qzp = qzp.elixir();
#endif

      if (ppm_type == 0) {

        dq.resize(obx, QVAR);
        Elixir elix_dq = dq.elixir();

#pragma gpu
        ctu_plm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]),
                       BL_TO_FORTRAN_ANYD(shk),
                       BL_TO_FORTRAN_ANYD(dq),
                       BL_TO_FORTRAN_ANYD(qxm),
                       BL_TO_FORTRAN_ANYD(qxp),
#if AMREX_SPACEDIM >= 2
                       BL_TO_FORTRAN_ANYD(qym),
                       BL_TO_FORTRAN_ANYD(qyp),
#endif
#if AMREX_SPACEDIM == 3
                       BL_TO_FORTRAN_ANYD(qzm),
                       BL_TO_FORTRAN_ANYD(qzp),
#endif
                       AMREX_REAL_ANYD(dx), dt,
#if (AMREX_SPACEDIM < 3)
                       BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      } else {

        Ip.resize(obx, 3*QVAR);
        Elixir elix_Ip = Ip.elixir();

        Im.resize(obx, 3*QVAR);
        Elixir elix_Im = Im.elixir();

        Ip_src.resize(obx, 3*NQSRC);
        Elixir elix_Ip_src = Ip_src.elixir();

        Im_src.resize(obx, 3*NQSRC);
        Elixir elix_Im_src = Im_src.elixir();

        Ip_gc.resize(obx, 3);
        Elixir elix_Ip_gc = Ip_gc.elixir();

        Im_gc.resize(obx, 3);
        Elixir elix_Im_gc = Im_gc.elixir();

        sm.resize(obx, AMREX_SPACEDIM);
        Elixir elix_sm = sm.elixir();

        sp.resize(obx, AMREX_SPACEDIM);
        Elixir elix_sp = sp.elixir();

#pragma gpu
        ctu_ppm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]),
                       BL_TO_FORTRAN_ANYD(shk),
                       BL_TO_FORTRAN_ANYD(Ip),
                       BL_TO_FORTRAN_ANYD(Im),
                       BL_TO_FORTRAN_ANYD(Ip_src),
                       BL_TO_FORTRAN_ANYD(Im_src),
                       BL_TO_FORTRAN_ANYD(Ip_gc),
                       BL_TO_FORTRAN_ANYD(Im_gc),
                       BL_TO_FORTRAN_ANYD(sm),
                       BL_TO_FORTRAN_ANYD(sp),
                       BL_TO_FORTRAN_ANYD(qxm),
                       BL_TO_FORTRAN_ANYD(qxp),
#if AMREX_SPACEDIM >= 2
                       BL_TO_FORTRAN_ANYD(qym),
                       BL_TO_FORTRAN_ANYD(qyp),
#endif
#if AMREX_SPACEDIM == 3
                       BL_TO_FORTRAN_ANYD(qzm),
                       BL_TO_FORTRAN_ANYD(qzp),
#endif
                       AMREX_REAL_ANYD(dx), dt,
#if (AMREX_SPACEDIM < 3)
                       BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      }

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      // compute divu -- we'll use this later when doing the artifical viscosity
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div));

      q_int.resize(obx, QVAR);
      Elixir elix_q_int = q_int.elixir();

#ifdef RADIATION
      lambda_int.resize(obx, Radiation::nGroups);
      Elixir elix_lambda_int = lambda_int.elixir();
#endif

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();

#ifdef RADIATION
      rad_flux[0].resize(gxbx, Radiation::nGroups);
      Elixir elix_rad_flux_x = rad_flux[0].elixir();
#endif

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();

#ifdef RADIATION
      rad_flux[1].resize(gybx, Radiation::nGroups);
      Elixir elix_rad_flux_y = rad_flux[1].elixir();
#endif
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();

#ifdef RADIATION
      rad_flux[2].resize(gzbx, Radiation::nGroups);
      Elixir elix_rad_flux_z = rad_flux[2].elixir();
#endif
#endif

#if AMREX_SPACEDIM <= 2
      if (!Geometry::IsCartesian()) {
          pradial.resize(xbx, 1);
          Elixir elix_pradial = pradial.elixir();
      }
#endif

#if AMREX_SPACEDIM == 1
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

#endif // 1-d



#if AMREX_SPACEDIM >= 2
      ftmp1.resize(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1.elixir();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();

#ifdef RADIATION
      rftmp1.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp1 = rftmp1.elixir();

      rftmp2.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp2 = rftmp2.elixir();
#endif

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();

      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();

      ql.resize(obx, QVAR);
      Elixir elix_ql = ql.elixir();

      qr.resize(obx, QVAR);
      Elixir elix_qr = qr.elixir();
#endif



#if AMREX_SPACEDIM == 2

      const amrex::Real hdt = 0.5*dt;
      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];

      // compute F^x
      // [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^y
      // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fy
      // rftmp2 = rfy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cybx.loVect()), AMREX_INT_ANYD(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym),
                          BL_TO_FORTRAN_ANYD(qyp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // add the transverse flux difference in y to the x states
      // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]

      // ftmp2 = fy
      // rftmp2 = rfy
#pragma gpu
      transy_on_xstates(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(ql),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qr),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp2),
#endif
                        BL_TO_FORTRAN_ANYD(qe[1]),
                        hdtdy);

      // solve the final Riemann problem axross the x-interfaces

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // add the transverse flux difference in x to the y states
      // [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0]

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx

#pragma gpu
      transx_on_ystates(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(ql),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qr),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        BL_TO_FORTRAN_ANYD(area[0][mfi]),
                        BL_TO_FORTRAN_ANYD(volume[mfi]),
                        hdt, hdtdx);

      // solve the final Riemann problem axross the y-interfaces

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
#endif // 2-d



#if AMREX_SPACEDIM == 3

      const amrex::Real hdt = 0.5*dt;

      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];
      const amrex::Real hdtdz = 0.5*dt/dx[2];

      const amrex::Real cdtdx = dt/dx[0]/3.0;
      const amrex::Real cdtdy = dt/dx[1]/3.0;
      const amrex::Real cdtdz = dt/dx[2]/3.0;

      // compute F^x
      // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      qmyx.resize(tyxbx, QVAR);
      Elixir elix_qmyx = qmyx.elixir();

      qpyx.resize(tyxbx, QVAR);
      Elixir elix_qpyx = qpyx.elixir();

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
#pragma gpu
      transx_on_ystates(AMREX_INT_ANYD(tyxbx.loVect()), AMREX_INT_ANYD(tyxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(qmyx),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qpyx),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmzx.resize(tzxbx, QVAR);
      Elixir elix_qmzx = qmzx.elixir();

      qpzx.resize(tzxbx, QVAR);
      Elixir elix_qpzx = qpzx.elixir();

#pragma gpu
      transx_on_zstates(AMREX_INT_ANYD(tzxbx.loVect()), AMREX_INT_ANYD(tzxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm),
                        BL_TO_FORTRAN_ANYD(qmzx),
                        BL_TO_FORTRAN_ANYD(qzp),
                        BL_TO_FORTRAN_ANYD(qpzx),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cybx.loVect()), AMREX_INT_ANYD(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym),
                          BL_TO_FORTRAN_ANYD(qyp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      qmxy.resize(txybx, QVAR);
      Elixir elix_qmxy = qmxy.elixir();

      qpxy.resize(txybx, QVAR);
      Elixir elix_qpxy = qpxy.elixir();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      transy_on_xstates(AMREX_INT_ANYD(txybx.loVect()), AMREX_INT_ANYD(txybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(qmxy),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qpxy),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      qmzy.resize(tzybx, QVAR);
      Elixir elix_qmzy = qmzy.elixir();

      qpzy.resize(tzybx, QVAR);
      Elixir elix_qpzy = qpzy.elixir();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      transy_on_zstates(AMREX_INT_ANYD(tzybx.loVect()), AMREX_INT_ANYD(tzybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm),
                        BL_TO_FORTRAN_ANYD(qmzy),
                        BL_TO_FORTRAN_ANYD(qzp),
                        BL_TO_FORTRAN_ANYD(qpzy),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czbx.loVect()), AMREX_INT_ANYD(czbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qzm),
                          BL_TO_FORTRAN_ANYD(qzp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmxz.resize(txzbx, QVAR);
      Elixir elix_qmxz = qmxz.elixir();

      qpxz.resize(txzbx, QVAR);
      Elixir elix_qpxz = qpxz.elixir();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      transz_on_xstates(AMREX_INT_ANYD(txzbx.loVect()), AMREX_INT_ANYD(txzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(qmxz),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qpxz),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      qmyz.resize(tyzbx, QVAR);
      Elixir elix_qmyz = qmyz.elixir();

      qpyz.resize(tyzbx, QVAR);
      Elixir elix_qpyz = qpyz.elixir();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      transz_on_ystates(AMREX_INT_ANYD(tyzbx.loVect()), AMREX_INT_ANYD(tyzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(qmyz),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qpyz),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

      //
      // Use qx?, q?yz, q?zy to compute final x-flux
      //

      // compute F^{y|z}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fyz
      // rftmp1 = rfyz
      // qgdnvtmp1 = qgdnvyz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cyzbx.loVect()), AMREX_INT_ANYD(cyzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyz),
                          BL_TO_FORTRAN_ANYD(qpyz), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czybx.loVect()), AMREX_INT_ANYD(czybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzy),
                          BL_TO_FORTRAN_ANYD(qpzy), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

#pragma gpu
      transyz(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qxm),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qxp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdy, hdtdz);

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czxbx.loVect()), AMREX_INT_ANYD(czxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzx),
                          BL_TO_FORTRAN_ANYD(qpzx), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxzbx.loVect()), AMREX_INT_ANYD(cxzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxz),
                          BL_TO_FORTRAN_ANYD(qpxz), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

#pragma gpu
      transxz(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
              BL_TO_FORTRAN_ANYD(qym),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qyp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              hdt, hdtdx, hdtdz);

      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxybx.loVect()), AMREX_INT_ANYD(cxybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxy),
                          BL_TO_FORTRAN_ANYD(qpxy), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cyxbx.loVect()), AMREX_INT_ANYD(cyxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyx),
                          BL_TO_FORTRAN_ANYD(qpyx), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      transxy(AMREX_INT_ANYD(zbx.loVect()), AMREX_INT_ANYD(zbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qzm),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qzp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdx, hdtdy);

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(zbx.loVect()), AMREX_INT_ANYD(zbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[2]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[2]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[2]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

#endif // 3-d

      /*
      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
	amrex::Print()<<"max flux["<<idir<<"]="<<flux[idir].max()<<std::endl;
	amrex::Print()<<flux[idir]<<std::endl;
	}*/
      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          int idir_f = idir + 1;

	  /*          Array4<Real> const flux_arr = (flux[idir]).array();
          /*const int temp_comp = Temp_comp;
#ifdef SHOCK_VAR
          const int shk_comp = Shock;
#endif

          // Zero out shock and temp fluxes -- these are physically meaningless here
          AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
          {
              flux_arr(i,j,k,temp_comp) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,shk_comp) = 0.e0;
#endif
          });*/

#pragma gpu
          apply_av(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f, AMREX_REAL_ANYD(dx),
                   BL_TO_FORTRAN_ANYD(div),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]));

#ifdef RADIATION
#pragma gpu
          apply_av_rad(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                       idir_f, AMREX_REAL_ANYD(dx),
                       BL_TO_FORTRAN_ANYD(div),
                       BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                       BL_TO_FORTRAN_ANYD(rad_flux[idir]));
#endif

          if (0){//limit_fluxes_on_small_dens == 1) {
#pragma gpu
              limit_hydro_fluxes_on_small_dens
                  (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(volume[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                   dt, AMREX_REAL_ANYD(dx));
          }

#pragma gpu
          normalize_species_fluxes(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(flux[idir]));

	  //	  amrex::Print()<<"max flux["<<idir<<"]="<<flux[idir].max()<<std::endl;
      }


      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();

      ca_consup(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		BL_TO_FORTRAN(Sborder[mfi]),
		BL_TO_FORTRAN(hydro_source[mfi]),
		BL_TO_FORTRAN(flux[0]),
#if AMREX_SPACEDIM >= 2
                 BL_TO_FORTRAN(flux[1]),
#endif
#if AMREX_SPACEDIM == 3
                 BL_TO_FORTRAN(flux[2]),
#endif
		BL_TO_FORTRAN(div),
		AMREX_REAL_ANYD(dx),&dt,&a_old,&a_new);

      /*
      for(int idir=0;idir<=3;++idir)
	amrex::Print()<<"max flux["<<idir<<"]="<<flux[idir].max()<<std::endl;
      */
#ifdef RADIATION
      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif

#if AMREX_SPACEDIM <= 2
      Array4<Real> pradial_fab = pradial.array();
#endif

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

#pragma gpu
        scale_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
#if AMREX_SPACEDIM == 1
                   BL_TO_FORTRAN_ANYD(qe[idir]),
#endif
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

#ifdef RADIATION
#pragma gpu
        scale_rad_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                       BL_TO_FORTRAN_ANYD(rad_flux[idir]),
                       BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);
#endif

        if (idir == 0) {
            // get the scaled radial pressure -- we need to treat this specially
            Array4<Real> const qex_fab = qe[idir].array();
            const int prescomp = GDPRES;

#if AMREX_SPACEDIM == 1
            if (!Geometry::IsCartesian()) {
                AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                {
                    pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                });
            }
#endif

#if AMREX_SPACEDIM == 2
            if (!mom_flux_has_p[0][0]) {
                AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                {
                    pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                });
            }
#endif
        }

        // Store the fluxes from this advance.

        // For normal integration we want to add the fluxes from this advance
        // since we may be subcycling the timestep. But for simplified SDC integration
        // we want to copy the fluxes since we expect that there will not be
        // subcycling and we only want the last iteration's fluxes.

        Array4<Real> const flux_fab = (flux[idir]).array();
        Array4<Real> fluxes_fab = (fluxes[idir]).array(mfi);
        const int numcomp = NUM_STATE;

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
            {
                fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
            });

#ifdef RADIATION
        Array4<Real> const rad_flux_fab = (rad_flux[idir]).array();
        Array4<Real> rad_fluxes_fab = (*rad_fluxes[idir]).array(mfi);
        const int radcomp = Radiation::nGroups;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), radcomp, i, j, k, n,
            {
                rad_fluxes_fab(i,j,k,n) = rad_flux_fab(i,j,k,n);
            });

        } else {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), radcomp, i, j, k, n,
            {
                rad_fluxes_fab(i,j,k,n) += rad_flux_fab(i,j,k,n);
            });

        }
#endif

      } // idir loop

#if AMREX_SPACEDIM <= 2
      if (!Geometry::IsCartesian()) {

          Array4<Real> P_radial_fab = P_radial.array(mfi);

          if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
              {
                  P_radial_fab(i,j,k,0) = pradial_fab(i,j,k,0);
              });

          } else {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
              {
                  P_radial_fab(i,j,k,0) += pradial_fab(i,j,k,0);
              });

          }

      }
#endif

      for (int i = 0; i < BL_SPACEDIM; ++i) 
	fluxes[i][mfi].copy(flux[i], mfi.nodaltilebox(i));

      //took out track_grid_losses

    } // MFIter loop

  } // OMP loop


  BL_PROFILE_VAR_STOP(CA_UMDRV);

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Leaving construct_ctu_hydro_sources()" << std::endl << std::endl;

  if (verbose > 0)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Nyx::construct_ctu_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

    if (add_to_flux_register)
    {
       if (do_reflux) {
         if (current) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
             current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1);
           }
         }
         if (fine) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
	         fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
           }
         }
       }
    }

}
