#include "Nyx.H"
#include "Nyx_F.H"

#define BL_ARR4_TO_FORTRAN_3D(a) a.p,&((a).begin.x),amrex::GpuArray<int,3>{(a).end.x-1,(a).end.y-1,(a).end.z-1}.data()
#define BL_ARR4_TO_FORTRAN(a) (a).p, AMREX_ARLIM(&((a).begin.x)), (a).end.x-1,(a).end.y-1,(a).end.z-1

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
    
    hydro_source.setVal(0.0);

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

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
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
  {

    amrex::Cuda::setLaunchRegion(true);

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).

    FArrayBox flatn;
    FArrayBox dq;
    FArrayBox Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc;
    FArrayBox sm, sp;
    FArrayBox shk;
    FArrayBox qxm, qxp;
    FArrayBox qym, qyp;
    FArrayBox qzm, qzp;
    FArrayBox div;
    FArrayBox q_int;
    FArrayBox ftmp1, ftmp2;
    FArrayBox qgdnvtmp1, qgdnvtmp2;
    FArrayBox ql, qr;
    FArrayBox flux[AMREX_SPACEDIM];
    FArrayBox qe[AMREX_SPACEDIM];
    FArrayBox qmyx, qpyx;
    FArrayBox qmzx, qpzx;
    FArrayBox qmxy, qpxy;
    FArrayBox qmzy, qpzy;
    FArrayBox qmxz, qpxz;
    FArrayBox qmyz, qpyz;
    FArrayBox pdivu;

    /*    //    FArrayBox* fab_flatn = &flatn;
    const auto fab_flatn = flatn.array();
    const auto fab_dq = dq.array();
    const auto fab_Ip = Ip.array();
    const auto fab_Im = Im.array();
    const auto fab_Ip_src = Ip_src.array();
    const auto fab_Im_src = Im_src.array();
    const auto fab_Ip_gc = Ip_gc.array();
    const auto fab_Im_gc = Im_gc.array();
    const auto fab_sm = sm.array();
    const auto fab_sp = sp.array();
    const auto fab_shk = shk.array();
    const auto fab_qxm = qxm.array();
    const auto fab_qxp = qxp.array();
    const auto fab_qym = qym.array();
    const auto fab_qyp = qyp.array();
    const auto fab_qzm = qzm.array();
    const auto fab_qzp = qzp.array();
    const auto fab_div = div.array();
    const auto fab_q_int = q_int.array();
    const auto fab_ftmp1 = ftmp1.array();
    const auto fab_ftmp2 = ftmp2.array();
    const auto fab_qgdnvtmp1 = qgdnvtmp1.array();
    const auto fab_qgdnvtmp2 = qgdnvtmp2.array();
    const auto fab_ql = ql.array();
    const auto fab_qr = qr.array();

    const auto fab_qyx = qmyx.array();
    const auto fab_qpyx = qpyx.array();
    const auto fab_qzx = qmzx.array();
    const auto fab_qpzx = qpzx.array();
    const auto fab_qxy = qmxy.array();
    const auto fab_qpxy = qpxy.array();
    const auto fab_qzy = qmzy.array();
    const auto fab_qpzy = qpzy.array();
    const auto fab_qxz = qmxz.array();
    const auto fab_qpxz = qpxz.array();
    const auto fab_qyz = qmyz.array();
    const auto fab_qpyz = qpyz.array();
    const auto fab_pdivu = pdivu.array();*/

    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_flux{
      AMREX_D_DECL(flux[0].array(), flux[1].array(),flux[2].array())};
    /*    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_flux;
    AMREX_D_TERM(fab_flux[0] = flux[0].array();,
		 fab_flux[1] = flux[1].array();,
		 fab_flux[2] = flux[2].array(););*/

    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_qe{
      AMREX_D_DECL(qe[0].array(), qe[1].array(),qe[2].array())};
    /*    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_qe;
    AMREX_D_TERM(fab_qe[0] = qe[0].array();,
                 fab_qe[1] = qe[1].array();,
                 fab_qe[2] = qe[2].array(););*/


    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      //      for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      FArrayBox* fab_Sborder = Sborder.fabPtr(mfi);
      FArrayBox* fab_hydro_source = hydro_source.fabPtr(mfi);

      FArrayBox* fab_q = q.fabPtr(mfi);
      FArrayBox* fab_qaux = qaux.fabPtr(mfi);
      FArrayBox* fab_src_q = src_q.fabPtr(mfi);

      GpuArray<FArrayBox*, AMREX_SPACEDIM> fab_area;
      AMREX_D_TERM(fab_area[0] = area[0].fabPtr(mfi);,
		   fab_area[1] = area[1].fabPtr(mfi);,
		   fab_area[2] = area[2].fabPtr(mfi););

      GpuArray<FArrayBox*, AMREX_SPACEDIM> fab_fluxes;
      AMREX_D_TERM(fab_fluxes[0] = fluxes[0].fabPtr(mfi);,
		   fab_fluxes[1] = fluxes[1].fabPtr(mfi);,
		   fab_fluxes[2] = fluxes[2].fabPtr(mfi););
      
      //      q.resize(obx, 1);
      //      Elixir elix_q = q.elixir();
      
      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();
      Array4<Real> fab_flatn = flatn.array();
      amrex::Print()<<bx<<obx<<std::endl;
      // compute the flattening coefficient
      // compute the flattening coefficient

      Array4<Real> flatn_arr = flatn.array();
      ////      amrex::Print()<<"created array"<<std::endl;
      int pres_comp = QPRES;

      //amrex::Cuda::setLaunchRegion(true);
      if (use_flattening == 1) {
      AMREX_LAUNCH_DEVICE_LAMBDA(obx, tobx,
		{
        ca_uflatten(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(*fab_q),
		    BL_ARR4_TO_FORTRAN_3D(fab_flatn),
		    /*                    flatn_arr.p,(int[AMREX_SPACEDIM]){flatn_arr.begin.x,flatn_arr.begin.y,flatn_arr.begin.z},
					  (int[AMREX_SPACEDIM]){flatn_arr.end.x-1,flatn_arr.end.y-1,flatn_arr.end.z-1}, */
		    pres_comp);
		});
      } else {
        AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 1.0; });
      }

      //      amrex::Print()<<"flatn"<<std::endl;
      //amrex::Cuda::setLaunchRegion(false);

      const Box& xbx = amrex::surroundingNodes(bx, 0);
      const Box& gxbx = amrex::grow(xbx, 1);
      const Box& ybx = amrex::surroundingNodes(bx, 1);
      const Box& gybx = amrex::grow(ybx, 1);
      const Box& zbx = amrex::surroundingNodes(bx, 2);
      const Box& gzbx = amrex::grow(zbx, 1);

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();

      qxm.resize(obx, QVAR);
      Elixir elix_qxm = qxm.elixir();

      qxp.resize(obx, QVAR);
      Elixir elix_qxp = qxp.elixir();

      qym.resize(obx, QVAR);
      Elixir elix_qym = qym.elixir();

      qyp.resize(obx, QVAR);
      Elixir elix_qyp = qyp.elixir();

      qzm.resize(obx, QVAR);
      Elixir elix_qzm = qzm.elixir();

      qzp.resize(obx, QVAR);
      Elixir elix_qzp = qzp.elixir();

    const auto fab_shk = shk.array();
    const auto fab_qxm = qxm.array();
    const auto fab_qxp = qxp.array();
    const auto fab_qym = qym.array();
    const auto fab_qyp = qyp.array();
    const auto fab_qzm = qzm.array();
    const auto fab_qzp = qzp.array();

    dq.resize(obx, QVAR);
    Elixir elix_dq = dq.elixir();
    const auto fab_dq = dq.array();

      if (ppm_type == 0) {

#pragma gpu
	      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(obx, tobx,
      {
        ctu_plm_states(AMREX_INT_ANYD(tobx.loVect()), AMREX_INT_ANYD(tobx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(*fab_q),
		       BL_ARR4_TO_FORTRAN_3D(fab_flatn),
                       BL_TO_FORTRAN_ANYD(*fab_qaux),
                       BL_TO_FORTRAN_ANYD(*fab_src_q),
                       BL_ARR4_TO_FORTRAN_3D(fab_shk),
                       BL_ARR4_TO_FORTRAN_3D(fab_dq),
                       BL_ARR4_TO_FORTRAN_3D(fab_qxm),
                       BL_ARR4_TO_FORTRAN_3D(fab_qxp),
                       BL_ARR4_TO_FORTRAN_3D(fab_qym),
                       BL_ARR4_TO_FORTRAN_3D(fab_qyp),
                       BL_ARR4_TO_FORTRAN_3D(fab_qzm),
                       BL_ARR4_TO_FORTRAN_3D(fab_qzp),
                       AMREX_REAL_ANYD(dx), dt,
		       a_old, a_new,
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //      amrex::Print()<<"1"<<std::endl;
      amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

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
    const auto fab_Ip = Ip.array();
    const auto fab_Im = Im.array();
    const auto fab_Ip_src = Ip_src.array();
    const auto fab_Im_src = Im_src.array();
    const auto fab_Ip_gc = Ip_gc.array();
    const auto fab_Im_gc = Im_gc.array();
    const auto fab_sm = sm.array();
    const auto fab_sp = sp.array();
#pragma gpu
        ctu_ppm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(*fab_q),
                       BL_ARR4_TO_FORTRAN_3D(fab_flatn),
                       BL_TO_FORTRAN_ANYD(*fab_qaux),
                       BL_TO_FORTRAN_ANYD(*fab_src_q),
                       BL_ARR4_TO_FORTRAN_3D(fab_shk),
                       BL_ARR4_TO_FORTRAN_3D(fab_Ip),
                       BL_ARR4_TO_FORTRAN_3D(fab_Im),
                       BL_ARR4_TO_FORTRAN_3D(fab_Ip_src),
                       BL_ARR4_TO_FORTRAN_3D(fab_Im_src),
                       BL_ARR4_TO_FORTRAN_3D(fab_Ip_gc),
                       BL_ARR4_TO_FORTRAN_3D(fab_Im_gc),
                       BL_ARR4_TO_FORTRAN_3D(fab_sm),
                       BL_ARR4_TO_FORTRAN_3D(fab_sp),
                       BL_ARR4_TO_FORTRAN_3D(fab_qxm),
                       BL_ARR4_TO_FORTRAN_3D(fab_qxp),
                       BL_ARR4_TO_FORTRAN_3D(fab_qym),
                       BL_ARR4_TO_FORTRAN_3D(fab_qyp),
                       BL_ARR4_TO_FORTRAN_3D(fab_qzm),
                       BL_ARR4_TO_FORTRAN_3D(fab_qzp),
                       AMREX_REAL_ANYD(dx), dt,
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      }

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();
    const auto fab_div = div.array();
      // compute divu -- we'll use this later when doing the artifical viscosity
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(obx, tobx,
      {
      divu(AMREX_INT_ANYD(tobx.loVect()), AMREX_INT_ANYD(tobx.hiVect()),
           BL_TO_FORTRAN_ANYD(*fab_q),
           AMREX_REAL_ANYD(dx),
           BL_ARR4_TO_FORTRAN_3D(fab_div));
      });
      ///      amrex::Print()<<"1"<<std::endl;
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

      q_int.resize(obx, QVAR);
      Elixir elix_q_int = q_int.elixir();

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();

      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();

      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();

      ftmp1.resize(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1.elixir();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();

      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();

      ql.resize(obx, QVAR);
      Elixir elix_ql = ql.elixir();

      qr.resize(obx, QVAR);
      Elixir elix_qr = qr.elixir();

    const auto fab_q_int = q_int.array();
    const auto fab_ftmp1 = ftmp1.array();
    const auto fab_ftmp2 = ftmp2.array();
    const auto fab_qgdnvtmp1 = qgdnvtmp1.array();
    const auto fab_qgdnvtmp2 = qgdnvtmp2.array();
    const auto fab_ql = ql.array();
    const auto fab_qr = qr.array();

    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_flux{
      AMREX_D_DECL(flux[0].array(), flux[1].array(),flux[2].array())};
    /*    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_flux;
    AMREX_D_TERM(fab_flux[0] = flux[0].array();,
		 fab_flux[1] = flux[1].array();,
		 fab_flux[2] = flux[2].array(););*/

    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_qe{
      AMREX_D_DECL(qe[0].array(), qe[1].array(),qe[2].array())};
    /*    GpuArray<Array4<Real>, AMREX_SPACEDIM> fab_qe;
    AMREX_D_TERM(fab_qe[0] = qe[0].array();,
                 fab_qe[1] = qe[1].array();,
                 fab_qe[2] = qe[2].array(););*/

      const amrex::Real hdt = 0.5*dt;
      const amrex::Real a_half = 0.5 * (a_old + a_new);

      const amrex::Real hdtdx = 0.5*dt/dx[0]/a_half;
      const amrex::Real hdtdy = 0.5*dt/dx[1]/a_half;
      const amrex::Real hdtdz = 0.5*dt/dx[2]/a_half;

      const amrex::Real cdtdx = dt/dx[0]/3.0/a_half;
      const amrex::Real cdtdy = dt/dx[1]/3.0/a_half;
      const amrex::Real cdtdz = dt/dx[2]/3.0/a_half;

      // compute F^x
      // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cxbx, tcxbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcxbx.loVect()), AMREX_INT_ANYD(tcxbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qxm),
                          BL_ARR4_TO_FORTRAN_3D(fab_qxp), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      ///      amrex::Print()<<"1"<<std::endl;
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      qmyx.resize(tyxbx, QVAR);
      Elixir elix_qmyx = qmyx.elixir();

      qpyx.resize(tyxbx, QVAR);
      Elixir elix_qpyx = qpyx.elixir();

    const auto fab_qyx = qmyx.array();
    const auto fab_qpyx = qpyx.array();

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(tyxbx, ttyxbx,
      {
      transx_on_ystates(AMREX_INT_ANYD(ttyxbx.loVect()), AMREX_INT_ANYD(ttyxbx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qym),
                        BL_ARR4_TO_FORTRAN_3D(fab_qyx),
                        BL_ARR4_TO_FORTRAN_3D(fab_qyp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpyx),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        hdt, cdtdx);
      });
      ///      amrex::Print()<<"1xy"<<std::endl;
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmzx.resize(tzxbx, QVAR);
      Elixir elix_qmzx = qmzx.elixir();

      qpzx.resize(tzxbx, QVAR);
      Elixir elix_qpzx = qpzx.elixir();

    const auto fab_qzx = qmzx.array();
    const auto fab_qpzx = qpzx.array();

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(tzxbx, ttzxbx,
      {
      transx_on_zstates(AMREX_INT_ANYD(ttzxbx.loVect()), AMREX_INT_ANYD(ttzxbx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzm),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzx),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpzx),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        hdt, cdtdx);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cybx, tcybx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcybx.loVect()), AMREX_INT_ANYD(tcybx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qym),
                          BL_ARR4_TO_FORTRAN_3D(fab_qyp), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      qmxy.resize(txybx, QVAR);
      Elixir elix_qmxy = qmxy.elixir();

      qpxy.resize(txybx, QVAR);
      Elixir elix_qpxy = qpxy.elixir();
    const auto fab_qxy = qmxy.array();
    const auto fab_qpxy = qpxy.array();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(txybx, ttxybx,
      {
      transy_on_xstates(AMREX_INT_ANYD(ttxybx.loVect()), AMREX_INT_ANYD(ttxybx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxm),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxy),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpxy),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        cdtdy);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      qmzy.resize(tzybx, QVAR);
      Elixir elix_qmzy = qmzy.elixir();

      qpzy.resize(tzybx, QVAR);
      Elixir elix_qpzy = qpzy.elixir();

    const auto fab_qzy = qmzy.array();
    const auto fab_qpzy = qpzy.array();
      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(tzybx, ttzybx,
      {
      transy_on_zstates(AMREX_INT_ANYD(ttzybx.loVect()), AMREX_INT_ANYD(ttzybx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzm),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzy),
                        BL_ARR4_TO_FORTRAN_3D(fab_qzp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpzy),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        cdtdy);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(czbx, tczbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tczbx.loVect()), AMREX_INT_ANYD(tczbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qzm),
                          BL_ARR4_TO_FORTRAN_3D(fab_qzp), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmxz.resize(txzbx, QVAR);
      Elixir elix_qmxz = qmxz.elixir();

      qpxz.resize(txzbx, QVAR);
      Elixir elix_qpxz = qpxz.elixir();

    const auto fab_qxz = qmxz.array();
    const auto fab_qpxz = qpxz.array();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(txzbx, ttxzbx,
      {
      transz_on_xstates(AMREX_INT_ANYD(ttxzbx.loVect()), AMREX_INT_ANYD(ttxzbx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxm),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxz),
                        BL_ARR4_TO_FORTRAN_3D(fab_qxp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpxz),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        cdtdz);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      qmyz.resize(tyzbx, QVAR);
      Elixir elix_qmyz = qmyz.elixir();

      qpyz.resize(tyzbx, QVAR);
      Elixir elix_qpyz = qpyz.elixir();

    const auto fab_qyz = qmyz.array();
    const auto fab_qpyz = qpyz.array();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(tyzbx, ttyzbx,
      {
      transz_on_ystates(AMREX_INT_ANYD(ttyzbx.loVect()), AMREX_INT_ANYD(ttyzbx.hiVect()),
                        BL_ARR4_TO_FORTRAN_3D(fab_qym),
                        BL_ARR4_TO_FORTRAN_3D(fab_qyz),
                        BL_ARR4_TO_FORTRAN_3D(fab_qyp),
                        BL_ARR4_TO_FORTRAN_3D(fab_qpyz),
                        BL_TO_FORTRAN_ANYD(*fab_qaux),
                        BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                        BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                        cdtdz);
      });
      ////amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
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
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cyzbx, tcyzbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcyzbx.loVect()), AMREX_INT_ANYD(tcyzbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qyz),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpyz), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(czybx, tczybx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tczybx.loVect()), AMREX_INT_ANYD(tczybx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qzy),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpzy), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(xbx, txbx,
      {
      transyz(AMREX_INT_ANYD(txbx.loVect()), AMREX_INT_ANYD(txbx.hiVect()),
              BL_ARR4_TO_FORTRAN_3D(fab_qxm),
              BL_ARR4_TO_FORTRAN_3D(fab_ql),
              BL_ARR4_TO_FORTRAN_3D(fab_qxp),
              BL_ARR4_TO_FORTRAN_3D(fab_qr),
              BL_TO_FORTRAN_ANYD(*fab_qaux),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
              BL_TO_FORTRAN_ANYD(*fab_src_q),
              hdt, hdtdy, hdtdz, a_old, a_new);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(xbx, txbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(txbx.loVect()), AMREX_INT_ANYD(txbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_ql),
                          BL_ARR4_TO_FORTRAN_3D(fab_qr), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_flux[0]),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qe[0]),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
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
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(czxbx, tczxbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tczxbx.loVect()), AMREX_INT_ANYD(tczxbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qzx),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpzx), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cxzbx, tcxzbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcxzbx.loVect()), AMREX_INT_ANYD(tcxzbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qxz),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpxz), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      ////amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(ybx, tybx,
      {
      transxz(AMREX_INT_ANYD(tybx.loVect()), AMREX_INT_ANYD(tybx.hiVect()),
              BL_ARR4_TO_FORTRAN_3D(fab_qym),
              BL_ARR4_TO_FORTRAN_3D(fab_ql),
              BL_ARR4_TO_FORTRAN_3D(fab_qyp),
              BL_ARR4_TO_FORTRAN_3D(fab_qr),
              BL_TO_FORTRAN_ANYD(*fab_qaux),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
	      BL_TO_FORTRAN_ANYD(*fab_src_q),
              hdt, hdtdx, hdtdz, a_old, a_new);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(ybx, tybx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tybx.loVect()), AMREX_INT_ANYD(tybx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_ql),
                          BL_ARR4_TO_FORTRAN_3D(fab_qr), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_flux[1]),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qe[1]),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
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
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cxybx, tcxybx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcxybx.loVect()), AMREX_INT_ANYD(tcxybx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qxy),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpxy), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(cyxbx, tcyxbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tcyxbx.loVect()), AMREX_INT_ANYD(tcyxbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_qyx),
                          BL_ARR4_TO_FORTRAN_3D(fab_qpyx), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(zbx, tzbx,
      {
      transxy(AMREX_INT_ANYD(tzbx.loVect()), AMREX_INT_ANYD(tzbx.hiVect()),
              BL_ARR4_TO_FORTRAN_3D(fab_qzm),
              BL_ARR4_TO_FORTRAN_3D(fab_ql),
              BL_ARR4_TO_FORTRAN_3D(fab_qzp),
              BL_ARR4_TO_FORTRAN_3D(fab_qr),
              BL_TO_FORTRAN_ANYD(*fab_qaux),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp1),
              BL_ARR4_TO_FORTRAN_3D(fab_ftmp2),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp1),
              BL_ARR4_TO_FORTRAN_3D(fab_qgdnvtmp2),
	      BL_TO_FORTRAN_ANYD(*fab_src_q),
              hdt, hdtdx, hdtdy, a_old, a_new);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(zbx, tzbx,
      {
      cmpflx_plus_godunov(AMREX_INT_ANYD(tzbx.loVect()), AMREX_INT_ANYD(tzbx.hiVect()),
                          BL_ARR4_TO_FORTRAN_3D(fab_ql),
                          BL_ARR4_TO_FORTRAN_3D(fab_qr), 1, 1,
                          BL_ARR4_TO_FORTRAN_3D(fab_flux[2]),
                          BL_ARR4_TO_FORTRAN_3D(fab_q_int),
                          BL_ARR4_TO_FORTRAN_3D(fab_qe[2]),
                          BL_TO_FORTRAN_ANYD(*fab_qaux),
                          BL_ARR4_TO_FORTRAN_3D(fab_shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          int idir_f = idir + 1;

	  /*          Array4<Real> const flux_arr = (flux[idir]).array();
          /*const int temp_comp = Temp_comp;

          // Zero out shock and temp fluxes -- these are physically meaningless here
          AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
          {
              flux_arr(i,j,k,temp_comp) = 0.e0;
          });*/

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(nbx, tnbx,
      {
          apply_av(AMREX_INT_ANYD(tnbx.loVect()), AMREX_INT_ANYD(tnbx.hiVect()),
                   idir_f, AMREX_REAL_ANYD(dx),
                   BL_ARR4_TO_FORTRAN_3D(fab_div),
                   BL_TO_FORTRAN_ANYD(*fab_Sborder),
                   BL_ARR4_TO_FORTRAN_3D(fab_flux[idir]),&dt);
	  /*          if (0){//limit_fluxes_on_small_dens == 1) {
#pragma gpu
              limit_hydro_fluxes_on_small_dens
                  (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(*fab_Sborder),
                   BL_TO_FORTRAN_ANYD(*fab_q),
                   BL_TO_FORTRAN_ANYD(*volume[mfi]),
                   BL_ARR4_TO_FORTRAN_3D(fab_flux[idir]),
                   BL_TO_FORTRAN_ANYD(*area[idir][mfi]),
                   dt, AMREX_REAL_ANYD(dx));
		   }*/
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(nbx, tnbx,
      {
	  normalize_species_fluxes(AMREX_INT_ANYD(tnbx.loVect()), AMREX_INT_ANYD(tnbx.hiVect()),
                                   BL_ARR4_TO_FORTRAN_3D(fab_flux[idir]));
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

	  //	  amrex::Print()<<"max flux["<<idir<<"]="<<flux[idir].max()<<std::endl;
      }


      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();
    const auto fab_pdivu = pdivu.array();

      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(bx, tbx,
      {
      ca_consup(AMREX_INT_ANYD(tbx.loVect()), AMREX_INT_ANYD(tbx.hiVect()),
                BL_TO_FORTRAN(*fab_Sborder),
                BL_TO_FORTRAN(*fab_hydro_source),
                BL_ARR4_TO_FORTRAN(fab_flux[0]),
                BL_ARR4_TO_FORTRAN(fab_flux[1]),
                BL_ARR4_TO_FORTRAN(fab_flux[2]),
                BL_ARR4_TO_FORTRAN_3D(fab_qe[0]),
                BL_ARR4_TO_FORTRAN_3D(fab_qe[1]),
                BL_ARR4_TO_FORTRAN_3D(fab_qe[2]),
                BL_ARR4_TO_FORTRAN_3D(fab_pdivu),
                AMREX_REAL_ANYD(dx),&dt,&a_old,&a_new);
            });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

#pragma gpu
      //amrex::Cuda::setLaunchRegion(true);
      AMREX_LAUNCH_DEVICE_LAMBDA(nbx, tnbx,
      {
        scale_flux(AMREX_INT_ANYD(tnbx.loVect()), AMREX_INT_ANYD(tnbx.hiVect()),
                   BL_ARR4_TO_FORTRAN_3D(fab_flux[idir]),
                   BL_TO_FORTRAN_ANYD(*fab_area[idir]), dt,&a_old,&a_new);
      });
      //amrex::Gpu::Device::streamSynchronize();
      //amrex::Cuda::setLaunchRegion(false);
        if (idir == 0) {
            // get the scaled radial pressure -- we need to treat this specially
            Array4<Real> const qex_fab = qe[idir].array();
            const int prescomp = GDPRES;
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
      } // idir loop

      //took out track_grid_losses
      amrex::Gpu::Device::streamSynchronize();
    } // MFIter loop

    amrex::Cuda::setLaunchRegion(false);

  } // OMP loop

  amrex::Gpu::Device::streamSynchronize();
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
