#include <Godunov.H>
#include <PLM.H>
#include <PPM.H>

using namespace amrex;

extern void
trace_ppm(const Box& bx,
          const int idir,
          Array4<Real const> const& q_arr,
          Array4<Real const> const& srcQ,
          Array4<Real> const& qm,
          Array4<Real> const& qp,
          const Box& vbx,
          const Real dt, const Real* del,
          const Real gamma,
          const Real small_dens, const Real small_pres, 
#ifndef CONST_SPECIES
          const int NumSpec,
#endif
          const Real a_old);

// Host function to call gpu hydro functions
void
pc_umeth_3D(
  Box const& bx,
  Array4<const Real> const& q,
  Array4<const Real> const& srcQ, 
  Array4<Real> const& flx1,
  Array4<Real> const& flx2,
  Array4<Real> const& flx3, 
  Array4<Real> const& q1,
  Array4<Real> const& q2,
  Array4<Real> const& q3,
  Array4<Real> const& pdivu,
  const Real* del,
  const Real dt,
  const Real a_old,
  const Real a_new, 
  const Real gamma, const Real gamma_minus_1,
  const Real small_dens, const Real small_pres, 
  const Real small,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const int ppm_type)
{
#ifndef CONST_SPECIES
  const int FirstSpec_comp_loc = FirstSpec_comp;
  const int NumSpec_loc   = NumSpec;
#endif

  const Real a_half = 0.5 * (a_old + a_new);
  Real const dx = del[0];
  Real const dy = del[1];
  Real const dz = del[2];

  Real const hdtdx = 0.5 * dt / dx / a_half;
  Real const hdtdy = 0.5 * dt / dy / a_half;
  Real const hdtdz = 0.5 * dt / dz / a_half;
  Real const cdtdx = 1.0 / 3.0 * dt / dx / a_half;
  Real const cdtdy = 1.0 / 3.0 * dt / dy / a_half;
  Real const cdtdz = 1.0 / 3.0 * dt / dz / a_half;

  const amrex::Box& bxg1 = grow(bx, 1);
  const amrex::Box& bxg2 = grow(bx, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qxm(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qym(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  // Z data
  cdir = 2;
  const amrex::Box& zmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& zflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qzm(zmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qzmarr = qzm.array();
  auto const& qzparr = qzp.array();

  // Put the PLM and slopes in the same kernel launch to avoid unnecessary
  // launch overhead Pelec_Slope_* are SIMD as well as PeleC_plm_* which loop
  // over the same box
  
  if(ppm_type == 0 )
  {
      amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real slope[8];

      const amrex::Real c = std::sqrt((gamma_minus_1+1.0) * q(i,j,k,QPRES)/q(i,j,k,QRHO));

      // X slopes and interp
      for (int n = 0; n < QVAR; ++n)
          slope[n] = plm_slope(i, j, k, n, 0, q, small_pres);
      pc_plm_x(i, j, k, qxmarr, qxparr, srcQ, slope, q, c, a_old, dx, dt, 
#ifndef CONST_SPECIES
               NumSpec_loc,
#endif
               gamma_minus_1, small_dens, small_pres);

      // Y slopes and interp
      for (int n = 0; n < QVAR; n++)
          slope[n] = plm_slope(i, j, k, n, 1, q, small_pres);
      pc_plm_y(i, j, k, qymarr, qyparr, srcQ, slope, q, c, a_old, dy, dt,
#ifndef CONST_SPECIES
               NumSpec_loc,
#endif
               gamma_minus_1, small_dens, small_pres);

      // Z slopes and interp
      for (int n = 0; n < QVAR; ++n)
          slope[n] = plm_slope(i, j, k, n, 2, q, small_pres);
      pc_plm_z(i, j, k, qzmarr, qzparr, srcQ, slope, q, c, a_old, dz, dt,
#ifndef CONST_SPECIES
               NumSpec_loc,
#endif
               gamma_minus_1, small_dens, small_pres);

     });

  } else {

      // Compute the normal interface states by reconstructing
      // the primitive variables using the piecewise parabolic method
      // and doing characteristic tracing.  We do not apply the
      // transverse terms here.

      int idir = 0;
      trace_ppm(bxg2,
                idir,
                q, srcQ, 
                qxmarr, qxparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
#ifndef CONST_SPECIES
                NumSpec_loc,
#endif
                a_old);

      idir = 1;
      trace_ppm(bxg2,
                idir,
                q, srcQ,
                qymarr, qyparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
#ifndef CONST_SPECIES
                NumSpec_loc,
#endif
                a_old);

      idir = 2;
      trace_ppm(bxg2,
                idir,
                q, srcQ,
                qzmarr, qzparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
#ifndef CONST_SPECIES
                NumSpec_loc,
#endif
                a_old);
  }


  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  cdir = 0;
  amrex::FArrayBox fx(xflxbx, flx1.nComp(), amrex::The_Async_Arena());
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(xflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempx = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, qxmarr, qxparr, fxarr, gdtempx, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // Y initial fluxes
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, flx2.nComp(), amrex::The_Async_Arena());
  auto const& fyarr = fy.array();
  amrex::FArrayBox qgdy(yflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempy = qgdy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, qymarr, qyparr, fyarr, gdtempy, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // Z initial fluxes
  cdir = 2;
  amrex::FArrayBox fz(zflxbx, flx3.nComp(), amrex::The_Async_Arena());
  auto const& fzarr = fz.array();
  amrex::FArrayBox qgdz(zflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempz = qgdz.array();
  amrex::ParallelFor(
    zflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, qzmarr, qzparr, fzarr, gdtempz, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // X interface corrections
  cdir = 0;
  const amrex::Box& txbx = grow(bxg1, cdir, 1);
  const amrex::Box& txbxm = growHi(txbx, cdir, 1);
  amrex::FArrayBox qxym(txbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxyp(txbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxy = qxym.array();
  auto const& qpxy = qxyp.array();

  amrex::FArrayBox qxzm(txbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxzp(txbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxz = qxzm.array();
  auto const& qpxz = qxzp.array();

  amrex::ParallelFor(txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // X|Y
    pc_transy1(
      i, j, k, qmxy, qpxy, qxmarr, qxparr, fyarr, gdtempy, cdtdy,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
    // X|Z
    pc_transz1(
      i, j, k, qmxz, qpxz, qxmarr, qxparr, fzarr, gdtempz, cdtdz,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  const amrex::Box& txfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxxy(txfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox fluxxz(txfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox gdvxyfab(txfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvxzfab(txfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flxy = fluxxy.array();
  auto const& flxz = fluxxz.array();
  auto const& qxy = gdvxyfab.array();
  auto const& qxz = gdvxzfab.array();

  // Riemann problem X|Y X|Z
  amrex::ParallelFor(
    txfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // X|Y
      pc_cmpflx(
        i, j, k, qmxy, qpxy, flxy, qxy, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
      // X|Z
      pc_cmpflx(
        i, j, k, qmxz, qpxz, flxz, qxz, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // Y interface corrections
  cdir = 1;
  const amrex::Box& tybx = grow(bxg1, cdir, 1);
  const amrex::Box& tybxm = growHi(tybx, cdir, 1);
  amrex::FArrayBox qyxm(tybxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyxp(tybx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyzm(tybxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyzp(tybx, QVAR, amrex::The_Async_Arena());
  auto const& qmyx = qyxm.array();
  auto const& qpyx = qyxp.array();
  auto const& qmyz = qyzm.array();
  auto const& qpyz = qyzp.array();

  amrex::ParallelFor(tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|X
    pc_transx1(
      i, j, k, qmyx, qpyx, qymarr, qyparr, fxarr, gdtempx, cdtdx,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
    // Y|Z
    pc_transz2(
      i, j, k, qmyz, qpyz, qymarr, qyparr, fzarr, gdtempz, cdtdz,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  // Riemann problem Y|X Y|Z
  const amrex::Box& tyfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxyx(tyfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox fluxyz(tyfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox gdvyxfab(tyfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvyzfab(tyfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flyx = fluxyx.array();
  auto const& flyz = fluxyz.array();
  auto const& qyx = gdvyxfab.array();
  auto const& qyz = gdvyzfab.array();

  amrex::ParallelFor(
    tyfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Y|X
      pc_cmpflx(
        i, j, k, qmyx, qpyx, flyx, qyx, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
      // Y|Z
      pc_cmpflx(
        i, j, k, qmyz, qpyz, flyz, qyz, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // Z interface corrections
  cdir = 2;
  const amrex::Box& tzbx = grow(bxg1, cdir, 1);
  const amrex::Box& tzbxm = growHi(tzbx, cdir, 1);
  amrex::FArrayBox qzxm(tzbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzxp(tzbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzym(tzbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzyp(tzbx, QVAR, amrex::The_Async_Arena());

  auto const& qmzx = qzxm.array();
  auto const& qpzx = qzxp.array();
  auto const& qmzy = qzym.array();
  auto const& qpzy = qzyp.array();

  amrex::ParallelFor(tzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|X
    pc_transx2(
      i, j, k, qmzx, qpzx, qzmarr, qzparr, fxarr, gdtempx, cdtdx,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
    // Z|Y
    pc_transy2(
      i, j, k, qmzy, qpzy, qzmarr, qzparr, fyarr, gdtempy, cdtdy,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  // Riemann problem Z|X Z|Y
  const amrex::Box& tzfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxzx(tzfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox fluxzy(tzfxbx, flx1.nComp(), amrex::The_Async_Arena());
  amrex::FArrayBox gdvzxfab(tzfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvzyfab(tzfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flzx = fluxzx.array();
  auto const& flzy = fluxzy.array();
  auto const& qzx = gdvzxfab.array();
  auto const& qzy = gdvzyfab.array();

  amrex::ParallelFor(
    tzfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Z|X
      pc_cmpflx(
        i, j, k, qmzx, qpzx, flzx, qzx, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
      // Z|Y
      pc_cmpflx(
        i, j, k, qmzy, qpzy, flzy, qzy, 
        q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
        FirstSpec_comp_loc, NumSpec_loc,
#endif
        cdir);
    });

  // Temp Fabs for Final Fluxes
  amrex::FArrayBox qmfab(bxg2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qpfab(bxg1, QVAR, amrex::The_Async_Arena());
  auto const& qm = qmfab.array();
  auto const& qp = qpfab.array();

  // X | Y&Z
  cdir = 0;
  const amrex::Box& xfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& tyzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(tyzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transyz(
      i, j, k, qm, qp, qxmarr, qxparr, flyz, flzy, qyz, qzy,
      hdtdy, hdtdz,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  // Final X flux
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, qm, qp, flx1, q1, 
              q, small_dens, small_pres, small, gamma, 
#ifndef CONST_SPECIES
              FirstSpec_comp_loc, NumSpec_loc,
#endif
              cdir);
  });

  // Y | X&Z
  cdir = 1;
  const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(txzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxz(
      i, j, k, qm, qp, qymarr, qyparr, flxz, flzx, qxz, qzx,
      hdtdx, hdtdz,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  // Final Y flux
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, qm, qp, flx2, q2, 
              q, small_dens, small_pres, small, gamma,
#ifndef CONST_SPECIES
              FirstSpec_comp_loc, NumSpec_loc,
#endif
              cdir);
  });

  // Z | X&Y
  cdir = 2;
  const amrex::Box& zfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txybx = grow(bx, cdir, 1);
  amrex::ParallelFor(txybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxy(
      i, j, k, qm, qp, qzmarr, qzparr, flxy, flyx, qxy, qyx, hdtdx, hdtdy,
#ifndef CONST_SPECIES
      NumSpec,
#endif
      gamma, small_pres);
  });

  // Final Z flux
  amrex::ParallelFor(zfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, qm, qp, flx3, q3, 
              q, small_dens, small_pres, small, gamma,
#ifndef CONST_SPECIES
              FirstSpec_comp_loc, NumSpec_loc,
#endif
              cdir);
  });

  // Construct p div{U}
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_pdivu(i, j, k, pdivu, q1, q2, q3, dx, dy, dz);
  });
}

