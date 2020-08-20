#include "Godunov.H"
#include "PLM.H"
#include "PPM.H"

using namespace amrex;

extern void
trace_ppm(const Box& bx,
          const int idir,
          Array4<Real const> const& q_arr,
          Array4<Real const> const& srcQ,
          const int nq,
          Array4<Real> const& qm,
          Array4<Real> const& qp,
          const Box& vbx,
          const Real dt, const Real* del,
          const Real gamma,
          const Real small_dens, const Real small_pres, 
          const Real small_vel , const Real small,
          const int FirstSpec, const int NumSpec,
          const int use_flattening,
          const Real a_old);

// Host function to call gpu hydro functions
void
pc_umeth_3D(
  Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  Array4<const Real> const& q,
  const int nq,
  Array4<const Real> const& srcQ, 
  Array4<Real> const& flx1,
  Array4<Real> const& flx2,
  Array4<Real> const& flx3, 
  Array4<Real> const& q1,
  Array4<Real> const& q2,
  Array4<Real> const& q3,
  Array4<Real> const& pdivu,
  Array4<const Real> const& vol,
  const Real* del,
  const Real dt,
  const Real a_old,
  const Real a_new, 
  const int NumSpec,
  const Real gamma, const Real gamma_minus_1,
  const Real small_dens, const Real small_pres, 
  const Real small_vel , const Real small,
  const int ppm_type, const int use_flattening)
{
  const int FirstSpec_loc = FirstSpec;
  const int NumSpec_loc   = NumSpec;
  const int use_flattening_loc   = use_flattening;

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
  Real const hdt = 0.5 * dt / a_old;

  const int bclx = bclo[0];
  const int bcly = bclo[1];
  const int bclz = bclo[2];
  const int bchx = bchi[0];
  const int bchy = bchi[1];
  const int bchz = bchi[2];
  const int dlx = domlo[0];
  const int dly = domlo[1];
  const int dlz = domlo[2];
  const int dhx = domhi[0];
  const int dhy = domhi[1];
  const int dhz = domhi[2];

  // auto const& bcMaskarr = bcMask.array();
  const amrex::Box& bxg1 = grow(bx, 1);
  const amrex::Box& bxg2 = grow(bx, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qxm(xmbx, nq);
  amrex::FArrayBox qxp(bxg2, nq);
  amrex::Elixir qxmeli = qxm.elixir();
  amrex::Elixir qxpeli = qxp.elixir();
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qym(ymbx, nq);
  amrex::FArrayBox qyp(bxg2, nq);
  amrex::Elixir qymeli = qym.elixir();
  amrex::Elixir qypeli = qyp.elixir();
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  // Z data
  cdir = 2;
  const amrex::Box& zmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& zflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qzm(zmbx, nq);
  amrex::FArrayBox qzp(bxg2, nq);
  amrex::Elixir qzmeli = qzm.elixir();
  amrex::Elixir qzpeli = qzp.elixir();
  auto const& qzmarr = qzm.array();
  auto const& qzparr = qzp.array();

  // Put the PLM and slopes in the same kernel launch to avoid unnecessary
  // launch overhead Pelec_Slope_* are SIMD as well as PeleC_plm_* which loop
  // over the same box
  
  if(ppm_type == 0 )
  {
      amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real slope[8];

      const amrex::Real c = std::sqrt((gamma_minus_1+1.0) * q(i,j,k,QREINT)/q(i,j,k,QRHO)*gamma_minus_1);

      // X slopes and interp
      for (int n = 0; n < nq; ++n)
        slope[n] = plm_slope(i, j, k, n, 0, use_flattening_loc, q);
      pc_plm_x(i, j, k, qxmarr, qxparr, srcQ, slope, q, c, a_old, dx, dt, NumSpec, 
               gamma_minus_1, small_dens, small_pres);

      // Y slopes and interp
      for (int n = 0; n < nq; n++)
        slope[n] = plm_slope(i, j, k, n, 1, use_flattening_loc, q);
      pc_plm_y(i, j, k, qymarr, qyparr, srcQ, slope, q, c, a_old, dy, dt, NumSpec,
               gamma_minus_1, small_dens, small_pres);

      // Z slopes and interp
      for (int n = 0; n < nq; ++n)
        slope[n] = plm_slope(i, j, k, n, 2, use_flattening_loc, q);
      pc_plm_z(i, j, k, qzmarr, qzparr, srcQ, slope, q, c, a_old, dz, dt, NumSpec,
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
                q, srcQ, nq,
                qxmarr, qxparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
                small_vel, small,
                FirstSpec, NumSpec,
                use_flattening, a_old);

      idir = 1;
      trace_ppm(bxg2,
                idir,
                q, srcQ, nq,
                qymarr, qyparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
                small_vel, small,
                FirstSpec, NumSpec,
                use_flattening, a_old);

      idir = 2;
      trace_ppm(bxg2,
                idir,
                q, srcQ, nq,
                qzmarr, qzparr,
                bxg2, dt, del, gamma,
                small_dens, small_pres,
                small_vel, small,
                FirstSpec, NumSpec,
                use_flattening, a_old);

  }


  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  cdir = 0;
  amrex::FArrayBox fx(xflxbx, flx1.nComp());
  amrex::Elixir fxeli = fx.elixir();
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(xflxbx, NGDNV);
  amrex::Elixir qgdxeli = qgdx.elixir();
  auto const& gdtempx = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtempx, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  // Y initial fluxes
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, flx2.nComp());
  amrex::Elixir fyeli = fy.elixir();
  auto const& fyarr = fy.array();
  amrex::FArrayBox qgdy(yflxbx, NGDNV);
  amrex::Elixir qgdyeli = qgdy.elixir();
  auto const& gdtempy = qgdy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, gdtempy, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  // Z initial fluxes
  cdir = 2;
  amrex::FArrayBox fz(zflxbx, flx3.nComp());
  amrex::Elixir fzeli = fz.elixir();
  auto const& fzarr = fz.array();
  amrex::FArrayBox qgdz(zflxbx, NGDNV);
  amrex::Elixir qgdzeli = qgdz.elixir();
  auto const& gdtempz = qgdz.array();
  amrex::ParallelFor(
    zflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qzmarr, qzparr, fzarr, gdtempz, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  // X interface corrections
  cdir = 0;
  const amrex::Box& txbx = grow(bxg1, cdir, 1);
  const amrex::Box& txbxm = growHi(txbx, cdir, 1);
  amrex::FArrayBox qxym(txbxm, nq);
  amrex::Elixir qxymeli = qxym.elixir();
  amrex::FArrayBox qxyp(txbx, nq);
  amrex::Elixir qxypeli = qxyp.elixir();
  auto const& qmxy = qxym.array();
  auto const& qpxy = qxyp.array();

  amrex::FArrayBox qxzm(txbxm, nq);
  amrex::Elixir qxzmeli = qxzm.elixir();
  amrex::FArrayBox qxzp(txbx, nq);
  amrex::Elixir qxzpeli = qxzp.elixir();
  auto const& qmxz = qxzm.array();
  auto const& qpxz = qxzp.array();

  amrex::ParallelFor(txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // X|Y
    pc_transy1(
      i, j, k, qmxy, qpxy, qxmarr, qxparr, fyarr, gdtempy, cdtdy, NumSpec, gamma, small_pres);
    // X|Z
    pc_transz1(
      i, j, k, qmxz, qpxz, qxmarr, qxparr, fzarr, gdtempz, cdtdz, NumSpec, gamma, small_pres);
  });

  const amrex::Box& txfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxxy(txfxbx, flx1.nComp());
  amrex::FArrayBox fluxxz(txfxbx, flx1.nComp());
  amrex::FArrayBox gdvxyfab(txfxbx, NGDNV);
  amrex::FArrayBox gdvxzfab(txfxbx, NGDNV);
  amrex::Elixir fluxxyeli = fluxxy.elixir(), gdvxyeli = gdvxyfab.elixir();
  amrex::Elixir fluxxzeli = fluxxz.elixir(), gdvxzeli = gdvxzfab.elixir();

  auto const& flxy = fluxxy.array();
  auto const& flxz = fluxxz.array();
  auto const& qxy = gdvxyfab.array();
  auto const& qxz = gdvxzfab.array();

  // Riemann problem X|Y X|Z
  amrex::ParallelFor(
    txfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // X|Y
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qmxy, qpxy, flxy, qxy, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
      // X|Z
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qmxz, qpxz, flxz, qxz, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  qxymeli.clear();
  qxypeli.clear();
  qxzmeli.clear();
  qxzpeli.clear();

  // Y interface corrections
  cdir = 1;
  const amrex::Box& tybx = grow(bxg1, cdir, 1);
  const amrex::Box& tybxm = growHi(tybx, cdir, 1);
  amrex::FArrayBox qyxm(tybxm, nq);
  amrex::FArrayBox qyxp(tybx, nq);
  amrex::FArrayBox qyzm(tybxm, nq);
  amrex::FArrayBox qyzp(tybx, nq);
  amrex::Elixir qyxmeli = qyxm.elixir(), qyxpeli = qyxp.elixir();
  amrex::Elixir qyzmeli = qyzm.elixir(), qyzpeli = qyzp.elixir();
  auto const& qmyx = qyxm.array();
  auto const& qpyx = qyxp.array();
  auto const& qmyz = qyzm.array();
  auto const& qpyz = qyzp.array();

  amrex::ParallelFor(tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|X
    pc_transx1(
      i, j, k, qmyx, qpyx, qymarr, qyparr, fxarr, gdtempx, cdtdx, NumSpec, gamma, small_pres);
    // Y|Z
    pc_transz2(
      i, j, k, qmyz, qpyz, qymarr, qyparr, fzarr, gdtempz, cdtdz, NumSpec, gamma, small_pres);
  });

  fzeli.clear();
  qgdzeli.clear();

  // Riemann problem Y|X Y|Z
  const amrex::Box& tyfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxyx(tyfxbx, flx1.nComp());
  amrex::FArrayBox fluxyz(tyfxbx, flx1.nComp());
  amrex::FArrayBox gdvyxfab(tyfxbx, NGDNV);
  amrex::FArrayBox gdvyzfab(tyfxbx, NGDNV);
  amrex::Elixir fluxyxeli = fluxyx.elixir(), gdvyxeli = gdvyxfab.elixir();
  amrex::Elixir fluxyzeli = fluxyz.elixir(), gdvyzeli = gdvyzfab.elixir();

  auto const& flyx = fluxyx.array();
  auto const& flyz = fluxyz.array();
  auto const& qyx = gdvyxfab.array();
  auto const& qyz = gdvyzfab.array();

  amrex::ParallelFor(
    tyfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Y|X
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qmyx, qpyx, flyx, qyx, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
      // Y|Z
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qmyz, qpyz, flyz, qyz, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  qyxmeli.clear();
  qyxpeli.clear();
  qyzmeli.clear();
  qyzpeli.clear();

  // Z interface corrections
  cdir = 2;
  const amrex::Box& tzbx = grow(bxg1, cdir, 1);
  const amrex::Box& tzbxm = growHi(tzbx, cdir, 1);
  amrex::FArrayBox qzxm(tzbxm, nq);
  amrex::FArrayBox qzxp(tzbx, nq);
  amrex::FArrayBox qzym(tzbxm, nq);
  amrex::FArrayBox qzyp(tzbx, nq);
  amrex::Elixir qzxmeli = qzxm.elixir(), qzxpeli = qzxp.elixir();
  amrex::Elixir qzymeli = qzym.elixir(), qzypeli = qzyp.elixir();

  auto const& qmzx = qzxm.array();
  auto const& qpzx = qzxp.array();
  auto const& qmzy = qzym.array();
  auto const& qpzy = qzyp.array();

  amrex::ParallelFor(tzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|X
    pc_transx2(
      i, j, k, qmzx, qpzx, qzmarr, qzparr, fxarr, gdtempx, cdtdx, NumSpec, gamma, small_pres);
    // Z|Y
    pc_transy2(
      i, j, k, qmzy, qpzy, qzmarr, qzparr, fyarr, gdtempy, cdtdy, NumSpec, gamma, small_pres);
  });

  fxeli.clear();
  fyeli.clear();
  qgdxeli.clear();
  qgdyeli.clear();

  // Riemann problem Z|X Z|Y
  const amrex::Box& tzfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxzx(tzfxbx, flx1.nComp());
  amrex::FArrayBox fluxzy(tzfxbx, flx1.nComp());
  amrex::FArrayBox gdvzxfab(tzfxbx, NGDNV);
  amrex::FArrayBox gdvzyfab(tzfxbx, NGDNV);
  amrex::Elixir fluxzxeli = fluxzx.elixir(), gdvzxeli = gdvzxfab.elixir();
  amrex::Elixir fluxzyeli = fluxzy.elixir(), gdvzyeli = gdvzyfab.elixir();

  auto const& flzx = fluxzx.array();
  auto const& flzy = fluxzy.array();
  auto const& qzx = gdvzxfab.array();
  auto const& qzy = gdvzyfab.array();

  amrex::ParallelFor(
    tzfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Z|X
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qmzx, qpzx, flzx, qzx, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
      // Z|Y
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qmzy, qpzy, flzy, qzy, 
        q, small_dens, small_pres, small_vel, small, gamma, 
        FirstSpec_loc, NumSpec_loc, cdir);
    });

  qzxmeli.clear();
  qzxpeli.clear();
  qzymeli.clear();
  qzypeli.clear();

  // Temp Fabs for Final Fluxes
  amrex::FArrayBox qmfab(bxg2, nq);
  amrex::FArrayBox qpfab(bxg1, nq);
  amrex::Elixir qmeli = qmfab.elixir();
  amrex::Elixir qpeli = qpfab.elixir();
  auto const& qm = qmfab.array();
  auto const& qp = qpfab.array();

  // X | Y&Z
  cdir = 0;
  const amrex::Box& xfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& tyzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(tyzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transyz(
      i, j, k, qm, qp, qxmarr, qxparr, flyz, flzy, qyz, qzy, hdt,
      hdtdy, hdtdz, NumSpec, gamma, small_pres);
  });

  fluxzyeli.clear();
  gdvzyeli.clear();
  gdvyzeli.clear();
  fluxyzeli.clear();
  qxmeli.clear();
  qxpeli.clear();
  // Final X flux
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bclx, bchx, dlx, dhx, qm, qp, flx1, q1, 
              q, small_dens, small_pres, small_vel, small, gamma, 
              FirstSpec_loc, NumSpec_loc, cdir);
  });

  // Y | X&Z
  cdir = 1;
  const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(txzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxz(
      i, j, k, qm, qp, qymarr, qyparr, flxz, flzx, qxz, qzx, hdt,
      hdtdx, hdtdz, NumSpec, gamma, small_pres);
  });

  fluxzxeli.clear();
  gdvzxeli.clear();
  gdvxzeli.clear();
  fluxxzeli.clear();
  qymeli.clear();
  qypeli.clear();
  // Final Y flux
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bcly, bchy, dly, dhy, qm, qp, flx2, q2, 
              q, small_dens, small_pres, small_vel, small, gamma,
              FirstSpec_loc, NumSpec_loc, cdir);
  });

  // Z | X&Y
  cdir = 2;
  const amrex::Box& zfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txybx = grow(bx, cdir, 1);
  amrex::ParallelFor(txybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxy(
      i, j, k, qm, qp, qzmarr, qzparr, flxy, flyx, qxy, qyx, hdt,
      hdtdx, hdtdy, NumSpec, gamma, small_pres);
  });

  gdvyxeli.clear();
  fluxyxeli.clear();
  gdvxyeli.clear();
  fluxxyeli.clear();
  qzmeli.clear();
  qzpeli.clear();
  // Final Z flux
  amrex::ParallelFor(zfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bclz, bchz, dlz, dhz, qm, qp, flx3, q3, 
              q, small_dens, small_pres, small_vel, small, gamma,
              FirstSpec_loc, NumSpec_loc, cdir);
  });

  qmeli.clear();
  qpeli.clear();
  // Construct p div{U}
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_pdivu(i, j, k, pdivu, q1, q2, q3, dx, dy, dz);
  });
}

