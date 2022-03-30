#ifndef _Godunov_H_
#define _Godunov_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include <IndexDefines.H>
#include <Riemann.H>

/*This header file contains functions and declarations for the 3D hydro portion
  of PeleC cuda. As per the convention of AMReX, inlined functions are defined
  here. Where as non-inline functions are declared here. */

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
flatten(
  const int i,
  const int j,
  const int k,
  const int dir,
  amrex::Array4<const amrex::Real> const& q,
  const amrex::Real small_pres)
{
  const int bdim[3] = {dir == 0, dir == 1, dir == 2};
  if((!q.contains(i+3 * bdim[0], j+3 * bdim[1], k+3 * bdim[2])||
      !q.contains(i-3 * bdim[0], j-3 * bdim[1], k-3 * bdim[2])))
  {
      return 1.0;
  }
  const int n = QPRES;
  //Parameters from uflatten
  const amrex::Real shktst = 0.33e0;
  const amrex::Real zcut1 = 0.75e0;
  const amrex::Real zcut2 = 0.85e0;
  const amrex::Real dzcut = 1.0/(zcut2-zcut1);

  const amrex::Real qm =
    q(i - 1 * bdim[0], j - 1 * bdim[1], k - 1 * bdim[2], n);
  const amrex::Real qp =
    q(i + 1 * bdim[0], j + 1 * bdim[1], k + 1 * bdim[2], n);
  //assumes QU,QV,QW ordering
  const amrex::Real qmun =
    q(i - 1 * bdim[0], j - 1 * bdim[1], k - 1 * bdim[2], QU+dir);
  const amrex::Real qpun =
    q(i + 1 * bdim[0], j + 1 * bdim[1], k + 1 * bdim[2], QU+dir);
  const amrex::Real qmm =
    q(i - 2 * bdim[0], j - 2 * bdim[1], k - 2 * bdim[2], n);
  const amrex::Real qpp =
    q(i + 2 * bdim[0], j + 2 * bdim[1], k + 2 * bdim[2], n);

  amrex::Real dp = qp - qm;

  const int ishft = (dp > 0.0) ? 1 : -1;

  amrex::Real denom = amrex::max(small_pres, amrex::Math::abs(qpp-qmm));
  amrex::Real zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real    z = amrex::min(amrex::Real(1.0), amrex::max(amrex::Real(0.0), dzcut * (zeta - zcut1)));

  amrex::Real tst = (qmun-qpun >= 0.0) ? 1.0 : 0.0;

  amrex::Real tmp = amrex::min(qp,qm);

  const amrex::Real chi = (amrex::Math::abs(dp)/tmp) > shktst ? tst : 0.0;

  const amrex::Real qpshft =
          q(i + (1-ishft) * bdim[0], j + (1-ishft) * bdim[1], k + (1-ishft) * bdim[2], n);
  const amrex::Real qmshft =
          q(i - (1+ishft) * bdim[0], j - (1+ishft) * bdim[1], k - (1+ishft) * bdim[2], n);
  const amrex::Real qpunshft =
          q(i + (1-ishft) * bdim[0], j + (1-ishft) * bdim[1], k + (1-ishft) * bdim[2], QU+dir);
  const amrex::Real qmunshft =
          q(i - (1+ishft) * bdim[0], j - (1+ishft) * bdim[1], k - (1+ishft) * bdim[2], QU+dir);
  const amrex::Real qppshft =
          q(i + (2-ishft) * bdim[0], j + (2-ishft) * bdim[1], k + (2-ishft) * bdim[2], n);
  const amrex::Real qmmshft =
          q(i - (2+ishft) * bdim[0], j - (2+ishft) * bdim[1], k - (2+ishft) * bdim[2], n);

  dp = qpshft - qmshft;

  denom = amrex::max(small_pres, amrex::Math::abs(qppshft-qmmshft));
  zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real z2 = amrex::min(amrex::Real(1.0), amrex::max(amrex::Real(0.0), dzcut * (zeta - zcut1)));

  tst =  (qmunshft - qpunshft >= 0.0) ? 1.0 : 0.0;

  tmp = amrex::min(qpshft, qmshft);

  const amrex::Real chi2 = (amrex::Math::abs(dp)/tmp) > shktst ? tst : 0.0;

  return 1.0 - amrex::max(chi2 * z2, chi * z);

}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpflx(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& ql,
  amrex::Array4<const amrex::Real> const& qr,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<amrex::Real> const& q,
  // amrex::Array4<const int> const& bcMask,
  amrex::Array4<const amrex::Real> const& qprim,
  const amrex::Real small_dens,
  const amrex::Real small_pres,
  const amrex::Real small,
  const amrex::Real gamma,
#ifndef CONST_SPECIES
  const int FirstSpec_comp, const int NumSpec,
#endif
  const int dir)
{
  amrex::Real csmall, cav, ustar;
  amrex::Real ul, ur, vl, vr, v2l, v2r, rel, rer;
  amrex::Real rl, rr, pl, pr;
  int IU, IV, IV2;
  int GU, GV, GV2;
  int f_idx[3];

  if (dir == 0) {
    IU = QU;
    IV = QV;
    IV2 = QW;
    GU = QU;
    GV = QV;
    GV2 = QW;
    csmall = amrex::max(small, small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i-1,j,k,QREINT) / qprim(i-1,j,k,QRHO))), small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k  ,QREINT) / qprim(i,j,k  ,QRHO))) );
    cav = 0.5 * ( std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i-1,j,k,QREINT) / qprim(i-1,j,k,QRHO)))
                  +std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i  ,j,k,QREINT) / qprim(i  ,j,k,QRHO))) );
    f_idx[0] = Xmom_comp;
    f_idx[1] = Ymom_comp;
    f_idx[2] = Zmom_comp;
  } else if (dir == 1) {
    IU = QV;
    IV = QU;
    IV2 = QW;
    GU = QV;
    GV = QU;
    GV2 = QW;
    csmall = amrex::max(small, small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j-1,k,QREINT) / qprim(i,j-1,k,QRHO))), small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k  ,QREINT) / qprim(i,j,k  ,QRHO))) );
    cav = 0.5 * ( std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j-1,k,QREINT) / qprim(i,j-1,k,QRHO)))
                                  +std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j  ,k,QREINT) / qprim(i,j  ,k,QRHO))) );
    f_idx[0] = Ymom_comp;
    f_idx[1] = Xmom_comp;
    f_idx[2] = Zmom_comp;
  } else {
    IU = QW;
    IV = QU;
    IV2 = QV;
    GU = QW;
    GV = QU;
    GV2 = QV;
    csmall = amrex::max(small, small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k-1,QREINT) / qprim(i,j,k-1,QRHO))), small * std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k  ,QREINT) / qprim(i,j,k  ,QRHO))) );
    cav = 0.5 * ( std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k-1,QREINT) / qprim(i,j,k-1,QRHO)))
                                  +std::sqrt(amrex::Math::abs(gamma * (gamma - amrex::Real(1.0)) * qprim(i,j,k  ,QREINT) / qprim(i,j,k  ,QRHO))) );
    f_idx[0] = Zmom_comp;
    f_idx[1] = Xmom_comp;
    f_idx[2] = Ymom_comp;
  }

  rl = ql(i,j,k,QRHO);
  ul = ql(i, j, k, IU);
  vl = ql(i, j, k, IV);
  v2l = ql(i, j, k, IV2);
  rel = ql(i, j, k, QREINT);
  pl = ql(i,j,k,QPRES);

  rr = qr(i,j,k,QRHO);
  ur = qr(i, j, k, IU);
  vr = qr(i, j, k, IV);
  v2r = qr(i, j, k, IV2);
  rer = qr(i, j, k, QREINT);
  pr = qr(i,j,k,QPRES);

  riemann(
    rl, ul, vl, v2l, pl, rel,
    rr, ur, vr, v2r, pr, rer,
    cav, ustar, flx(i, j, k, Density_comp),
    flx(i, j, k, f_idx[0]), flx(i, j, k, f_idx[1]), flx(i, j, k, f_idx[2]),
    flx(i, j, k, Eden_comp), flx(i, j, k, Eint_comp), q(i, j, k, GU), q(i, j, k, GV),
    q(i, j, k, GV2), q(i, j, k, QPRES), q(i, j, k, QREINT), 
    small_dens, small_pres, small, csmall, gamma);

#ifndef CONST_SPECIES
  amrex::Real flxrho = flx(i, j, k, Density_comp);
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    if (ustar > 0.0)
      flx(i, j, k, n) = flxrho * ql(i, j, k, n);
    else if (ustar < 0.0)
      flx(i, j, k, n) = flxrho * qr(i, j, k, n);
    else
      flx(i, j, k, n) = flxrho * 0.50 * (ql(i, j, k, n) + qr(i, j, k, n));
  }
#endif
}

// First Transverse Correction for Predicted Y-states, using X-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transx1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdx, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real dup, pav, du;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, Density_comp) - flxx(i, j, k, Density_comp));
  flxu = cdtdx * (flxx(i + 1, j, k, Xmom_comp) - flxx(i, j, k, Xmom_comp));
  flxv = cdtdx * (flxx(i + 1, j, k, Ymom_comp) - flxx(i, j, k, Ymom_comp));
  flxw = cdtdx * (flxx(i + 1, j, k, Zmom_comp) - flxx(i, j, k, Zmom_comp));
  flxe = cdtdx * (flxx(i + 1, j, k, Eden_comp) - flxx(i, j, k, Eden_comp));

  // Update passive variables
#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, n) - compn;
    qm(i, j + 1, k, n) = compo / rrnew;
  }
#endif

  // Update hydro vars
  pggp = qint(i + 1, j, k, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i + 1, j, k, QU);
  ugm = qint(i, j, k, QU);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qyp(i, j, k, QRHO);
  rur = qyp(i, j, k, QU);
  rvr = qyp(i, j, k, QV);
  rwr = qyp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qyp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  qp(i, j, k, QREINT) = renewr - rhoekinr;

  amrex::Real pnewr = qp(i, j, k, QREINT) * (gamma - amrex::Real(1.0));

  if (pnewr < small_pres)
  {
       pnewr = qyp(i,j,k,QPRES) - cdtdx * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qym(i, j + 1, k, QRHO);
  rul = qym(i, j + 1, k, QU);
  rvl = qym(i, j + 1, k, QV);
  rwl = qym(i, j + 1, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl = qm(i, j + 1, k, QREINT) * (gamma - amrex::Real(1.0));
  if(pnewl < small_pres)
  {
       pnewl = qym(i,j+1,k,QPRES) - cdtdx * ( dup + pav*du*(gamma - 1.));
  }
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, small_pres);
  qm(i, j + 1, k, QREINT) = qm(i, j + 1, k, QPRES) / (gamma - amrex::Real(1.0));
}

// Corrects the left and right x states with y flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transy1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdy, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{

  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, Density_comp) - flxy(i, j, k, Density_comp));
  flxu = cdtdy * (flxy(i, j + 1, k, Xmom_comp) - flxy(i, j, k, Xmom_comp));
  flxv = cdtdy * (flxy(i, j + 1, k, Ymom_comp) - flxy(i, j, k, Ymom_comp));
  flxw = cdtdy * (flxy(i, j + 1, k, Zmom_comp) - flxy(i, j, k, Zmom_comp));
  flxe = cdtdy * (flxy(i, j + 1, k, Eden_comp) - flxy(i, j, k, Eden_comp));

  // -  Update Passive Variables
#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, n) - compn;
    qm(i + 1, j, k, n) = compo / rrnew;
  }
#endif

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j + 1, k, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i, j + 1, k, QV);
  ugm = qint(i, j, k, QV);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qxp(i, j, k, QRHO);
  rur = qxp(i, j, k, QU);
  rvr = qxp(i, j, k, QV);
  rwr = qxp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qxp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  qp(i, j, k, QREINT) = renewr - rhoekinr;
  pnewr = qp(i, j, k, QREINT) * (gamma-1.0);
  if(pnewr < small_pres)
  {
       pnewr = qxp(i,j,k,QPRES) - cdtdy * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qxm(i + 1, j, k, QRHO);
  rul = qxm(i + 1, j, k, QU);
  rvl = qxm(i + 1, j, k, QV);
  rwl = qxm(i + 1, j, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
  pnewl = qm(i + 1, j, k, QREINT) * (gamma-1.0);
  if(pnewl < small_pres)
  {
       pnewl = qxm(i+1,j,k,QPRES) - cdtdy * ( dup + pav*du*(gamma - 1.));
  }
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, small_pres);
  qm(i + 1, j, k, QREINT) = qm(i + 1, j, k, QPRES) / (gamma - amrex::Real(1.0));

}

/* Transverse Correction for Predicted Z-states, using X-Flux */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transx2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdx, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dup, pav, du;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, Density_comp) - flxx(i, j, k, Density_comp));
  flxu = cdtdx * (flxx(i + 1, j, k, Xmom_comp) - flxx(i, j, k, Xmom_comp));
  flxv = cdtdx * (flxx(i + 1, j, k, Ymom_comp) - flxx(i, j, k, Ymom_comp));
  flxw = cdtdx * (flxx(i + 1, j, k, Zmom_comp) - flxx(i, j, k, Zmom_comp));
  flxe = cdtdx * (flxx(i + 1, j, k, Eden_comp) - flxx(i, j, k, Eden_comp));

  // -  Update Passive Variables
#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, n) - compn;
    qm(i, j, k + 1, n) = compo / rrnew;
  }
#endif
  //- Update hydro vars -

  pggp = qint(i + 1, j, k, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i + 1, j, k, QU);
  ugm = qint(i, j, k, QU);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qzp(i, j, k, QRHO);
  rur = qzp(i, j, k, QU);
  rvr = qzp(i, j, k, QV);
  rwr = qzp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qzp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real pnewr = qp(i, j, k, QREINT) * (gamma-1.0);
  if(pnewr < small_pres)
  {
       pnewr = qzp(i,j,k,QPRES) - cdtdx * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qzm(i, j, k + 1, QRHO);
  rul = qzm(i, j, k + 1, QU);
  rvl = qzm(i, j, k + 1, QV);
  rwl = qzm(i, j, k + 1, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl = qm(i, j, k + 1, QREINT) * (gamma-1.0);
  if(pnewl < small_pres)
  {
       pnewl = qzm(i,j,k+1,QPRES) - cdtdx * ( dup + pav*du*(gamma - 1.));
  }
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, small_pres);
  qm(i, j, k + 1, QREINT) = qm(i, j, k + 1, QPRES) / (gamma - amrex::Real(1.0));

}

// Corrects the left and right z states with y flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transy2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdy,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dup, pav, du;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, Density_comp) - flxy(i, j, k, Density_comp));
  flxu = cdtdy * (flxy(i, j + 1, k, Xmom_comp) - flxy(i, j, k, Xmom_comp));
  flxv = cdtdy * (flxy(i, j + 1, k, Ymom_comp) - flxy(i, j, k, Ymom_comp));
  flxw = cdtdy * (flxy(i, j + 1, k, Zmom_comp) - flxy(i, j, k, Zmom_comp));
  flxe = cdtdy * (flxy(i, j + 1, k, Eden_comp) - flxy(i, j, k, Eden_comp));

  // -  Update Passive Variables
#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, n) - compn;
    qm(i, j, k + 1, n) = compo / rrnew;
  }
#endif

  //- Update hydro vars -
  pggp = qint(i, j + 1, k, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i, j + 1, k, QV);
  ugm = qint(i, j, k, QV);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qzp(i, j, k, QRHO);
  rur = qzp(i, j, k, QU);
  rvr = qzp(i, j, k, QV);
  rwr = qzp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qzp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real pnewr =
    qp(i, j, k, QREINT) * (gamma-1.0);
  if(pnewr < small_pres)
  {
       pnewr = qzp(i,j,k,QPRES) - cdtdy * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qzm(i, j, k + 1, QRHO);
  rul = qzm(i, j, k + 1, QU);
  rvl = qzm(i, j, k + 1, QV);
  rwl = qzm(i, j, k + 1, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qm(i, j, k + 1, QREINT) * (gamma-1.0);
  if(pnewl < small_pres)
  {
       pnewl = qzm(i,j,k+1,QPRES) - cdtdy * ( dup + pav*du*(gamma - 1.));
  }
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, small_pres);
  qm(i, j, k + 1, QREINT) = qm(i, j, k + 1, QPRES) / (gamma - amrex::Real(1.0));

}

// Corrects the left and right x states with z flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transz1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdz,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, Density_comp) - flxz(i, j, k, Density_comp));
  flxu = cdtdz * (flxz(i, j, k + 1, Xmom_comp) - flxz(i, j, k, Xmom_comp));
  flxv = cdtdz * (flxz(i, j, k + 1, Ymom_comp) - flxz(i, j, k, Ymom_comp));
  flxw = cdtdz * (flxz(i, j, k + 1, Zmom_comp) - flxz(i, j, k, Zmom_comp));
  flxe = cdtdz * (flxz(i, j, k + 1, Eden_comp) - flxz(i, j, k, Eden_comp));

#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, n) - compn;
    qm(i + 1, j, k, n) = compo / rrnew;
  }
#endif

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j, k + 1, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i, j, k + 1, QW);
  ugm = qint(i, j, k, QW);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qxp(i, j, k, QRHO);
  rur = qxp(i, j, k, QU);
  rvr = qxp(i, j, k, QV);
  rwr = qxp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qxp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  qp(i, j, k, QREINT) = renewr - rhoekinr;
  pnewr = qp(i, j, k, QREINT) * (gamma-1.0);
  if(pnewr < small_pres)
  {
       pnewr = qxp(i,j,k,QPRES) - cdtdz * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qxm(i + 1, j, k, QRHO);
  rul = qxm(i + 1, j, k, QU);
  rvl = qxm(i + 1, j, k, QV);
  rwl = qxm(i + 1, j, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
  pnewl = qm(i + 1, j, k, QREINT) * (gamma-1.0);
  if(pnewl < small_pres)
  {
       pnewl = qxm(i+1,j,k,QPRES) - cdtdz * ( dup + pav*du*(gamma - 1.));
  }
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, small_pres);
  qm(i + 1, j, k, QREINT) = qm(i + 1, j, k, QPRES) / (gamma - amrex::Real(1.0));

}

// Transverse Correction for Predicted Y-states, using Z-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transz2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdz,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dup, pav, du;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, Density_comp) - flxz(i, j, k, Density_comp));
  flxu = cdtdz * (flxz(i, j, k + 1, Xmom_comp) - flxz(i, j, k, Xmom_comp));
  flxv = cdtdz * (flxz(i, j, k + 1, Ymom_comp) - flxz(i, j, k, Ymom_comp));
  flxw = cdtdz * (flxz(i, j, k + 1, Zmom_comp) - flxz(i, j, k, Zmom_comp));
  flxe = cdtdz * (flxz(i, j, k + 1, Eden_comp) - flxz(i, j, k, Eden_comp));

#ifndef CONST_SPECIES
  amrex::Real compn, rr, rrnew, compo;
  // Update passive variables
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, n) - compn;
    qm(i, j + 1, k, n) = compo / rrnew;
  }
#endif

  // Update hydro vars
  pggp = qint(i, j, k + 1, QPRES);
  pggm = qint(i, j, k, QPRES);
  ugp = qint(i, j, k + 1, QW);
  ugm = qint(i, j, k, QW);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qyp(i, j, k, QRHO);
  rur = qyp(i, j, k, QU);
  rvr = qyp(i, j, k, QV);
  rwr = qyp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;
  rer = qyp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real pnewr =
    qp(i, j, k, QREINT) * (gamma-1.0);
  if(pnewr < small_pres)
  {
       pnewr = qyp(i,j,k,QPRES) - cdtdz * ( dup + pav*du*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(pnewr, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM

  // Conversion to Conservative
  rrl = qym(i, j + 1, k, QRHO);
  rul = qym(i, j + 1, k, QU);
  rvl = qym(i, j + 1, k, QV);
  rwl = qym(i, j + 1, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qm(i, j + 1, k, QREINT) * (gamma-1.0);
  if(pnewl < small_pres)
  {
       pnewl = qym(i,j+1,k,QPRES) - cdtdz * ( dup + pav*du*(gamma - 1.));
  }
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, small_pres);
  qm(i, j + 1, k, QREINT) = qm(i, j + 1, k, QPRES) / (gamma - amrex::Real(1.0));
}

// Z corrected from xy
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transxy(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qx,
  amrex::Array4<const amrex::Real> const& qy,
  const amrex::Real cdtdx,
  const amrex::Real cdtdy, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, Density_comp) - flxx(i, j, k, Density_comp)) +
           cdtdy * (flxy(i, j + 1, k, Density_comp) - flxy(i, j, k, Density_comp));
  flxu = cdtdx * (flxx(i + 1, j, k, Xmom_comp) - flxx(i, j, k, Xmom_comp)) +
         cdtdy * (flxy(i, j + 1, k, Xmom_comp) - flxy(i, j, k, Xmom_comp));
  flxv = cdtdx * (flxx(i + 1, j, k, Ymom_comp) - flxx(i, j, k, Ymom_comp)) +
         cdtdy * (flxy(i, j + 1, k, Ymom_comp) - flxy(i, j, k, Ymom_comp));
  flxw = cdtdx * (flxx(i + 1, j, k, Zmom_comp) - flxx(i, j, k, Zmom_comp)) +
         cdtdy * (flxy(i, j + 1, k, Zmom_comp) - flxy(i, j, k, Zmom_comp));
  flxe = cdtdx * (flxx(i + 1, j, k, Eden_comp) - flxx(i, j, k, Eden_comp)) +
         cdtdy * (flxy(i, j + 1, k, Eden_comp) - flxy(i, j, k, Eden_comp));

  // Update passive variables
  rrr = qzp(i, j, k, QRHO);
  rrl = qzm(i, j, k + 1, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;

#ifndef CONST_SPECIES
  amrex::Real compn, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    // qp
    compo = rrr * qzp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnewr;
    // qm
    compo = rrl * qzm(i, j, k + 1, n) - compn;
    qm(i, j, k + 1, n) = compo / rrnewl;
  }
#endif

  // Update hydro vars
  pggpx = qx(i + 1, j, k, QPRES);
  pggmx = qx(i, j, k, QPRES);
  ugpx = qx(i + 1, j, k, QU);
  ugmx = qx(i, j, k, QU);

  pggpy = qy(i, j + 1, k, QPRES);
  pggmy = qy(i, j, k, QPRES);
  ugpy = qy(i, j + 1, k, QV);
  ugmy = qy(i, j, k, QV);

  // qp state
  rur = rrr * qzp(i, j, k, QU);
  rvr = rrr * qzp(i, j, k, QV);
  rwr = rrr * qzp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qzp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupx = pggpx * ugpx - pggmx * ugmx;
       amrex::Real pavx = 0.5 * (pggpx + pggmx);
       amrex::Real dux = ugpx - ugmx;

       amrex::Real dupy = pggpy * ugpy - pggmy * ugmy;
       amrex::Real pavy = 0.5 * (pggpy + pggmy);
       amrex::Real duy = ugpy - ugmy;

       temppres = qzp(i,j,k,QPRES) - cdtdx * ( dupx + pavx*dux*(gamma - 1.))
                                   - cdtdy * ( dupy + pavy*duy*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(temppres, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // qm state
  rul = rrl * qzm(i, j, k + 1, QU);
  rvl = rrl * qzm(i, j, k + 1, QV);
  rwl = rrl * qzm(i, j, k + 1, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  temppres = qm(i, j, k + 1, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupx = pggpx * ugpx - pggmx * ugmx;
       amrex::Real pavx = 0.5 * (pggpx + pggmx);
       amrex::Real dux = ugpx - ugmx;

       amrex::Real dupy = pggpy * ugpy - pggmy * ugmy;
       amrex::Real pavy = 0.5 * (pggpy + pggmy);
       amrex::Real duy = ugpy - ugmy;

       temppres = qzm(i,j,k+1,QPRES) - cdtdx * ( dupx + pavx*dux*(gamma - 1.))
                                     - cdtdy * ( dupy + pavy*duy*(gamma - 1.));
  }
  qm(i, j, k + 1, QPRES) = amrex::max(temppres, small_pres);
  qm(i, j, k + 1, QREINT) = qm(i, j, k + 1, QPRES) / (gamma - amrex::Real(1.0));

}

// Y corrected by x and z
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transxz(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qx,
  amrex::Array4<const amrex::Real> const& qz,
  const amrex::Real cdtdx,
  const amrex::Real cdtdz, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, Density_comp) - flxx(i, j, k, Density_comp)) +
           cdtdz * (flxz(i, j, k + 1, Density_comp) - flxz(i, j, k, Density_comp));
  flxu = cdtdx * (flxx(i + 1, j, k, Xmom_comp) - flxx(i, j, k, Xmom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Xmom_comp) - flxz(i, j, k, Xmom_comp));
  flxv = cdtdx * (flxx(i + 1, j, k, Ymom_comp) - flxx(i, j, k, Ymom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Ymom_comp) - flxz(i, j, k, Ymom_comp));
  flxw = cdtdx * (flxx(i + 1, j, k, Zmom_comp) - flxx(i, j, k, Zmom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Zmom_comp) - flxz(i, j, k, Zmom_comp));
  flxe = cdtdx * (flxx(i + 1, j, k, Eden_comp) - flxx(i, j, k, Eden_comp)) +
         cdtdz * (flxz(i, j, k + 1, Eden_comp) - flxz(i, j, k, Eden_comp));

  // Update passive variables
  rrr = qyp(i, j, k, QRHO);
  rrl = qym(i, j + 1, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;

#ifndef CONST_SPECIES
  amrex::Real compn, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qyp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnewr;
    // qm
    compo = rrl * qym(i, j + 1, k, n) - compn;
    qm(i, j + 1, k, n) = compo / rrnewl;
  }
#endif

  // Update hydro vars
  pggpx = qx(i + 1, j, k, QPRES);
  pggmx = qx(i, j, k, QPRES);
  ugpx = qx(i + 1, j, k, QU);
  ugmx = qx(i, j, k, QU);

  pggpz = qz(i, j, k + 1, QPRES);
  pggmz = qz(i, j, k, QPRES);
  ugpz = qz(i, j, k + 1, QW);
  ugmz = qz(i, j, k, QW);

  // QP
  rur = rrr * qyp(i, j, k, QU);
  rvr = rrr * qyp(i, j, k, QV);
  rwr = rrr * qyp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qyp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupx = pggpx * ugpx - pggmx * ugmx;
       amrex::Real pavx = 0.5 * (pggpx + pggmx);
       amrex::Real dux = ugpx - ugmx;

       amrex::Real dupz = pggpz * ugpz - pggmz * ugmz;
       amrex::Real pavz = 0.5 * (pggpz + pggmz);
       amrex::Real duz = ugpz - ugmz;

       temppres = qyp(i,j,k,QPRES) - cdtdx * ( dupx + pavx*dux*(gamma - 1.))
                                   - cdtdz * ( dupz + pavz*duz*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(temppres, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM
  rul = rrl * qym(i, j + 1, k, QU);
  rvl = rrl * qym(i, j + 1, k, QV);
  rwl = rrl * qym(i, j + 1, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  temppres = qm(i, j + 1, k, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupx = pggpx * ugpx - pggmx * ugmx;
       amrex::Real pavx = 0.5 * (pggpx + pggmx);
       amrex::Real dux = ugpx - ugmx;

       amrex::Real dupz = pggpz * ugpz - pggmz * ugmz;
       amrex::Real pavz = 0.5 * (pggpz + pggmz);
       amrex::Real duz = ugpz - ugmz;

       temppres = qym(i,j+1,k,QPRES) - cdtdx * ( dupx + pavx*dux*(gamma - 1.))
                                     - cdtdz * ( dupz + pavz*duz*(gamma - 1.));
  }
  qm(i, j + 1, k, QPRES) = amrex::max(temppres, small_pres);
  qm(i, j + 1, k, QREINT) = qm(i, j + 1, k, QPRES) / (gamma - amrex::Real(1.0));

}

// x corrected by y and z
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transyz(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qy,
  amrex::Array4<const amrex::Real> const& qz,
  const amrex::Real cdtdy,
  const amrex::Real cdtdz, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const amrex::Real gamma, 
  const amrex::Real small_pres)
{
  amrex::Real rrr, rur, rvr, rwr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, Density_comp) - flxy(i, j, k, Density_comp)) +
           cdtdz * (flxz(i, j, k + 1, Density_comp) - flxz(i, j, k, Density_comp));
  flxu = cdtdy * (flxy(i, j + 1, k, Xmom_comp) - flxy(i, j, k, Xmom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Xmom_comp) - flxz(i, j, k, Xmom_comp));
  flxv = cdtdy * (flxy(i, j + 1, k, Ymom_comp) - flxy(i, j, k, Ymom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Ymom_comp) - flxz(i, j, k, Ymom_comp));
  flxw = cdtdy * (flxy(i, j + 1, k, Zmom_comp) - flxy(i, j, k, Zmom_comp)) +
         cdtdz * (flxz(i, j, k + 1, Zmom_comp) - flxz(i, j, k, Zmom_comp));
  flxe = cdtdy * (flxy(i, j + 1, k, Eden_comp) - flxy(i, j, k, Eden_comp)) +
         cdtdz * (flxz(i, j, k + 1, Eden_comp) - flxz(i, j, k, Eden_comp));

  // Update passive variables
  rrr = qxp(i, j, k, QRHO);
  rrl = qxm(i + 1, j, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;

#ifndef CONST_SPECIES
  amrex::Real compn, compo;
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qxp(i, j, k, n) - compn;
    qp(i, j, k, n) = compo / rrnewr;
    // qm
    compo = rrl * qxm(i + 1, j, k, n) - compn;
    qm(i + 1, j, k, n) = compo / rrnewl;
  }
#endif

  // Update hydro vars
  pggpz = qz(i, j, k + 1, QPRES);
  pggmz = qz(i, j, k, QPRES);
  ugpz = qz(i, j, k + 1, QW);
  ugmz = qz(i, j, k, QW);

  pggpy = qy(i, j + 1, k, QPRES);
  pggmy = qy(i, j, k, QPRES);
  ugpy = qy(i, j + 1, k, QV);
  ugmy = qy(i, j, k, QV);

  // QP
  rur = rrr * qxp(i, j, k, QU);
  rvr = rrr * qxp(i, j, k, QV);
  rwr = rrr * qxp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qxp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Reset to original value if adding transverse terms made density negative
  bool mask = rrnewr < 0.0;
  rrnewr = mask ? rrr : rrnewr;
  runewr = mask ? rur : runewr;
  rvnewr = mask ? rvr : rvnewr;
  rwnewr = mask ? rwr : rwnewr;
  renewr = mask ? rer : renewr;

  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupy = pggpy * ugpy - pggmy * ugmy;
       amrex::Real pavy = 0.5 * (pggpy + pggmy);
       amrex::Real duy = ugpy - ugmy;

       amrex::Real dupz = pggpz * ugpz - pggmz * ugmz;
       amrex::Real pavz = 0.5 * (pggpz + pggmz);
       amrex::Real duz = ugpz - ugmz;

       temppres = qxp(i,j,k,QPRES) - cdtdy * ( dupy + pavy*duy*(gamma - 1.))
                                   - cdtdz * ( dupz + pavz*duz*(gamma - 1.));
  }
  qp(i, j, k, QPRES) = amrex::max(temppres, small_pres);
  qp(i, j, k, QREINT) = qp(i, j, k, QPRES) / (gamma - amrex::Real(1.0));

  // QM
  rul = rrl * qxm(i + 1, j, k, QU);
  rvl = rrl * qxm(i + 1, j, k, QV);
  rwl = rrl * qxm(i + 1, j, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  // Reset to original value if adding transverse terms made density negative
  mask = rrnewl < 0.0;
  rrnewl = mask ? rrl : rrnewl;
  runewl = mask ? rul : runewl;
  rvnewl = mask ? rvl : rvnewl;
  rwnewl = mask ? rwl : rwnewl;
  renewl = mask ? rel : renewl;

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
  temppres = qm(i + 1, j, k, QREINT) * (gamma - amrex::Real(1.0));
  if(temppres < small_pres)
  {

       amrex::Real dupy = pggpy * ugpy - pggmy * ugmy;
       amrex::Real pavy = 0.5 * (pggpy + pggmy);
       amrex::Real duy = ugpy - ugmy;

       amrex::Real dupz = pggpz * ugpz - pggmz * ugmz;
       amrex::Real pavz = 0.5 * (pggpz + pggmz);
       amrex::Real duz = ugpz - ugmz;

       temppres = qxm(i+1,j,k,QPRES) - cdtdy * ( dupy + pavy*duy*(gamma - 1.))
                                     - cdtdz * ( dupz + pavz*duz*(gamma - 1.));
  }
  qm(i + 1, j, k, QPRES) = amrex::max(temppres, small_pres);
  qm(i + 1, j, k, QREINT) = qm(i + 1, j, k, QPRES) / (gamma - amrex::Real(1.0));
}

// Use interface states from Riemann solver for pdivu.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_pdivu(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& q1,
  amrex::Array4<const amrex::Real> const& q2,
  amrex::Array4<const amrex::Real> const& q3,
  const amrex::Real dx, const amrex::Real dy, const amrex::Real dz)
{
        pdivu(i,j,k) = -(
                                  (q1(i+1,j,k,QU) - q1(i,j,k,QU))/dx +
                                  (q2(i,j+1,k,QV) - q2(i,j,k,QV))/dy +
                                  (q3(i,j,k+1,QW) - q3(i,j,k,QW))/dz);
}

// Consup Support Functions
// Designed for CPU or GPU
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_artif_visc(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& div,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Real const dx,
  amrex::Real const difmag,
  const int dir)
{
  amrex::Real div1;
  if (dir == 0) { // X direction
    AMREX_D_PICK(div1 = div(i, j, k);
                 , div1 = 0.5 * (div(i, j, k) + div(i, j + 1, k));
                 , div1 = 0.25 * (div(i, j, k) + div(i, j + 1, k) +
                                  div(i, j, k + 1) + div(i, j + 1, k + 1)););
    div1 = difmag * amrex::min(amrex::Real(0.0), div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        // Passing in difmag TODO not pass in.
        flx(i, j, k, n) += dx * div1 * (u(i, j, k, n) - u(i - 1, j, k, n));
    }
  } else if (dir == 1) { // Y direction
    AMREX_D_PICK(, div1 = 0.5 * (div(i, j, k) + div(i + 1, j, k));
                 , div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) +
                                  div(i, j, k + 1) + div(i + 1, j, k + 1)););
    div1 = difmag * amrex::min(amrex::Real(0.0), div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j - 1, k, n)); // Here dx is dy
    }
  } else { // Z direction
    div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) + div(i, j + 1, k) +
                   div(i + 1, j + 1, k));
    div1 = difmag * amrex::min(amrex::Real(0.0), div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j, k - 1, n)); // Here dx is dz
    }
  }
}

// Host Functions
void pc_umeth_3D(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<amrex::Real> const& flx1,
  amrex::Array4<amrex::Real> const& flx2,
  amrex::Array4<amrex::Real> const& flx3,
  amrex::Array4<amrex::Real> const& q1,
  amrex::Array4<amrex::Real> const& q2,
  amrex::Array4<amrex::Real> const& q3,
  amrex::Array4<amrex::Real> const& pdivu,
  const amrex::Real* del,
  const amrex::Real dt,
  const amrex::Real a_old,
  const amrex::Real a_new, 
  const amrex::Real gamma, const amrex::Real gamma_minus_1,
  const amrex::Real small_dens, const amrex::Real small_pres, 
  const amrex::Real small,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  const int ppm_type);
#endif
