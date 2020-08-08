#ifndef _PLM_H_
#define _PLM_H_

#include <cmath>

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include "Constants.H"
#include "IndexDefines.H"
#include "Riemann.H"

/*This header file contains functions and declarations for the 3D hydro portion
  of PeleC cuda. As per the convention of AMReX, inlined functions are defined
  here. Where as non-inline functions are declared here. */

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
plm_slope(
  const int i,
  const int j,
  const int k,
  const int n,
  const int dir,
  const int use_flattening,
  amrex::Array4<const amrex::Real> const& q)
{
  const int bdim[3] = {dir == 0, dir == 1, dir == 2};

  amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, drgtm1, dlftp1, slop;
  const amrex::Real qc = q(i, j, k, n);
  const amrex::Real qm =
    q(i - 1 * bdim[0], j - 1 * bdim[1], k - 1 * bdim[2], n);
  const amrex::Real qp =
    q(i + 1 * bdim[0], j + 1 * bdim[1], k + 1 * bdim[2], n);

  dlft = (qm - q(i - 2 * bdim[0], j - 2 * bdim[1], k - 2 * bdim[2], n))*2.0;
  drgt = (qc - qm)*2.0;

  dcen = 0.25 * (dlft + drgt);
  dsgn = amrex::Math::copysign(1.0, dcen);
  dlim = (dlft * drgt >= 0.0)
           ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt))
           : 0.0;
  dfm = dsgn * amrex::min(dlim, amrex::Math::abs(dcen));
  drgtm1=drgt;

  dlft = (qp - qc)*2.0;
  drgt = (q(i + 2 * bdim[0], j + 2 * bdim[1], k + 2 * bdim[2], n) - qp)*2.0;
  dcen = 0.25 * (dlft + drgt);
  dsgn = amrex::Math::copysign(1.0, dcen);
  dlim = (dlft * drgt >= 0.0)
           ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt))
           : 0.0;
  dfp = dsgn * amrex::min(dlim, amrex::Math::abs(dcen));
  dlftp1=dlft;

  //Replace with limited fourth order slopes
  dlft = drgtm1;
  drgt = dlftp1;
  dcen = 0.25 * (dlft + drgt);
  dsgn = amrex::Math::copysign(1.0, dcen);
  dlim = (dlft * drgt >= 0.0)
           ? amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt))
           : 0.0;
  //Use limiter on dq, not dcen

  dtemp = 4.0 / 3.0 * dcen - 1.0 / 6.0 * (dfp + dfm);

  amrex::Real flatn = 1.0;
  //Calculate flattening in-place
  if (use_flattening == 1)
  {
    for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++)
      flatn = amrex::min(flatn,flatten(i, j, k, dir_flat, q));
  }

  return flatn * dsgn * amrex::min(dlim, amrex::Math::abs(dtemp));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_x(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  // amrex::Array4<const amrex::Real> const& dloga,
  amrex::Real const a_old,
  amrex::Real const dx,
  amrex::Real const dt, const int NumSpec)
{
  amrex::Real e[3];
  amrex::Real dtdx = dt / (a_old * dx);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  amrex::Real rhoe = q(i, j, k, QREINT);
  amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - du) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + du) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0v = dv;
  amrex::Real alpha0w = dw;

  e[0] = u - cc;
  e[1] = u;
  e[2] = u + cc;
  //  Construct the right state on the i-1/2 interface
  amrex::Real rho_ref = rho;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * drho;
  amrex::Real u_ref = u;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * du;
  amrex::Real v_ref = v;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * dv;
  amrex::Real w_ref = w;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * dw;
  amrex::Real p_ref = p;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * dp;
  amrex::Real rhoe_ref =
      rhoe;// - 0.5 * (1.0 + dtdx * amrex::min(e[0], 0.0)) * drhoe;
  amrex::Real spminus;
  amrex::Real spplus;
  amrex::Real spzero;

   if (u-cc > 0.0)
           spminus = -1.0;
   else
           spminus = (u-cc)*dtdx;

   if (u+cc > 0.0)
           spplus = -1.0;
   else
           spplus = (u+cc)*dtdx;

   if (u > 0.0)
           spzero = -1.0;
   else
           spzero = u*dtdx;

  amrex::Real apright = 0.5*(-1.0 - spplus )*alphap;
  amrex::Real amright = 0.5*(-1.0 - spminus )*alpham;
  amrex::Real azrright = 0.5*(-1.0 - spzero )*alpha0r;
  amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azv1rght = 0.5*(-1.0 - spzero )*alpha0v;
  amrex::Real azw1rght = 0.5*(-1.0 - spzero )*alpha0w;
  
  qp(i, j, k, QRHO) = rho_ref + apright + amright + azrright;
  qp(i, j, k, QRHO) = amrex::max(SMALL_DENS, qp(i, j, k, QRHO));
  qp(i, j, k, QU) = u_ref + (apright - amright) * cc / rho;
  qp(i, j, k, QV) = 0.; qp(i, j, k, QW) = 0.;
  qp(i, j, k, QV) = v_ref + azv1rght;
  qp(i, j, k, QW) = w_ref + azw1rght;
  qp(i, j, k, QPRES) = p_ref + (apright + amright) * cs;
  qp(i, j, k, QPRES) = amrex::max(qp(i, j, k, QPRES), SMALL_PRES);
  qp(i, j, k, QREINT) = rhoe_ref + (apright + amright) * enth * cs + azeright;

  //  Construct the left state on the i+1/2 interface

  rho_ref = rho;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * drho;
  u_ref = u;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * du;
  v_ref = v;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * dv;
  w_ref = w;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * dw;
  p_ref = p;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * dp;
  rhoe_ref = rhoe;// + 0.5 * (1.0 - dtdx * amrex::max(e[2], 0.0)) * drhoe;

  //                  apleft = 0.25*dtdx*(e[2] - e[2])*(1.0 +
  //                  amrex::Math::copysign(1.0,e[2]))*alphap;

   if (u-cc >= 0.0)
           spminus = (u-cc)*dtdx;
   else
           spminus = 1.0;

   if (u+cc >= 0.0)
           spplus = (u+cc)*dtdx;
   else
           spplus = 1.0;
   if (u >= 0.0)
           spzero = u*dtdx;
   else
           spzero = 1.0;

  amrex::Real apleft = 0.5*(1.0 - spplus )*alphap;
  amrex::Real amleft = 0.5*(1.0 - spminus )*alpham;
  amrex::Real azrleft = 0.5*(1.0 - spzero )*alpha0r;
  amrex::Real azeleft = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azv1left = 0.5*(1.0 - spzero )*alpha0v;
  amrex::Real azw1left = 0.5*(1.0 - spzero )*alpha0w;

  qm(i + 1, j, k, QRHO) = rho_ref + apleft + amleft + azrleft;
  qm(i + 1, j, k, QRHO) = amrex::max(qm(i + 1, j, k, QRHO), SMALL_DENS);
  qm(i + 1, j, k, QU) = u_ref + (apleft - amleft) * cc / rho;
  qm(i + 1, j, k, QV) = v_ref + azv1left;
  qm(i + 1, j, k, QW) = w_ref + azw1left;
  qm(i + 1, j, k, QPRES) = p_ref + (apleft + amleft) * cs;
  qm(i + 1, j, k, QPRES) = amrex::max(qm(i + 1, j, k, QPRES), SMALL_PRES);
  qm(i + 1, j, k, QREINT) = rhoe_ref + (apleft + amleft) * enth * cs + azeleft;

  // Upwind the passive variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {

    // Right state
    u = q(i, j, k, QU);
    amrex::Real spzero;
    spzero = u > 0 ? -1.0 : u * dtdx;
    amrex::Real acmprght = 0.5 * (-1.0 - spzero) * slope[n];
    qp(i, j, k, n) = q(i, j, k, n) + acmprght;

    spzero = u >= 0 ? u * dtdx : 1.0;
    amrex::Real acmpleft = 0.5 * (1.0 - spzero) * slope[n];
    qm(i + 1, j, k, n) = q(i, j, k, n) + acmpleft;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_y(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  amrex::Real const a_old,
  amrex::Real const dy,
  amrex::Real const dt, const int NumSpec)
{
  amrex::Real e[3];
  amrex::Real dtdy = dt / (a_old * dy);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  amrex::Real rhoe = q(i, j, k, QREINT);
  amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - dv) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + dv) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0u = du;
  amrex::Real alpha0w = dw;

  e[0] = v - cc;
  e[1] = v;
  e[2] = v + cc;
  // Construct the "right" state on the j-1/2 interface
  amrex::Real rho_ref = rho;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * drho;
  amrex::Real u_ref = u;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * du;
  amrex::Real v_ref = v;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * dv;
  amrex::Real w_ref = w;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * dw;
  amrex::Real p_ref = p;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * dp;
  amrex::Real rhoe_ref =
    rhoe;// - 0.5 * (1.0 + dtdy * amrex::min(e[0], 0.0)) * drhoe;
  amrex::Real spminus;
  amrex::Real spplus;
  amrex::Real spzero;

   if (v-cc > 0.0)
           spminus = -1.0;
   else
           spminus = (v-cc)*dtdy;

   if (v+cc > 0.0)
           spplus = -1.0;
   else
           spplus = (v+cc)*dtdy;

   if (v > 0.0)
           spzero = -1.0;
   else
           spzero = v*dtdy;

  amrex::Real apright = 0.5*(-1.0 - spplus )*alphap;
  amrex::Real amright = 0.5*(-1.0 - spminus )*alpham;
  amrex::Real azrright = 0.5*(-1.0 - spzero )*alpha0r;
  amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azu1rght = 0.5*(-1.0 - spzero )*alpha0u;
  amrex::Real azw1rght = 0.5*(-1.0 - spzero )*alpha0w;

  qp(i, j, k, QRHO) = rho_ref + apright + amright + azrright;
  qp(i, j, k, QRHO) = amrex::max(SMALL_DENS, qp(i, j, k, QRHO));
  qp(i, j, k, QU) = u_ref + azu1rght;
  qp(i, j, k, QV) = v_ref + (apright - amright) * cc / rho;
  qp(i, j, k, QW) = w_ref + azw1rght;
  qp(i, j, k, QPRES) = p_ref + (apright + amright) * cs;
  qp(i, j, k, QPRES) = amrex::max(qp(i, j, k, QPRES), SMALL_PRES);
  qp(i, j, k, QREINT) = rhoe_ref + (apright + amright) * enth * cs + azeright;

  // Construct the "left" state on the j+1/2 interface

  rho_ref = rho;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * drho;
  u_ref = u;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * du;
  v_ref = v;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * dv;
  w_ref = w;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * dw;
  p_ref = p;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * dp;
  rhoe_ref = rhoe;// + 0.5 * (1.0 - dtdy * amrex::max(e[2], 0.0)) * drhoe;

  //                 apleft = 0.25*dtdy*(e[2] - e[2])*(1.0 +
  //                 amrex::Math::copysign(1.0,e[2]))*alphap;

   if (v-cc >= 0.0)
           spminus = (v-cc)*dtdy;
   else
           spminus = 1.0;

   if (v+cc >= 0.0)
           spplus = (v+cc)*dtdy;
   else
           spplus = 1.0;
   if (v >= 0.0)
           spzero = v*dtdy;
   else
           spzero = 1.0;

  amrex::Real apleft = 0.5*(1.0 - spplus )*alphap;
  amrex::Real amleft = 0.5*(1.0 - spminus )*alpham;
  amrex::Real azrleft = 0.5*(1.0 - spzero )*alpha0r;
  amrex::Real azeleft = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azu1left = 0.5*(1.0 - spzero )*alpha0u;
  amrex::Real azw1left = 0.5*(1.0 - spzero )*alpha0w;

  qm(i, j + 1, k, QRHO) = rho_ref + apleft + amleft + azrleft;
  qm(i, j + 1, k, QRHO) = amrex::max(qm(i, j + 1, k, QRHO), SMALL_DENS);
  qm(i, j + 1, k, QU) = u_ref + azu1left;
  qm(i, j + 1, k, QV) = v_ref + (apleft - amleft) * cc / rho;
  qm(i, j + 1, k, QW) = w_ref + azw1left;
  qm(i, j + 1, k, QPRES) = p_ref + (apleft + amleft) * cs;
  qm(i, j + 1, k, QPRES) = amrex::max(qm(i, j + 1, k, QPRES), SMALL_PRES);
  qm(i, j + 1, k, QREINT) = rhoe_ref + (apleft + amleft) * enth * cs + azeleft;

  // Upwind the passive variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {

    // Right state
    v = q(i, j, k, QV);
    amrex::Real spzero;
    spzero = v > 0 ? -1.0 : v * dtdy;
    amrex::Real acmptop = 0.5 * (-1.0 - spzero) * slope[n];
    qp(i, j, k, n) = q(i, j, k, n) + acmptop;

    spzero = v >= 0 ? v * dtdy : 1.0;
    amrex::Real acmpbot = 0.5 * (1.0 - spzero) * slope[n];
    qm(i, j + 1, k, n) = q(i, j, k, n) + acmpbot;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_z(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  amrex::Real const a_old,
  amrex::Real const dz,
  amrex::Real const dt, const int NumSpec)
{
  amrex::Real e[3];
  amrex::Real dtdz = dt / (a_old*dz);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  amrex::Real rhoe = q(i, j, k, QREINT);
  amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - dw) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + dw) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0u = du;
  amrex::Real alpha0v = dv;

  e[0] = w - cc;
  e[1] = w;
  e[2] = w + cc;

  // Construct the "right" state on the k-1/2 interface
  amrex::Real rho_ref = rho;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * drho;
  amrex::Real u_ref = u;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * du;
  amrex::Real v_ref = v;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dv;
  amrex::Real w_ref = w;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dw;
  amrex::Real p_ref = p;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dp;
  amrex::Real rhoe_ref =
    rhoe;// - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * drhoe;

  amrex::Real spminus;
  amrex::Real spplus;
  amrex::Real spzero;

   if (w-cc > 0.0)
           spminus = -1.0;
   else
           spminus = (w-cc)*dtdz;

   if (w+cc > 0.0)
           spplus = -1.0;
   else
           spplus = (w+cc)*dtdz;

   if (w > 0.0)
           spzero = -1.0;
   else
           spzero = w*dtdz;

  amrex::Real apright = 0.5*(-1.0 - spplus )*alphap;
  amrex::Real amright = 0.5*(-1.0 - spminus )*alpham;
  amrex::Real azrright = 0.5*(-1.0 - spzero )*alpha0r;
  amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azu1rght = 0.5*(-1.0 - spzero )*alpha0u;
  amrex::Real azv1rght = 0.5*(-1.0 - spzero )*alpha0v;

  qp(i, j, k, QRHO) = rho_ref + apright + amright + azrright;
  qp(i, j, k, QRHO) = amrex::max(SMALL_DENS, qp(i, j, k, QRHO));
  qp(i, j, k, QU) = u_ref + azu1rght;
  qp(i, j, k, QV) = v_ref + azv1rght;
  qp(i, j, k, QW) = w_ref + (apright - amright) * cc / rho;
  qp(i, j, k, QPRES) = p_ref + (apright + amright) * cs;
  qp(i, j, k, QPRES) = amrex::max(qp(i, j, k, QPRES), SMALL_PRES);
  qp(i, j, k, QREINT) = rhoe_ref + (apright + amright) * enth * cs + azeright;

  // Construct the "left" state on the j+1/2 interface

  rho_ref = rho;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * drho;
  u_ref = u;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * du;
  v_ref = v;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dv;
  w_ref = w;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dw;
  p_ref = p;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dp;
  rhoe_ref = rhoe;// + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * drhoe;
   if (w-cc >= 0.0)
           spminus = (w-cc)*dtdz;
   else
           spminus = 1.0;

   if (w+cc >= 0.0)
           spplus = (w+cc)*dtdz;
   else
           spplus = 1.0;
   if (w >= 0.0)
           spzero = w*dtdz;
   else
           spzero = 1.0;

  amrex::Real apleft = 0.5*(1.0 - spplus )*alphap;
  amrex::Real amleft = 0.5*(1.0 - spminus )*alpham;
  amrex::Real azrleft = 0.5*(1.0 - spzero )*alpha0r;
  amrex::Real azeleft = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azu1left = 0.5*(1.0 - spzero )*alpha0u;
  amrex::Real azv1left = 0.5*(1.0 - spzero )*alpha0v;

  qm(i, j, k + 1, QRHO) = rho_ref + apleft + amleft + azrleft;
  qm(i, j, k + 1, QRHO) = amrex::max(qm(i, j, k + 1, QRHO), SMALL_DENS);
  qm(i, j, k + 1, QU) = u_ref + azu1left;
  qm(i, j, k + 1, QV) = v_ref + azv1left;
  qm(i, j, k + 1, QW) = w_ref + (apleft - amleft) * cc / rho;
  qm(i, j, k + 1, QPRES) = p_ref + (apleft + amleft) * cs;
  qm(i, j, k + 1, QPRES) = amrex::max(qm(i, j, k + 1, QPRES), SMALL_PRES);
  qm(i, j, k + 1, QREINT) = rhoe_ref + (apleft + amleft) * enth * cs + azeleft;

  // Upwind the passive variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {

    // Right state
    w = q(i, j, k, QW);
    amrex::Real spzero;
    spzero = w > 0 ? -1.0 : w * dtdz;
    amrex::Real acmptop = 0.5 * (-1.0 - spzero) * slope[n];
    qp(i, j, k, n) = q(i, j, k, n) + acmptop;

    spzero = w >= 0 ? w * dtdz : 1.0;
    amrex::Real acmpbot = 0.5 * (1.0 - spzero) * slope[n];
    qm(i, j, k + 1, n) = q(i, j, k, n) + acmpbot;
  }
}
#endif
