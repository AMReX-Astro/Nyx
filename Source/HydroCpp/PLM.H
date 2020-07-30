#ifndef _PLM_H_
#define _PLM_H_

#include <cmath>

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include "Nyx.H"
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
  if(i==0&&j==2&&k==0&&dir==1)
  {
	  amrex::Print()<<qc<<"\t"<<qm<<"\t"<<dlft<<
		  "\t"<<dcen<<"\t"<<drgt<<"\t"<<dfm<<"\t"<<dlim<<"\t"<<n<<std::endl;
	  //	  if(n==FirstSpec+2)
		  //	  exit(0);
  }
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
  // flattening can go here.
  if(i==0&&j==2&&k==0&&dir==1)
  {
	  amrex::Print()<<dtemp<<"\t"<<dcen<<"\t"<<dfp<<"\t"<<dfm<<"\t"<<dlim<<std::endl;
	  //	  exit(0);
  }
  return dsgn * amrex::min(dlim, amrex::Math::abs(dtemp));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpflx(
  const int i,
  const int j,
  const int k,
  const int bclo,
  const int bchi,
  const int domlo,
  const int domhi,
  amrex::Array4<const amrex::Real> const& ql,
  amrex::Array4<const amrex::Real> const& qr,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<amrex::Real> const& q,
  // amrex::Array4<const int> const& bcMask,
  amrex::Array4<const amrex::Real> const& qprim,
  const amrex::Real small_dens,
  const amrex::Real small_pres,
  const amrex::Real small_vel,
  const amrex::Real small,
  const amrex::Real gamma,
  const int FirstSpec, const int NumSpec,
  const int dir)
{
  amrex::Real cav, ustar;
  amrex::Real ul, ur, vl, vr, v2l, v2r, rel, rer;
  int idx;
  int IU, IV, IV2;
  int GU, GV, GV2;
  int f_idx[3];

  if (dir == 0) {
    IU = QU;
    IV = QV;
    IV2 = QW;
    GU = GDU;
    GV = GDV;
    GV2 = GDW;
    cav = 0.5 * ( std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i-1,j,k,QREINT) / qprim(i-1,j,k,QRHO)))
                  +std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i  ,j,k,QREINT) / qprim(i  ,j,k,QRHO))) );
    f_idx[0] = UMX;
    f_idx[1] = UMY;
    f_idx[2] = UMZ;
  } else if (dir == 1) {
    IU = QV;
    IV = QU;
    IV2 = QW;
    GU = GDV;
    GV = GDU;
    GV2 = GDW;
    cav = 0.5 * ( std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i,j-1,k,QREINT) / qprim(i,j-1,k,QRHO)))
				  +std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i,j  ,k,QREINT) / qprim(i,j  ,k,QRHO))) );
    f_idx[0] = UMY;
    f_idx[1] = UMX;
    f_idx[2] = UMZ;
  } else {
    IU = QW;
    IV = QU;
    IV2 = QV;
    GU = GDW;
    GV = GDU;
    GV2 = GDV;
    cav = 0.5 * ( std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i,j,k-1,QREINT) / qprim(i,j,k-1,QRHO)))
				  +std::sqrt(std::abs(gamma * (gamma - 1.0) * qprim(i,j,k  ,QREINT) / qprim(i,j,k  ,QRHO))) );
    f_idx[0] = UMZ;
    f_idx[1] = UMX;
    f_idx[2] = UMY;
  }

  ul = ql(i, j, k, IU);
  vl = ql(i, j, k, IV);
  v2l = ql(i, j, k, IV2);
  rel = ql(i, j, k, QREINT);

  ur = qr(i, j, k, IU);
  vr = qr(i, j, k, IV);
  v2r = qr(i, j, k, IV2);
  rer = qr(i, j, k, QREINT);

  // Outflow Hack
  if (dir == 2)
    idx = k;
  else
    idx = (dir == 0) ? i : j;
  if (bclo == Outflow && idx == domlo) {
    ul = ur;
    vl = vr;
    v2l = v2r;
    rel = rer;
  }
  if (bchi == Outflow && idx == domhi + 1) {
    ur = ul;
    vr = vl;
    v2r = v2l;
    rer = rel;
  }

  const int bc_test_val = 1;

  //    ql(i, j, k, QRHO), ul, vl, v2l, ql(i, j, k, QPRES), rel, spl, 
  //    qr(i, j, k, QRHO), ur, vr, v2r, qr(i, j, k, QPRES), rer, spr, 

  amrex::Real c = std::sqrt((gamma) * qprim(i,j,k,QREINT)/qprim(i,j,k,QRHO)*(gamma - 1.0));//std::sqrt(std::abs(gamma * q(i,j,k,QPRES) / q(i,j,k,QRHO)));
  riemann(
    ql(i, j, k, QRHO), ul, vl, v2l, ql(i, j, k, QPRES), rel,
    qr(i, j, k, QRHO), ur, vr, v2r, qr(i, j, k, QPRES), rer,
    bc_test_val, c, cav, ustar, flx(i, j, k, URHO),
    flx(i, j, k, f_idx[0]), flx(i, j, k, f_idx[1]), flx(i, j, k, f_idx[2]),
    flx(i, j, k, UEDEN), flx(i, j, k, UEINT), q(i, j, k, GU), q(i, j, k, GV),
    q(i, j, k, GV2), q(i, j, k, GDPRES), q(i, j, k, GDGAME), 
    small_dens, small_pres, small_vel, small, gamma, NumSpec);

  amrex::Real flxrho = flx(i, j, k, URHO);
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    if (ustar > 0.0)
      flx(i, j, k, n) = flxrho * ql(i, j, k, nqp);
    else if (ustar < 0.0)
      flx(i, j, k, n) = flxrho * qr(i, j, k, nqp);
    else
      flx(i, j, k, n) = flxrho * 0.50 * (ql(i, j, k, nqp) + qr(i, j, k, nqp));
  }
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
  amrex::Real const dx,
  amrex::Real const dt, const int NumSpec, const amrex::Real gamma)
{
  amrex::Real e[3];
  amrex::Real dtdx = dt / dx;
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
  AMREX_D_TERM(qp(i, j, k, QU) = u_ref + (apright - amright) * cc / rho;
               qp(i, j, k, QV) = 0.; qp(i, j, k, QW) = 0.;
               , qp(i, j, k, QV) = v_ref + azv1rght;
               , qp(i, j, k, QW) = w_ref + azw1rght;);
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
  AMREX_D_TERM(qm(i + 1, j, k, QU) = u_ref + (apleft - amleft) * cc / rho;
               qm(i + 1, j, k, QV) = 0.; qm(i + 1, j, k, QW) = 0.;
               , qm(i + 1, j, k, QV) = v_ref + azv1left;
               , qm(i + 1, j, k, QW) = w_ref + azw1left;);
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
  amrex::Real const dy,
  amrex::Real const dt, const int NumSpec, const amrex::Real gamma)
{
  amrex::Real e[3];
  amrex::Real dtdy = dt / dy;
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
  AMREX_D_TERM(qp(i, j, k, QW) = 0.;, , qp(i, j, k, QW) = w_ref + azw1rght;);
  qp(i, j, k, QPRES) = p_ref + (apright + amright) * cs;
  qp(i, j, k, QPRES) = amrex::max(qp(i, j, k, QPRES), SMALL_PRES);
  qp(i, j, k, QREINT) = rhoe_ref + (apright + amright) * enth * cs + azeright;
  if(i==0&&j==2&&k==0)
	  {
		  amrex::Print()<<"plm_y\n"<<rho_ref<<"\t"<<apright<<"\t"<<amright<<"\t"<<azrright<<"\t"<<azu1rght<<"\t"<<azw1rght<< std::endl;
		  amrex::Print()<<u<<"\t"<<cc<<"\t"<<dp<<"\t"<<du<<std::endl;
		  //		  exit(0);
	  }
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
  AMREX_D_TERM(qm(i, j + 1, k, QW) = 0.;
               , , qm(i, j + 1, k, QW) = w_ref + azw1left;);
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
  amrex::Real const dz,
  amrex::Real const dt, const int NumSpec, const amrex::Real gamma)
{
  amrex::Real e[3];
  amrex::Real dtdz = dt / dz;
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
  amrex::Real rho_ref = rho - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * drho;
  amrex::Real u_ref = u - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * du;
  amrex::Real v_ref = v - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dv;
  amrex::Real w_ref = w - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dw;
  amrex::Real p_ref = p - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * dp;
  amrex::Real rhoe_ref =
    rhoe - 0.5 * (1.0 + dtdz * amrex::min(e[0], 0.0)) * drhoe;

  amrex::Real apright = 0.25 * dtdz * (e[0] - e[2]) *
                        (1.0 - amrex::Math::copysign(1.0, e[2])) * alphap;
  //     amright = 0.25*dtdz*(e[0] - e[0])*(1.0 - amrex::Math::copysign(1.0,
  //     e[0]))*alpham;
  amrex::Real amright = 0.0;

  amrex::Real azrright = 0.25 * dtdz * (e[0] - e[1]) *
                         (1.0 - amrex::Math::copysign(1.0, e[1])) * alpha0r;
  amrex::Real azeright = 0.25 * dtdz * (e[0] - e[1]) *
                         (1.0 - amrex::Math::copysign(1.0, e[1])) * alpha0e;
  amrex::Real azu1rght = 0.25 * dtdz * (e[0] - e[1]) *
                         (1.0 - amrex::Math::copysign(1.0, e[1])) * alpha0u;
  amrex::Real azv1rght = 0.25 * dtdz * (e[0] - e[1]) *
                         (1.0 - amrex::Math::copysign(1.0, e[1])) * alpha0v;

  qp(i, j, k, QRHO) = rho_ref + apright + amright + azrright;
  qp(i, j, k, QRHO) = amrex::max(SMALL_DENS, qp(i, j, k, QRHO));
  qp(i, j, k, QU) = u_ref + azu1rght;
  qp(i, j, k, QV) = v_ref + azv1rght;
  qp(i, j, k, QW) = w_ref + (apright - amright) * cc / rho;
  qp(i, j, k, QPRES) = p_ref + (apright + amright) * cs;
  qp(i, j, k, QPRES) = amrex::max(qp(i, j, k, QPRES), SMALL_PRES);
  qp(i, j, k, QREINT) = rhoe_ref + (apright + amright) * enth * cs + azeright;

  // Construct the "left" state on the j+1/2 interface

  rho_ref = rho + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * drho;
  u_ref = u + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * du;
  v_ref = v + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dv;
  w_ref = w + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dw;
  p_ref = p + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * dp;
  rhoe_ref = rhoe + 0.5 * (1.0 - dtdz * amrex::max(e[2], 0.0)) * drhoe;

  //                 apleft = 0.25*dtdz*(e[2] - e[2])*(1.0 +
  //                 amrex::Math::copysign(1.0,e[2]))*alphap;
  amrex::Real apleft = 0.0;
  amrex::Real amleft = 0.25 * dtdz * (e[2] - e[0]) *
                       (1.0 + amrex::Math::copysign(1.0, e[0])) * alpham;

  amrex::Real azrleft = 0.25 * dtdz * (e[2] - e[1]) *
                        (1.0 + amrex::Math::copysign(1.0, e[1])) * alpha0r;
  amrex::Real azeleft = 0.25 * dtdz * (e[2] - e[1]) *
                        (1.0 + amrex::Math::copysign(1.0, e[1])) * alpha0e;
  amrex::Real azu1left = 0.25 * dtdz * (e[2] - e[1]) *
                         (1.0 + amrex::Math::copysign(1.0, e[1])) * alpha0u;
  amrex::Real azv1left = 0.25 * dtdz * (e[2] - e[1]) *
                         (1.0 + amrex::Math::copysign(1.0, e[1])) * alpha0v;

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
  const amrex::Real cdtdx, const int NumSpec,
  const amrex::Real gamma) 
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN));

  // Update passive variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnew;
  }
  // Update hydro vars

  pggp = qint(i + 1, j, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i + 1, j, k, GDU);
  ugm = qint(i, j, k, GDU);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

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

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  amrex::Real pnewr = qp(i, j, k, QREINT) * (gamma - 1.0);
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);

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

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl = qm(i, j + 1, k, QREINT) * (gamma - 1.0);
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, SMALL_PRES);

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
  const amrex::Real cdtdy, const int NumSpec, 
  const amrex::Real gamma)
{

  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));

  // -  Update Passive Variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j + 1, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j + 1, k, GDV);
  ugm = qint(i, j, k, GDV);

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

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr; // + hdt*srcr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  qp(i, j, k, QREINT) = renewr - rhoekinr;
  pnewr = qp(i, j, k, QREINT) * (gamma-1.0);
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);

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

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
  pnewl = qm(i + 1, j, k, QREINT) * (gamma-1.0);
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
  //  amrex::Print()<<rrnewl<<"\t"<<rrl<<"\t"<<flxrho<<std::endl;
  //  exit(0);

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
  const amrex::Real cdtdx, const int NumSpec, 
  const amrex::Real gamma)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN));

  // -  Update Passive Variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  pggp = qint(i + 1, j, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i + 1, j, k, GDU);
  ugm = qint(i, j, k, GDU);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

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
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);

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

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qm(i, j, k + 1, QREINT) * (gamma-1.0);
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, SMALL_PRES);
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
  const amrex::Real cdtdy, const int NumSpec, 
  const amrex::Real gamma)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));

  // -  Update Passive Variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  pggp = qint(i, j + 1, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j + 1, k, GDV);
  ugm = qint(i, j, k, GDV);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

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
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);

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

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qm(i, j, k + 1, QREINT) * (gamma-1.0);
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, SMALL_PRES);
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
  const amrex::Real cdtdz, const int NumSpec, 
  const amrex::Real gamma)
{
  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));

  // -  Update Passive Variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j, k + 1, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j, k + 1, GDW);
  ugm = qint(i, j, k, GDW);

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

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr; // + hdt*srcr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rvnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  qp(i, j, k, QREINT) = renewr - rhoekinr;
  pnewr = qp(i, j, k, QREINT) * (gamma-1.0);
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);

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

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
  pnewl = qm(i + 1, j, k, QREINT) * (gamma-1.0);
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, SMALL_PRES);

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
  const amrex::Real cdtdz, const int NumSpec, 
  const amrex::Real gamma)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));

  // Update passive variables
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnew;
  }

  // Update hydro vars
  pggp = qint(i, j, k + 1, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j, k + 1, GDW);
  ugm = qint(i, j, k, GDW);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

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
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);


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

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qm(i, j + 1, k, QREINT) * (gamma-1.0);
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
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
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  const amrex::Real cdtdy, 
  const int NumSpec, const amrex::Real gamma)
{
  amrex::Real compn, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real dAupx, pavx, dAux;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real dAupy, pavy, dAuy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO)) +
           cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX)) +
         cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY)) +
         cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ)) +
         cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN)) +
         cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));

  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qzp(i, j, k, QRHO);
  rrl = qzm(i, j, k + 1, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    // qp
    compo = rrr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpx = qx(i + 1, j, k, GDPRES);
  pggmx = qx(i, j, k, GDPRES);
  ugpx = qx(i + 1, j, k, GDU);
  ugmx = qx(i, j, k, GDU);

  dAupx = pggpx * ugpx - pggmx * ugmx;
  pavx = 0.5 * (pggpx + pggmx);
  dAux = ugpx - ugmx;

  pggpy = qy(i, j + 1, k, GDPRES);
  pggmy = qy(i, j, k, GDPRES);
  ugpy = qy(i, j + 1, k, GDV);
  ugmy = qy(i, j, k, GDV);

  dAupy = pggpy * ugpy - pggmy * ugmy;
  pavy = 0.5 * (pggpy + pggmy);
  dAuy = ugpy - ugmy;
  amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (gamma-1.0));
  amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (gamma-1.0));

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

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - 1.0);
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

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

  qm(i, j, k + 1, QRHO) = rrnewl + hdt * srcrho;
  qm(i, j, k + 1, QU) = runewl / rrnewl + hdt * srcu;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j, k + 1, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qm(i, j, k + 1, QREINT) * (gamma - 1.0);
  qm(i, j, k + 1, QPRES) = amrex::max(temppres, SMALL_PRES);
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
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  const amrex::Real cdtdz, 
  const int NumSpec, const amrex::Real gamma)
{
  amrex::Real compn, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real dAupx, pavx, dAux;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real dAupz, pavz, dAuz;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO)) +
           cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX)) +
         cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY)) +
         cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ)) +
         cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN)) +
         cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));

  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qyp(i, j, k, QRHO);
  rrl = qym(i, j + 1, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;
  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpx = qx(i + 1, j, k, GDPRES);
  pggmx = qx(i, j, k, GDPRES);
  ugpx = qx(i + 1, j, k, GDU);
  ugmx = qx(i, j, k, GDU);

  dAupx = pggpx * ugpx - pggmx * ugmx;
  pavx = 0.5 * (pggpx + pggmx);
  dAux = ugpx - ugmx;

  pggpz = qz(i, j, k + 1, GDPRES);
  pggmz = qz(i, j, k, GDPRES);
  ugpz = qz(i, j, k + 1, GDW);
  ugmz = qz(i, j, k, GDW);

  dAupz = pggpz * ugpz - pggmz * ugmz;
  pavz = 0.5 * (pggpz + pggmz);
  dAuz = ugpz - ugmz;
  amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (gamma-1.0));
  amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (gamma-1.0));

  // qp state
  rur = rrr * qyp(i, j, k, QU);
  rvr = rrr * qyp(i, j, k, QV);
  rwr = rrr * qyp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qyp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - 1.0);
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

  // qm state
  rul = rrl * qym(i, j + 1, k, QU);
  rvl = rrl * qym(i, j + 1, k, QV);
  rwl = rrl * qym(i, j + 1, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j + 1, k, QRHO) = rrnewl + hdt * srcrho;
  qm(i, j + 1, k, QU) = runewl / rrnewl + hdt * srcu;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j + 1, k, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qm(i, j + 1, k, QREINT) * (gamma - 1.0);
  qm(i, j + 1, k, QPRES) = amrex::max(temppres, SMALL_PRES);
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
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdy,
  const amrex::Real cdtdz, 
  const int NumSpec, const amrex::Real gamma)
{
  amrex::Real compn, compo, rwr;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real dAupz, pavz, dAuz;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real dAupy, pavy, dAuy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO)) +
           cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX)) +
         cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY)) +
         cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ)) +
         cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN)) +
         cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));

  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qxp(i, j, k, QRHO);
  rrl = qxm(i + 1, j, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;

  for (int n = FirstSpec; n < FirstSpec + NumSpec; ++n) {
    int nqp = n;
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpz = qz(i, j, k + 1, GDPRES);
  pggmz = qz(i, j, k, GDPRES);
  ugpz = qz(i, j, k + 1, GDW);
  ugmz = qz(i, j, k, GDW);

  dAupz = pggpz * ugpz - pggmz * ugmz;
  pavz = 0.5 * (pggpz + pggmz);
  dAuz = ugpz - ugmz;

  pggpy = qy(i, j + 1, k, GDPRES);
  pggmy = qy(i, j, k, GDPRES);
  ugpy = qy(i, j + 1, k, GDV);
  ugmy = qy(i, j, k, GDV);

  dAupy = pggpy * ugpy - pggmy * ugmy;
  pavy = 0.5 * (pggpy + pggmy);
  dAuy = ugpy - ugmy;
  amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (gamma-1.0));
  amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (gamma-1.0));

  // qp state
  rur = rrr * qxp(i, j, k, QU);
  rvr = rrr * qxp(i, j, k, QV);
  rwr = rrr * qxp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qxp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qp(i, j, k, QREINT) * (gamma - 1.0);
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

  // qm state
  rul = rrl * qxm(i + 1, j, k, QU);
  rvl = rrl * qxm(i + 1, j, k, QV);
  rwl = rrl * qxm(i + 1, j, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i + 1, j, k, QRHO) = rrnewl + hdt * srcrho;
  qm(i + 1, j, k, QU) = runewl / rrnewl + hdt * srcu;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i + 1, j, k, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qm(i + 1, j, k, QREINT) * (gamma - 1.0);
  qm(i + 1, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);
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
  AMREX_D_DECL(
    amrex::Array4<const amrex::Real> const& q1,
    amrex::Array4<const amrex::Real> const& q2,
    amrex::Array4<const amrex::Real> const& q3),
  const amrex::Real dx, const amrex::Real dy, const amrex::Real dz)
{
	pdivu(i,j,k) = -(
				  (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx +
				  (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dy +
				  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dz);
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
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        // Passing in difmag TODO not pass in.
        flx(i, j, k, n) += dx * div1 * (u(i, j, k, n) - u(i - 1, j, k, n));
    }
  } else if (dir == 1) { // Y direction
    AMREX_D_PICK(, div1 = 0.5 * (div(i, j, k) + div(i + 1, j, k));
                 , div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) +
                                  div(i, j, k + 1) + div(i + 1, j, k + 1)););
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j - 1, k, n)); // Here dx is dy
    }
  } else { // Z direction
    div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) + div(i, j + 1, k) +
                   div(i + 1, j + 1, k));
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < flx.nComp(); ++n) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j, k - 1, n)); // Here dx is dz
    }
  }
}

// Host Functions
void pc_umeth_3D(
  amrex::Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  const int nq,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<amrex::Real> const& flx1,
  amrex::Array4<amrex::Real> const& flx2,
  amrex::Array4<amrex::Real> const& flx3,
  amrex::Array4<amrex::Real> const& q1,
  amrex::Array4<amrex::Real> const& q2,
  amrex::Array4<amrex::Real> const& q3,
  amrex::Array4<const amrex::Real> const& a1,
  amrex::Array4<const amrex::Real> const& a2,
  amrex::Array4<const amrex::Real> const& a3,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& vol,
  const int FirstSpec, const int NumSpec,
  const amrex::Real small_dens, const amrex::Real small_pres, 
  const amrex::Real small_vel , const amrex::Real small, 
  const amrex::Real gamma, const amrex::Real gamma_minus_1, 
  const amrex::Real* del,
  const amrex::Real dt);
#endif
