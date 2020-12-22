#ifndef _PLM_H_
#define _PLM_H_

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
plm_slope(
  const int i,
  const int j,
  const int k,
  const int n,
  const int dir,
  amrex::Array4<const amrex::Real> const& q,
  const amrex::Real small_pres)
{
  const int bdim[3] = {dir == 0, dir == 1, dir == 2};

  amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, drgtm1, dlftp1;
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
  for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++)
      flatn = amrex::min(flatn,flatten(i, j, k, dir_flat, q, small_pres));

  return flatn * dsgn * amrex::min(dlim, amrex::Math::abs(dtemp));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_x(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real      > const& qm,
  amrex::Array4<amrex::Real      > const& qp,
  amrex::Array4<amrex::Real const> const& srcQ,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  amrex::Real const a_old,
  amrex::Real const dx,
  amrex::Real const dt, 
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  amrex::Real const gamma_minus_1,
  amrex::Real const small_dens,
  amrex::Real const small_pres)
{
  amrex::Real dtdx = dt / (a_old * dx);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  // amrex::Real rhoe = q(i, j, k, QREINT);
  // amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  // amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - du) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + du) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  // amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0v = dv;
  amrex::Real alpha0w = dw;

  //  Construct the right state on the i-1/2 interface

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

  amrex::Real apright  = 0.5*(-1.0 - spplus )*alphap;
  amrex::Real amright  = 0.5*(-1.0 - spminus )*alpham;
  amrex::Real azrright = 0.5*(-1.0 - spzero )*alpha0r;
  // amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azv1rght = 0.5*(-1.0 - spzero )*alpha0v;
  amrex::Real azw1rght = 0.5*(-1.0 - spzero )*alpha0w;
  
  qp(i, j, k, QRHO  ) = rho + apright + amright + azrright;
  qp(i, j, k, QU    ) = u + (apright - amright) * cc / rho;
  qp(i, j, k, QV    ) = v + azv1rght;
  qp(i, j, k, QW    ) = w + azw1rght;
  // qp(i, j, k, QREINT) = rhoe + (apright + amright) * enth * cs + azeright;
  qp(i, j, k, QPRES ) = p + (apright + amright) * cs;

  // add the source terms
  qp(i, j, k, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qp(i, j, k, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qp(i, j, k, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qp(i, j, k, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qp(i, j, k, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qp(i, j, k, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qp(i, j, k, QRHO  ) = amrex::max(qp(i, j, k, QRHO ), small_dens);
  qp(i, j, k, QPRES ) = amrex::max(qp(i, j, k, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qp(i, j, k, QREINT ) = qp(i, j, k, QPRES) / gamma_minus_1;

  //  Construct the left state on the i+1/2 interface

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

  amrex::Real apleft   = 0.5*(1.0 - spplus )*alphap;
  amrex::Real amleft   = 0.5*(1.0 - spminus)*alpham;
  amrex::Real azrleft  = 0.5*(1.0 - spzero )*alpha0r;
  // amrex::Real azeleft  = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azv1left = 0.5*(1.0 - spzero )*alpha0v;
  amrex::Real azw1left = 0.5*(1.0 - spzero )*alpha0w;

  qm(i + 1, j, k, QRHO  ) = rho +  apleft + amleft + azrleft;
  qm(i + 1, j, k, QU    ) = u   + (apleft - amleft) * cc / rho;
  qm(i + 1, j, k, QV    ) = v   +  azv1left;
  qm(i + 1, j, k, QW    ) = w   +  azw1left;
  // qm(i + 1, j, k, QREINT) = rhoe + (apleft + amleft) * enth * cs + azeleft;
  qm(i + 1, j, k, QPRES ) = p   + (apleft + amleft) * cs;

  // add the source terms
  qm(i + 1, j, k, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qm(i + 1, j, k, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qm(i + 1, j, k, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qm(i + 1, j, k, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qm(i + 1, j, k, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qm(i + 1, j, k, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qm(i + 1, j, k, QRHO  ) = amrex::max(qm(i + 1, j, k, QRHO ), small_dens);
  qm(i + 1, j, k, QPRES ) = amrex::max(qm(i + 1, j, k, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qm(i + 1, j, k, QREINT ) = qm(i + 1, j, k, QPRES) / gamma_minus_1;

#ifndef CONST_SPECIES
  // Upwind the passive variables
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {

    u = q(i, j, k, QU);

    spzero = u > 0 ? -1.0 : u * dtdx;
    amrex::Real acmprght = 0.5 * (-1.0 - spzero) * slope[n];

    // Right state
    qp(i, j, k, n) = q(i, j, k, n) + acmprght + 0.5 * dt * srcQ(i, j, k, n) / a_old;

    spzero = u >= 0 ? u * dtdx : 1.0;
    amrex::Real acmpleft = 0.5 * (1.0 - spzero) * slope[n];

    // Left state
    qm(i + 1, j, k, n) = q(i, j, k, n) + acmpleft + 0.5 * dt * srcQ(i, j, k, n) / a_old;
  }
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_y(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real      > const& qm,
  amrex::Array4<amrex::Real      > const& qp,
  amrex::Array4<amrex::Real const> const& srcQ,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  amrex::Real const a_old,
  amrex::Real const dy,
  amrex::Real const dt,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  amrex::Real const gamma_minus_1,
  amrex::Real const small_dens,
  amrex::Real const small_pres)
{
  amrex::Real dtdy = dt / (a_old * dy);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  // amrex::Real rhoe = q(i, j, k, QREINT);
  // amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  // amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - dv) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + dv) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  // amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0u = du;
  amrex::Real alpha0w = dw;

  // Construct the "right" state on the j-1/2 interface
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
  // amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azu1rght = 0.5*(-1.0 - spzero )*alpha0u;
  amrex::Real azw1rght = 0.5*(-1.0 - spzero )*alpha0w;

  qp(i, j, k, QRHO  ) = rho +  apright + amright + azrright;
  qp(i, j, k, QU    ) = u   +  azu1rght;
  qp(i, j, k, QV    ) = v   + (apright - amright) * cc / rho;
  qp(i, j, k, QW    ) = w   +  azw1rght;
  // qp(i, j, k, QREINT) = rhoe + (apright + amright) * enth * cs + azeright;
  qp(i, j, k, QPRES ) = p   + (apright + amright) * cs;

  // add the source terms
  qp(i, j, k, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qp(i, j, k, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qp(i, j, k, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qp(i, j, k, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qp(i, j, k, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qp(i, j, k, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qp(i, j, k, QRHO  ) = amrex::max(qp(i, j, k, QRHO ), small_dens);
  qp(i, j, k, QPRES ) = amrex::max(qp(i, j, k, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qp(i, j, k, QREINT ) = qp(i, j, k, QPRES) / gamma_minus_1;

  // Construct the "left" state on the j+1/2 interface

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
  // amrex::Real azeleft = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azu1left = 0.5*(1.0 - spzero )*alpha0u;
  amrex::Real azw1left = 0.5*(1.0 - spzero )*alpha0w;

  qm(i, j + 1, k, QRHO  ) = rho + apleft + amleft + azrleft;
  qm(i, j + 1, k, QU    ) = u   + azu1left;
  qm(i, j + 1, k, QV    ) = v   + (apleft - amleft) * cc / rho;
  qm(i, j + 1, k, QW    ) = w   + azw1left;
  // qm(i, j + 1, k, QREINT) = rhoe + (apleft + amleft) * enth * cs + azeleft;
  qm(i, j + 1, k, QPRES ) = p   + (apleft + amleft) * cs;

  // add the source terms
  qm(i, j + 1, k, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qm(i, j + 1, k, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qm(i, j + 1, k, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qm(i, j + 1, k, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qm(i, j + 1, k, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qm(i, j + 1, k, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qm(i, j + 1, k, QRHO  ) = amrex::max(qm(i, j + 1, k, QRHO ), small_dens);
  qm(i, j + 1, k, QPRES ) = amrex::max(qm(i, j + 1, k, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qm(i, j + 1, k, QREINT ) = qm(i, j + 1, k, QPRES) / gamma_minus_1;

#ifndef CONST_SPECIES
  // Upwind the passive variables
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {

    // Right state
    v = q(i, j, k, QV);
    amrex::Real spzerot = v > 0 ? -1.0 : v * dtdy;
    amrex::Real acmptop = 0.5 * (-1.0 - spzerot) * slope[n];
    qp(i, j, k, n) = q(i, j, k, n) + acmptop + 0.5 * dt * srcQ(i, j, k, n) / a_old;

    amrex::Real spzerob = v >= 0 ? v * dtdy : 1.0;
    amrex::Real acmpbot = 0.5 * (1.0 - spzerob) * slope[n];
    qm(i, j + 1, k, n) = q(i, j, k, n) + acmpbot + 0.5 * dt * srcQ(i, j, k, n) / a_old;
  }
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_plm_z(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real      > const& qm,
  amrex::Array4<amrex::Real      > const& qp,
  amrex::Array4<amrex::Real const> const& srcQ,
  amrex::Real const slope[],
  amrex::Array4<const amrex::Real> const& q,
  amrex::Real const cc,
  amrex::Real const a_old,
  amrex::Real const dz,
  amrex::Real const dt,
#ifndef CONST_SPECIES
  const int NumSpec,
#endif
  amrex::Real const gamma_minus_1,
  amrex::Real const small_dens,
  amrex::Real const small_pres)
{
  amrex::Real dtdz = dt / (a_old*dz);
  amrex::Real cs = cc * cc;
  amrex::Real rho = q(i, j, k, QRHO);
  amrex::Real u = q(i, j, k, QU);
  amrex::Real v = q(i, j, k, QV);
  amrex::Real w = q(i, j, k, QW);
  amrex::Real p = q(i, j, k, QPRES);
  // amrex::Real rhoe = q(i, j, k, QREINT);
  // amrex::Real enth = ((rhoe + p) / rho) / cs;
  amrex::Real drho = slope[QRHO];
  amrex::Real du = slope[QU];
  amrex::Real dv = slope[QV];
  amrex::Real dw = slope[QW];
  amrex::Real dp = slope[QPRES];
  // amrex::Real drhoe = slope[QREINT];
  amrex::Real alpham = 0.5 * (dp / (rho * cc) - dw) * rho / cc;
  amrex::Real alphap = 0.5 * (dp / (rho * cc) + dw) * rho / cc;
  amrex::Real alpha0r = drho - dp / cs;
  // amrex::Real alpha0e = drhoe - dp * enth;
  amrex::Real alpha0u = du;
  amrex::Real alpha0v = dv;

  // Construct the "right" state on the k-1/2 interface

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
  // amrex::Real azeright = 0.5*(-1.0 - spzero )*alpha0e;
  amrex::Real azu1rght = 0.5*(-1.0 - spzero )*alpha0u;
  amrex::Real azv1rght = 0.5*(-1.0 - spzero )*alpha0v;

  qp(i, j, k, QRHO  ) = rho + apright + amright + azrright;
  qp(i, j, k, QU    ) = u   + azu1rght;
  qp(i, j, k, QV    ) = v   + azv1rght;
  qp(i, j, k, QW    ) = w   + (apright - amright) * cc / rho;
  // qp(i, j, k, QREINT) = rhoe + (apright + amright) * enth * cs + azeright;
  qp(i, j, k, QPRES ) = p   + (apright + amright) * cs;

  // add the source terms
  qp(i, j, k, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qp(i, j, k, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qp(i, j, k, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qp(i, j, k, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qp(i, j, k, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qp(i, j, k, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qp(i, j, k, QRHO  ) = amrex::max(qp(i, j, k, QRHO ), small_dens);
  qp(i, j, k, QPRES ) = amrex::max(qp(i, j, k, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qp(i, j, k, QREINT ) = qp(i, j, k, QPRES) / gamma_minus_1;

   // Construct the "left" state on the j+1/2 interface

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
  // amrex::Real azeleft = 0.5*(1.0 - spzero )*alpha0e;
  amrex::Real azu1left = 0.5*(1.0 - spzero )*alpha0u;
  amrex::Real azv1left = 0.5*(1.0 - spzero )*alpha0v;

  qm(i, j, k + 1, QRHO  ) = rho + apleft + amleft + azrleft;
  qm(i, j, k + 1, QU    ) = u   + azu1left;
  qm(i, j, k + 1, QV    ) = v   + azv1left;
  qm(i, j, k + 1, QW    ) = w   + (apleft - amleft) * cc / rho;
  // qm(i, j, k + 1, QREINT) = rhoe + (apleft + amleft) * enth * cs + azeleft;
  qm(i, j, k + 1, QPRES ) = p   + (apleft + amleft) * cs;

  // add the source terms
  qm(i, j, k + 1, QRHO  ) += 0.5*dt*srcQ(i, j, k, QRHO) / a_old;
  qm(i, j, k + 1, QU    ) += 0.5*dt*srcQ(i, j, k, QU) / a_old;
  qm(i, j, k + 1, QV    ) += 0.5*dt*srcQ(i, j, k, QV) / a_old;
  qm(i, j, k + 1, QW    ) += 0.5*dt*srcQ(i, j, k, QW) / a_old;
  // qm(i, j, k + 1, QREINT) += 0.5*dt*srcQ(i, j, k, QREINT) / a_old;
  qm(i, j, k + 1, QPRES ) += 0.5*dt*srcQ(i, j, k, QPRES) / a_old;

  qm(i, j, k + 1, QRHO  ) = amrex::max(qm(i, j, k + 1, QRHO ), small_dens);
  qm(i, j, k + 1, QPRES ) = amrex::max(qm(i, j, k + 1, QPRES), small_pres);

  // This allows the (rho e) to take advantage of (pressure > small_pres)
  qm(i, j, k + 1, QREINT ) = qm(i, j, k + 1, QPRES) / gamma_minus_1;

#ifndef CONST_SPECIES
  // Upwind the passive variables
  for (int n = FirstSpec_comp; n < FirstSpec_comp + NumSpec; ++n) {

    // Right state
    w = q(i, j, k, QW);
    amrex::Real spzerot = w > 0 ? -1.0 : w * dtdz;
    amrex::Real acmptop = 0.5 * (-1.0 - spzerot) * slope[n];
    qp(i, j, k, n) = q(i, j, k, n) + acmptop + 0.5 * dt * srcQ(i, j, k, n) / a_old;

    amrex::Real spzerob = w >= 0 ? w * dtdz : 1.0;
    amrex::Real acmpbot = 0.5 * (1.0 - spzerob) * slope[n];
    qm(i, j, k + 1, n) = q(i, j, k, n) + acmpbot + 0.5 * dt * srcQ(i, j, k, n) / a_old;
  }
#endif
}
#endif
