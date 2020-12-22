#ifndef _HYDRO_H_
#define _HYDRO_H_

#include <IndexDefines.H>
#include <Utilities.H>
#include <Godunov.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_srctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& grav,
  amrex::Array4<const amrex::Real> const& src,
  amrex::Array4<amrex::Real> const& srcQ,
  const amrex::Real a_dot,
  const int NumSpec, const amrex::Real gamma_minus_1)
{
  // Local loading Vars ==> Load once from global memory use multiple times
  const amrex::Real srcrho = src(i, j, k, Density_comp);
  const amrex::Real rhoinv = 1.0 / q(i, j, k, QRHO);
  const amrex::Real vx = q(i, j, k, QU);
  const amrex::Real vy = q(i, j, k, QV);
  const amrex::Real vz = q(i, j, k, QW);
  const amrex::Real umx = src(i, j, k, Xmom_comp);
  const amrex::Real umy = src(i, j, k, Ymom_comp);
  const amrex::Real umz = src(i, j, k, Zmom_comp);
  const amrex::Real gravx = grav(i, j, k, 0);
  const amrex::Real gravy = grav(i, j, k, 1);
  const amrex::Real gravz = grav(i, j, k, 2);

#if 1
  // srcrho != rhoinv
  // SrcQ Calc
  srcQ(i, j, k, QRHO) = srcrho;
  srcQ(i, j, k, QU) = umx * rhoinv - a_dot * vx + gravx;
  srcQ(i, j, k, QV) = umy * rhoinv - a_dot * vy + gravy;
  srcQ(i, j, k, QW) = umz * rhoinv - a_dot * vz + gravz;

  //QREINT src depends on total and internal energy relationship
  srcQ(i, j, k, QREINT) = src(i, j, k, Eden_comp) - vx * umx - vy * umy -
                          vz * umz +
                          - a_dot * 3.0 * gamma_minus_1 * q(i , j, k, QREINT);
  
  const amrex::Real dpde = gamma_minus_1 * q(i,j,k,QRHO);
  const amrex::Real dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO);

  srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoinv + dpdr * srcQ(i,j,k,QRHO);

  //Might not need correction from q
  for (int n = 0; n < NumSpec; ++n) {
    srcQ(i, j, k, FirstSpec_comp + n) = (src(i, j, k, FirstSpec_comp + n) - q(i, j, k, FirstSpec_comp + n) * srcrho) * rhoinv;
  }
#else
  srcQ(i, j, k, QRHO) = 0.;
  srcQ(i, j, k, QU)   = 0.;
  srcQ(i, j, k, QV)   = 0.;
  srcQ(i, j, k, QW)   = 0.;
  srcQ(i, j, k, QPRES)   = 0.;
  srcQ(i, j, k, QREINT)   = 0.;
  for (int n = 0; n < NumSpec; ++n) {
    srcQ(i, j, k, FirstSpec_comp + n) = 0.;
  }
#endif

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_divu(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  const amrex::Real dx, 
  const amrex::Real dy, 
  const amrex::Real dz,
  amrex::Array4<amrex::Real> const& divu)
{
  amrex::Real ux =
    0.25 *
    (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j, k - 1, QU) -
     q(i - 1, j, k - 1, QU) + q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU) +
     q(i, j - 1, k - 1, QU) - q(i - 1, j - 1, k - 1, QU)) /
    dx;
  amrex::Real vy =
    0.25 *
    (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i, j, k - 1, QV) -
     q(i, j - 1, k - 1, QV) + q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV) +
     q(i - 1, j, k - 1, QV) - q(i - 1, j - 1, k - 1, QV)) /
    dy;
  amrex::Real wz =
    0.25 *
    (q(i, j, k, QW) - q(i, j, k - 1, QW) + q(i - 1, j, k, QW) -
     q(i - 1, j, k - 1, QW) + q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW) +
     q(i - 1, j - 1, k, QW) - q(i - 1, j - 1, k - 1, QW)) /
    dz;

  divu(i, j, k) = AMREX_D_TERM(ux, +vy, +wz);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
normalize_species_fluxes(const int i, const int j, const int k, amrex::Array4<amrex::Real> const& flx, const int NumSpec)
{
  amrex::Real summ, fac;
  summ = 0.0;
  for (int n = 0; n < NumSpec; ++n)
    summ += flx(i, j, k, FirstSpec_comp + n);

  // We skip the normalization if the sum is zero or within epsilon.
  // There can be numerical problems here if the density flux is
  // approximately zero at the interface but not exactly, resulting in
  // division by a small number and/or resulting in one of the species
  // fluxes being negative because of roundoff error. There are also other
  // terms like artificial viscosity which can cause these problems.
  // So checking that sum is sufficiently large helps avoid this.


  fac = (std::abs(summ) > std::numeric_limits<amrex::Real>::epsilon() * std::abs(flx(i, j, k, Density_comp))) ? flx(i, j, k, Density_comp) / summ : 1.0;
  for (int n = 0; n < NumSpec; ++n)
    flx(i, j, k, FirstSpec_comp + n) *= fac;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ext_flx(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Real const area,
  amrex::Real const dt)
{
    for (int n = 0; n < flx.nComp(); ++n)
        flx(i,j,k,n) *= area*dt;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ext_flx_dt(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Real const a_old,
  amrex::Real const a_new)
{
  const amrex::Real a_half = 0.5*(a_old + a_new);
  const amrex::Real a_half_inv = 1.0 / a_half;
  const amrex::Real a_new_inv = 1.0 / a_new;
  const amrex::Real a_newsq_inv = 1.0 / (a_new * a_new);
  //Store all 8 factors
  const amrex::Real a_fact[8] = {a_half_inv,a_new_inv,a_new_inv,a_new_inv,a_half*a_newsq_inv,a_half*a_newsq_inv,a_half_inv,a_half_inv};

  for (int n = 0; n < flx.nComp(); ++n)
     flx(i,j,k,n)*=a_fact[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_update(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  amrex::Real const vol,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const a_old,
  amrex::Real const a_new,
  amrex::Real const dt,
  const amrex::Real gamma_minus_1)
{
  const amrex::Real volinv = 1 / vol;
  const amrex::Real a_half = 0.5*(a_old + a_new);
  const amrex::Real a_half_inv = 1.0 / a_half;
  //Store all 8 factors
  const amrex::Real a_fact[8] = {a_half_inv,1.0,1.0,1.0,a_half,a_half,a_half_inv,a_half_inv};
  for (int n = 0; n < update.nComp(); ++n) {
    update(i, j, k, n) += (  (flx[0](i, j, k, n) - flx[0](i + 1, j, k, n))
                            +(flx[1](i, j, k, n) - flx[1](i, j + 1, k, n))
                            +(flx[2](i, j, k, n) - flx[2](i, j, k + 1, n)) ) *
                          volinv * a_fact[n];
  }

  update(i, j, k, Eint_comp) += a_half * (a_new - a_old) * ( 2.0 - 3.0 * gamma_minus_1) * u(i,j,k,Eint_comp);
  update(i, j, k, Eden_comp) += a_half * (a_new - a_old) * ( 2.0 - 3.0 * gamma_minus_1) * u(i,j,k,Eint_comp);

  update(i, j, k, Eint_comp) += pdivu(i, j, k)*u(i,j,k,Eint_comp)*gamma_minus_1*dt*a_half;
}

// Host functions
void pc_umdrv(
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
  const int ppm_type);

void pc_consup(
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
  amrex::Real const difmag);
#endif
