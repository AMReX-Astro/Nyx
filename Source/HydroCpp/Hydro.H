#ifndef _HYDRO_H_
#define _HYDRO_H_

#include "IndexDefines.H"
#include "Utilities.H"
#include "Godunov.H"

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
  amrex::Array4<amrex::Real> const& srcq,
  const amrex::Real a_half, const amrex::Real a_dot,
  const int NumSpec, const amrex::Real gamma_minus_1)
{
  // Local loading Vars ==> Load once from global memory use multiple times
  const amrex::Real srcrho = src(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / q(i, j, k, QRHO);
  const amrex::Real vx = q(i, j, k, QU);
  const amrex::Real vy = q(i, j, k, QV);
  const amrex::Real vz = q(i, j, k, QW);
  const amrex::Real umx = src(i, j, k, UMX);
  const amrex::Real umy = src(i, j, k, UMY);
  const amrex::Real umz = src(i, j, k, UMZ);
  const amrex::Real gravx = grav(i, j, k, 0);
  const amrex::Real gravy = grav(i, j, k, 1);
  const amrex::Real gravz = grav(i, j, k, 2);

  // srcrho != rhoinv
  // SrcQ Calc
  srcq(i, j, k, QRHO) = srcrho;
  srcq(i, j, k, QU) = umx * rhoinv - a_dot * vx + gravx;
  srcq(i, j, k, QV) = umy * rhoinv - a_dot * vy + gravy;
  srcq(i, j, k, QW) = umz * rhoinv - a_dot * vz + gravz;

  //QREINT src depends on total and internal energy relationship
  /*
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy -
                          vz * umz +
                          srcrho * (vx * vx + vy * vy + vz * vz);
  */
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy -
                          vz * umz +
                          - a_dot * 3.0 * gamma_minus_1 * q(i , j, k, QREINT);
  
  const amrex::Real dpde = gamma_minus_1 * q(i,j,k,QRHO);
  const amrex::Real dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO);
  srcq(i,j,k,QPRES ) = dpde * srcq(i,j,k,QREINT) * rhoinv
          + dpdr * srcq(i,j,k,QRHO);

  //Might not need correction from q
  for (int n = 0; n < NumSpec; ++n) {
    srcq(i, j, k, FirstSpec + n) = (src(i, j, k, FirstSpec + n) - q(i, j, k, FirstSpec + n) * srcrho) * rhoinv;
  }

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
pc_norm_spec_flx(
                                 const int i, const int j, const int k, amrex::Array4<amrex::Real> const& flx, const int NumSpec)
{
  amrex::Real summ, fac;
  summ = 0.0;
  for (int n = 0; n < NumSpec; ++n)
    summ += flx(i, j, k, FirstSpec + n);

  fac = (summ != 0.0) ? flx(i, j, k, URHO) / summ : 1.0;
  for (int n = 0; n < NumSpec; ++n)
    flx(i, j, k, FirstSpec + n) *= fac;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ext_flx(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& a,
  amrex::Real const dt)
{
  amrex::Real area = a(i, j, k); // avoid reading from global memory flx1.nComp() times
  for (int n = 0; n < flx.nComp(); ++n)
          flx(i,j,k,n)*=area*dt;
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
  amrex::Real const a_new,
  amrex::Real const dt)
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
do_enforce_minimum_density(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& state,
  const int NumSpec,
  const amrex::Real small_dens)
{
    //Reset negative density to small_dens and zero out momenta:
    if(state(i,j,k,URHO) < small_dens)
    {
        for (int n = 0; n < NumSpec; ++n)
        {
            state(i,j,k,FirstSpec+n) *= small_dens / state(i,j,k,URHO);
        }

        state(i,j,k,URHO ) = small_dens;

        state(i,j,k,UMX) = 0.0;
        state(i,j,k,UMY) = 0.0;
        state(i,j,k,UMZ) = 0.0;

    /*
    Consider resetting UEINT here as well, like negative e in reset_internal_energy:
    Real rhoInv = 1.0 / u(i,j,k,URHO);
    Real Up     = u(i,j,k,UMX) * rhoInv;
    Real Vp     = u(i,j,k,UMY) * rhoInv;
    Real Wp     = u(i,j,k,UMZ) * rhoInv;
    Real ke     = 0.5 * u(i,j,k,URHO) * (Up*Up + Vp*Vp + Wp*Wp);

    Real rho_eint = u(i,j,k,UEDEN) - ke;
    Real dummy_pres = 1e200;
    Real eint_new = -1e200;
    nyx_eos_given_RT(gamma_minus_1, h_species, &eint_new, &dummy_pres, u(i,j,k,URHO), small_temp, 
                     d(i,j,k,Ne_comp),comoving_a);

    // Create reset source so u(i,j,k,UEINT) = u(i,j,k,UEINT) + r(i,j,k) = u(i,j,k,URHO) * eint_new
    r(i,j,k) = u(i,j,k,URHO) *  eint_new - u(i,j,k,UEINT);
    u(i,j,k,UEINT) = u(i,j,k,URHO) *  eint_new;

    u(i,j,k,UEDEN) = u(i,j,k,UEINT) + ke;
    //  state(i,j,k,UEINT) = state(i,j,k,UEINT);
    //  state(i,j,k,UEDEN) = state(i,j,k,UEINT);
        */
	}
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
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const a_old,
  amrex::Real const a_new,
  amrex::Real const dt,
  const amrex::Real gamma_minus_1)
{
  const amrex::Real volinv = 1 / vol(i, j, k);
  const amrex::Real a_half = 0.5*(a_old + a_new);
  const amrex::Real a_half_inv = 1.0 / a_half;
  //Store all 8 factors
  const amrex::Real a_fact[8] = {a_half_inv,1.0,1.0,1.0,a_half,a_half,a_half_inv,a_half_inv};
  for (int n = 0; n < update.nComp(); ++n) {
    update(i, j, k, n) += (AMREX_D_TERM(
                            flx[0](i, j, k, n) - flx[0](i + 1, j, k, n),
                            +flx[1](i, j, k, n) - flx[1](i, j + 1, k, n),
                            +flx[2](i, j, k, n) - flx[2](i, j, k + 1, n))) *
                          volinv * a_fact[n];
  }

  update(i, j, k, UEINT) += pdivu(i, j, k)*u(i,j,k,UEINT)*gamma_minus_1*dt*a_half;
}

// Host functions
void pc_umdrv(
  const int is_finest_level,
  const amrex::Real time,
  amrex::Box const& bx,
  const int* dlo,
  const int* dhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::Real* dx,
  const amrex::Real dt,
  const amrex::Real a_old,
  const amrex::Real a_new,
  const amrex::Real gamma, const amrex::Real gamma_minus_1, 
  const int NumSpec,
  const amrex::Real small_dens, const amrex::Real small_pres, 
  const amrex::Real small_vel , const amrex::Real small, 
  amrex::Real cflLoc,
  const int ppm_type,
  const int use_flattening,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<amrex::Real> const& vol);

void pc_consup(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<const amrex::Real> const& vol,
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
