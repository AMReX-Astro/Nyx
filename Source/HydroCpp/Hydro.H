#ifndef _HYDRO_H_
#define _HYDRO_H_

#include "Nyx.H"
#include "IndexDefines.H"
#include "Utilities.H"
#include "PLM.H"

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

  // Assume srcrho == rhoinv
  // SrcQ Calc
  srcq(i, j, k, QRHO) = srcrho;
  srcq(i, j, k, QU) = (umx - a_dot * vx * srcrho) * rhoinv + gravx;
  srcq(i, j, k, QV) = (umy - a_dot * vy * srcrho) * rhoinv + gravy;
  srcq(i, j, k, QW) = (umz - a_dot * vz * srcrho) * rhoinv + gravz;

  //QREINT src depends on total and internal energy relationship
  /*
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy -
                          vz * umz +
                          srcrho * (vx * vx + vy * vy + vz * vz);
  */
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy -
                          vz * umz +
                          - a_dot * 3.0 * gamma_minus_1 * q(i , j, k, QREINT);
  /*
  if(i==1&&j==2&&k==2)
	  {
		  amrex::Print()<<srcq(i,j,k,QU)<<"\t"<<umx<<"\t"<<vx<<"\t"<<srcrho<<"\t"<<rhoinv<<std::endl;
		  exit(0			   );
	  }*/
  /*
             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
                  grav(i,j,k,1)
             srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
                  grav(i,j,k,2)
             srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
                  grav(i,j,k,3)
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
                  q(i,j,k,QV)*src(i,j,k,UMY) - &
                  q(i,j,k,QW)*src(i,j,k,UMZ) - &
                  a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)
  */
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
  AMREX_D_DECL(
    const amrex::Real dx, const amrex::Real dy, const amrex::Real dz),
  amrex::Array4<amrex::Real> const& divu)
{
#if AMREX_SPACEDIM == 1
  amrex::Real ux = (q(i, j, k, QU) - q(i - 1, j, k, QU)) / dx;
#elif AMREX_SPACEDIM == 2
  amrex::Real ux = 0.5 *
                   (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j - 1, k, QU) -
                    q(i - 1, j - 1, k, QU)) /
                   dx;
  amrex::Real vy = 0.5 *
                   (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i - 1, j, k, QV) -
                    q(i - 1, j - 1, k, QV)) /
                   dy;
#elif AMREX_SPACEDIM == 3
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
#endif
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
pc_update(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const dt,
  const amrex::Real gamma_minus_1)
{
  amrex::Real v = vol(i, j, k);
  for (int n = 0; n < update.nComp(); ++n) {
    update(i, j, k, n) += (AMREX_D_TERM(
                            flx[0](i, j, k, n) - flx[0](i + 1, j, k, n),
                            +flx[1](i, j, k, n) - flx[1](i, j + 1, k, n),
                            +flx[2](i, j, k, n) - flx[2](i, j, k + 1, n))) /
                          v;
  }

  update(i, j, k, UEINT) += pdivu(i, j, k)*u(i,j,k,UEINT)*gamma_minus_1*dt;
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
  amrex::Real const* del,
  amrex::Real const difmag,
  amrex::Real const dt,
  const int NumSpec,
  amrex::Real const gamma_minus_1);
#endif
