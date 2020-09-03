#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <AMReX_FArrayBox.H>
#include "IndexDefines.H"

/*
AMREX_GPU_DEVICE
void pc_cmpTemp(
  const int i, const int j, const int k, 
  const int FirstSpec, const int NumSpec,
  amrex::Array4<amrex::Real> const& S);
*/

AMREX_GPU_DEVICE
void pc_rst_int_e(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& S);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& q,
  const int NumSpec, const amrex::Real gamma_minus_1)
{
  const amrex::Real rho = u(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / rho;
  const amrex::Real vx = u(i, j, k, UMX) * rhoinv;
  const amrex::Real vy = u(i, j, k, UMY) * rhoinv;
  const amrex::Real vz = u(i, j, k, UMZ) * rhoinv;
  const amrex::Real e = u(i, j, k, UEINT) * rhoinv;
  //  const amrex::Real kineng = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
  q(i, j, k, QRHO) = rho;
  q(i, j, k, QU) = vx;
  q(i, j, k, QV) = vy;
  q(i, j, k, QW) = vz;
  q(i,j,k,QREINT) = e;

  // Convert "e" back to "rho e"
  q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO);
  
  // Pressure = (gamma - 1) * rho * e
  q(i,j,k,QPRES) = gamma_minus_1 * q(i,j,k,QREINT);

  for (int n = 0; n < NumSpec; ++n) {
    q(i, j, k, FirstSpec + n) = u(i, j, k, FirstSpec + n) / rho;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_move_transcoefs_to_ec(
  const int i,
  const int j,
  const int k,
  const int n,
  const amrex::Array4<const amrex::Real>& carr,
  amrex::Real* earr,
  const int dir,
  const int do_harmonic)

{
  if (do_harmonic == 0) {
    if (dir == 0)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i - 1, j, k, n));
    else if (dir == 1)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j - 1, k, n));
    else
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j, k - 1, n));
  } else {
    amrex::Real a = carr(i, j, k, n);
    if (dir == 0) {
      amrex::Real b = carr(i - 1, j, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else if (dir == 1) {
      amrex::Real b = carr(i, j - 1, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else {
      amrex::Real b = carr(i, j, k - 1, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    }
  }
}

AMREX_FORCE_INLINE
void
setV(
  const amrex::Box box,
  const int Ncomp,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      array(i, j, k, n) = val;
    });
}

AMREX_FORCE_INLINE
void
setC(
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      array(i, j, k, n) = val;
  });
}

AMREX_FORCE_INLINE
void
copy_array4(
  const amrex::Box box,
  const int Ncomp,
  const amrex::Array4<const amrex::Real>& in,
  const amrex::Array4<amrex::Real>& out)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      out(i, j, k, n) = in(i, j, k, n);
    });
}

AMREX_FORCE_INLINE
void
lincomb_array4(
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real alpha,
  const amrex::Real beta,
  const amrex::Array4<amrex::Real>& d)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      d(i, j, k, n) = alpha * f1(i, j, k, n) + beta * f2(i, j, k, n);
  });
}

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void read_binary(
  const std::string iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<amrex::Real>& data);

void read_csv(
  const std::string iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data);

AMREX_GPU_HOST_DEVICE
void locate(const amrex::Real* xtable, const int n, amrex::Real& x, int& idxlo);

#endif
