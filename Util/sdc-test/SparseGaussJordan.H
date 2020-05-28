#ifndef _SPARSE_GAUSS_JORDAN_H
#define _SPARSE_GAUSS_JORDAN_H
#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Extension.H"

class SparseGaussJordan {
public:
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  static void solve(Real* A, Real* x, Real* b) {

    Real scratch_0 = b[0]/A[0];

    x[0] = scratch_0;
    x[1] = scratch_0;
    x[2] = scratch_0;

  }
};

#endif
