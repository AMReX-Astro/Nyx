#ifndef _Gravity_F_H_
#define _Gravity_F_H_
#include <AMReX_BLFort.H>

#ifdef CGRAV
BL_FORT_PROC_DECL(FORT_PRESCRIBE_GRAV,fort_prescribe_grav)
    (const int lo[], const int hi[],
     const amrex::Real* dx,
     const BL_FORT_FAB_ARG(S),
     const amrex::Real* problo,
     const int* addToSelfGrav);
#endif
#endif

