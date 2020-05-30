#ifndef _Derive_H_
#define _Derive_H_
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void derpres
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dereint1
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dereint2
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dervel
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dermagvel
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dermaggrav
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dermagmom
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void derkineng (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                     amrex::Real time, const int* bcrec, int level);

  void dernull
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void derspec
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void derstate
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);


#ifdef __cplusplus
}
#endif

#endif
