#ifndef _Derive_nohydro_H_
#define _Derive_nohydro_H_
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void dernull
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void dermaggrav
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void derdenvol
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

  void deroverden
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void deroverdenzoom
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

#ifdef __cplusplus
}
#endif

#endif
