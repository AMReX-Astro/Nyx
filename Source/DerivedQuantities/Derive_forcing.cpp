#include <AMReX_REAL.H>
#include <Derive.H>
#include <Nyx.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif


    void derforcex(const Box& /*bx*/, FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
                   const FArrayBox& /*datfab*/, const Geometry& /*geomdata*/,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
    }

    void derforcey(const Box& /*bx*/, FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
                   const FArrayBox& /*datfab*/, const Geometry& /*geomdata*/,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
    }

    void derforcez(const Box& /*bx*/, FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
                   const FArrayBox& /*datfab*/, const Geometry& /*geomdata*/,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
    }

#ifdef __cplusplus
}
#endif
