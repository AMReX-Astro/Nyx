#include "IndexDefines.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
         const amrex::Real x[AMREX_SPACEDIM],
         const amrex::Real s_int[QVAR],
         amrex::Real s_ext[QVAR],
         const int idir,
         const int sgn,
         const amrex::Real time,
         amrex::GeometryData const& geomdata)
{
    for (int n = 0; n < QVAR; n++) {
        s_ext[n] = s_int[n];
    }
}

void nyx_bcfill (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);
