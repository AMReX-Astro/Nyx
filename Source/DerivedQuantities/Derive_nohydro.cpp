#include <AMReX_REAL.H>

#include <Nyx.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

  void dernull(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datfab, const Geometry& geomdata,
               Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    // This routine is used by particle_count.  Yes it does nothing.

  }

    void dermaggrav(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = std::sqrt(dat(i,j,k,0)*dat(i,j,k,0) +
                                 dat(i,j,k,1)*dat(i,j,k,1) +
                                 dat(i,j,k,2)*dat(i,j,k,2));

      });
    }

    void derdenvol(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      // Here dat contains (Density, Xmom, Ymom, Zmom)
      const Real V_cell = dx[0] * dx[1] * dx[2];
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = V_cell * dat(i,j,k,0);

      });
    }

    void deroverden(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int level)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      // Here dat contains (Density, Xmom, Ymom, Zmom)
      const Real over_den = Nyx::average_total_density * std::pow(8,level+1);

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = dat(i,j,k,0) / over_den;

      });
    }

    void deroverdenzoom(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const FArrayBox& datfab, const Geometry& geomdata,
                        Real /*time*/, const int* /*bcrec*/, int level)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      //Assume Domain is a cube
      int idim = 0;
      int domlo = geomdata.Domain().smallEnd(idim);
      int domhi = geomdata.Domain().bigEnd(idim);

      int ref_size = domhi / (2*std::pow(2,(level+1)));
      int center   = (domhi-domlo+1) / 2;

      int ilo      = amrex::max(center-ref_size+1, bx.smallEnd(idim));
      int ihi      = amrex::min(center+ref_size,   bx.bigEnd(idim));
      auto const bx_ref = Box(IntVect(AMREX_D_DECL(amrex::max(center-ref_size+1, bx.smallEnd(0)),
                                                   amrex::max(center-ref_size+1, bx.smallEnd(1)),
                                                   amrex::max(center-ref_size+1, bx.smallEnd(2)))),
                              IntVect(AMREX_D_DECL(amrex::min(center+ref_size,   bx.bigEnd(0)),
                                                   amrex::min(center+ref_size,   bx.bigEnd(1)),
                                                   amrex::min(center+ref_size,   bx.bigEnd(2))) ));
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = 0.0;

      });
      amrex::ParallelFor(bx_ref,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = 1.0;

      });
      
    }

#ifdef __cplusplus
}
#endif
