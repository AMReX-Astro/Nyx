#include "AMReX_REAL.H"

#include "Derive.H"
#include <Nyx.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

  void derstate(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        // density
        der(i,j,k,0) = dat(i,j,k,0);

        // temperature
        der(i,j,k,1) = dat(i,j,k,1);

        if(der.nComp()>=3 && dat.nComp()>=3)
        {
            // (rho X)_1 = X_1
            der(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,0);
        }

      });

    }

    void dervel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
             der(i,j,k,0) = dat(i,j,k,1) / dat(i,j,k,0);
      });
    }



    void dermagvel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
        	   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = std::sqrt( (dat(i,j,k,1) / dat(i,j,k,0))*(dat(i,j,k,1) / dat(i,j,k,0)) +
				  (dat(i,j,k,2) / dat(i,j,k,0))*(dat(i,j,k,2) / dat(i,j,k,0)) +
				  (dat(i,j,k,3) / dat(i,j,k,0))*(dat(i,j,k,3) / dat(i,j,k,0)));
      });
    }

    /*    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real deninv = 1.0_rt/dat(i,j,k,0);

        der(i,j,k,0) = std::sqrt( (dat(i,j,k,1) * dat(i,j,k,1) +
                                   dat(i,j,k,2) * dat(i,j,k,2) +
                                   dat(i,j,k,3) * dat(i,j,k,3)) ) * deninv;
      });
    }*/


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

  void dermagmom(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
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

    void derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::Real gamma_minus_1_local = Nyx::gamma - 1.0_rt;
      //No abort call for non-positive energy
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = gamma_minus_1_local * dat(i,j,k,UEINT);

      });
    }


    void dereint1(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real rhoInv = 1.0_rt/dat(i,j,k,URHO);
        Real ux = dat(i,j,k,UMX)*rhoInv;
        Real uy = dat(i,j,k,UMY)*rhoInv;
        Real uz = dat(i,j,k,UMZ)*rhoInv;

        der(i,j,k,0) = dat(i,j,k,UEDEN)*rhoInv -
          0.5_rt * (ux*ux + uy*uy + uz*uz);
      });
    }

    void dereint2(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = dat(i,j,k,UEINT) / dat(i,j,k,URHO);
      });
    }

  void derkineng (const Box& bx, FArrayBox& kinengfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& /*geomdata*/,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    auto const dat = datfab.array();
    auto const kineng = kinengfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                         kineng(i,j,k,0) = 0.5_rt / dat(i,j,k,0) * ( dat(i,j,k,1)*dat(i,j,k,1) +
                                                                     dat(i,j,k,2)*dat(i,j,k,2) +
                                                                     dat(i,j,k,3)*dat(i,j,k,3) );
                       });

  }

  void dernull(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                  const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                  const int* lo, const int* hi,
                  const int* domain_lo, const int* domain_hi,
                  const Real* delta, const Real* xlo,
                  const Real* time, const Real* dt, const int* bcrec, 
                  const int* level, const int* grid_no)
  {

    // This routine is used by particle_count.  Yes it does nothing.

  }

  void derspec(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
             der(i,j,k,0) = dat(i,j,k,1) / dat(i,j,k,0);
      });
    }

#ifdef __cplusplus
}
#endif
