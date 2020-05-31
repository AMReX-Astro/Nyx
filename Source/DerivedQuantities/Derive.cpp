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

    void dermagvort(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
        	    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSize();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real uy = (dat(i,j+1,k,1)/dat(i,j+1,k,0) - dat(i,j-1,k,1)/dat(i,j-1,k,0)) / dx[1];
        Real uz = (dat(i,j,k+1,1)/dat(i,j,k+1,0) - dat(i,j,k-1,1)/dat(i,j,k-1,0)) / dx[2];
        Real vx = (dat(i+1,j,k,2)/dat(i+1,j,k,0) - dat(i-1,j,k,2)/dat(i-1,j,k,0)) / dx[0];
        Real vz = (dat(i,j,k+1,2)/dat(i,j,k+1,0) - dat(i,j,k-1,2)/dat(i,j,k-1,0)) / dx[2];
        Real wx = (dat(i+1,j,k,3)/dat(i+1,j,k,0) - dat(i-1,j,k,3)/dat(i-1,j,k,0)) / dx[0];
        Real wy = (dat(i,j+1,k,3)/dat(i,j+1,k,0) - dat(i,j-1,k,3)/dat(i,j-1,k,0)) / dx[1];
 
        Real v1 = 0.5 * (wy - vz);  
        Real v2 = 0.5 * (uz - wx);  
        Real v3 = 0.5 * (vx - uy);  

        der(i,j,k,0) = std::sqrt(v1*v1 + v2*v2 + v3*v3); 
      });
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

        der(i,j,k,0) = gamma_minus_1_local * dat(i,j,k,Eint);

      });
    }

    void derlogden(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
             der(i,j,k,0) = std::log(dat(i,j,k,0));
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

        Real rhoInv = 1.0_rt/dat(i,j,k,Density);
        Real ux = dat(i,j,k,Xmom)*rhoInv;
        Real uy = dat(i,j,k,Ymom)*rhoInv;
        Real uz = dat(i,j,k,Zmom)*rhoInv;

        der(i,j,k,0) = dat(i,j,k,Eden)*rhoInv -
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

        der(i,j,k,0) = dat(i,j,k,Eint) / dat(i,j,k,Density);
      });
    }


    void dersoundspeed(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

	  Real sound_speed_factor=sqrt(Nyx::gamma*(Nyx::gamma-1.0_rt));

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real rhoInv = 1.0_rt/dat(i,j,k,Density);
        Real ux = dat(i,j,k,Xmom)*rhoInv;
        Real uy = dat(i,j,k,Ymom)*rhoInv;
        Real uz = dat(i,j,k,Zmom)*rhoInv;

		// Use internal energy for calculating dt 
		Real e  = dat(i,j,k,Eint)*rhoInv;

		// Protect against negative e
#ifdef HEATCOOL
		if (e > 0.0)
		    der(i,j,k,0)=sound_speed_factor*std::sqrt(e);
#else
		if (e > 0.0)
		    der(i,j,k,0)=sound_speed_factor*std::sqrt(dat(i,j,k,Density)*e/dat(i,j,k,Density));
#endif
		else
			der(i,j,k,0) = 0.0;

      });
    }


    void dermachnumber(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

	  Real sound_speed_factor=sqrt(Nyx::gamma*(Nyx::gamma-1.0_rt));

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real rhoInv = 1.0_rt/dat(i,j,k,Density);
        Real ux = dat(i,j,k,Xmom)*rhoInv;
        Real uy = dat(i,j,k,Ymom)*rhoInv;
        Real uz = dat(i,j,k,Zmom)*rhoInv;

		// Use internal energy for calculating dt 
		Real e  = dat(i,j,k,Eint)*rhoInv;

		Real c;
		// Protect against negative e
#ifdef HEATCOOL
		if (e > 0.0)
			c=sound_speed_factor*std::sqrt(e);
#else
		if (e > 0.0)
		    c=sound_speed_factor*std::sqrt(dat(i,j,k,Density)*e/dat(i,j,k,Density));
#endif
		else
			c = 0.0;

		der(i,j,k,0) = std::sqrt((ux*ux + uy*uy + uz*uz)) / c;

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

  void dernull(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datfab, const Geometry& geomdata,
               Real /*time*/, const int* /*bcrec*/, int /*level*/)
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

    void derdivu(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      // Here dat contains (Density, Xmom, Ymom, Zmom)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real rho = dat(i,j,k,0);

        Real uhi = dat(i+1,j,k,1) / rho;
        Real ulo = dat(i-1,j,k,1) / rho;
        Real vhi = dat(i,j+1,k,2) / rho;
        Real vlo = dat(i,j-1,k,2) / rho;
        Real whi = dat(i,j,k+1,3) / rho;
        Real wlo = dat(i,j,k-1,3) / rho;

        der(i,j,k,0) = 0.5 * ( (uhi-ulo) / dx[0] + (vhi-vlo) / dx[1] + (whi-wlo) / dx[2] );

      });
    }

#ifdef __cplusplus
}
#endif
