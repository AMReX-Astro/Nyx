#include "AMReX_REAL.H"

#include "Derive.H"
#include <Nyx.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    void derforcex(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      // We only pass in the density here

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      Real twicePi = 2.0*M_PI;

#if 0
      Real Lx = prob_hi(1)-prob_lo(1);
      Real Ly = prob_hi(2)-prob_lo(2);
      Real Lz = prob_hi(3)-prob_lo(3);
    
      Real Lmin = amrex::min(Lx,Ly,Lz);
      Real kappaMax = static_cast<Real>(modes) / Lmin + 1.0e-8;
      int  nxmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lx/Lmin));
      int  nymodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Ly/Lmin));
      int  nzmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lz/Lmin));
    
      int xstep = int(Lx/Lmin+0.5);
      int ystep = int(Ly/Lmin+0.5);
      int zstep = int(Lz/Lmin+0.5);
    
      Real HLx = Lx;
      Real HLy = Ly;
      Real HLz = Lz;
#endif

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
#if 0
        Real x = (i+0.5) * dx[0];
        Real y = (j+0.5) * dx[1];
        Real z = (k+0.5) * dx[2];
               
        Real f1 = 0.0;
        for (int kz = mode_start*zstep; kz < nzmodes; kz += zstep)
        {
           Real kzd = static_cast<Real>(kz);
           Real freqz = twicePi*kzd*HLz;

           for (int ky = mode_start*ystep; ky < nymodes; ky += ystep)
           {
              Real kyd = static_cast<Real>(ky);
              Real freqy=twicePi*kyd/HLy;

              for (int kx = mode_start*ystep; kx < nxmodes; kx += xstep)
              {
                 Real kxd = static_cast<Real>(kx);
                 Real kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) );
                 Real freqx = twicePi*kxd/HLx;

                 if (kappa <= kappaMax) 
                 {
                    Real xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz))
                    f1 += xT * ( FAZ(kx,ky,kz)*freqy*sin(freqx*x+FPZX(kx,ky,kz)) * cos(freqy*y+FPZY(kx,ky,kz)) * 
                                                     sin(freqz*z+FPZZ(kx,ky,kz))
                                -FAY(kx,ky,kz)*freqz*sin(freqx*x+FPYX(kx,ky,kz)) * sin(freqy*y+FPYY(kx,ky,kz)) * 
                                                     cos(freqz*z+FPYZ(kx,ky,kz)) );
                 }
              }
           }
        }

        der(i,j,k,0) = dat(i,j,k,0)*f1;
#endif

      });
    }

    void derforcey(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      // We only pass in the density here

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      Real twicePi = 2.0*M_PI;

#if 0
      Real Lx = prob_hi(1)-prob_lo(1);
      Real Ly = prob_hi(2)-prob_lo(2);
      Real Lz = prob_hi(3)-prob_lo(3);
    
      Real Lmin = amrex::min(Lx,Ly,Lz);
      Real kappaMax = static_cast<Real>(modes) / Lmin + 1.0e-8;
      int  nxmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lx/Lmin));
      int  nymodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Ly/Lmin));
      int  nzmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lz/Lmin));
    
      int xstep = int(Lx/Lmin+0.5);
      int ystep = int(Ly/Lmin+0.5);
      int zstep = int(Lz/Lmin+0.5);
    
      Real HLx = Lx;
      Real HLy = Ly;
      Real HLz = Lz;
#endif

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real f2 = 0.0;
#if 0
        for (int kz = mode_start*zstep; kz < nzmodes; kz += zstep)
        {
           Real kzd = static_cast<Real>(kz);
           Real freqz = twicePi*kzd*HLz;

           for (int ky = mode_start*ystep; ky < nymodes; ky += ystep)
           {
              Real kyd = static_cast<Real>(ky);
              Real freqy=twicePi*kyd/HLy;

              for (int kx = mode_start*ystep; kx < nxmodes; kx += xstep)
              {
                 Real kxd = static_cast<Real>(kx);
                 Real kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) );
                 Real freqx = twicePi*kxd/HLx;

                 if (kappa <= kappaMax) 
                 {

                     if (kappa <= kappaMax) 
                     {
                        xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz));
                        f2 = f2 + xT * ( FAX(kx,ky,kz)*freqz*sin(freqx*x+FPXX(kx,ky,kz)) * sin(freqy*y+FPXY(kx,ky,kz)) *
                                                             cos(freqz*z+FPXZ(kx,ky,kz))
                                        -FAZ(kx,ky,kz)*freqx*cos(freqx*x+FPZX(kx,ky,kz)) * sin(freqy*y+FPZY(kx,ky,kz)) *
                                                             sin(freqz*z+FPZZ(kx,ky,kz)) );
                     }
                 }
               
               der(i,j,k,0) = dat(i,j,k,0)*f2;
              }
           }
        }
#endif
      });
}

    void derforcez(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      // We only pass in the density here

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      Real twicePi = 2.0*M_PI;

#if 0
      Real Lx = prob_hi(1)-prob_lo(1);
      Real Ly = prob_hi(2)-prob_lo(2);
      Real Lz = prob_hi(3)-prob_lo(3);
    
      Real Lmin = amrex::min(Lx,Ly,Lz);
      Real kappaMax = static_cast<Real>(modes) / Lmin + 1.0e-8;
      int  nxmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lx/Lmin));
      int  nymodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Ly/Lmin));
      int  nzmodes = nmodes* static_cast<int>(amrex::Math::floor(0.5+Lz/Lmin));
    
      int xstep = int(Lx/Lmin+0.5);
      int ystep = int(Ly/Lmin+0.5);
      int zstep = int(Lz/Lmin+0.5);
    
      Real HLx = Lx;
      Real HLy = Ly;
      Real HLz = Lz;
#endif

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real f3 = 0.0;
#if 0
        for (int kz = mode_start*zstep; kz < nzmodes; kz += zstep)
        {
           Real kzd = static_cast<Real>(kz);
           Real freqz = twicePi*kzd*HLz;

           for (int ky = mode_start*ystep; ky < nymodes; ky += ystep)
           {
              Real kyd = static_cast<Real>(ky);
              Real freqy=twicePi*kyd/HLy;

              for (int kx = mode_start*ystep; kx < nxmodes; kx += xstep)
              {
                 Real kxd = static_cast<Real>(kx);
                 Real kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) );
                 Real freqx = twicePi*kxd/HLx;

                 if (kappa <= kappaMax) 
                 {

                     if (kappa <= kappaMax) 
                     {
                         xT = cos(FTX(kx,ky,kz)*time+TAT(kx,ky,kz))
                         f3 = f3 + xT * ( FAY(kx,ky,kz)*freqx*cos(freqx*x+FPYX(kx,ky,kz)) * sin(freqy*y+FPYY(kx,ky,kz)) *
                                                              sin(freqz*z+FPYZ(kx,ky,kz))
                                         -FAX(kx,ky,kz)*freqy*sin(freqx*x+FPXX(kx,ky,kz)) * cos(freqy*y+FPXY(kx,ky,kz)) *
                                                              sin(freqz*z+FPXZ(kx,ky,kz)) );
                     }
                 }
               
               der(i,j,k,0) = dat(i,j,k,0)*f3;
              }
           }
        }
#endif
      });
}

#ifdef __cplusplus
}
#endif
