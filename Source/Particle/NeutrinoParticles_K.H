#ifndef NEUTRINO_PARTICLES_K_H_
#define NEUTRINO_PARTICLES_K_H_

#include <NeutrinoParticleContainer.H>
#include <AMReX_REAL.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void neutrino_deposit_relativistic_cic (NeutrinoParticleContainer::ParticleType const& p,
                                        amrex::Array4<amrex::Real> const& rho,
                                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    using namespace amrex::literals;
    constexpr amrex::Real csq = 1.0;  // placeholder for now

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;

    amrex::Real sx[] = {1.0_rt - xint, xint};
    amrex::Real sy[] = {1.0_rt - yint, yint};
    amrex::Real sz[] = {1.0_rt - zint, zint};

    amrex::Real vel_x = p.rdata(1);
    amrex::Real vel_y = p.rdata(2);
    amrex::Real vel_z = p.rdata(3);
    amrex::Real vsq = vel_x*vel_x + vel_y*vel_y + vel_z*vel_z;
    amrex::Real gamma = 1.0 / std::sqrt ( 1.0 - vsq/csq);

    for (int kk = 0; kk <= 1; ++kk) {
        for (int jj = 0; jj <= 1; ++jj) {
            for (int ii = 0; ii <= 1; ++ii) {
                amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, 0),
                                        static_cast<amrex::Real>(sx[ii]*sy[jj]*sz[kk]*p.rdata(0))*gamma);
            }
        }
    }

    for (int comp=1; comp < AMREX_SPACEDIM; ++comp) {
        for (int kk = 0; kk <= 1; ++kk) {
            for (int jj = 0; jj <= 1; ++jj) {
                for (int ii = 0; ii <= 1; ++ii) {
                    amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, comp),
                                            static_cast<amrex::Real>(sx[ii]*sy[jj]*sz[kk]*p.rdata(0)*p.rdata(comp)));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void neutrino_deposit_particle_dx_relativistic_cic (
                                     NeutrinoParticleContainer::ParticleType const& p,
                                     amrex::Array4<amrex::Real> const& rho,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& pdxi)
{
    using namespace amrex::literals;
    amrex::Real factor = (pdxi[0]/dxi[0])*(pdxi[1]/dxi[1])*(pdxi[2]/dxi[2]);

    amrex::Real lx = (p.pos(0) - plo[0] - 0.5/pdxi[0]) * dxi[0];
    amrex::Real ly = (p.pos(1) - plo[1] - 0.5/pdxi[1]) * dxi[1];
    amrex::Real lz = (p.pos(2) - plo[2] - 0.5/pdxi[2]) * dxi[2];

    amrex::Real hx = (p.pos(0) - plo[0] + 0.5/pdxi[0]) * dxi[0];
    amrex::Real hy = (p.pos(1) - plo[1] + 0.5/pdxi[1]) * dxi[1];
    amrex::Real hz = (p.pos(2) - plo[2] + 0.5/pdxi[2]) * dxi[2];

    int lo_x = static_cast<int>(amrex::Math::floor(lx));
    int lo_y = static_cast<int>(amrex::Math::floor(ly));
    int lo_z = static_cast<int>(amrex::Math::floor(lz));

    int hi_x = static_cast<int>(amrex::Math::floor(hx));
    int hi_y = static_cast<int>(amrex::Math::floor(hy));
    int hi_z = static_cast<int>(amrex::Math::floor(hz));

    for (int k = lo_z; k <= hi_z; ++k) {
        if (k < rho.begin.z || k >= rho.end.z) continue;
        amrex::Real wz = amrex::min(hz - k, amrex::Real(1.0)) - amrex::max(lz - k, amrex::Real(0.0));
        for (int j = lo_y; j <= hi_y; ++j) {
            if (j < rho.begin.y || j >= rho.end.y) continue;
            amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
            for (int i = lo_x; i <= hi_x; ++i) {
                if (i < rho.begin.x || i >= rho.end.x) continue;
                amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));

                amrex::Real weight = wx*wy*wz*factor;

                rho(i, j, k, 0) += weight*p.rdata(0);
            }
        }
    }

    for (int comp = 1; comp < AMREX_SPACEDIM; ++comp) {
        for (int k = lo_z; k <= hi_z; ++k) {
            if (k < rho.begin.z || k >= rho.end.z) continue;
            amrex::Real wz = amrex::min(hz - k, amrex::Real(1.0)) - amrex::max(lz - k, amrex::Real(0.0));
            for (int j = lo_y; j <= hi_y; ++j) {
                if (j < rho.begin.y || j >= rho.end.y) continue;
                amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
                for (int i = lo_x; i <= hi_x; ++i) {
                    if (i < rho.begin.x || i >= rho.end.x) continue;
                    amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
                    amrex::Real weight = wx*wy*wz*factor;
                    rho(i, j, k, comp) += weight*p.rdata(0)*p.rdata(comp);
                }
            }
        }
    }
}

#endif
