
#ifndef _DarkMatterParticleContainer_H_
#define _DarkMatterParticleContainer_H_

#include <NyxParticleContainer.H>

class DarkMatterParticleContainer
    : public NyxParticleContainer<1+AMREX_SPACEDIM>
{
public:
    DarkMatterParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<1+AMREX_SPACEDIM>(amr)
    {
      real_comp_names.clear();
      real_comp_names.push_back("mass");
      real_comp_names.push_back("xvel");
      real_comp_names.push_back("yvel");
      real_comp_names.push_back("zvel");
    }

    using MyParIter = amrex::ParIter<1+AMREX_SPACEDIM>;
    using MyConstParIter = amrex::ParConstIter<1+AMREX_SPACEDIM>;

    virtual ~DarkMatterParticleContainer () {}

    void InitCosmo1ppcMultiLevel(amrex::MultiFab& mf, const amrex::Real disp_fac[], const amrex::Real vel_fac[],
                                 const amrex::Real particleMass, int disp_idx, int vel_idx,
                                 amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void AssignDensityAndVels (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0) const;

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_old = 1.0, amrex::Real a_half = 1.0, int where_width = 0);
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_new = 1.0, amrex::Real a_half = 1.0);

    void InitFromBinaryMortonFile(const std::string& particle_directory, int nextra, int skip_factor);

};

AMREX_GPU_HOST_DEVICE AMREX_INLINE void update_dm_particle_single (amrex::ParticleContainer<4, 0>::SuperParticleType&  p,
                                                                   const int nc,
                                                                   amrex::Array4<amrex::Real const> const& acc,
                                                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                                                   const amrex::Real& dt, const amrex::Real& a_prev,
                                                                   const amrex::Real& a_cur, const int& do_move);

#endif /* _DarkMatterParticleContainer_H_ */
