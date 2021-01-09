#ifdef NEUTRINO_PARTICLES
#ifndef _NeutrinoParticleContainer_H_
#define _NeutrinoParticleContainer_H_

#include <NyxParticleContainer.H>

class NeutrinoParticleContainer
    : public NyxParticleContainer<2+AMREX_SPACEDIM>
{

private:
    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    amrex::Real        m_csq;   // the square of the speed of light -- used to compute relativistic effects

public:
    NeutrinoParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<2+AMREX_SPACEDIM>(amr)
    { }

    using MyParIter = amrex::ParIter<2+AMREX_SPACEDIM>;
    using MyConstParIter = amrex::ParConstIter<2+AMREX_SPACEDIM>;

    virtual ~NeutrinoParticleContainer () {}

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (amrex::Real csq) { m_csq = csq; }

#ifndef NEUTRINO_DARK_PARTICLES
    void AssignDensity (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, 
                        int ncomp = 1, int finest_level = -1, int ngrow = 2) const
        {  AssignRelativisticDensity (mf,lev_min,ncomp,finest_level,ngrow); }

    void AssignRelativisticDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignRelativisticDensity (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, 
                                    int lev_min = 0, int ncomp = 1, int finest_level = -1, int ngrow = 2) const;
#endif

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_old = 1.0, amrex::Real a_half = 1.0, int where_width = 0);
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep,
                                amrex::Real a_new = 1.0, amrex::Real a_half = 1.0);
};

#endif /*_NeutrinoParticleContainer_H_*/
#endif /*NEUTRINO_PARTICLES*/
