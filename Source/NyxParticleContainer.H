
#ifndef _NyxParticleContainer_H_
#define _NyxParticleContainer_H_

#include "Amr.H"
#include "AmrLevel.H"
#include "AmrParticles.H"

// We make DarkMatterParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class DarkMatterParticleContainer
    : public AmrParticleContainer<1+BL_SPACEDIM>
{
public:
    DarkMatterParticleContainer (Amr* amr)
        : AmrParticleContainer<1+BL_SPACEDIM>(amr)
    { }

    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass);
    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part,
                     const Real particleMass, const Real shift[]);
    void InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass);
    void InitCosmo1ppcMultiLevel(MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                                 const Real particleMass, int disp_idx, int vel_idx, 
                                 BoxArray &baWhereNot, int lev);

    void AssignDensityAndVels (PArray<MultiFab>& mf, int lev_min = 0) const;

};

// We make NeutrinoParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class NeutrinoParticleContainer
    : public AmrParticleContainer<2+BL_SPACEDIM>
{

    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    Real        m_csq;          // the square of the speed of light -- used to compute relativistic effects

public:
    NeutrinoParticleContainer (Amr* amr)
        : AmrParticleContainer<2+BL_SPACEDIM>(amr)
    { }

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (Real csq) { m_csq = csq; }

    void AssignDensity (PArray<MultiFab>& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const BL_OVERRIDE;

    void AssignRelativisticDensitySingleLevel (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignRelativisticDensity (PArray<MultiFab>& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

};

#endif /*_NyxParticleContainer_H_*/
