
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

};
#endif /*_NyxParticleContainer_H_*/
