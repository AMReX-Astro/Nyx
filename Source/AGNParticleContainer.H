
#ifndef _AGNParticleContainer_H_
#define _AGNParticleContainer_H_

#include "NyxParticleContainer.H"

// We make AGNParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class AGNParticleContainer
    : public NyxParticleContainer<1+BL_SPACEDIM+3>
{
public:
    AGNParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<1+BL_SPACEDIM+3>(amr)
    { }

    virtual ~AGNParticleContainer () {}
};

#endif /* _AGNParticleContainer_H_ */
