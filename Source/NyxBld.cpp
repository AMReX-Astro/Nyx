
#include "AMReX_LevelBld.H"
#include "Nyx.H"

using namespace amrex;

class NyxBld
    :
    public LevelBld
{
    virtual void variable_setup();
    virtual void variable_setup_for_new_comp_procs();
    virtual void variable_cleanup();

    // hack copies for BoxLib overriding
    virtual void variableSetUp();
    virtual void variableSetUpForNewCompProcs();
    virtual void variableCleanUp();

    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr& papa, int lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba, 
                                  const DistributionMapping& dm, 
				  Real time);
};

NyxBld Nyx_bld;

LevelBld*
get_level_bld()
{
    return &Nyx_bld;
}

void
NyxBld::variable_setup()
{
    Nyx::variable_setup();
}

void
NyxBld::variable_setup_for_new_comp_procs()
{
    Nyx::variable_setup_for_new_comp_procs();
}

void
NyxBld::variable_cleanup()
{
    Nyx::variable_cleanup();
}

AmrLevel*
NyxBld::operator() ()
{
    return new Nyx;
}

AmrLevel*
NyxBld::operator() (Amr& papa, int lev, const Geometry& level_geom,
                    const BoxArray& ba, const DistributionMapping& dm, Real time)
{
    return new Nyx(papa, lev, level_geom, ba, dm, time);
}

// override hacks, copies of above
LevelBld*
getLevelBld()
{
    return &Nyx_bld;
}

void NyxBld::variableSetUp()
{
    Nyx::variable_setup();
}
void NyxBld::variableSetUpForNewCompProcs()
{
    Nyx::variable_setup_for_new_comp_procs();
}
void NyxBld::variableCleanUp()
{
    Nyx::variable_cleanup();
}
