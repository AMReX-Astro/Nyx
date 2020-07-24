#include <iostream>

#include <Nyx.H>

#include <henson/context.h>
#include <henson/data.h>

int main()
{
    using namespace amrex;

    Amr* amr;
    DarkMatterParticleContainer* DMPC;
    henson_load_pointer("amr",  (void**) &amr);
    henson_load_pointer("dmpc", (void**) &DMPC);

    Vector<Real> locations;
    for (unsigned lev = 0; lev < DMPC->GetParticles().size(); lev++) {
        const auto& pmap = DMPC->GetParticles(lev);
        for (auto& kv : pmap) {
            const auto& pbx = kv.second.GetArrayOfStructs();
            for (const auto& p : pbx) {
                if (p.id() > 0) {
                    // Load positions
                    for (int d=0; d < BL_SPACEDIM; d++)
                        locations.push_back(p.m_rdata.pos[d]);
                }
            }
        }
    }

    //size_t total = ...;                 // FIXME: could get this via MPI_Reduce, but amr probably has this already
    //henson_save_size_t("total", total);

    size_t count = locations.size()/BL_SPACEDIM;    // local count
    henson_save_size_t("count", count);
    henson_save_array("x", &locations[0], sizeof(Real), count, BL_SPACEDIM*sizeof(Real));
    henson_save_array("y", &locations[1], sizeof(Real), count, BL_SPACEDIM*sizeof(Real));
    henson_save_array("z", &locations[2], sizeof(Real), count, BL_SPACEDIM*sizeof(Real));

    henson_yield();
}
