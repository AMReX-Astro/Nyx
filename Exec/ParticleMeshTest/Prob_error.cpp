#include "Nyx.H"
#include "Prob.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
    amrex::AMRErrorTagInfo info;
    errtags.push_back(amrex::AMRErrorTag(1,amrex::AMRErrorTag::GREATER,"total_particle_count",info));
}
