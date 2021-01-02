#include "Nyx.H"
#include "Prob.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
    AMRErrorTagInfo info;
    info.SetMaxLevel(0);
    errtags.push_back(AMRErrorTag(1.e18,AMRErrorTag::GREATER,"denerr",info));
    errtags.push_back(AMRErrorTag(2.e8,AMRErrorTag::GREATER,"dengrad",info));
}
