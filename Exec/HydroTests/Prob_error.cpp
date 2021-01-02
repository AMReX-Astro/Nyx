#include "Nyx.H"
#include "Prob.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
    AMRErrorTagInfo info;
    info.SetMaxLevel(3);
    errtags.push_back(AMRErrorTag(3,AMRErrorTag::GREATER,"density",info));
    errtags.push_back(AMRErrorTag(0.01,AMRErrorTag::GRAD,"density",info));
    errtags.push_back(AMRErrorTag(3,AMRErrorTag::GREATER,"pressure",info));
    errtags.push_back(AMRErrorTag(0.01,AMRErrorTag::GRAD,"pressure",info));
}
