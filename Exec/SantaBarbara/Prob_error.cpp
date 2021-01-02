#include "Nyx.H"
#include "Prob.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
#ifndef NO_HYDRO
    AMRErrorTagInfo info;
    errtags.push_back(AMRErrorTag(1,AMRErrorTag::GREATER,"overden",info));
#endif
}
