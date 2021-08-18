#include "Nyx.H"
#include "Prob.H"

void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
    //Only include default tagging if NO_HYDRO=FALSE
#ifndef NO_HYDRO
    AMRErrorTagInfo info;
    errtags.push_back(AMRErrorTag(3.5e9,AMRErrorTag::GREATER,"denvol",info));
#endif
}
