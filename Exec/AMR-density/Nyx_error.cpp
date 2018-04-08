
#include "Nyx.H"
#include "Nyx_error_F.H"

using namespace amrex;

void
Nyx::error_setup()
{
    err_list.add("total_density",1,ErrorRec::UseAverage,
                 BL_FORT_PROC_CALL(TAG_OVERDENSITY, tag_overdensity));
}

void
Nyx::manual_tags_placement (TagBoxArray&    tags,
                            const Vector<IntVect>& bf_lev)
{
}
