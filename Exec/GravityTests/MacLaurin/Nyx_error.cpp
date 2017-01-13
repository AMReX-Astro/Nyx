#include <AMReX_winstd.H>

#include "Nyx.H"
#include "Nyx_error_F.H"

void
Nyx::error_setup()
{
    // The lines below define routines to be called to tag cells for error
    // estimation -- the arguments of each "add" call are:
    //   1. Name of variable (state variable or derived quantity) which will be
    //      passed into the Fortran subroutine.
    //   2. Number of ghost cells each array needs in each call to the Fortran
    //      subroutine
    //   3. Type of Fortran subroutine -- this determines the argument list of
    //      the Fortran subroutine. These types are pre-defined and are
    //      currently restricted to `ErrorRec::Standard` and
    //      `ErrorRec::UseAverage`.
    //   4. Name of Fortran subroutine.

//  err_list.add("density", 0, ErrorRec::UseAverage,
//               BL_FORT_PROC_CALL(TAG_OVERDENSITY, tag_overdensity));
    err_list.add("density", 0, ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_REGION, tag_region));
}
