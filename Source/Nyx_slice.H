#ifndef _slice_util_H_
#define _slice_util_H_

#include <AMReX_MultiFab.H>

namespace slice_util
{
    std::unique_ptr<amrex::MultiFab> getSliceData(int dir, const amrex::MultiFab& cell_centered_data, 
                                                  int fstart, int ncomp, const amrex::Geometry& geom, amrex::Real dir_coord);
}
#endif
