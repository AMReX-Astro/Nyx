#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>
#include <math.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <AMReX_CONSTANTS.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <Nyx_slice.H>

using namespace amrex;

namespace slice_util
{

    static Box 
    getIndexBox(const RealBox& real_box, const Geometry& geom) {
        IntVect slice_lo, slice_hi;
    
        D_TERM(slice_lo[0]=floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_lo[1]=floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_lo[2]=floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
        D_TERM(slice_hi[0]=floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_hi[1]=floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_hi[2]=floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
        return Box(slice_lo, slice_hi) & geom.Domain();
    }    

    static 
    std::unique_ptr<MultiFab> allocateSlice(int dir, const MultiFab& cell_centered_data, 
                                            int ncomp, const Geometry& geom, Real dir_coord,
                                            Array<int>& slice_to_full_ba_map) {
 
        // Get our slice and convert to index space
        RealBox real_slice = geom.ProbDomain();
        real_slice.setLo(dir, dir_coord);
        real_slice.setHi(dir, dir_coord);
        Box slice_box = getIndexBox(real_slice, geom);
        
        // define the multifab that stores slice
        BoxArray ba = cell_centered_data.boxArray();
        const DistributionMapping& dm = cell_centered_data.DistributionMap();
        std::vector< std::pair<int, Box> > isects;
        ba.intersections(slice_box, isects, false, 0);
        Array<Box> boxes;
        Array<int> procs;
        for (int i = 0; i < isects.size(); ++i) {
            procs.push_back(dm[isects[i].first]);
            boxes.push_back(isects[i].second);
            slice_to_full_ba_map.push_back(isects[i].first);
        }
        procs.push_back(ParallelDescriptor::MyProc());
        BoxArray slice_ba(&boxes[0], boxes.size());
        DistributionMapping slice_dmap(procs);
        std::unique_ptr<MultiFab> slice(new MultiFab(slice_ba, slice_dmap, ncomp, 0));
        return slice;
    }
    
    std::unique_ptr<MultiFab> getSliceData(int dir, const MultiFab& cell_centered_data, 
                                           int fstart, int ncomp, const Geometry& geom, Real dir_coord) {
        
        Array<int> slice_to_full_ba_map;
        std::unique_ptr<MultiFab> slice = allocateSlice(dir, cell_centered_data, ncomp, geom, dir_coord,
                                                        slice_to_full_ba_map);

        // Fill the slice with sampled data
        int nf = cell_centered_data.nComp();
        const BoxArray& ba = cell_centered_data.boxArray();
        for (MFIter mfi(*slice); mfi.isValid(); ++mfi) {
            int slice_gid = mfi.index();
            int full_gid = slice_to_full_ba_map[slice_gid];
            
            const Box& slice_box = mfi.validbox();
            const Box& full_box  = cell_centered_data[full_gid].box();
            const Box& tile_box  = mfi.tilebox();

            fill_slice(cell_centered_data[full_gid].dataPtr(),
                       full_box.loVect(), full_box.hiVect(),
                       &fstart, &nf, 
                       (*slice)[slice_gid].dataPtr(),
                       slice_box.loVect(), slice_box.hiVect(),
                       tile_box.loVect(), tile_box.hiVect(), &ncomp);
        }
        
        return slice;
    }
}
