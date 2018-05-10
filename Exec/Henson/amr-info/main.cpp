#include <iostream>

#include <Nyx.H>

#include <henson/context.h>
#include <henson/data.h>

int main()
{
    using namespace amrex;

    Amr* amr;
    henson_load_pointer("amr", (void**) &amr);

    const Box&      domain     = amr->getLevel(0).Domain();
    const IntVect&  refinement = amr->getLevel(0).fineRatio();

    auto finest_level = amr->finestLevel();
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const MultiFab&                 mf = amr->getLevel(lev).get_old_data(PhiGrav_Type);       // TODO: might want a different way to specify the data type we want
        const BoxArray&                 ba = mf.boxArray();
        std::vector<std::pair<int,Box>> isects;
        const std::vector<IntVect>&     pshifts = amr->Geom(lev).periodicity().shiftIntVect();
        int                             ng = mf.nGrow();

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) // Loop over grids
        {
            // This is the valid Box of the current FArrayBox.
            // By "valid", we mean the original ungrown Box in BoxArray.
            const Box& box = mfi.validbox();

            // A reference to the current FArrayBox in this loop iteration.
            const FArrayBox& fab = mf[mfi];

            // Pointer to the floating point data of this FArrayBox.
            const Real* a = fab.dataPtr();

            // This is the Box on which the FArrayBox is defined.
            // Note that "abox" includes ghost cells (if there are any),
            // and is thus larger than or equal to "box".
            const Box& abox = fab.box();

            // We can now pass the information to a function that does
            // work on the region (specified by box) of the data pointed to
            // by Real* a.  The data should be viewed as a multidimensional
            // with bounds specified by abox.
            // Function f1 has the signature of
            // void f1(const int*, const int*, Real*, const int*, const int*);
            //f1(box.loVect(), box.hiVect(), a, abox.loVect(), abox.hiVect());

            std::cout << mfi.index() << ": " << box << " / " << abox << std::endl;

            // TODO: this only compute neighbors at the current level
            std::cout << "Neighbors:" << std::endl;
            Box gbx = grow(box,1);
            for (const auto& piv : pshifts)
            {
                ba.intersections(gbx + piv, isects);

                for (const auto& is : isects)
                {
                    // is.first is the index of neighbor box
                    // ba[is.first] is the neighbor box
                    const Box&  nbr_box         = ba[is.first];
                    Box         nbr_ghost_box   = grow(nbr_box,ng);
                    std::cout << "  " << is.first << ": " << nbr_box << " / " << nbr_ghost_box << std::endl;
                    //const Box& bx = is.second; // intersection of ungrown unshifted neighbor and
                    //                           // current box grown by 1 and shifted by piv
                    //const Box& bx2 = bx - piv;  // "unshift" intersection
                }
            }
        }
    }
}


