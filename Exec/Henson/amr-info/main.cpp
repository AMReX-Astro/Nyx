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
        const MultiFab& mf = amr->getLevel(lev).get_old_data(PhiGrav_Type);       // TODO: might want a different way to specify the data type we want

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

            std::cout << "amr-info: " << box << " " << abox << std::endl;
        }
    }
}
