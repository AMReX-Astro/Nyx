#include <iostream>

#include <Nyx.H>

#include <henson/context.h>
#include <henson/data.h>

int main()
{
    amrex::Amr* amr;
    henson_load_pointer("amr", (void**) &amr);

    std::cout << "sample-analysis: amr->cumTime() = " << amr->cumTime() << std::endl;
}
