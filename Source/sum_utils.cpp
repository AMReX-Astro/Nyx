#include <iostream>
#include <iomanip>

#include <Nyx.H>
#include <Nyx_F.H>

using namespace amrex;

Real
Nyx::vol_weight_sum (const std::string& name,
                     Real               time,
                     bool               masked)
{
    Real        sum = 0;
    const Real* dx  = geom.CellSize();
    auto        mf  = derive(name, time, 0);

    BL_PROFILE("vol_weight_sum(name)");
    BL_ASSERT(mf != 0);

    if (masked)
    {
        int flev = parent->finestLevel();
        while (parent->getAmrLevels()[flev] == nullptr)
            flev--;

        if (level < flev)
        {
            Nyx* fine_level = dynamic_cast<Nyx*>(&(parent->getLevel(level+1)));
            const MultiFab* mask = fine_level->build_fine_mask();
            MultiFab::Multiply(*mf, *mask, 0, 0, 1, 0);
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sum)
#endif
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real       s;
        const Box& box = mfi.tilebox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

         sum_over_level
            (BL_TO_FORTRAN(fab), lo, hi, dx, &s);
        sum += s;
    }

    ParallelDescriptor::ReduceRealSum(sum);

    if (!masked) 
        sum /= geom.ProbSize();

    return sum;
}

Real
Nyx::vol_weight_sum (MultiFab& mf, bool masked)
{
    BL_PROFILE("vol_weight_sum");
    Real        sum = 0;
    const Real* dx  = geom.CellSize();

    MultiFab* mask = 0;
    BL_PROFILE_VAR("Nyx::vol_weight_sum()::mask",volmask);
    if (masked)
    {
        int flev = parent->finestLevel();
        while (parent->getAmrLevels()[flev] == nullptr) flev--;

        if (level < flev)
        {
            Nyx* fine_level = dynamic_cast<Nyx*>(&(parent->getLevel(level+1)));
            mask = fine_level->build_fine_mask();
        }
    }

    BL_PROFILE_VAR_STOP(volmask);

#ifdef AMREX_USE_CUDA
        if ( !masked || (mask == 0) )
        {
	  //	  sum = mf.sum(0,true)*dx[0]*dx[1]*dx[2]
	  BL_PROFILE("vol_weight_nomask");
	  sum = mf.sum(0)*dx[0]*dx[1]*dx[2];
        }
        else
        {
	  BL_PROFILE("vol_weight_mask");
	  //	  sum = amrex::MultiFab::Dot(mf,0,(*mask),0,1,0,true)*dx[0]*dx[1]*dx[2];
	  sum = amrex::MultiFab::Dot(mf,0,(*mask),0,1,0)*dx[0]*dx[1]*dx[2];
        }

#else
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sum)
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (mf)[mfi];

        Real       s;
        const Box& box = mfi.tilebox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        if ( !masked || (mask == 0) )
        {
            sum_over_level
                (BL_TO_FORTRAN(fab), lo, hi, dx, &s);
        }
        else
        {
            FArrayBox& fab2 = (*mask)[mfi];
            sum_product
                (BL_TO_FORTRAN(fab), BL_TO_FORTRAN(fab2), lo, hi, dx, &s);
        }

        sum += s;
    }

    ParallelDescriptor::ReduceRealSum(sum);
#endif
    if (!masked) 
        sum /= geom.ProbSize();

    return sum;
}

Real
Nyx::vol_weight_squared_sum_level (const std::string& name,
                                   Real               time)
{
    Real        sum = 0;
    const Real* dx  = geom.CellSize();
    auto        mf  = derive(name, time, 0);
    BL_PROFILE("vol_weight_squared_sum_level");
    BL_ASSERT(mf != 0);

    Real lev_vol = parent->boxArray(level).d_numPts() * dx[0] * dx[1] * dx[2];

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sum)
#endif
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real       s;
        const Box& box = mfi.tilebox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        sum_product
            (BL_TO_FORTRAN(fab), BL_TO_FORTRAN(fab), lo, hi, dx, &s);
        sum += s;
    }

    ParallelDescriptor::ReduceRealSum(sum);

    sum /= lev_vol;

    return sum;
}

Real
Nyx::vol_weight_squared_sum (const std::string& name,
                             Real               time)
{

    BL_PROFILE("vol_weight_squared_sum");
    Real        sum = 0;
    const Real* dx  = geom.CellSize();
    auto        mf  = derive(name, time, 0);

    BL_ASSERT(mf != 0);

    int flev = parent->finestLevel();
    while (parent->getAmrLevels()[flev] == nullptr) flev--;

    MultiFab* mask = 0;
    if (level < flev)
    {
        Nyx* fine_level = dynamic_cast<Nyx*>(&(parent->getLevel(level+1)));
        mask = fine_level->build_fine_mask();
    }

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sum)
#endif
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real       s;
        const Box& box = mfi.tilebox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        if (mask == 0)
        {
            sum_product
                (BL_TO_FORTRAN(fab), BL_TO_FORTRAN(fab), 
                 lo, hi, dx, &s);
        }
        else
        {
            FArrayBox& fab2 = (*mask)[mfi];
            sum_prod_prod
                (BL_TO_FORTRAN(fab), BL_TO_FORTRAN(fab), BL_TO_FORTRAN(fab2), 
                 lo, hi, dx, &s);
        }

        sum += s;
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

MultiFab*
Nyx::build_fine_mask()
{
    BL_ASSERT(level > 0); // because we are building a mask for the coarser level

    if (fine_mask != 0) return fine_mask;

    BoxArray baf = parent->boxArray(level);
    baf.coarsen(crse_ratio);

    const BoxArray& bac = parent->boxArray(level-1);
    fine_mask = new MultiFab(bac,parent->DistributionMap(level-1), 1,0);
    fine_mask->setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*fine_mask); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*fine_mask)[mfi];

        std::vector< std::pair<int,Box> > isects = baf.intersections(fab.box());

        for (int ii = 0; ii < isects.size(); ii++)
        {
            fab.setVal<RunOn::Host>(0.0,isects[ii].second,0);
        }
    }
/*
    if (fine_mask.empty()) {
        fine_mask = makeFineMask(parent->boxArray(level-1),
                                 parent->DistributionMap(level-1),
                                 parent->boxArray(level), crse_ratio,
                                 1.0,  // coarse
                                 0.0); // fine
    }
*/
    return fine_mask;
}
