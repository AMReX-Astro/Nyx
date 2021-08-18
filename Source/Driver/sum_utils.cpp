#include <iostream>
#include <iomanip>

#include <Nyx.H>

using namespace amrex;

Real
Nyx::vol_weight_sum (const std::string& name,
                     Real               time,
                     bool               masked)
{
    BL_PROFILE("vol_weight_sum(name)");

    auto        mf  = derive(name, time, 0);

    Real sum = vol_weight_sum(*mf, masked);
    return sum;
}

Real
Nyx::vol_weight_sum (MultiFab& mf, bool masked)
{
    BL_PROFILE("vol_weight_sum");

    const auto dx  = geom.CellSizeArray();

    MultiFab* mask = 0;
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

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    if ( !masked || (mask == 0) )
    {
        BL_PROFILE("Nyx::vol_weight_sum()::ReduceOpsOnDevice");
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto fab = mf.array(mfi);
            const Box& tbx = mfi.tilebox();
  
            reduce_op.eval(tbx, reduce_data,
            [fab]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                 Real x = fab(i,j,k);
                 return x;
            });
        }

    } else {
        BL_PROFILE("Nyx::vol_weight_sum()::ReduceOpsOnDevice");
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto fab = mf.array(mfi);
            const auto msk = mask->array(mfi);
            const Box& tbx = mfi.tilebox();
  
            reduce_op.eval(tbx, reduce_data,
            [fab,msk]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                 Real x = fab(i,j,k)*msk(i,j,k);
                 return x;
            });
        }
    }

    ReduceTuple hv = reduce_data.value();
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(hv));

    Real sum = get<0>(hv) * (dx[0] * dx[1] * dx[2]);

    if (!masked) 
        sum /= geom.ProbSize();

    return sum;
}

Real
Nyx::vol_weight_squared_sum_level (const std::string& name,
                                   Real               time)
{
    BL_PROFILE("vol_weight_squared_sum_level");

    auto        mf  = derive(name, time, 0);
    AMREX_ASSERT(mf != 0);

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
#endif
    for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto fab = mf->array(mfi);
        const Box& tbx = mfi.tilebox();

        reduce_op.eval(tbx, reduce_data,
        [fab]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
             Real x = fab(i,j,k)*fab(i,j,k);
             return x;
        });
    }
    
    ReduceTuple hv = reduce_data.value();
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(hv));

    Real sum = get<0>(hv) / grids.d_numPts();

    return sum;
}

Real
Nyx::vol_weight_squared_sum (const std::string& name,
                             Real               time)
{

    BL_PROFILE("vol_weight_squared_sum");

    auto        mf  = derive(name, time, 0);
    AMREX_ASSERT(mf != 0);

    const auto dx = geom.CellSizeArray();

    bool masked = true;

    MultiFab* mask = 0;
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

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    if ( !masked || (mask == 0) )
    {
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
#endif
        for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto fab = mf->array(mfi);
            const Box& tbx = mfi.tilebox();
  
            reduce_op.eval(tbx, reduce_data,
            [fab]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                 Real x = fab(i,j,k)*fab(i,j,k);
                 return x;
            });
        }

    } else {
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
#endif
        for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto fab = mf->array(mfi);
            const auto msk = mask->array(mfi);
            const Box& tbx = mfi.tilebox();
  
            reduce_op.eval(tbx, reduce_data,
            [fab,msk]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                 Real x = fab(i,j,k)*fab(i,j,k)*msk(i,j,k);
                 return x;
            });
        }
    }

    ReduceTuple hv = reduce_data.value();
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(hv));

    Real sum = get<0>(hv) * (dx[0] * dx[1] * dx[2]);

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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*fine_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*fine_mask)[mfi];

        std::vector< std::pair<int,Box> > isects = baf.intersections(fab.box());

        for (int ii = 0; ii < isects.size(); ii++)
        {
            fab.setVal<RunOn::Device>(0.0,isects[ii].second,0);
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
