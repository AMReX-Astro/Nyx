#ifndef NO_HYDRO
#include <Nyx.H>
#include <constants_cosmo.H>

#include <Gravity.H>

using namespace amrex;

void
Nyx::correct_gsrc(int lev, Real time, Real prev_time, Real cur_time, Real dt)
{
    const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
    const auto& dm = get_level(lev).get_new_data(State_Type).DistributionMap();

    // These vectors are only used for the call to correct_gsrc so they 
    //    don't need any ghost cells
    MultiFab grav_vec_old(ba, dm, AMREX_SPACEDIM, 0);
    MultiFab grav_vec_new(ba, dm, AMREX_SPACEDIM, 0);

    get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
    get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);

    MultiFab& S_old = get_level(lev).get_old_data(State_Type);
    MultiFab& S_new = get_level(lev).get_new_data(State_Type);

    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        grav_vec_old[mfi].prefetchToDevice();
        grav_vec_new[mfi].prefetchToDevice();

        S_old[mfi].prefetchToDevice();
        S_new[mfi].prefetchToDevice();

        const auto g_old = grav_vec_old.array(mfi);
        const auto g_new = grav_vec_new.array(mfi);
        const auto s_old = S_old.array(mfi);
        const auto s_new = S_new.array(mfi);

        int  iden = Density_comp; 
        int ieden = Eden_comp; 

        amrex::ParallelFor(bx, [g_old,g_new,s_old,s_new,a_old,a_new,iden,ieden,dt]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real a_half    = 0.50 * (a_old + a_new);
            Real a_newsq   = a_new*a_new;
            Real a_new_inv = 1.00 / a_new;

            int ixmom = iden+1; 
            int iymom = iden+2; 
            int izmom = iden+3; 

            Real rhoo    = s_old(i,j,k,iden);
            Real rhooinv = 1.0 / rhoo;
            Real Upo     = s_old(i,j,k,ixmom) * rhooinv;
            Real Vpo     = s_old(i,j,k,iymom) * rhooinv;
            Real Wpo     = s_old(i,j,k,izmom) * rhooinv;

            // Define old source terms
            Real SrU_old = rhoo * g_old(i,j,k,0);
            Real SrV_old = rhoo * g_old(i,j,k,1);
            Real SrW_old = rhoo * g_old(i,j,k,2);

            Real rhon    = s_new(i,j,k,iden);
            Real rhoninv = 1.0 / rhon;
            Real Upn     = s_new(i,j,k,ixmom) * rhoninv;
            Real Vpn     = s_new(i,j,k,iymom) * rhoninv;
            Real Wpn     = s_new(i,j,k,izmom) * rhoninv;

            // Define new source terms
            Real SrU_new = rhon * g_new(i,j,k,0);
            Real SrV_new = rhon * g_new(i,j,k,1);
            Real SrW_new = rhon * g_new(i,j,k,2);

            // Define corrections to source terms
            Real SrUcorr = 0.5*(SrU_new - SrU_old);
            Real SrVcorr = 0.5*(SrV_new - SrV_old);
            Real SrWcorr = 0.5*(SrW_new - SrW_old);
            
            // Correct state with correction terms
            s_new(i,j,k,ixmom) += SrUcorr*dt * a_new_inv;
            s_new(i,j,k,iymom) += SrVcorr*dt * a_new_inv;
            s_new(i,j,k,izmom) += SrWcorr*dt * a_new_inv;

            Real SrEcorr =  0.5 * ( (SrU_new * Upn - SrU_old * Upo) +
                                    (SrV_new * Vpn - SrV_old * Vpo) +
                                    (SrW_new * Wpn - SrW_old * Wpo) );
            s_new(i,j,k,ieden) += SrEcorr * dt * (a_half / a_newsq);
         });
    }
}
#endif
