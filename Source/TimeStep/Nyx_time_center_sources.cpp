#include <Nyx.H>

using namespace amrex;

#ifndef NO_HYDRO

void
Nyx::time_center_source_terms (MultiFab& S_new,
                               MultiFab& ext_src_old,
                               MultiFab& ext_src_new,
                               Real      dt)
{
    BL_PROFILE("Nyx::time_center_source_terms()");

    // Subtract off half of the old source term, and add half of the new.
    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    Real a_old = get_comoving_a(prev_time);
    Real a_new = get_comoving_a(cur_time);

    int iden  = Density_comp;
    int ieint = Eint_comp;
    int ieden = Eden_comp;

    int ncomp = S_new.nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto s_arr = S_new.array(mfi);

        const auto src_old = ext_src_old.array(mfi);
        const auto src_new = ext_src_new.array(mfi);

        amrex::ParallelFor(bx, ncomp,
          [s_arr,src_old,src_new,a_old,a_new,iden,ieint,ieden,dt]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real a_half    = 0.5 * (a_old + a_new);
              Real a_newsq_inv   = 1.0 / (a_new*a_new);
              Real a_half_inv = 1.0 / a_half;
              Real a_new_inv  = 1.0 / a_new;

              int ixmom = iden+1;
              int izmom = iden+3;

              // Density
              if (n == iden)
              {
                  s_arr(i,j,k,n) += 0.5 * dt * (src_new(i,j,k,n) - src_old(i,j,k,n)) * a_half_inv;
              }
              // Momentum
              else if (n >= ixmom && n <= izmom)
              {
                  s_arr(i,j,k,n) += 0.5 * dt * (src_new(i,j,k,n) - src_old(i,j,k,n)) * a_new_inv;
              }

              // (rho e) and (rho E)
              else if (n == ieint)
              {
                  s_arr(i,j,k,ieint) += 0.5 * dt * a_half * (src_new(i,j,k,ieint) - src_old(i,j,k,ieint)) * a_newsq_inv;
                  s_arr(i,j,k,ieden) += 0.5 * dt * a_half * (src_new(i,j,k,ieden) - src_old(i,j,k,ieden)) * a_newsq_inv;
              }

              // Don't do anything here because we did this when we did eint
              else if (n == ieden)
              {
              }

              // (rho X_i) and (rho adv_i) and (rho aux_i)
              else
              {
                s_arr(i,j,k,n) += 0.50 * dt * (src_new(i,j,k,n) - src_old(i,j,k,n)) * a_half_inv;
              }
          });
    }

}
#endif
