#include "Nyx.H"
#include "Prob.H"

void prob_param_special_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{}

#ifndef NO_HYDRO
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
  // This is the case where we have compiled with states defined
  //  but they have only one component each so we fill them this way.
  if (state.nComp() == 1)
  {
    //Could be replaced with setVal
    state(i,j,k,0)    = 0.00;
  }
#ifndef NO_HYDRO
  // This is the regular case with NO_HYDRO = FALSE
  else if (state.nComp() > 1)
  {
    state(i,j,k,Density_comp)    = 0.00; //1.5d0 * small_dens
    state(i,j,k,Xmom_comp) = 0.00;
    state(i,j,k,Ymom_comp) = 0.00;
    state(i,j,k,Zmom_comp) = 0.00;

    // These will both be set later in the call to init_e.
    state(i,j,k,Eint_comp) = 0.0;
    state(i,j,k,Eden_comp) = 0.0;

#ifndef CONST_SPECIES
      state(i,j,k,FirstSpec_comp  ) = prob_param[ h_species_comp];
      state(i,j,k,FirstSpec_comp+1) = prob_param[he_species_comp];
#endif
  }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(const int i,
                   const int j,
                   const int k,
                   amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<amrex::Real> const& diag_eos,
                   amrex::GeometryData const& geomdata,
                   const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
  // This is the case where we have compiled with states defined
  //  but they have only one component each so we fill them this way.
  if (state.nComp() == 1 && diag_eos.nComp() == 1)
  {
    //Could be replaced with setVal
    diag_eos(i,j,k,0)    = 0.00;
  }
#ifndef NO_HYDRO
  // This is the regular case with NO_HYDRO = FALSE
  else if (state.nComp() > 1 && diag_eos.nComp() >= 2)
  {
    diag_eos(i,j,k,Temp_comp) = 0.0210*(1.00 + prob_param[z_in_comp])*(1.00 + prob_param[z_in_comp]);
    diag_eos(i,j,k,  Ne_comp) = 0.0;

    //Should be equivalent to inhomo_reion > 0 Nyx_setup.cpp
    if (diag_eos.nComp() > 2)
      diag_eos(i,j,k, Zhi_comp) = 7.5;
  }
#endif
}

void prob_initdata_on_box(const Box& bx,
                          Array4<amrex::Real> const& state,
                          Array4<amrex::Real> const& diag_eos,
                          GeometryData const& geomdata,
                          const GpuArray<Real,max_prob_param>& prob_param)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        prob_initdata_state(i, j ,k, state, geomdata, prob_param);
        prob_initdata      (i, j ,k, state, diag_eos, geomdata, prob_param);
    });
}

void prob_initdata_state_on_box(const Box& bx,
                                Array4<amrex::Real> const& state,
                                GeometryData const& geomdata,
                                const GpuArray<Real,max_prob_param>& prob_param)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        prob_initdata_state(i, j ,k, state, geomdata, prob_param);
    });
}
#endif
