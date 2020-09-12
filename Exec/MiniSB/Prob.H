#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

enum Prob_Type_Index {
  z_in_comp = 0,
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param, const amrex::Real initial_z)
{
    prob_param[z_in_comp] = initial_z;
}

static void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag> errtags)
{
    AMRErrorTagInfo info;
    errtags.push_back(AMRErrorTag(3.5e9,AMRErrorTag::GREATER,"denvol",info));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                           amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
    // HACK HACK HACK -- we need to get this from the ParmParse ....
    amrex::Real h_species = 0.75;

    state(i,j,k,Density)    = 0.00; //1.5d0 * small_dens
    state(i,j,k,Xmom) = 0.00;
    state(i,j,k,Ymom) = 0.00;
    state(i,j,k,Zmom) = 0.00;

    // These will both be set later in the call to init_e.
    state(i,j,k,Eint) = 0.0;
    state(i,j,k,Eden) = 0.0;

    amrex::Real FirstSpec = Eint + 1;   
    if (state.nComp() > FirstSpec)
    {
      state(i,j,k,FirstSpec  ) = h_species;
      state(i,j,k,FirstSpec+1) = (1.0 - h_species);
    }
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
  prob_initdata_state(i, j ,k, state, geomdata, prob_param);

  // This is the regular case with NO_HYDRO = FALSE
  {
    diag_eos(i,j,k,Temp_comp) = 1000.0;
    diag_eos(i,j,k,  Ne_comp) = 1.0;
  }
}

#endif
