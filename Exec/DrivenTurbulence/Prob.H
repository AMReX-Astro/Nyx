#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

enum Prob_Type_Index {
  z_in_comp = 0,
  gamma_comp = 1
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param, const amrex::Real initial_z)
{
    prob_param[z_in_comp] = initial_z;
    prob_param[gamma_comp] = 5.0/3.0;

    amrex::ParmParse pp_nyx("nyx");
    pp_nyx.query("gamma", prob_param[gamma_comp]);
}

static void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag>& errtags)
{
    AMRErrorTagInfo info;
    info.SetMaxLevel(0);
    errtags.push_back(AMRErrorTag(1.e18,AMRErrorTag::GREATER,"denerr",info));
    errtags.push_back(AMRErrorTag(2.e8,AMRErrorTag::GREATER,"dengrad",info));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{

  amrex::Real alpha=100.0;
  amrex::Real rho0=1.0;
  amrex::Real temp0=10.0;
  amrex::Real ne0=1.0;

  const amrex::Real h_species_in=1;
  const amrex::Real gamma_minus_1_in=prob_param[gamma_comp] - 1.0;
  const amrex::Real a=1.0/(prob_param[z_in_comp]+1.0);
  amrex::Real eint, dummy_pres;
  
  nyx_eos_given_RT(NULL,gamma_minus_1_in, h_species_in, &eint, &dummy_pres, state(i,j,k,Density), temp0, ne0, a);
  amrex::Real rhoe0 = rho0 * eint;

    state(i,j,k,Density)    = rho0; //1.5d0 * small_dens
    state(i,j,k,Xmom) = 0.00;
    state(i,j,k,Ymom) = 0.00;
    state(i,j,k,Zmom) = 0.00;

    // These will both be set later in the call to init_e.
    state(i,j,k,Eint) = rhoe0;
    state(i,j,k,Eden) = rhoe0 + 0.5 *
                                   (state(i, j, k, UMX) * state(i, j, k, UMX) +
                                    state(i, j, k, UMY) * state(i, j, k, UMY) +
                                    state(i, j, k, UMZ) * state(i, j, k, UMZ)) /
                                   state(i, j, k, URHO);

    int FirstSpec = Eint + 1;   
    if (state.nComp() > FirstSpec)
    {
      //XHYDROGEN hardcoded to .76 init_managed
      state(i,j,k,FirstSpec  ) = .76;
      state(i,j,k,FirstSpec+1) = (1.0 - .76);
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

  amrex::Real alpha=100.0;
  amrex::Real rho0=1.0;
  amrex::Real temp0=10.0;
  amrex::Real ne0=1.0;
  
  diag_eos(i,j,k,Temp_comp) = temp0;
  diag_eos(i,j,k,  Ne_comp) = ne0;
  
  prob_initdata_state(i, j ,k, state, geomdata, prob_param);

}

#endif
