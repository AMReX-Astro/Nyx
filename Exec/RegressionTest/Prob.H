#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx,
                         amrex::Real z_in)
{
  // This is the case where we have compiled with states defined 
  //  but they have only one component each so we fill them this way.
  if (state.nComp() == 1)
  {
    //Could be replaced with setVal
    state(i,j,k,URHO)    = 1.00;
  }
  // This is the regular case with NO_HYDRO = FALSE
  else if (state.nComp() > 1)
  {
    state(i,j,k,Density)    = 1.00; //1.5d0 * small_dens
    state(i,j,k,Xmom) = 0.00;
    state(i,j,k,Ymom) = 0.00;
    state(i,j,k,Zmom) = 0.00;

    // These will both be set later in the call to init_e.
    state(i,j,k,Eint) = 10.0;
    state(i,j,k,Eden) = state(i,j,k,UEINT) + 0.5 / state(i,j,k,URHO) *
      (state(i,j,k,UMX)*state(i,j,k,UMX) + state(i,j,k,UMY)*state(i,j,k,UMY) + state(i,j,k,UMZ)*state(i,j,k,UMZ));

#ifndef CONST_SPECIES
      //XHYDROGEN hardcoded to .76
      state(i,j,k,FirstSpec  ) = .76;
      state(i,j,k,FirstSpec+1) = (1.0 - .76);
#endif
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(const int i,
                   const int j,
                   const int k,
                   amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<amrex::Real> const& diag_eos,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx,
                   amrex::Real z_in)
{

  prob_initdata_state(i, j ,k, state, dx,z_in);

  // This is the case where we have compiled with states defined 
  //  but they have only one component each so we fill them this way.
  if (state.nComp() == 1 && diag_eos.nComp() == 1)
  {
    //Could be replaced with setVal
    diag_eos(i,j,k,0)    = 0.00;
  }
  // This is the regular case with NO_HYDRO = FALSE
  else if (state.nComp() > 1 && diag_eos.nComp() >= 2)
  {
    diag_eos(i,j,k,Temp_comp) = 1000.0;
    diag_eos(i,j,k,  Ne_comp) = 0.0;

    //Should be equivalent to inhomo_reion > 0 Nyx_setup.cpp
    if (diag_eos.nComp() > 2)
      diag_eos(i,j,k, Zhi_comp) = 7.5;
  }
}

#endif
