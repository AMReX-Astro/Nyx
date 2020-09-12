#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

enum Prob_Type_Index {
  z_in_comp = 0,
  gamma_comp,
  p_l_comp,
  u_l_comp,
  rho_l_comp,
  rhoe_l_comp,
  T_l_comp,
  p_r_comp,
  u_r_comp,
  rho_r_comp,
  rhoe_r_comp,
  T_r_comp,
  frac_comp,
  use_Tinit_comp,
  idir_comp
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param, const amrex::Real initial_z)
{
    prob_param[z_in_comp] = initial_z;
    //Duplicate gamma-filling logic from Nyx.cpp
    prob_param[gamma_comp] = 5.0/3.0;;

    prob_param[p_l_comp] = 1.0;  // left pressure (erg/cc)
    prob_param[u_l_comp] = 0.0; // left velocity (cm/s)
    prob_param[rho_l_comp] = 1.0; // left density (g/cc)

    prob_param[p_r_comp] = 0.1; // right pressure (erg/cc)
    prob_param[u_r_comp] = 0.0; // right velocity (cm/s)
    prob_param[rho_r_comp] = 0.125; // right density (g/cc)

    prob_param[idir_comp] = 1; // direction across which to jump
    prob_param[frac_comp] = 0.5;  // fraction of the domain for the interface

    // Parse params
    amrex::ParmParse pp("prob");
    pp.query("p_l",prob_param[p_l_comp]);
    pp.query("u_l",prob_param[u_l_comp]);
    pp.query("rho_l",prob_param[rho_l_comp]);
    pp.query("T_l",prob_param[T_l_comp]);
    pp.query("p_r",prob_param[p_r_comp]);
    pp.query("u_r",prob_param[u_r_comp]);
    pp.query("rho_r",prob_param[rho_r_comp]);
    pp.query("T_r",prob_param[T_r_comp]);
    pp.query("frac",prob_param[frac_comp]);
    pp.query("idir",prob_param[idir_comp]);
    pp.query("use_Tinit",prob_param[use_Tinit_comp]);
    amrex::ParmParse pp_nyx("nyx");
    pp_nyx.query("gamma", prob_param[gamma_comp]);
    // compute the internal energy (erg/cc) for the left and right state
    prob_param[rhoe_l_comp] = prob_param[p_l_comp]/(prob_param[gamma_comp] - 1.0);
    prob_param[rhoe_r_comp] = prob_param[p_r_comp]/(prob_param[gamma_comp] - 1.0);

}

static void prob_errtags_default(amrex::Vector<amrex::AMRErrorTag> errtags)
{
    AMRErrorTagInfo info;
    info.SetMaxLevel(3);
    errtags.push_back(AMRErrorTag(3,AMRErrorTag::GREATER,"density",info));
    errtags.push_back(AMRErrorTag(0.01,AMRErrorTag::GRAD,"density",info));
    errtags.push_back(AMRErrorTag(3,AMRErrorTag::GREATER,"pressure",info));
    errtags.push_back(AMRErrorTag(0.01,AMRErrorTag::GRAD,"pressure",info));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  //Use middle for comparison, and assume lo[0]=0:
  //      prob_lo[0] + (i+0.5) * dx[0] != prob_lo[0] + (float(i - lo[0])+0.5) * dx[0]
  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];)

  AMREX_D_TERM(
    amrex::Real ctr_x = prob_param[frac_comp] * (prob_lo[0] + prob_hi[0]);
    , amrex::Real ctr_y = prob_param[frac_comp] * (prob_lo[1] + prob_hi[1]);
    , amrex::Real ctr_z = prob_param[frac_comp] * (prob_lo[2] + prob_hi[2]);)

  // Set the states
  if (int(prob_param[idir_comp]) == 1) {
    if (x <= ctr_x) {
      state(i, j, k, URHO) = prob_param[rho_l_comp];
      state(i, j, k, UMX) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, URHO) = prob_param[rho_r_comp];
      state(i, j, k, UMX) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_r_comp];
    }
  } else if (int(prob_param[idir_comp]) == 2) {
    if (y <= ctr_y) {
      state(i, j, k, URHO) = prob_param[rho_l_comp];
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, URHO) = prob_param[rho_r_comp];
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_r_comp];
    }
      } else if (int(prob_param[idir_comp]) == 3) {
    if (z <= ctr_z) {
      state(i, j, k, URHO) = prob_param[rho_l_comp];
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, UEDEN) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, URHO) = prob_param[rho_r_comp];
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, UEDEN) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, UEINT) = prob_param[rhoe_r_comp];
    }
  } else {
    amrex::Abort("invalid idir");
  }

  //This uses 2 species, rather than just 1 (Doesn't set state(i, j, k, UFS) = state(i, j, k, URHO);)
    amrex::Real FirstSpec = Eint + 1;
    if (state.nComp() > FirstSpec)
    {
      //XHYDROGEN hardcoded to .76 init_managed
        state(i,j,k,FirstSpec  ) = .76 * state(i,j,k,URHO);
        state(i,j,k,FirstSpec+1) = (1.0 - .76) * state(i,j,k,URHO);
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

  if(int(prob_param[use_Tinit_comp])==0)
  {

  const int JH = 1;
  const int JHe = 1;

  const amrex::Real h_species_in=.76;
  const amrex::Real gamma_minus_1_in=prob_param[gamma_comp] - 1.0;
  const amrex::Real a=1.0/(prob_param[z_in_comp]+1.0);
  const amrex::Real rhoInv = 1.e0 / state(i,j,k,Density);
  const amrex::Real eint = state(i,j,k,Eint) * rhoInv;
  nyx_eos_T_given_Re_device(gamma_minus_1_in, h_species_in, JH, JHe, &diag_eos(i,j,k,Temp_comp), &diag_eos(i,j,k,Ne_comp), state(i,j,k,Density), eint, a);
  }
  else
  {
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  //Use middle for comparison, and assume lo[0]=0:
  //      prob_lo[0] + (i+0.5) * dx[0] != prob_lo[0] + (float(i - lo[0])+0.5) * dx[0]
  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];)

  AMREX_D_TERM(
    amrex::Real ctr_x = prob_param[frac_comp] * (prob_lo[0] + prob_hi[0]);
    , amrex::Real ctr_y = prob_param[frac_comp] * (prob_lo[1] + prob_hi[2]);
    , amrex::Real ctr_z = prob_param[frac_comp] * (prob_lo[1] + prob_hi[2]);)

  // Set the states
  if (int(prob_param[idir_comp]) == 1) {
    if (x <= ctr_x) {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_l_comp];
    } else {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_r_comp];
    }
  } else if (int(prob_param[idir_comp]) == 2) {
    if (y <= ctr_y) {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_l_comp];
    } else {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_r_comp];
    }
      } else if (int(prob_param[idir_comp]) == 3) {
    if (z <= ctr_z) {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_l_comp];
    } else {
      diag_eos(i,j,k,Temp_comp) = prob_param[T_r_comp];
    }
  } else {
    amrex::Abort("invalid idir");
  }
      diag_eos(i,j,k,Ne_comp) = 0.0;
  }
    //Should be equivalent to inhomo_reion > 0 Nyx_setup.cpp
    if (diag_eos.nComp() > 2)
      diag_eos(i,j,k, Zhi_comp) = 7.5;
}

#endif
