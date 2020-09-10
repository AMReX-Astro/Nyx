#ifndef _Prob_H_
#define _Prob_H_

#include "Nyx.H"

enum Prob_Type_Index {
  z_in_comp = 0,
  gamma_comp,
  p_ambient_comp,
  dens_ambient_comp,
  exp_energy_comp,
  r_init_comp,
  nsub_comp
};

static void prob_param_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param, const amrex::Real initial_z)
{
    prob_param[z_in_comp] = initial_z;
    //Duplicate gamma-filling logic from Nyx.cpp
    prob_param[gamma_comp] = 5.0/3.0;;
    amrex::ParmParse pp("prob");
    pp.query("p_ambient", prob_param[p_ambient_comp]);
    pp.query("dens_ambient", prob_param[dens_ambient_comp]);
    pp.query("exp_energy", prob_param[exp_energy_comp]);
    pp.query("r_init", prob_param[r_init_comp]);
    pp.query("nsub", prob_param[nsub_comp]);
    amrex::ParmParse pp_nyx("nyx");
    pp_nyx.query("gamma", prob_param[gamma_comp]);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& geomdata,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
  // Set explosion pressure -- we will convert the point-explosion energy into
  // a corresponding pressure distributed throughout the perturbed volume
  amrex::Real vctr =
    4.0 / 3.0 * M_PI * (prob_param[r_init_comp] * prob_param[r_init_comp] * prob_param[r_init_comp]);
  amrex::Real p_exp = (prob_param[gamma_comp] - 1.0) * prob_param[exp_energy_comp] / vctr;

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  //Use left edge for min, and assume lo[0]=0:
  //      prob_lo[0] + (i) * dx[0] != prob_lo[0] + (i + 0.5) * dx[0]
  //      prob_lo[0] + (i) * dx[0] != prob_lo[0] + (i - lo[0]) * dx[0]
  AMREX_D_TERM(const amrex::Real xmin = prob_lo[0] + (i) * dx[0];
               , const amrex::Real ymin = prob_lo[1] + (j) * dx[1];
               , const amrex::Real zmin = prob_lo[2] + (k) * dx[2];)

  AMREX_D_TERM(
    amrex::Real dx_sub = dx[0] / static_cast<amrex::Real>(prob_param[nsub_comp]);
    amrex::Real ctr_x = 0.5 * (prob_lo[0] + prob_hi[0]);
    , amrex::Real dy_sub = dx[1] / static_cast<amrex::Real>(prob_param[nsub_comp]);
    amrex::Real ctr_y = 0.5 * (prob_lo[1] + prob_hi[2]);
    , amrex::Real dz_sub = dx[2] / static_cast<amrex::Real>(prob_param[nsub_comp]);
    amrex::Real ctr_z = 0.5 * (prob_lo[1] + prob_hi[2]);)

  // We initialize by summing over subvolumes of each cell
  int npert = 0;
  int n_amb = 0;

#if AMREX_SPACEDIM > 2
  for (int kk = 0; kk < prob_param[nsub_comp]; kk++) {
    const amrex::Real z = zmin + (kk + 0.5) * dz_sub - ctr_z;
#endif
#if AMREX_SPACEDIM > 1
    for (int jj = 0; jj < prob_param[nsub_comp]; jj++) {
      const amrex::Real y = ymin + (jj + 0.5) * dy_sub - ctr_y;
#endif
      for (int ii = 0; ii < prob_param[nsub_comp]; ii++) {
        const amrex::Real x = xmin + (ii + 0.5) * dx_sub - ctr_x;
        amrex::Real dist = AMREX_D_TERM(x * x, +y * y, +z * z);

        if (dist <= prob_param[r_init_comp] * prob_param[r_init_comp]) {
          npert += 1;
        } else {
          n_amb += 1;
        }
      }
#if AMREX_SPACEDIM > 1
    }
#endif
#if AMREX_SPACEDIM > 2
  }
#endif

  amrex::Real p_zone =
    (static_cast<amrex::Real>(npert) * p_exp +
     static_cast<amrex::Real>(n_amb) * prob_param[p_ambient_comp]) /
    static_cast<amrex::Real>(prob_param[nsub_comp] * prob_param[nsub_comp] * prob_param[nsub_comp]);

  amrex::Real eint = p_zone / (prob_param[gamma_comp] - 1.0);

  state(i, j, k, URHO) = prob_param[dens_ambient_comp];
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;

  state(i, j, k, UEDEN) = eint + 0.5 *
                                   (state(i, j, k, UMX) * state(i, j, k, UMX) +
                                    state(i, j, k, UMY) * state(i, j, k, UMY) +
                                    state(i, j, k, UMZ) * state(i, j, k, UMZ)) /
                                   state(i, j, k, URHO);

  state(i, j, k, UEINT) = eint;

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

  const int JH = 1;
  const int JHe = 1;

  const amrex::Real h_species_in=.76;
  const amrex::Real gamma_minus_1_in=prob_param[gamma_comp] - 1.0;
  const amrex::Real a=1.0/(prob_param[z_in_comp]+1.0);
  const amrex::Real rhoInv = 1.e0 / state(i,j,k,Density);
  const amrex::Real eint = state(i,j,k,Eint) * rhoInv;
  nyx_eos_T_given_Re_device(gamma_minus_1_in, h_species_in, JH, JHe, &diag_eos(i,j,k,Temp_comp), &diag_eos(i,j,k,Ne_comp), state(i,j,k,Density), eint, a);

    //Should be equivalent to inhomo_reion > 0 Nyx_setup.cpp
    if (diag_eos.nComp() > 2)
      diag_eos(i,j,k, Zhi_comp) = 7.5;
}

#endif
