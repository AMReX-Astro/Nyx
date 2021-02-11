#include "Nyx.H"
#include "Prob.H"

enum Prob_Type_Index_Sod {
  p_l_comp = 6, // Note this must be one greater than the final index in Prob_param.H
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

enum Prob_Type_Index_Sedov {
  p_ambient_comp = 6, // Note this must be one greater than the final index in Prob_param.H
  dens_ambient_comp,
  exp_energy_comp,
  r_init_comp,
  nsub_comp
};

void prob_param_special_fill(amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{
    if(amrex::Math::round(prob_param[prob_type_comp]) == 0)
    {
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
    // compute the internal energy (erg/cc) for the left and right state
    prob_param[rhoe_l_comp] = prob_param[p_l_comp]/(prob_param[gamma_comp] - 1.0);
    prob_param[rhoe_r_comp] = prob_param[p_r_comp]/(prob_param[gamma_comp] - 1.0);
    }
    else
    {
    prob_param[gamma_comp] = 5.0/3.0;;
    prob_param[r_init_comp] = 0.01;
    prob_param[p_ambient_comp] = 1.e-5;
    prob_param[dens_ambient_comp] = 1.0;
    prob_param[exp_energy_comp] = 1.0;
    prob_param[nsub_comp] = 10;
    amrex::ParmParse pp("prob");
    pp.query("p_ambient", prob_param[p_ambient_comp]);
    pp.query("dens_ambient", prob_param[dens_ambient_comp]);
    pp.query("exp_energy", prob_param[exp_energy_comp]);
    pp.query("r_init", prob_param[r_init_comp]);
    pp.query("nsub", prob_param[nsub_comp]);
    }
}

#ifndef NO_HYDRO
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

  if(amrex::Math::round(prob_param[prob_type_comp]) == 0)
  {
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
  if (amrex::Math::round(prob_param[idir_comp]) == 1) {
    if (x <= ctr_x) {
      state(i, j, k, Density_comp) = prob_param[rho_l_comp];
      state(i, j, k, Xmom_comp) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, Ymom_comp) = 0.0;
      state(i, j, k, Zmom_comp) = 0.0;
      state(i, j, k, Eden_comp) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, Density_comp) = prob_param[rho_r_comp];
      state(i, j, k, Xmom_comp) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, Ymom_comp) = 0.0;
      state(i, j, k, Zmom_comp) = 0.0;
      state(i, j, k, Eden_comp) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_r_comp];
    }
  } else if (amrex::Math::round(prob_param[idir_comp]) == 2) {
    if (y <= ctr_y) {
      state(i, j, k, Density_comp) = prob_param[rho_l_comp];
      state(i, j, k, Xmom_comp) = 0.0;
      state(i, j, k, Ymom_comp) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, Zmom_comp) = 0.0;
      state(i, j, k, Eden_comp) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, Density_comp) = prob_param[rho_r_comp];
      state(i, j, k, Xmom_comp) = 0.0;
      state(i, j, k, Ymom_comp) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, Zmom_comp) = 0.0;
      state(i, j, k, Eden_comp) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_r_comp];
    }
      } else if (amrex::Math::round(prob_param[idir_comp]) == 3) {
    if (z <= ctr_z) {
      state(i, j, k, Density_comp) = prob_param[rho_l_comp];
      state(i, j, k, Xmom_comp) = 0.0;
      state(i, j, k, Ymom_comp) = 0.0;
      state(i, j, k, Zmom_comp) = prob_param[rho_l_comp] * prob_param[u_l_comp];
      state(i, j, k, Eden_comp) = prob_param[rhoe_l_comp] + 0.5 * prob_param[rho_l_comp] *
                                                   prob_param[u_l_comp] *
                                                   prob_param[u_l_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_l_comp];
    } else {
      state(i, j, k, Density_comp) = prob_param[rho_r_comp];
      state(i, j, k, Xmom_comp) = 0.0;
      state(i, j, k, Ymom_comp) = 0.0;
      state(i, j, k, Zmom_comp) = prob_param[rho_r_comp] * prob_param[u_r_comp];
      state(i, j, k, Eden_comp) = prob_param[rhoe_r_comp] + 0.5 * prob_param[rho_r_comp] *
                                                   prob_param[u_r_comp] *
                                                   prob_param[u_r_comp];
      state(i, j, k, Eint_comp) = prob_param[rhoe_r_comp];
    }
  } else {
    amrex::Abort("invalid idir");
  }
  }
  else
  {
  // Set explosion pressure -- we will convert the point-explosion energy into
  // a corresponding pressure distributed throughout the perturbed volume
  amrex::Real vctr =
    4.0 / 3.0 * M_PI * (prob_param[r_init_comp] * prob_param[r_init_comp] * prob_param[r_init_comp]);
  amrex::Real p_exp = (prob_param[gamma_comp] - 1.0) * prob_param[exp_energy_comp] / vctr;

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
    amrex::Real ctr_y = 0.5 * (prob_lo[1] + prob_hi[1]);
    , amrex::Real dz_sub = dx[2] / static_cast<amrex::Real>(prob_param[nsub_comp]);
    amrex::Real ctr_z = 0.5 * (prob_lo[2] + prob_hi[2]);)

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

  state(i, j, k, Density_comp) = prob_param[dens_ambient_comp];
  state(i, j, k, Xmom_comp) = 0.0;
  state(i, j, k, Ymom_comp) = 0.0;
  state(i, j, k, Zmom_comp) = 0.0;

  state(i, j, k, Eden_comp) = eint + 0.5 *
                                   (state(i, j, k, Xmom_comp) * state(i, j, k, Xmom_comp) +
                                    state(i, j, k, Ymom_comp) * state(i, j, k, Ymom_comp) +
                                    state(i, j, k, Zmom_comp) * state(i, j, k, Zmom_comp)) /
                                    state(i, j, k, Density_comp);

  state(i, j, k, Eint_comp) = eint;

  }

#ifndef CONST_SPECIES
     state(i,j,k,FirstSpec_comp  ) = prob_param[ h_species_comp];
     state(i,j,k,FirstSpec_comp+1) = prob_param[he_species_comp];
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

  prob_initdata_state(i, j ,k, state, geomdata, prob_param);

  if ( amrex::Math::round(prob_param[prob_type_comp]) != 0 ||
      (amrex::Math::round(prob_param[prob_type_comp]) == 0 &&
       amrex::Math::round(prob_param[use_Tinit_comp]) == 0) )
  {

      const int JH = 1;
      const int JHe = 1;

      const amrex::Real h_species_in = prob_param[h_species_comp];
      const amrex::Real gamma_minus_1_in = prob_param[gamma_comp] - 1.0;
      const amrex::Real a = 1.0/(prob_param[z_in_comp]+1.0);
      const amrex::Real rhoInv = 1.e0 / state(i,j,k,Density_comp);
      const amrex::Real eint = state(i,j,k,Eint_comp) * rhoInv;

      nyx_eos_T_given_Re_device(NULL,gamma_minus_1_in, h_species_in, JH, JHe,
                               &diag_eos(i,j,k,Temp_comp), &diag_eos(i,j,k,Ne_comp),
                               state(i,j,k,Density_comp), eint, a);
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
                 , amrex::Real ctr_y = prob_param[frac_comp] * (prob_lo[1] + prob_hi[1]);
                 , amrex::Real ctr_z = prob_param[frac_comp] * (prob_lo[2] + prob_hi[2]);)

      // Set the states
      if (amrex::Math::round(prob_param[idir_comp]) == 1) {
        if (x <= ctr_x) {
          diag_eos(i,j,k,Temp_comp) = prob_param[T_l_comp];
        } else {
          diag_eos(i,j,k,Temp_comp) = prob_param[T_r_comp];
        }
      } else if (amrex::Math::round(prob_param[idir_comp]) == 2) {
        if (y <= ctr_y) {
          diag_eos(i,j,k,Temp_comp) = prob_param[T_l_comp];
        } else {
          diag_eos(i,j,k,Temp_comp) = prob_param[T_r_comp];
        }
          } else if (amrex::Math::round(prob_param[idir_comp]) == 3) {
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
