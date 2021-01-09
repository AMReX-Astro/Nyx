#include "Nyx.H"
#include "Prob.H"

void prob_param_special_fill(amrex::GpuArray<amrex::Real,max_prob_param>& /*prob_param*/)
{}

#ifndef NO_HYDRO
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         amrex::Array4<amrex::Real> const& state,
                         amrex::GeometryData const& /*geomdata*/,
                         const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{

  amrex::Real rho0=1.0;
  amrex::Real temp0=10.0;
  amrex::Real ne0=1.0;

  const amrex::Real h_species_in=1;
  const amrex::Real gamma_minus_1_in=prob_param[gamma_comp] - 1.0;
  const amrex::Real a=1.0/(prob_param[z_in_comp]+1.0);
  amrex::Real eint, dummy_pres;

  nyx_eos_given_RT(NULL,gamma_minus_1_in, h_species_in, &eint, &dummy_pres, state(i,j,k,Density_comp), temp0, ne0, a);
  amrex::Real rhoe0 = rho0 * eint;

    state(i,j,k,Density_comp)    = rho0; //1.5d0 * small_dens
    state(i,j,k,Xmom_comp) = 0.00;
    state(i,j,k,Ymom_comp) = 0.00;
    state(i,j,k,Zmom_comp) = 0.00;

    // These will both be set later in the call to init_e.
    state(i,j,k,Eint_comp) = rhoe0;
    state(i,j,k,Eden_comp) = rhoe0 + 0.5 *
                                   (state(i, j, k, Xmom_comp) * state(i, j, k, Xmom_comp) +
                                    state(i, j, k, Ymom_comp) * state(i, j, k, Ymom_comp) +
                                    state(i, j, k, Zmom_comp) * state(i, j, k, Zmom_comp)) /
                                   state(i, j, k, Density_comp);

#ifndef CONST_SPECIES
      state(i,j,k,FirstSpec_comp  ) = prob_param[h_species_comp];
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
  amrex::Real temp0=10.0;
  amrex::Real ne0=1.0;

  diag_eos(i,j,k,Temp_comp) = temp0;
  diag_eos(i,j,k,  Ne_comp) = ne0;

  prob_initdata_state(i, j ,k, state, geomdata, prob_param);

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
