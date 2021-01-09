
#include <Nyx.H>
#include <atomic_rates_data.H>
#include <constants_cosmo.H>

using namespace amrex;

extern AtomicRates* atomic_rates_glob;
/*
#ifdef __cplusplus
extern "C"
{
#endif

  AMREX_GPU_DEVICE void fort_nyx_eos_T_given_Re_device(int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a);
  AMREX_GPU_DEVICE void fort_nyx_eos_given_RT(Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a);

#ifdef __cplusplus
}
#endif
*/

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_T_given_Re_device(
    AtomicRates* /*atomic_rates*/, Real gamma_minus_1, Real h_species, 
    int /*JH*/, int /*JHe*/, Real* T, Real* Ne, Real /*R*/, Real e, Real /*comoving_a*/)
{
  Real xn_eos[2];

  xn_eos[0] = h_species;
  xn_eos[1] = (1.e0 - h_species);

  const int nspecd = 2;

  Real aiond[nspecd];
  
  aiond[0] = 1.0;
  aiond[1] = 4.0;

  Real ziond[nspecd];
  
  ziond[0] = 1.0;
  ziond[1] = 2.0;
    
  //    !-------------------------------------------------------------------------
  //    ! compute mu -- the mean molecular weight
  //    !-------------------------------------------------------------------------

  Real sum_y  = 0.0;
  Real ymass[nspecd];

  // assume completely ionized species
  for(int n=0;n< nspecd; n++)
    {
      ymass[n] = xn_eos[n]*(1.e0 + ziond[n])/aiond[n];
      sum_y = sum_y + ymass[n];
    }

  Real mu = 1.0/sum_y;

  //    !-------------------------------------------------------------------------
  //    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
  //    ! and temp as needed from the inputs
  //    !-------------------------------------------------------------------------

  //    ! These are both very large numbers so we pre-divide
  Real  m_nucleon_over_kB = m_nucleon / k_B;

  // Solve e = k T / [(mu m_nucleon)*(gamma-1)] for T
  *T = e * (mu * m_nucleon_over_kB * gamma_minus_1);

  //    !-------------------------------------------------------------------------
  //    ! now we have the density and temperature (and mass fractions / mu)
  //    ! regardless of the inputs.
  //    !-------------------------------------------------------------------------

   *Ne = 1.0;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_given_RT(
     AtomicRates* /*atomic_rates*/, Real gamma_minus_1, Real h_species, 
     Real* e, Real* P, Real R, Real T, Real /*Ne*/, Real comoving_a)
{
  
  Real comoving_a_cubed= (comoving_a*comoving_a*comoving_a);

  // Density is the only variable we convert from comoving to proper coordinates
  Real den_eos = R / comoving_a_cubed;

  Real temp_eos = T;

  Real xn_eos[2];

  xn_eos[0] = h_species;
  xn_eos[1] = (1.0 - h_species);

  const int nspecd = 2;

  Real aiond[nspecd];
  
  aiond[0] = 1.0;
  aiond[1] = 4.0;

  Real ziond[nspecd];
  
  ziond[0] = 1.0;
  ziond[1] = 2.0;
    
  //    !-------------------------------------------------------------------------
  //    ! compute mu -- the mean molecular weight
  //    !-------------------------------------------------------------------------

  Real sum_y  = 0.0;
  Real ymass[nspecd];

  // assume completely ionized species
  for(int n=0;n< nspecd; n++)
    {
      ymass[n] = xn_eos[n]*(1.e0 + ziond[n])/aiond[n];
      sum_y = sum_y + ymass[n];
    }

  Real mu = 1.0/sum_y;

  //    !-------------------------------------------------------------------------
  //    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
  //    ! and temp as needed from the inputs
  //    !-------------------------------------------------------------------------

  //    ! These are both very large numbers so we pre-divide
  Real  m_nucleon_over_kB = m_nucleon / k_B;

  // Compute e = k T / [(mu m_nucleon)*(gamma-1)] from T
  *e = temp_eos / (mu * m_nucleon_over_kB * gamma_minus_1);

  // -------------------------------------------------------------------------
  //  now we have the density and temperature (and mass fractions / mu)
  //  regardless of the inputs.
  // -------------------------------------------------------------------------

  // Compute the pressure from the ideal gas law -- 
  //    note it must be converted from proper to comoving coordinates
  *P = ( gamma_minus_1 * den_eos * (*e) ) * comoving_a_cubed;;

}
