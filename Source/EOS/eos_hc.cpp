
#include "Nyx.H"
#include "Nyx_F.H"
#include <atomic_rates_data.H>
#include <constants_cosmo.H>
#include <eos_hc.H>

using namespace amrex;

AMREX_GPU_DEVICE void Nyx::nyx_eos_T_given_Re_device(Real gamma_minus_1, Real h_species, int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a)
{

    amrex::Real z, rho, U, nh, nh0, nhep, nhp, nhe0, nhepp;
    // This converts from code units to CGS
    rho = R * density_to_cgs / (comoving_a*comoving_a*comoving_a);
    U = e * e_to_cgs;
    nh  = rho*XHYDROGEN/MPROTON;

    z   = 1.e0/comoving_a - 1.e0;

    amrex::Real T_in = *T;
    amrex::Real Ne_in = *Ne;
    iterate_ne_device(JH, JHe, z, U, T, nh, Ne, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1);

}

AMREX_GPU_DEVICE void Nyx::nyx_eos_given_RT(Real gamma_minus_1, Real h_species, Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a)
{
  
  Real mu = (1.0+4.0*YHELIUM) / (1.0+YHELIUM+Ne);
  *e = T / (gamma_minus_1 * mp_over_kb * mu);
  *P  = gamma_minus_1 * (R) * (*e);

}
