#ifndef _eos_hc_H_
#define _eos_hc_H_
#include <Nyx.H>
#include <atomic_rates_data.H>
#include <constants_cosmo.H>

using namespace amrex;

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void interp_to_this_z(AtomicRates* atomic_rates, const Real z, Real & ggh0, Real& gghe0, Real& gghep,
                      Real& eh0, Real& ehe0, Real& ehep)
{
    Real lopz, fact;
    int i, j;
    lopz   = std::log10(1.0e0 + z);

    if (lopz >= atomic_rates->lzr(NCOOLFILE))
    {
        ggh0  = 0.0e0;
        gghe0 = 0.0e0;
        gghep = 0.0e0;
        eh0   = 0.0e0;
        ehe0  = 0.0e0;
        ehep  = 0.0e0;
        return;
    }

    if (lopz <= atomic_rates->lzr(1))
        j = 1;
    else
        for(i = 2; i<= NCOOLFILE; i++)
        {
            if (lopz < atomic_rates->lzr(i))
            {
                j = i-1;
                break;
            }
        }

    fact  = (lopz-atomic_rates->lzr(j))/(atomic_rates->lzr(j+1)-atomic_rates->lzr(j));

    ggh0  = atomic_rates->rggh0(j)  + (atomic_rates->rggh0(j+1)-atomic_rates->rggh0(j))*fact;
    gghe0 = atomic_rates->rgghe0(j) + (atomic_rates->rgghe0(j+1)-atomic_rates->rgghe0(j))*fact;
    gghep = atomic_rates->rgghep(j) + (atomic_rates->rgghep(j+1)-atomic_rates->rgghep(j))*fact;
    eh0   = atomic_rates->reh0(j)   + (atomic_rates->reh0(j+1)-atomic_rates->reh0(j))*fact;
    ehe0  = atomic_rates->rehe0(j)  + (atomic_rates->rehe0(j+1)-atomic_rates->rehe0(j))*fact;
    ehep  = atomic_rates->rehep(j)  + (atomic_rates->rehep(j+1)-atomic_rates->rehep(j))*fact;
    return;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void ion_n_device(AtomicRates* atomic_rates, const int JH, const int JHe,
                                        const Real U, const Real nh, const Real& ne,
                                        Real& nhp, Real& nhep, Real& nhepp,
                                        Real& t, const Real gamma_minus_1, const Real h_species,
                                        const Real z)
{

    Real ggh0, gghe0, gghep, eh0, ehe0, ehep;
    Real ahp, ahep, ahepp, ad, geh0, gehe0, gehep;
    Real ggh0ne, gghe0ne, gghepne;
    Real mu, tmp, logT, flo, fhi;
    const Real smallest_val=std::numeric_limits<Real>::min();
    int j;

    Real deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB;
    Real YHELIUM = (Real(1.0)-h_species)/(Real(4.0)*h_species);
    mu = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+(ne));
    /*
    Real MPROTON = m_proton * M_unit;
    Real BOLTZMANN = k_B / (T_unit*T_unit/(M_unit*L_unit*L_unit));*/
    //    *t  = gamma_minus_1*mp_over_kb/e_to_cgs * U * mu;
    //    t = U / e_to_cgs * (gamma_minus_1 * mp_over_kb * mu);
    t  = gamma_minus_1*MPROTON/BOLTZMANN * U * mu;

    logT = std::log10(t);

    if (logT >= TCOOLMAX) // Fully ionized plasma
    {
        nhp   = 1.0e0;
        nhep  = 0.0e0;
        nhepp = YHELIUM;
        return;
    }

    // Temperature floor
    if (logT <= TCOOLMIN) logT = TCOOLMIN + 0.5e0*deltaT;

    // Interpolate rates
    tmp = (logT-TCOOLMIN)/deltaT;
    j = amrex::Math::floor(tmp); // replacing int(tmp)
    fhi = tmp - j;
    flo = 1.0e0 - fhi;
    j = j + 1; // F90 arrays start with 1, using Array1D<Real,1,size>

    ahp   = flo*atomic_rates->AlphaHp  (j) + fhi*atomic_rates->AlphaHp  (j+1);
    ahep  = flo*atomic_rates->AlphaHep (j) + fhi*atomic_rates->AlphaHep (j+1);
    ahepp = flo*atomic_rates->AlphaHepp(j) + fhi*atomic_rates->AlphaHepp(j+1);
    ad    = flo*atomic_rates->Alphad   (j) + fhi*atomic_rates->Alphad   (j+1);
    geh0  = flo*atomic_rates->GammaeH0 (j) + fhi*atomic_rates->GammaeH0 (j+1);
    gehe0 = flo*atomic_rates->GammaeHe0(j) + fhi*atomic_rates->GammaeHe0(j+1);
    gehep = flo*atomic_rates->GammaeHep(j) + fhi*atomic_rates->GammaeHep(j+1);

    interp_to_this_z(atomic_rates, z, ggh0, gghe0, gghep, eh0, ehe0, ehep);

    if ((ne) > 0.0)
    {
        ggh0ne   = JH  * ggh0  / (ne*nh);
        gghe0ne  = JH  * gghe0 / (ne*nh);
        gghepne  = JHe * gghep / (ne*nh);
    }
    else
    {
        ggh0ne   = 0.0e0;
        gghe0ne  = 0.0e0;
        gghepne  = 0.0e0;
    }

    // H+
    nhp = 1.0e0 - ahp/(ahp + geh0 + ggh0ne);

    // He+
    if ((gehe0 + gghe0ne) > smallest_val)
    {
        nhep  = YHELIUM/(1.0e0 + (ahep  + ad     )/(gehe0 + gghe0ne) 
                         + (gehep + gghepne)/ahepp);
    }
    else
        nhep  = 0.0e0;

    // He++
    if (nhep > 0.0e0)
        nhepp = nhep*(gehep + gghepne)/ahepp;
    else
        nhepp = 0.0e0;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void iterate_ne_device(AtomicRates* atomic_rates, const int JH, const int JHe, const Real z,
                       const Real U, Real& t, const Real nh,
                       Real& ne, Real& nh0, Real& nhp,
                       Real& nhe0, Real& nhep, Real& nhepp,
                       const Real gamma_minus_1, const Real h_species)
{
    int i;
    Real f, df, eps, ne2;
    Real nhp_plus, nhep_plus, nhepp_plus;
    Real dnhp_dne, dnhep_dne, dnhepp_dne, dne;
    Real YHELIUM = (Real(1.0)-h_species)/(Real(4.0)*h_species);

    i = 0;
    ne = 1.0e0; // 0 is a bad guess

    for(i = i+1;i<=15;i++)  // Newton-Raphson solver
    {
        // Ion number densities
        ion_n_device(atomic_rates, JH, JHe, U, nh, ne, nhp, nhep, nhepp, t, gamma_minus_1, h_species, z);

        // Forward difference derivatives
        if (ne > 0.0e0)
            eps = xacc*(ne);
        else
            eps = 1.0e-24;

        ne2 = ne+eps;
        ion_n_device(atomic_rates, JH, JHe, U, nh, ne2, nhp_plus, nhep_plus, nhepp_plus, t, gamma_minus_1, h_species, z);

        dnhp_dne   = (nhp_plus   - nhp)   / eps;
        dnhep_dne  = (nhep_plus  - nhep)  / eps;
        dnhepp_dne = (nhepp_plus - nhepp) / eps;

        f   = ne - nhp - nhep - 2.0e0*nhepp;
        df  = 1.0e0 - dnhp_dne - dnhep_dne - 2.0e0*dnhepp_dne;
        dne = f/df;

        ne = amrex::max((ne-dne), 0.0e0_rt);

        if (amrex::Math::abs(dne) < xacc)
            break;

    }

    // Get rates for the final ne
    ion_n_device(atomic_rates, JH, JHe, U, nh, ne, nhp, nhep, nhepp, t, gamma_minus_1, h_species, z);

    // Neutral fractions:
    nh0   = 1.0e0 - nhp;
    nhe0  = YHELIUM - (nhep + nhepp);
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_T_given_Re_device(
    AtomicRates* atomic_rates, Real gamma_minus_1,
    Real h_species, int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a,
    Real& nh0, Real& nhp, Real& nhe0, Real& nhep, Real& nhepp)
{
    Real z, rho, U, nh;
    // This converts from code units to CGS
    rho = R * density_to_cgs / (comoving_a*comoving_a*comoving_a);
    U = e * e_to_cgs;

    nh  = rho*h_species/MPROTON;

    z   = 1.e0/comoving_a - 1.e0;

    iterate_ne_device(atomic_rates, JH, JHe, z, U, *T, nh, *Ne, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1, h_species);
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_T_given_Re_device(
    AtomicRates* atomic_rates, Real gamma_minus_1,
    Real h_species, int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a)
{
    Real z, rho, U, nh, nh0, nhep, nhp, nhe0, nhepp;
    // This converts from code units to CGS
    rho = R * density_to_cgs / (comoving_a*comoving_a*comoving_a);
    U = e * e_to_cgs;

    nh  = rho*h_species/MPROTON;

    z   = 1.e0/comoving_a - 1.e0;

    iterate_ne_device(atomic_rates, JH, JHe, z, U, *T, nh, *Ne, nh0, nhp, nhe0, nhep, nhepp, gamma_minus_1, h_species);
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_given_RT(
    AtomicRates* atomic_rates, Real gamma_minus_1,
    Real h_species, Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a)
{
    Real YHELIUM = (1.0-h_species)/(4.0*h_species);
    Real mu = (1.0+4.0*YHELIUM) / (1.0+YHELIUM+Ne);
    *e = T / (gamma_minus_1 * mp_over_kb * mu);
    *P  = gamma_minus_1 * (R) * (*e);
}
#endif
