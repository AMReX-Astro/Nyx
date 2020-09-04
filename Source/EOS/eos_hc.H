
#include "Nyx.H"
#include "Nyx_F.H"
#include <atomic_rates_data.H>
#include <constants_cosmo.H>

using namespace amrex;

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void interp_to_this_z(const amrex::Real z, amrex::Real & ggh0, amrex::Real& gghe0, amrex::Real& gghep,
                      amrex::Real& eh0, amrex::Real& ehe0, amrex::Real& ehep)
{
    amrex::Real lopz, fact;
    int i, j;
    lopz   = std::log10(1.0e0 + z);

    if (lopz >= lzr(NCOOLFILE))
    {
        ggh0  = 0.0e0;
        gghe0 = 0.0e0;
        gghep = 0.0e0;
        eh0   = 0.0e0;
        ehe0  = 0.0e0;
        ehep  = 0.0e0;
        return;
    }

    if (lopz <= lzr(1))
        j = 1;
    else
        for(i = 2; i<= NCOOLFILE; i++)
        {
            if (lopz < lzr(i))
            {
                j = i-1;
                break;
            }
        }

    fact  = (lopz-lzr(j))/(lzr(j+1)-lzr(j));

    ggh0  = rggh0(j)  + (rggh0(j+1)-rggh0(j))*fact;
    gghe0 = rgghe0(j) + (rgghe0(j+1)-rgghe0(j))*fact;
    gghep = rgghep(j) + (rgghep(j+1)-rgghep(j))*fact;
    eh0   = reh0(j)   + (reh0(j+1)-reh0(j))*fact;
    ehe0  = rehe0(j)  + (rehe0(j+1)-rehe0(j))*fact;
    ehep  = rehep(j)  + (rehep(j+1)-rehep(j))*fact;
    return;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void ion_n_device(const int JH, const int JHe, const amrex::Real U,
                                        const amrex::Real nh, const amrex::Real* ne,
                                        amrex::Real& nhp, amrex::Real& nhep, amrex::Real& nhepp,
                                        amrex::Real* t, const amrex::Real gamma_minus_1,
                                        const amrex::Real z)
{

    amrex::Real ggh0, gghe0, gghep, eh0, ehe0, ehep;
    amrex::Real ahp, ahep, ahepp, ad, geh0, gehe0, gehep;
    amrex::Real ggh0ne, gghe0ne, gghepne;
    amrex::Real mu, tmp, logT, flo, fhi;
    const amrex::Real smallest_val=std::numeric_limits<amrex::Real>::min();
    int j;

    amrex::Real deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB;
    mu = (1.0e0+4.0e0*YHELIUM) / (1.0e0+YHELIUM+*ne);
    /*
    amrex::Real MPROTON = m_proton * M_unit;
    amrex::Real BOLTZMANN = k_B / (T_unit*T_unit/(M_unit*L_unit*L_unit));*/
    //    t  = gamma_minus_1*mp_over_kb * U * mu;
    *t  = gamma_minus_1*MPROTON/BOLTZMANN * U * mu;

    logT = std::log10(*t);

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

    ahp   = flo*AlphaHp  (j) + fhi*AlphaHp  (j+1);
    ahep  = flo*AlphaHep (j) + fhi*AlphaHep (j+1);
    ahepp = flo*AlphaHepp(j) + fhi*AlphaHepp(j+1);
    ad    = flo*Alphad   (j) + fhi*Alphad   (j+1);
    geh0  = flo*GammaeH0 (j) + fhi*GammaeH0 (j+1);
    gehe0 = flo*GammaeHe0(j) + fhi*GammaeHe0(j+1);
    gehep = flo*GammaeHep(j) + fhi*GammaeHep(j+1);

    interp_to_this_z(z, ggh0, gghe0, gghep, eh0, ehe0, ehep);

    if ((*ne) > 0.0)
    {
        ggh0ne   = JH  * ggh0  / (*ne*nh);
        gghe0ne  = JH  * gghe0 / (*ne*nh);
        gghepne  = JHe * gghep / (*ne*nh);
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
void iterate_ne_device(const int JH, const int JHe, const amrex::Real z,
                       const amrex::Real U, amrex::Real* t, const amrex::Real nh,
                       amrex::Real* ne, amrex::Real& nh0, amrex::Real& nhp,
                       amrex::Real& nhe0, amrex::Real& nhep, amrex::Real& nhepp,
                       const amrex::Real gamma_minus_1)
{
    int i;
    amrex::Real f, df, eps, ne2;
    amrex::Real nhp_plus, nhep_plus, nhepp_plus;
    amrex::Real dnhp_dne, dnhep_dne, dnhepp_dne, dne;

    i = 0;
    *ne = 1.0e0; // 0 is a bad guess

    for(i = i+1;i<=15;i++)  // Newton-Raphson solver
    {
        // Ion number densities
        ion_n_device(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t, gamma_minus_1, z);

        // Forward difference derivatives
        if (*ne > 0.0e0)
            eps = xacc*(*ne);
        else
            eps = 1.0e-24;

        ne2 = *ne+eps;
        ion_n_device(JH, JHe, U, nh, &ne2, nhp_plus, nhep_plus, nhepp_plus, t, gamma_minus_1, z);

        dnhp_dne   = (nhp_plus   - nhp)   / eps;
        dnhep_dne  = (nhep_plus  - nhep)  / eps;
        dnhepp_dne = (nhepp_plus - nhepp) / eps;

        f   = *ne - nhp - nhep - 2.0e0*nhepp;
        df  = 1.0e0 - dnhp_dne - dnhep_dne - 2.0e0*dnhepp_dne;
        dne = f/df;

        *ne = amrex::max((*ne-dne), 0.0e0);

        if (amrex::Math::abs(dne) < xacc)
            break;

    }

    // Get rates for the final ne
    ion_n_device(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t, gamma_minus_1, z);

    // Neutral fractions:
    nh0   = 1.0e0 - nhp;
    nhe0  = YHELIUM - (nhep + nhepp);
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_T_given_Re_device(Real gamma_minus_1, Real h_species, int JH, int JHe, Real* T, Real* Ne, Real R,Real e,Real comoving_a)
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

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void nyx_eos_given_RT(Real gamma_minus_1, Real h_species, Real* e, Real* P, Real R, Real T, Real Ne,Real comoving_a)
{
  
  Real mu = (1.0+4.0*YHELIUM) / (1.0+YHELIUM+Ne);
  *e = T / (gamma_minus_1 * mp_over_kb * mu);
  *P  = gamma_minus_1 * (R) * (*e);

}
