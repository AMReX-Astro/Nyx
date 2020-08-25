#ifndef __f_rhs__
#define __f_rhs__

#include <atomic_rates_data.H>
#include <constants_cosmo.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
f_rhs_rpar(amrex::Real const time, amrex::Real e_in, amrex::Real energy, amrex::Real* rpar)
{

    const amrex::Real compt_c = 1.01765467e-37;
    const amrex::Real T_cmb = 2.725e0;

    amrex::Real logT, tmp, fhi, flo;
    amrex::Real ggh0, gghe0, gghep, eh0, ehe0, ehep;
    amrex::Real ahp, ahep, ahepp, ad, geh0, gehe0, gehep;
    amrex::Real bh0, bhe0, bhep, bff1, bff4, rhp, rhep, rhepp;
    amrex::Real lambda_c, lambda_ff, lambda, heat;
    amrex::Real rho, U, a, rho_heat;
    amrex::Real nh, nh0, nhp, nhe0, nhep, nhepp;
    int j;
    amrex::Real z_vode, rho_vode, T_vode, ne_vode;

    T_vode=rpar[0];
    ne_vode=rpar[1];
    rho_vode=rpar[2];
    z_vode=rpar[3];

    if(e_in<=0) e_in=std::numeric_limits<amrex::Real>::lowest();

    // Converts from code units to CGS
    rho = rho_vode * density_to_cgs * (1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode))*(1.0e0+amrex::Math::abs(z_vode));
    U = e_in * e_to_cgs;
    nh  = rho*XHYDROGEN/MPROTON;

    // Get gas temperature and individual ionization species
    // testing different memory structures
    iterate_ne_device(JH_vode, JHe_vode, z_vode, U, &T_vode, nh, &ne_vode, nh0, nhp, nhe0, nhep, nhepp);

    // Convert species to CGS units: 
    ne_vode = nh * ne_vode;
    nh0   = nh * nh0;
    nhp   = nh * nhp;
    nhe0  = nh * nhe0;
    nhep  = nh * nhep;
    nhepp = nh * nhepp;

    logT = std::log10(T_vode);
    if (logT >= TCOOLMAX) // Only free-free and Compton cooling are relevant
    {
        lambda_ff = 1.42e-27 * std::sqrt(T_vode) * (1.1e0 + 0.34e0*std::exp(-(5.5e0 - logT)*(5.5e0 - logT) / 3.0e0))
            * (nhp + 4.0e0*nhepp)*ne_vode;
        lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));

        energy  = (-lambda_ff -lambda_c) * heat_from_cgs/((1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode)));

        // Convert to the actual term to be used in e_out = e_in + dt*energy
        energy  = energy / rho_vode * (1.0e0+amrex::Math::abs(z_vode));
        ne_vode = ne_vode / nh;
        rpar[0]=T_vode;
        rpar[1]=ne_vode;
        rpar[2]=rho_vode;
        return;
    }

    // Temperature floor
    if (logT <= TCOOLMIN)  logT = TCOOLMIN + 0.5e0*deltaT;

    // Interpolate rates
    // Interpolate rates
    tmp = (logT-TCOOLMIN)/deltaT;
    j = std::floor(tmp); // replacing int(tmp)
    fhi = tmp - j;
    flo = 1.0e0 - fhi;
    j = j + 1; // F90 arrays start with 1, using Array1D<Real,1,size>

    bh0   = flo*BetaH0   (j) + fhi*BetaH0   (j+1);
    bhe0  = flo*BetaHe0  (j) + fhi*BetaHe0  (j+1);
    bhep  = flo*BetaHep  (j) + fhi*BetaHep  (j+1);
    bff1  = flo*Betaff1  (j) + fhi*Betaff1  (j+1);
    bff4  = flo*Betaff4  (j) + fhi*Betaff4  (j+1);
    rhp   = flo*RecHp    (j) + fhi*RecHp    (j+1);
    rhep  = flo*RecHep   (j) + fhi*RecHep   (j+1);
    rhepp = flo*RecHepp  (j) + fhi*RecHepp  (j+1);

    // Cooling: 
    lambda = ( bh0*nh0 + bhe0*nhe0 + bhep*nhep + 
               rhp*nhp + rhep*nhep + rhepp*nhepp + 
               bff1*(nhp+nhep) + bff4*nhepp ) * ne_vode;

    lambda_c  = compt_c*T_cmb*T_cmb*T_cmb*T_cmb * ne_vode * (T_vode - T_cmb*(1.0e0+amrex::Math::abs(z_vode)))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode))*(1.0e0 + amrex::Math::abs(z_vode));
    lambda = lambda + lambda_c;

    interp_to_this_z(z_vode, ggh0, gghe0, gghep, eh0, ehe0, ehep);
    // Heating terms
    heat = JH_vode*nh0*eh0 + JH_vode*nhe0*ehe0 + JHe_vode*nhep*ehep;
    
    rho_heat = uvb_density_A * (rho_vode/mean_rhob)**uvb_density_B;
    heat = rho_heat*heat;

    // Convert back to code units
    ne_vode     = ne_vode / nh;
    energy  = (heat - lambda)*heat_from_cgs* (1.0e0+amrex::Math::abs(z_vode));

    // Convert to the actual term to be used in e_out = e_in + dt*energy
    a = 1.e0 / (1.e0 + amrex::Math::abs(z_vode));
    energy = energy / rho_vode / a;
    rpar[0]=T_vode;
    rpar[1]=ne_vode;
    rpar[2]=rho_vode;
}
#endif
