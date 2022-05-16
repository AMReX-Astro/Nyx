#ifndef _atomic_rates_H_
#define _atomic_rates_H_

// This is a constant gamma equation of state, using an ideal gas.
//
// This a simplified version of the more general eos_gamma_general.

#include <atomic_rates_data.H>

void tabulate_rates(std::string file_in, amrex::Real mean_rhob)
{
    int Katz96 = 0;
    amrex::Real t3=1.0e3;
    amrex::Real t5=1.0e5;
    amrex::Real t6=1.0e6;
    amrex::Real deltaT = (TCOOLMAX - TCOOLMIN)/NCOOLTAB;
    amrex::Real t, U, E, y, sqrt_t, corr_term, tmp;
    AtomicRates atomic_rates_host;
    atomic_rates_host.mean_rhob=mean_rhob;
    int i = 1;

    if(!amrex::FileSystem::Exists(file_in))
        amrex::Abort("Treecool file does not exist at "+file_in);

    std::ifstream input_stream(file_in);
        //    while 
    for(i = 1; i <= NCOOLFILE; i++)
    {
        input_stream >> atomic_rates_host.lzr(i);
                input_stream >> atomic_rates_host.rggh0(i);
                input_stream >> atomic_rates_host.rgghe0(i);
                input_stream >> atomic_rates_host.rgghep(i);
                input_stream >> atomic_rates_host.reh0(i);
                input_stream >> atomic_rates_host.rehe0(i);
                input_stream >> atomic_rates_host.rehep(i);
    }
    Real check_length[7];
    auto count_length=0;
    while (!input_stream.eof())
    {
        for(i=0;i<7;i++)
	    input_stream>>check_length[i];
        count_length++;
    }

    if(count_length>1)
    {
        std::string abort_message="Using TREECOOL file longer than NCOOLFILE, change NCOOLFILE in atomic_rates_data.H and recompile. ";
	abort_message+="Use NCOOLFILE="+std::to_string(NCOOLFILE+count_length-1);
        amrex::Abort(abort_message);
    }

    // Initialize cooling tables
    t = std::pow(10,TCOOLMIN);
    if (Katz96)
    {
        for(i = 1; i <= NCOOLTAB+1; i++)
        {
            // Rates are as in Katz et al. 1996
            sqrt_t = std::sqrt(t);
            corr_term    = 1.e0 / (1.0e0 + sqrt_t/std::sqrt(t5));

            // Recombination rates
            // Alphad: dielectronic recombination rate of singly ioniozed helium
            atomic_rates_host.Alphad(i)    = 1.90e-03/(t*sqrt_t) * std::exp(-4.7e5/t) * (1.0e0+0.3e0*std::exp(-9.4e4/t));
            atomic_rates_host.AlphaHp(i)   = 8.40e-11/sqrt_t * std::pow((t/t3),(-0.2e0)) / (1.0e0 + std::pow((t/t6),0.7e0));
            atomic_rates_host.AlphaHep(i)  = 1.50e-10 * std::pow(t,(-0.6353e0));
            atomic_rates_host.AlphaHepp(i) = 3.36e-10/sqrt_t * std::pow((t/t3),(-0.2e0)) / (1.0e0 + std::pow((t/t6),0.7e0));

            // Collisional ionization rates
            atomic_rates_host.GammaeH0(i)  = 5.85e-11*sqrt_t * std::exp(-157809.1e0/t) * corr_term;
            atomic_rates_host.GammaeHe0(i) = 2.38e-11*sqrt_t * std::exp(-285335.4e0/t) * corr_term;
            atomic_rates_host.GammaeHep(i) = 5.68e-12*sqrt_t * std::exp(-631515.0e0/t) * corr_term;

            // Collisional ionization & excitation cooling rates
            atomic_rates_host.BetaH0(i)  = 7.5e-19 * std::exp(-118348.0e0/t) * corr_term + 2.171e-11*atomic_rates_host.GammaeH0(i);
            atomic_rates_host.BetaHe0(i) = 3.941e-11 * atomic_rates_host.GammaeHe0(i);
            atomic_rates_host.BetaHep(i) = 5.54e-17 * std::pow(t,(-0.397e0)) * std::exp(-473638.0e0/t) * corr_term + 0.715e-11 * atomic_rates_host.GammaeHep(i);

            // Recombination cooling rates
            atomic_rates_host.RecHp(i)   = 1.036e-16 * t * atomic_rates_host.AlphaHp(i);
            atomic_rates_host.RecHep(i)  = 1.036e-16 * t * atomic_rates_host.AlphaHep(i) + 6.526e-11 * atomic_rates_host.Alphad(i);
            atomic_rates_host.RecHepp(i) = 1.036e-16 * t * atomic_rates_host.AlphaHepp(i);

            // Free-free cooling rate
            atomic_rates_host.Betaff1(i) = 1.42e-27 * sqrt_t * (1.1e0 + 0.34e0*std::exp(-(5.5e0 - std::log10(t)) * (5.5e0 - std::log10(t)) / 3.0e0));
            atomic_rates_host.Betaff4(i) = atomic_rates_host.Betaff1(i);

            t = t*std::pow(10,deltaT);
        }
    }
    else
    {
        for(i = 1; i <= NCOOLTAB+1; i++)
        {
            // Rates are as in Lukic et al.
            sqrt_t = std::sqrt(t);

            // Recombination rates
            // Alphad: dielectronic recombination rate of singly ioniozed helium
            atomic_rates_host.Alphad(i)    = 1.90e-03/(t*sqrt_t) * std::exp(-4.7e5/t) * (1.0e0+0.3e0*std::exp(-9.4e4/t));
            atomic_rates_host.AlphaHp(i)   = 7.982e-11 / (std::sqrt(t/3.148e0)* std::pow((1.0e0+std::sqrt(t/3.148e0)),0.252) * 
                                        std::pow((1.0e0+std::sqrt(t/7.036e5)),1.748));
               if (t <= 1.0e6)
                   atomic_rates_host.AlphaHep(i)  = 3.294e-11 / (std::sqrt(t/15.54e0)* std::pow((1.0e0+std::sqrt(t/15.54e0)),0.309) * 
                                               std::pow((1.0e0+std::sqrt(t/3.676e7)),1.691));
               else
                   atomic_rates_host.AlphaHep(i)  = 9.356e-10 / (std::sqrt(t/4.266e-2)* std::pow((1.0e0+std::sqrt(t/4.266e-2)),0.2108) * 
                                               std::pow((1.0e0+std::sqrt(t/4.677e6)),1.7892));

               atomic_rates_host.AlphaHepp(i) = 1.891e-10 / (std::sqrt(t/9.37e0)* std::pow((1.0e0+std::sqrt(t/9.37e0)),0.2476) * 
                                           std::pow((1.0e0+std::sqrt(t/2.774e6)),1.7524));

               // Collisional ionization rates
               E = 13.6e0;
               U = 1.16045e4*E/t;
               atomic_rates_host.GammaeH0(i)  = 2.91e-8*std::pow(U,0.39)*std::exp(-U) / (0.232e0+U);
               E = 24.6e0;
               U = 1.16045e4*E/t;
               atomic_rates_host.GammaeHe0(i) = 1.75e-8*std::pow(U,0.35)*std::exp(-U) / (0.18e0+U);
               E = 54.4e0;
               U = 1.16045e4*E/t;
               atomic_rates_host.GammaeHep(i) = 2.05e-9*(1.0e0+std::sqrt(U))*std::pow(U,0.25)*std::exp(-U) / (0.265e0+U);

               // Collisional ionization & excitation cooling rates
               corr_term = 1.e0 / (1.0e0 + sqrt_t/std::sqrt(5.0e7));
               y = std::log(t);
               if (t <= 1.0e5)
                  atomic_rates_host.BetaH0(i)  = 1.0e-20 * std::exp( 2.137913e2 - 1.139492e2*y + 2.506062e1*y*y - 
                                               2.762755e0*y*y*y + 1.515352e-1*y*y*y*y - 
                                               3.290382e-3*y*y*y*y*y - 1.18415e5/t );
               else
                  atomic_rates_host.BetaH0(i)  = 1.0e-20 * std::exp( 2.7125446e2 - 9.8019455e1*y + 1.400728e1*y*y - 
                                               9.780842e-1*y*y*y + 3.356289e-2*y*y*y*y - 
                                               4.553323e-4*y*y*y*y*y - 1.18415e5/t );

               atomic_rates_host.BetaHe0(i) = 9.38e-22 * sqrt_t * std::exp(-285335.4e0/t) * corr_term;
               atomic_rates_host.BetaHep(i) = (5.54e-17 * std::pow(t,(-0.397e0)) * std::exp(-473638.0e0/t) +
                             4.95e-22 * sqrt_t * std::exp(-631515.0e0/t) )*corr_term;

               // Recombination cooling rates
               atomic_rates_host.RecHp(i)   = 2.851e-27 * sqrt_t * (5.914e0-0.5e0*std::log(t)+1.184e-2*std::pow(t,(1.0e0/3.0e0)));
               atomic_rates_host.RecHep(i)  = 1.55e-26 * std::pow(t,0.3647) + 1.24e-13/(t*sqrt_t) * std::exp(-4.7e5/t) * (1.0e0+0.3e0*std::exp(-9.4e4/t));
               atomic_rates_host.RecHepp(i) = 1.14e-26 * sqrt_t * (6.607e0-0.5e0*std::log(t)+7.459e-3*std::pow(t,(1.0e0/3.0e0)));

               // Free-free cooling rate
               if (t <= 3.2e5)
                   atomic_rates_host.Betaff1(i) = 1.426e-27 * sqrt_t * (0.79464e0 + 0.1243e0*std::log10(t));
               else
                   atomic_rates_host.Betaff1(i) = 1.426e-27 * sqrt_t * (2.13164e0 - 0.1240e0*std::log10(t));

               if (t/4.0e0 <= 3.2e5)
                   atomic_rates_host.Betaff4(i) = 1.426e-27 * sqrt_t * 4.0e0*(0.79464e0 + 0.1243e0*std::log10(t/4.0e0));
               else
                   atomic_rates_host.Betaff4(i) = 1.426e-27 * sqrt_t * 4.0e0*(2.13164e0 - 0.1240e0*std::log10(t/4.0e0));

               t = t*std::pow(10,deltaT);
        }
    }
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(atomic_rates_glob,&atomic_rates_host,sizeof(AtomicRates));
#else
    std::memcpy(atomic_rates_glob,&atomic_rates_host,sizeof(AtomicRates));
#endif
    return;
}

/*
//Used to store ggh0, gghe0, gghep, eh0, ehe0, ehep for iterate_ne
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
interp_to_this_z()
{

}*/

#endif
