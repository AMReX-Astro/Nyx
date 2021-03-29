#include <Nyx.H>

#include <Forcing.H>

using namespace amrex;

unsigned long int mt_random();

int StochasticForcing::verbose      = 0;
int StochasticForcing::SpectralRank = 3;

//
//  Default constructor
//
StochasticForcing::StochasticForcing() 
{
    i1 = i2 = 0;
    j1 = j2 = 0;
    k1 = k2 = 0;
    NumModes = 0;
    NumNonZeroModes = 0;
    decay = 0; 
    seed = 27011974;

    SpectProfile = Parabolic;

    AmpltThresh = 1.051250e-1;
    SolenoidalWeight = 1.0;
    DecayInitTime = 0.0;

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        alpha[dim]         = 2;         
        BandWidth[dim]     = 1.0;
        IntgrVelocity[dim] = 0.0; 
        IntgrLength[dim]   = 0.0;
        WaveNumber[dim]    = 0.0;
        IntgrTime[dim]     = 0.0;
        AutoCorrlTime[dim] = 1.0;

        Amplitude[dim]     = NULL;
        InjectionEven[dim] = NULL;
        InjectionOdd[dim]  = NULL;
        SpectrumEven[dim]  = NULL;
        SpectrumOdd[dim]   = NULL;
        wavevectors[dim]   = NULL;
        modes_even[dim]    = NULL;
        modes_odd[dim]     = NULL;
    }
    mask = NULL;
}

//
//  Default destructor
//
StochasticForcing::~StochasticForcing() 
{
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        delete [] Amplitude[dim];
        delete [] InjectionEven[dim];
        delete [] InjectionOdd[dim];
        delete [] SpectrumEven[dim];
        delete [] SpectrumOdd[dim];
        delete [] wavevectors[dim];
        delete [] modes_even[dim];
        delete [] modes_odd[dim];
    }
    delete [] mask;
}

/***********************************************************************
 *
 *  STOCHASTIC FORCING CLASS: evolve
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:  Oct, 2014: updated to support Enzo 2.4 // P. Grete
 *  modified2:  May, 2017: ported to Nyx
 *
 *  PURPOSE: evolves the random forcing spectrum in the fashion of
 *           a multi-dimensional Ornstein-Uhlenbeck process
 *           
 *           Parameters:
 *           dt -- time step (small compared to AutoCorrlTime)
 *
 *  AUXILIARIES: inject, gauss_deviate, distribute, rms
 *
 ***********************************************************************/

void StochasticForcing::evolve(Real dt)
{
    if (ParallelDescriptor::IOProcessor()) {

        Real DriftCoeff[MAX_DIMENSION], DiffCoeff[MAX_DIMENSION];

        if (decay == 0) {

            inject();
                    
            /* Increment forcing spectrum (drift and random diffusion) 
             * For general properties of Ornstein-Uhlenbeck process, see e.g.
             * Turbulent Flows by Pope (2000) Appendix J with 
             * drift and diffusion coefficients given eq (J.41)
             */

            for (int dim = 0; dim < SpectralRank; dim++) {
                DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
                DiffCoeff[dim]  = sqrt(1 - DriftCoeff[dim]*DriftCoeff[dim]);
                for (int n = 0, m = 0; n < NumModes; n++)
                    if (mask[n]) {
                        SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m] +
                                               DiffCoeff [dim] * InjectionEven[dim][n];
                        SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m] + 
                                               DiffCoeff [dim] * InjectionOdd [dim][n];
                        ++m;
                    }
            }

        } else {

            /* increment forcing spectrum (drift only) */

            for (int dim = 0; dim < SpectralRank; dim++) {
                DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
                for (int m = 0; m < NumNonZeroModes; m++) {
                    SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m];
                    SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m];
                }
            }
        }
    }

    /* communicate spectrum among processors */

    distribute();
}

//
// Compute new random injection
//
void StochasticForcing::inject(void)
{
    if (ParallelDescriptor::IOProcessor()) {

        int i, j, k, n, dim;
        Real a, b, contr;

        /* compute Gaussian deviates */

        for (dim = 0; dim < SpectralRank; dim++)
            for (n = 0; n < NumModes; n++) {
                if (mask[n]) {
                    gauss_deviate(Amplitude[dim][n], &a, &b);
                } else {
                    a = 0.0; b = 0.0;
                }
                InjectionEven[dim][n] = a;
                InjectionOdd[dim][n]  = b;
            }

        /* project modes 
         * see eq (8) in Schmidt et al., A&A (2009)
         * http://dx.doi.org/10.1051/0004-6361:200809967 */

        for (i = 0; i < i2; i++) { // wave vectors in positive x-direction
            InjectionEven[0][i] = (1.0 - SolenoidalWeight) * InjectionEven[0][i];
            InjectionOdd[0][i]  = (1.0 - SolenoidalWeight) * InjectionOdd[0][i];
        }

        if (SpectralRank > 1) {
            
            for (n = 0; n < i2; n++) { // wave vectors in positive x-direction
                InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n];
                InjectionOdd [1][n] = SolenoidalWeight * InjectionOdd [1][n];
            }
            
            n = i2;
            for (j = 1; j <= j2; j++) { // wave vectors in xy-plane
                for (i = i1; i <= i2; i++) {
                    contr = (1.0 - 2.0 * SolenoidalWeight) * 
                        (i*InjectionEven[0][n] + 
                         j*InjectionEven[1][n]) / Real(i*i + j*j);
                    InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
                    InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
                    contr = (1.0 - 2.0 * SolenoidalWeight) * 
                        (i*InjectionOdd[0][n] + 
                         j*InjectionOdd[1][n]) / Real(i*i + j*j);
                    InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
                    InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
                    ++n;
                }
            }
            
            if (SpectralRank > 2) {
                
                for (n = 0; n < i2 + j2*(i2-i1+1); n++) { // wave vectors in xy-plane
                    InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n];
                    InjectionOdd[2][n]  = SolenoidalWeight * InjectionOdd [2][n];
                }
                
                for (k = 1; k <= k2; k++) { // wave vectors not aligned to xy-plane
                    for (j = j1; j <= j2; j++) {
                        for (i = i1; i <= i2; i++) {
                            contr = (1.0 - 2.0 * SolenoidalWeight) * 
                                (i*InjectionEven[0][n] + 
                                 j*InjectionEven[1][n] + 
                                 k*InjectionEven[2][n] ) / Real(i*i + j*j + k*k);
                            InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
                            InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
                            InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n] + k*contr;
                            contr = (1.0 - 2.0 * SolenoidalWeight) * 
                                (i*InjectionOdd[0][n] + 
                                 j*InjectionOdd[1][n] +
                                 k*InjectionOdd[2][n]) / Real(i*i + j*j + k*k);
                            InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
                            InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
                            InjectionOdd[2][n] = SolenoidalWeight * InjectionOdd[2][n] + k*contr;
                            ++n;

                        }
                    }
                }
            }
        }

    }
}

//
// Generate couple of normally distributed random deviates (Box-Muller-Algorithm)
//
void StochasticForcing::gauss_deviate(Real amplt, Real *x, Real *y)
{
        Real v_sqr, v1, v2;
        Real norm;

        do {
            v1 = 2.0* (Real)(mt_random()%2147483563)/(2147483563.0) - 1.0;
            v2 = 2.0* (Real)(mt_random()%2147483563)/(2147483563.0) - 1.0;
            v_sqr = v1*v1+v2*v2;
        } while (v_sqr >= 1.0 || v_sqr == 0.0);
        
        norm = amplt * sqrt(-2.0*log(v_sqr)/v_sqr);

        *x = norm * v1; *y = norm * v2;
}

//
// Distribute the spectrum
//
void StochasticForcing::distribute(void)
{
    /* communicate spectrum among processors */

    for (int dim = 0; dim < SpectralRank; dim++) {
        ParallelDescriptor::Bcast(SpectrumEven[dim], NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(SpectrumOdd[dim],  NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
    }

    /* copy sepctrum to forcing_spect_module */

    for (int dim = 0; dim < SpectralRank; dim++)
        for (int l = 0; l < NumNonZeroModes; l++) {
            modes_even[dim][l]=SpectrumEven[dim][l];
            modes_odd[dim][l]=SpectrumOdd[dim][l];
        }
}

//
// Compute RMS magnitude
// 
Real StochasticForcing::rms(void)
{
    int m;
    Real sum_even = 0.0, sum_odd = 0.0;

    for (int dim = 0; dim < SpectralRank; dim++) {
        for (m = 0; m < NumNonZeroModes; m++)
            sum_even += SpectrumEven[dim][m] * SpectrumEven[dim][m];
        for (m = 0; m < NumNonZeroModes; m++)
            sum_odd  += SpectrumOdd[dim][m]  * SpectrumOdd[dim][m];
    }

    return sqrt(sum_even + sum_odd);
}

void StochasticForcing::integrate_state_force(
  amrex::Box const& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& diag_eos,
  amrex::GeometryData const& geomdata, amrex::Real a,
  amrex::Real half_dt,
  amrex::Real small_eint, amrex::Real small_temp)
{
    int mi, mj, mk;
    int num_phases[3];
    Real delta_phase[3];
    Real phase_lo[3];
    Real accel[3];

    int num_modes = NumNonZeroModes;

    Vector<Real> buf(num_modes);
    Vector<Real> phasefct_init_even(num_modes);
    Vector<Real> phasefct_init_odd(num_modes);

    Vector<Real> phasefct_mult_even_x(num_modes);
    Vector<Real> phasefct_mult_even_y(num_modes);
    Vector<Real> phasefct_mult_even_z(num_modes);

    Vector<Real> phasefct_mult_odd_x(num_modes);
    Vector<Real> phasefct_mult_odd_y(num_modes);
    Vector<Real> phasefct_mult_odd_z(num_modes);
    Vector<Real> phasefct_yz0(num_modes);
    Vector<Real> phasefct_yz1(num_modes);

    Real *phasefct_even_x, *phasefct_even_y, *phasefct_even_z;
    Real *phasefct_odd_x, *phasefct_odd_y, *phasefct_odd_z;

    amrex::Real alpha_const =100.0;
    amrex::Real temp0_const =10.0;

    // Note that (lo,hi) define the region of the box containing the grow cells
    // Do *not* assume this is just the valid region
    // apply heating-cooling to Eden_comp and Eint_comp
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            // Original values
            Real rho        = state(i,j,k,Density_comp);
            Real rho_e_orig = state(i,j,k,Eint_comp);
            Real rho_K_res  = state(i,j,k,Eden_comp) - state(i,j,k,Eint_comp);
            Real T_orig     = diag_eos(i,j,k,Temp_comp);
                //ne         = diag_eos(i,j,k,  NE_COMP)

                //                if (rho_e_orig < 0.0) neg_e_count = neg_e_count + 1

             // Compute temperature increment and ensure that new temperature is positive
             Real delta = half_dt * alpha_const * (temp0_const - T_orig) / a;
            diag_eos(i,j,k,Temp_comp) = amrex::max(T_orig + delta, small_temp);

            // Call EOS to get internal energy for constant equilibrium temperature
            Real eint0 = temp0_const / (small_temp / small_eint);
            Real delta_re = half_dt * alpha_const * (rho*eint0 - rho_e_orig) / a;

            // Call EOS to get the internal energy floor

            // Update cell quantities
            state(i,j,k,Eint_comp) = amrex::max(rho_e_orig + delta_re, rho*small_eint);
            state(i,j,k,Eden_comp) = state(i,j,k,Eint_comp) + rho_K_res;
        });

    const auto prob_hi = geomdata.ProbHi();
    const auto prob_lo = geomdata.ProbLo();
    const auto dx = geomdata.CellSize();
    for (int dim = 0; dim < SpectralRank; dim++) {
        delta_phase[dim] = 2.0*M_PI * dx[dim] / (prob_hi[dim] - prob_lo[dim]); // phase increment per cell
        phase_lo[dim] = (double(bx.smallEnd(dim)) + 0.5) * delta_phase[dim];           // phase of low corner
        num_phases[dim] = (bx.bigEnd(dim)-bx.smallEnd(dim)+1)*num_modes;
    }

    phasefct_even_x = new Real[num_phases[0]];
    phasefct_even_y = new Real[num_phases[1]];
    phasefct_even_z = new Real[num_phases[2]];
    phasefct_odd_x  = new Real[num_phases[0]];
    phasefct_odd_y  = new Real[num_phases[1]];
    phasefct_odd_z  = new Real[num_phases[2]];

    for (int m = 0; m < num_modes; m++) {
        int i = wavevectors[0][m];
        int j = wavevectors[1][m];
        int k = wavevectors[2][m];
        phasefct_init_even[m] =
            (cos(i*phase_lo[0]) * cos(j*phase_lo[1]) -
             sin(i*phase_lo[0]) * sin(j*phase_lo[1])) * cos(k*phase_lo[2]) -
            (cos(i*phase_lo[0]) * sin(j*phase_lo[1]) +
             sin(i*phase_lo[0]) * cos(j*phase_lo[1])) * sin(k*phase_lo[2]);

        phasefct_init_odd[m] =
            (cos(i*phase_lo[0]) * cos(j*phase_lo[1]) -
             sin(i*phase_lo[0]) * sin(j*phase_lo[1])) * sin(k*phase_lo[2]) +
            (cos(i*phase_lo[0]) * sin(j*phase_lo[1]) +
             sin(i*phase_lo[0]) * cos(j*phase_lo[1])) * cos(k*phase_lo[2]);

        phasefct_mult_even_x[m] = cos(i*delta_phase[0]);
        phasefct_mult_odd_x[m] = sin(i*delta_phase[0]);

        phasefct_mult_even_y[m] = cos(j*delta_phase[1]);
        phasefct_mult_odd_y[m] = sin(j*delta_phase[1]);

        phasefct_mult_even_z[m] = cos(k*delta_phase[2]);
        phasefct_mult_odd_z[m] = sin(k*delta_phase[2]);
    }

    // initialize phase factors for each coordinate axis:
    // since phase factors for inverse FT are given by
    // exp(i*(k1*x + k2*y + k3*z)) = exp(i*k1*x) * exp(i*k2*y)*...,
    // we iteratively multiply with exp(i*k1*delta_x), etc.
    for (int m = 0; m < num_modes; m++) {
       phasefct_even_x[m] = 1.0;
       phasefct_odd_x[m]  = 0.0;
    }
    for (int i = bx.smallEnd(0)+1; i <= bx.bigEnd(0); i++) {
       mi = (i-bx.smallEnd(0))*num_modes;
       for (int m = 0; m < num_modes; m++) {
           buf[m] = phasefct_even_x[mi-num_modes];
           phasefct_even_x[mi] = phasefct_mult_even_x[m] * phasefct_even_x[mi-num_modes] -
                                 phasefct_mult_odd_x[m]  * phasefct_odd_x[mi-num_modes];
           phasefct_odd_x[mi]  = phasefct_mult_even_x[m] * phasefct_odd_x[mi-num_modes] +
                                 phasefct_mult_odd_x[m]  * buf[m];
           mi = mi + 1;
       }
    }

    for (int m = 0; m < num_modes; m++) {
       phasefct_even_y[m] = 1.0;
       phasefct_odd_y[m]  = 0.0;
    }
    for (int j = bx.smallEnd(1)+1; j <= bx.bigEnd(1); j++) {
       mj = (j-bx.smallEnd(1))*num_modes;
       for (int m = 0; m < num_modes; m++) {
           buf[m] = phasefct_even_y[mj-num_modes];
           phasefct_even_y[mj] = phasefct_mult_even_y[m] * phasefct_even_y[mj-num_modes] -
                                 phasefct_mult_odd_y[m]  * phasefct_odd_y[mj-num_modes];
           phasefct_odd_y[mj]  = phasefct_mult_even_y[m] * phasefct_odd_y[mj-num_modes] +
                                 phasefct_mult_odd_y[m]  * buf[m];
           mj = mj + 1;
       }
    }

    for (int m = 0; m < num_modes; m++) {
       phasefct_even_z[m] = phasefct_init_even[m];
       phasefct_odd_z[m]  = phasefct_init_odd[m];
    }
    for (int k = bx.smallEnd(2)+1; k <= bx.bigEnd(2); k++) {
       mk = (k-bx.smallEnd(2))*num_modes;
       for (int m = 0; m < num_modes; m++) {
           buf[m] = phasefct_even_z[mk-num_modes];
           phasefct_even_z[mk] = phasefct_mult_even_z[m] * phasefct_even_z[mk-num_modes] -
                                 phasefct_mult_odd_z[m]  * phasefct_odd_z[mk-num_modes];
           phasefct_odd_z[mk]  = phasefct_mult_even_z[m] * phasefct_odd_z[mk-num_modes] +
                                 phasefct_mult_odd_z[m]  * buf[m];
           mk = mk + 1;
       }
    }

    // apply forcing in physical space
    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
            mj = (j-bx.smallEnd(1))*num_modes;
            mk = (k-bx.smallEnd(2))*num_modes;

            // pre-compute products of phase factors depending on y- and z-coordinates
            for (int m = 0; m < num_modes; m++) {
                phasefct_yz0[m] = phasefct_even_y[mj] * phasefct_even_z[mk] - phasefct_odd_y[mj]  * phasefct_odd_z[mk];
                phasefct_yz1[m] = phasefct_odd_y[mj]  * phasefct_even_z[mk] + phasefct_even_y[mj] * phasefct_odd_z[mk];
                mj = mj + 1;
                mk = mk + 1;
            }

            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {

                accel[0] = 0.0;
                accel[1] = 0.0;
                accel[2] = 0.0;

                // compute components of acceleration via inverse FT
                for (int n = 0; n < SpectralRank; n++) {
                    mi = (i-bx.smallEnd(0))*num_modes;

                    for (int m = 0; m < num_modes; m++) {
                        // sum up even modes
                        accel[n] = accel[n] + (phasefct_even_x[mi] * phasefct_yz0[m] -
                                               phasefct_odd_x[mi]  * phasefct_yz1[m]) * modes_even[n][m];
                        // sum up odd modes
                        accel[n] = accel[n] - (phasefct_even_x[mi] * phasefct_yz1[m] +
                                               phasefct_odd_x[mi]  * phasefct_yz0[m]) * modes_odd[n][m];
                        mi = mi + 1;
                    }

                    accel[n] = M_SQRT2 * accel[n];
                }
                // add forcing to state     
                state(i,j,k,Xmom_comp) = state(i,j,k,Xmom_comp) + half_dt * state(i,j,k,Density_comp)*accel[0] / a;
                state(i,j,k,Ymom_comp) = state(i,j,k,Ymom_comp) + half_dt * state(i,j,k,Density_comp)*accel[1] / a;
                state(i,j,k,Zmom_comp) = state(i,j,k,Zmom_comp) + half_dt * state(i,j,k,Density_comp)*accel[2] / a;
            }
        }
    }

    // update total energy
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
             state(i,j,k,Eden_comp) = state(i,j,k,Eint_comp) + 0.5  *(state(i,j,k,Xmom_comp)*state(i,j,k,Xmom_comp) +
                                                              state(i,j,k,Ymom_comp)*state(i,j,k,Ymom_comp) +
                                                              state(i,j,k,Zmom_comp)*state(i,j,k,Zmom_comp))/state(i,j,k,Density_comp);
        });

}
