// The Lyman-alpha forest calculation needs the gas density and internal energy
// (and the atomic rates), so it is only built when hydro is present.
#ifndef NO_HYDRO

// The H I fraction comes from the heating/cooling network (eos_hc.H), which is only
// available with HEATCOOL. Without it, the routines below are stubs that abort.
#ifdef HEATCOOL

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <string>

#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include <Nyx.H>
#include <IndexDefines.H>        // Density_comp, Eint_comp
#include <atomic_rates_data.H>   // AtomicRates, atomic_rates_glob, MPROTON
#include <constants_cosmo.H>     // L_unit, density_to_cgs
// nyx_eos_T_given_Re_device comes from eos_hc.H, included through Nyx.H under HEATCOOL.

using namespace amrex;

static void ComputeRaySphereIntersections(
    Real R,
    Real center_x,
    Real center_y,
    Real center_z,
    Vector<Real>& ray_intersect_x,
    Vector<Real>& ray_intersect_y,
    Vector<Real>& ray_intersect_z)
{
    const int ntheta = 100;
    const int nphi   = 100;

    // Equal-area sampling of the full sphere: cos(theta) uniform in (-1,1),
    // phi uniform in [0,2pi). Cell-centered in cos(theta), so no ring sits on a pole.
    Vector<Real> cos_theta(ntheta);
    Vector<Real> phi(nphi);

    for (int i = 0; i < ntheta; ++i) {
        cos_theta[i] = 1.0 - 2.0 * (i + 0.5) / ntheta;
    }

    for (int j = 0; j < nphi; ++j) {
        phi[j] = (2.0 * M_PI) * (j + 0.5) / nphi;
    }

    // Full ntheta x nphi grid of rays; ray index n = i*nphi + j
    const int n_rays = ntheta * nphi;

    ray_intersect_x.resize(n_rays);
    ray_intersect_y.resize(n_rays);
    ray_intersect_z.resize(n_rays);

    for (int i = 0; i < ntheta; ++i) {
        const Real cth = cos_theta[i];
        const Real sth = std::sqrt(1.0 - cth*cth);

        for (int j = 0; j < nphi; ++j) {
            const int n = i*nphi + j;

            ray_intersect_x[n] = center_x + R * sth * std::cos(phi[j]);
            ray_intersect_y[n] = center_y + R * sth * std::sin(phi[j]);
            ray_intersect_z[n] = center_z + R * cth;
        }
    }
}

void
Nyx::InitializeLyAOpticalDepth()
{
    ParmParse pp("nyx");
    pp.query("lya_num_lambda", lya_num_lambda);
    pp.query("lya_lambda_min", lya_lambda_min);
    pp.query("lya_lambda_max", lya_lambda_max);
    pp.query("lya_source_z",   lya_source_z);

    // Validate the inputs needed for the Lyman-alpha forest calculation.
    {
        std::string err;
        if (lya_num_lambda < 2) {
            err += "  nyx.lya_num_lambda must be >= 2 (got " + std::to_string(lya_num_lambda) + ")\n";
        }
        if (!(lya_lambda_min > 0.0)) {
            err += "  nyx.lya_lambda_min must be > 0 (got " + std::to_string(lya_lambda_min) + ")\n";
        }
        if (!(lya_lambda_max > lya_lambda_min)) {
            err += "  nyx.lya_lambda_max (" + std::to_string(lya_lambda_max)
                 + ") must be > nyx.lya_lambda_min (" + std::to_string(lya_lambda_min) + ")\n";
        }
        if (!(lya_source_z > 0.0)) {
            err += "  nyx.lya_source_z must be > 0 (got " + std::to_string(lya_source_z) + ")\n";
        }
        // The light cone only covers final_z <= z <= initial_z.
        if (!(lya_source_z > final_z)) {
            err += "  nyx.lya_source_z (" + std::to_string(lya_source_z)
                 + ") must be > nyx.final_z (" + std::to_string(final_z)
                 + "), otherwise no absorption is ever accumulated\n";
        }
        if (!(lya_source_z <= initial_z)) {
            err += "  nyx.lya_source_z (" + std::to_string(lya_source_z)
                 + ") must be <= nyx.initial_z (" + std::to_string(initial_z) + ")\n";
        }
        if (!(comoving_h > 0.0)) {
            err += "  comoving_h must be > 0 (the distance integral uses H0 = 100*h); "
                   "is this a comoving run?\n";
        }
        if (!do_hydro) {
            err += "  nyx.do_hydro must be on (the gas internal energy and density are needed)\n";
        }
        if (!err.empty()) {
            amrex::Abort("nyx.write_lya_forest_file = true but the Lyman-alpha inputs are invalid:\n" + err);
        }
        if (parent->maxLevel() > 0) {
            amrex::Warning("nyx.write_lya_forest_file: only level-0 data is used for the optical depth; "
                           "refined regions are ignored.");
        }
    }

    // The number of rays does not depend on the radius.
    Vector<Real> rx, ry, rz;
    ComputeRaySphereIntersections(0.0, 0.0, 0.0, 0.0, rx, ry, rz);

    // tau lives on the device for the whole run.
    const std::size_t ntau = static_cast<std::size_t>(rx.size()) * lya_num_lambda;
    lya_tau_of_lambda.resize(ntau);
    Real* tau_ptr = lya_tau_of_lambda.data();
    ParallelFor(static_cast<Long>(ntau), [=] AMREX_GPU_DEVICE (Long n) noexcept
    {
        tau_ptr[n] = 0.0;
    });
    Gpu::streamSynchronize();
}

void
Nyx::AccumulateLyAOpticalDepth(
    const Vector<Real>& ray_intersect_x,
    const Vector<Real>& ray_intersect_y,
    const Vector<Real>& ray_intersect_z)
{

    const Real sigma_alpha = 0.01105;
    const Real nu_0        = 2.4661e15;   // Lyman-alpha line center [Hz]

    // CGS constants for the Doppler width (code-unit versions in constants_cosmo.H)
    const Real c_cgs  = 2.99792458e10;    // [cm/s]
    const Real kB_cgs = 1.38064852e-16;   // [erg/K]
    const Real mH_cgs = 1.6735575e-24;    // [g]
    const Real sqrt_pi = 1.7724538509055159;

    // Emitted (source rest-frame) wavelength grid [Angstrom], num_lambda points, inclusive.
    // lya_source_z is the redshift of the background source (e.g. quasar).
    // Local copies (the lambda below runs on the device and cannot capture class members).
    const Real lam_min = lya_lambda_min;
    const Real lam_max = lya_lambda_max;
    const Real src_z = lya_source_z;

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();
    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

    Real radius_old, radius_new;
    integrate_distance_given_a(a_old, 1.0, radius_old);
    integrate_distance_given_a(a_new, 1.0, radius_new);

    const Real a_mid = 0.5*(a_old + a_new);

    // Gas behind the background source (z > lya_source_z) cannot absorb its light.
    if (1.0/a_mid - 1.0 > src_z) {
        return;
    }

    // radius_* are comoving distances in Mpc (H0 = 100*h km/s/Mpc, c in km/s);
    // multiply by L_unit [cm/Mpc] to get physical path length in cm.

    const int num_lambda = lya_num_lambda;
    Gpu::DeviceVector<Real>& tau_of_lambda = lya_tau_of_lambda;

    const MultiFab& S_new = get_new_data(State_Type);

    // Length of the (periodic) domain, used for the duplicate-box shifts.
    const Real lenx = geom.ProbLength(0);
    const Real leny = geom.ProbLength(1);
    const Real lenz = geom.ProbLength(2);

    const int num_rays = static_cast<int>(ray_intersect_x.size());

    AMREX_ALWAYS_ASSERT(ray_intersect_y.size() == num_rays);
    AMREX_ALWAYS_ASSERT(ray_intersect_z.size() == num_rays);
    AMREX_ALWAYS_ASSERT(num_lambda > 0);
    AMREX_ALWAYS_ASSERT(tau_of_lambda.size() ==
                        static_cast<std::size_t>(num_rays) * num_lambda);

    const auto dx = geom.CellSizeArray();
    const auto dom = geom.Domain();
    const auto dom_lo = dom.smallEnd();
    const auto plo = geom.ProbLoArray();

    // The step spans z_new <= z <= z_old; the profile is averaged over this range.
    const Real z_old = 1.0/a_old - 1.0;
    const Real z_new = 1.0/a_new - 1.0;
    const Real dlambda = (num_lambda > 1)
        ? (lam_max - lam_min)/(num_lambda - 1) : 0.0;

    // Copy the ray positions to the device (tau itself stays resident on the device).
    Gpu::DeviceVector<Real> d_rx(num_rays), d_ry(num_rays), d_rz(num_rays);
    Gpu::copyAsync(Gpu::hostToDevice, ray_intersect_x.begin(), ray_intersect_x.end(), d_rx.begin());
    Gpu::copyAsync(Gpu::hostToDevice, ray_intersect_y.begin(), ray_intersect_y.end(), d_ry.begin());
    Gpu::copyAsync(Gpu::hostToDevice, ray_intersect_z.begin(), ray_intersect_z.end(), d_rz.begin());
    Gpu::streamSynchronize();

    const Real* rx_ptr = d_rx.data();
    const Real* ry_ptr = d_ry.data();
    const Real* rz_ptr = d_rz.data();
    Real* tau_ptr = tau_of_lambda.data();

    // Quantities needed to get the H I number density in each cell.
    const Real h_species_in       = Nyx::h_species;
    const Real gamma_minus_1_in   = Nyx::gamma - 1.0;
    AtomicRates* atomic_rates     = atomic_rates_glob;

    const Real plo0 = plo[0], plo1 = plo[1], plo2 = plo[2];
    const Real dx0 = dx[0], dx1 = dx[1], dx2 = dx[2];
    const int dlo0 = dom_lo[0], dlo1 = dom_lo[1], dlo2 = dom_lo[2];

    // Domain center, and the line-of-sight walk from radius_new to radius_old with a step
    // of half the smallest cell width.
    const Real cen0 = plo0 + 0.5*lenx;
    const Real cen1 = plo1 + 0.5*leny;
    const Real cen2 = plo2 + 0.5*lenz;

    const Real r_near = radius_new;
    const Real seg_len = radius_old - radius_new;
    const Real dl_target = 0.5*amrex::min(dx0, amrex::min(dx1, dx2));
    const int nsamp = amrex::max(1, static_cast<int>(std::ceil(seg_len/dl_target)));
    const Real dl = seg_len/nsamp;

    // One GPU thread per ray per box. Different boxes may run concurrently on different
    // streams and share a ray's tau entries, so the accumulation below is atomic.
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const auto arr = S_new.const_array(mfi);

        ParallelFor(num_rays, [=] AMREX_GPU_DEVICE (int r) noexcept
        {
            // Ray direction from the domain center through the stored point.
            const Real ux0 = rx_ptr[r] - cen0;
            const Real uy0 = ry_ptr[r] - cen1;
            const Real uz0 = rz_ptr[r] - cen2;
            const Real rnorm = std::sqrt(ux0*ux0 + uy0*uy0 + uz0*uz0);
            if (rnorm <= 0.0) { return; }
            const Real ux = ux0/rnorm, uy = uy0/rnorm, uz = uz0/rnorm;

            // Walk the line of sight from radius_new to radius_old in sub-cell steps.
            for (int m = 0; m < nsamp; ++m)
            {
                const Real frac0 = static_cast<Real>(m)   / nsamp;
                const Real frac1 = static_cast<Real>(m+1) / nsamp;

                // Redshift range spanned by this piece (linear in distance over a step).
                const Real za = z_new + (z_old - z_new)*frac0;
                const Real zb = z_new + (z_old - z_new)*frac1;
                const Real z_s = 0.5*(za + zb);

                // Gas behind the background source cannot absorb its light.
                if (z_s > src_z) { continue; }

                // Point at the middle of this piece, wrapped into the periodic domain.
                const Real sdist = r_near + (m + 0.5)*dl;
                const Real px = cen0 + ux*sdist;
                const Real py = cen1 + uy*sdist;
                const Real pz = cen2 + uz*sdist;

                const int dup_i = static_cast<int>(std::floor((px - plo0) / lenx));
                const int dup_j = static_cast<int>(std::floor((py - plo1) / leny));
                const int dup_k = static_cast<int>(std::floor((pz - plo2) / lenz));

                const int i = static_cast<int>(std::floor((px - dup_i*lenx - plo0) / dx0)) + dlo0;
                const int j = static_cast<int>(std::floor((py - dup_j*leny - plo1) / dx1)) + dlo1;
                const int k = static_cast<int>(std::floor((pz - dup_k*lenz - plo2) / dx2)) + dlo2;

                if (!bx.contains(IntVect(AMREX_D_DECL(i, j, k)))) { continue; }

                // Proper H I number density [cm^-3] in this cell.
                const Real rho_comoving = arr(i, j, k, Density_comp);
                const Real eint_per_mass = arr(i, j, k, Eint_comp) / rho_comoving;

                Real T = 0.0, Ne = 0.0;   // outputs of the EOS solve
                Real nh0, nhp, nhe0, nhep, nhepp;
                nyx_eos_T_given_Re_device(
                    atomic_rates, gamma_minus_1_in, h_species_in, 1, 1,
                    &T, &Ne, rho_comoving, eint_per_mass, a_mid,
                    nh0, nhp, nhe0, nhep, nhepp);

                // nh0 is a fraction of n_H; convert to proper density.
                const Real rho_cgs = rho_comoving * density_to_cgs / (a_mid * a_mid * a_mid);
                const Real nh      = rho_cgs * h_species_in / MPROTON;
                const Real n_HI    = nh * nh0;

                // Physical path length of this piece [cm].
                const Real ds_phys = dl * L_unit / (1.0 + z_s);

                // Thermal Doppler width [Hz].
                const Real dnu_D = (nu_0/c_cgs) * std::sqrt(2.0*kB_cgs*T/mH_cgs);

                // Only wavelengths within +-xmax Doppler widths of the line contribute,
                // so restrict to that window of the lambda grid. The resonance sits at
                // lambda_em = (c/nu)*(1+z)/(1+z_src) for nu in [nu_0 - w, nu_0 + w].
                const Real xmax = 6.0;
                const Real w = xmax*dnu_D;
                const Real lam_lo = (c_cgs/(nu_0 + w))*(1.0 + za)/(1.0 + src_z)*1.0e8;
                const Real lam_hi = (c_cgs/(nu_0 - w))*(1.0 + zb)/(1.0 + src_z)*1.0e8;

                int il_lo = static_cast<int>(std::floor((lam_lo - lam_min)/dlambda));
                int il_hi = static_cast<int>(std::ceil ((lam_hi - lam_min)/dlambda));
                il_lo = amrex::max(il_lo, 0);
                il_hi = amrex::min(il_hi, num_lambda - 1);

                for (int ilambda = il_lo; ilambda <= il_hi; ++ilambda)
                {
                    // nu(z) = (c/lambda_em)*(1+z)/(1+z_src): frequency of a photon emitted
                    // at lambda_em that reaches gas at redshift z. Peculiar velocity not
                    // included. Average the profile over the piece's frequency range:
                    //   <phi> = [erf(x2) - erf(x1)] / (2 (nu2 - nu1)),  x = (nu-nu_0)/dnu_D.
                    const Real lambda_em = lam_min + ilambda*dlambda;   // [A]
                    const Real nu_per_opz = c_cgs/(lambda_em*1.0e-8) / (1.0 + src_z);
                    const Real nu1 = nu_per_opz * (1.0 + za);
                    const Real nu2 = nu_per_opz * (1.0 + zb);

                    const Real x1 = (nu1 - nu_0)/dnu_D;
                    const Real x2 = (nu2 - nu_0)/dnu_D;

                    Real phi;
                    if (std::abs(nu2 - nu1) > 1.0e-12*nu1) {
                        // erfc form avoids cancellation in the far wings.
                        const Real d_erf = (x1 + x2 > 0.0)
                            ? (std::erfc(x1) - std::erfc(x2))
                            : (std::erf(x2) - std::erf(x1));
                        phi = 0.5*d_erf/(nu2 - nu1);
                    } else {
                        const Real xm = 0.5*(x1 + x2);
                        phi = std::exp(-xm*xm)/(sqrt_pi*dnu_D);
                    }

                    // Atomic: kernels for different boxes run on different GPU streams
                    // and can update the same (ray, lambda) entry concurrently.
                    Gpu::Atomic::AddNoRet(
                        &tau_ptr[static_cast<std::size_t>(r)*num_lambda + ilambda],
                        sigma_alpha * n_HI * phi * ds_phys);
                }
            }
        });
    }

    Gpu::streamSynchronize();
}

void
Nyx::FinalizeLyAOpticalDepth()
{
    // Each rank only accumulated the cells it owns; sum over ranks so that
    // every rank holds the full optical depth. Call once, at the end of the run.
    // Bring the device-resident tau back to the host first.
    Vector<Real> h_tau(lya_tau_of_lambda.size());
    if (!h_tau.empty()) {
        Gpu::copy(Gpu::deviceToHost, lya_tau_of_lambda.begin(), lya_tau_of_lambda.end(), h_tau.begin());
        ParallelDescriptor::ReduceRealSum(h_tau.data(), static_cast<int>(h_tau.size()));
    }
    // Release the device memory before amrex::Finalize (this is a static member).
    lya_tau_of_lambda.clear();
    lya_tau_of_lambda.shrink_to_fit();

    // Write emitted wavelength and the transmitted flux exp(-tau) of a single ray
    // (nyx.lya_ray_index, default 0), not the ray average.
    if (ParallelDescriptor::IOProcessor())
    {
        int ray_index = 0;
        ParmParse pp("nyx");
        pp.query("lya_ray_index", ray_index);

        const int num_lambda = lya_num_lambda;
        AMREX_ALWAYS_ASSERT(num_lambda > 0);
        const std::size_t num_rays = h_tau.size() / num_lambda;
        AMREX_ALWAYS_ASSERT(num_rays > 0);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ray_index >= 0 && static_cast<std::size_t>(ray_index) < num_rays,
                                         "nyx.lya_ray_index is out of range");
        const Real dlambda = (num_lambda > 1)
            ? (lya_lambda_max - lya_lambda_min)/(num_lambda - 1) : 0.0;

        std::ofstream ofs("LyA_forest.txt");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ofs.good(), "Could not open LyA_forest.txt");

        ofs << "# columns: lambda_emitted[Angstrom]  flux=exp(-tau) for a single ray\n";
        ofs << "# ray_index = " << ray_index << "  num_rays = " << num_rays
            << "  num_lambda = " << num_lambda << "\n";
        ofs << std::scientific << std::setprecision(8);

        const std::size_t offset = static_cast<std::size_t>(ray_index) * num_lambda;
        for (int ilambda = 0; ilambda < num_lambda; ++ilambda)
        {
            ofs << (lya_lambda_min + ilambda*dlambda)
                << " " << std::exp(-h_tau[offset + ilambda]) << "\n";
        }
    }
}

void
Nyx::ComputeLymanAlphaForestFractionvsWavelength()
{
    // Find the comoving radii of the shell radius_inner and radius_outer
    // Loop over the rays and find their intersection with the
    // shell center sphere as (x,y,z) and store it in
    // ray_intersect_x, ray_intersect_y, ray_intersect_z

    // Every rank gets allocated with num_rays*num_lambda
    // Loop from finest level to 0
    // If not the finest level create the fine_mask
    // Loop over the boxes Mfiter
    // Loop over the rays and find if their intersection (x,y,z) is the in box. If yes, find the cell.
    // If not the finest level, check if fine_mask is not 0, then do the integrand
    // Loop over just that cell with ParallelFor
    // Find the integrand in that cell and add to the integral

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();
    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

    Real radius_old, radius_new;
    integrate_distance_given_a(a_old, 1.0, radius_old);
    integrate_distance_given_a(a_new, 1.0, radius_new);

    Print() << "The values of radius_old and radius_new are " << radius_old << " " << radius_new << std::endl;

    const Real comoving_radius_intersect_sphere = 0.5 * (radius_old + radius_new);

    Vector<Real> ray_intersect_x;
    Vector<Real> ray_intersect_y;
    Vector<Real> ray_intersect_z;

    // Rays emerge from the center of the domain.
    const Real center_x = 0.5 * (geom.ProbLo(0) + geom.ProbHi(0));
    const Real center_y = 0.5 * (geom.ProbLo(1) + geom.ProbHi(1));
    const Real center_z = 0.5 * (geom.ProbLo(2) + geom.ProbHi(2));

    ComputeRaySphereIntersections(comoving_radius_intersect_sphere,
                                  center_x, center_y, center_z,
                                  ray_intersect_x,
                                  ray_intersect_y,
                                  ray_intersect_z);

    AccumulateLyAOpticalDepth(ray_intersect_x,
                              ray_intersect_y,
                              ray_intersect_z);
}

#else // !HEATCOOL

#include <AMReX_Abort.H>
#include <Nyx.H>

using namespace amrex;

namespace {
    constexpr const char* lya_heatcool_msg =
        "The Lyman-alpha forest (nyx.write_lya_forest_file) requires a build with USE_HEATCOOL = TRUE";
}

void Nyx::InitializeLyAOpticalDepth()  { amrex::Abort(lya_heatcool_msg); }

void Nyx::AccumulateLyAOpticalDepth(const Vector<Real>&,
                                    const Vector<Real>&,
                                    const Vector<Real>&)
{
    amrex::Abort(lya_heatcool_msg);
}

void Nyx::FinalizeLyAOpticalDepth()    { amrex::Abort(lya_heatcool_msg); }

void Nyx::ComputeLymanAlphaForestFractionvsWavelength() { amrex::Abort(lya_heatcool_msg); }

#endif // HEATCOOL

#endif // NO_HYDRO
