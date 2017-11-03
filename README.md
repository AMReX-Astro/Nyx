# Nyx
*an adaptive mesh, massively-parallel, cosmological simulation code*

******************************************************

## About

Nyx code solves the equations of compressible hydrodynamics on an adaptive grid
hierarchy coupled with an N-body treatment of dark matter. The gasdynamics in
Nyx uses a finite volume methodology on an adaptive set of 3-D Eulerian grids;
dark matter is represented as discrete particles moving under the influence of
gravity. Particles are evolved via a particle-mesh method, using Cloud-in-Cell
deposition/interpolation scheme. Both baryonic and dark matter contribute to
the gravitational field. In addition, Nyx currently includes physics needed to
accurately model the intergalactic medium: in optically thin limit and ionization
equilibrium, the code solves for the abundance of six primordial species: neutral and ionized
hydrogen, neutral, once and twice ionized helium, and free electrons. For these
species, all relevant atomic processes -- ionization, recombination, and free-free
transition are included, as well as different prescriptions for the ionizing ultra-violet
background (UVB). Additional physics capabilities are currently under development.

Nyx is parallelized with MPI + OpenMP, and has been run at parallel concurrency
of up to 2,097,152 (on NERSC's Cori).

More information on Nyx can be found here:
http://amrex-astro.github.io/Nyx/

If you prefer to run depreciated BoxLib-based version of the code, then 
you can use the `boxlib` branch (which will no longer be updated).

******************************************************

## Getting Started

To compile the code, we require Fortran 2003 and C++11 compliant compilers that
support (if parallelism is desired) OpenMP 4.5 or better, and MPI-2 or higher implementation.

To use Nyx, you also need to download AMReX from
https://github.com/AMReX-codes/amrex
which is the only required software dependency.

There is a User's Guide in `Nyx/UsersGuide/` (type `make` to build
from LaTeX source) that will guide you through running your first
problem.

## Development Model

New features are committed to the `development` branch.  We use nightly
regression testing to ensure that no answers change (or if they do, that
the changes were expected).  No changes should be pushed directly into
`master`. Approximately once a month, we perform a merge of `development`
into `master`.

## Physics References

For the description of the N-body and adiabatic hydro algorithms in Nyx, see
Almgren, Bell, Lijewski, Lukic & Van Andel (2013), ApJ, 765, 39:
http://adsabs.harvard.edu/abs/2013ApJ...765...39A

For the reaction and thermal rates of the primordial chemical composition gas 
(and convergence tests in the context of the Lyman-alpha forest), see
Lukic, Stark, Nugent, White, Meiksin & Almgren (2015), MNRAS, 446, 3697:
http://adsabs.harvard.edu/abs/2015MNRAS.446.3697L

For considerations regarding the synthesis model of the UV background, 
which provides the photo-ionization and photo-heating rates, see Onorbe,
Hennawi & Lukic (2017), ApJ, 847, 63:
http://adsabs.harvard.edu/abs/2017ApJ...847...63O

## Output

Nyx outputs certain global diagnostics at each timestep and plot files at regular
intervals, or at user-specified redshifts. The visualization packages VisIt, Paraview
and yt have built-in support for the AMReX file format used by Nyx. For more
details, see the Nyx User Guide.

In addition, Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber
uses topological methods to construct merge trees of scalar fields, which Nyx in
turn uses to find halos. Gimlet computes a variety of quantities
related to the Lyman-alpha forest science. These suites are fully MPI-parallel and can
be run either "in situ" or "in-transit", or with a combination of both. A detailed
description of their usage is provided in the Nyx User Guide.

## Contact

For questions on how to use Nyx, email Ann Almgren at ASAlmgren@lbl.gov

For more information about the science with Nyx, email Zarija Lukic at zarija@lbl.gov
