# Nyx
*an adaptive mesh, cosmological simulation code*

*******************************************************

PLEASE NOTE:  Nyx has been migrated from BoxLib to AMReX.

If you prefer to run with the BoxLib-based version, then 
you can use the "boxlib" branch (which will no longer be updated).

If you would like to use the latest version of Nyx, you will
need to download AMReX.  Please email ASAlmgren@lbl.gov to get
access to AMReX (it will become public as soon as it makes
it through the licensing process).


`Nyx` is an N-body hydro cosmological simulation code with 
adaptive mesh refinement (AMR).  `Nyx` uses an unsplit PPM approach 
for the hydrodynamics, a particle representation of dark matter, 
and multigrid to solve the Poisson equation for self-gravity.

`Nyx` is parallelized with MPI + OpenMP.

More information on Nyx can be found here:

http://boxlib-codes.github.io/Nyx/

*******************************************************

## Getting Started

To build `Nyx`, you need a copy of the `AMReX` library  --
please email ASAlmgren@lbl.gov to get access to AMReX 
(it will become public as soon as it makes it through the licensing process).

There is a User's Guide in `Nyx/UsersGuide/` (type `make` to build
from LaTeX source) that will guide you through running your first
problem.  

## Development Model

New features are committed to the `development` branch.  Nightly
regression testing is used to ensure that no answers change (or if
they do, that the changes were expected).  No changes should ever
be pushed directly into `master`.

On the first workday of each month, we perform a merge of
`development` into `master`, in coordination with `amrex`. 
For this merge to take place, we need to be passing the regression tests.  
To accommodate this need, we close the
merge window into `development` a few days before the merge day.
While the merge window is closed, only bug fixes should be pushed into
`development`.  Once the merge from `development` -> `master` is done,
the merge window reopens.

## Physics Included
For the description of the N-body and adiabatic hydro algorithms in Nyx, see
Almgren, Bell, Lijewski, Lukic & Van Andel (2013), ApJ, 765, 39:
http://adsabs.harvard.edu/abs/2013ApJ...765...39A

For the reaction and thermal rates of the primordial chemical composition gas 
(and convergence tests in the context of the Lyman-alpha forest), see
Lukic, Stark, Nugent, White, Meiksin & Almgren (2015), MNRAS, 446, 3697:
http://adsabs.harvard.edu/abs/2015MNRAS.446.3697L

For the synthesis model of the UV background, 
which provides the photo-ionization and photo-heating rates, 
we include the prescription of Onorbe, Hennawi & Lukic (2016):
http://adsabs.harvard.edu/abs/2016arXiv160704218O

## Post-processing

Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber uses
topological methods to construct merge trees of scalar fields, which Nyx in
turn uses to find dark matter halos. Gimlet computes a variety of quantities
about the simulation, including optical depths, Lyman-alpha fluxes, power
spectra (both 1-D "line-of-sight" as well as fully 3-D), and probability
distribution functions. These suites are fully MPI-parallel and can be run
either "in situ" or "in-transit", or with a combination of both. A detailed
description of their usage is provided in the Nyx User Guide.

## Mailing List

For questions on how to use `Nyx`, email Ann Almgren at ASAlmgren@lbl.gov

For more information about the science being done with `Nyx`, 
email Zarija Lukic at Zarija@lbl.gov

You can also subscribe to the nyx-help mailing list at google groups:

https://groups.google.com/forum/#!forum/nyx-help
