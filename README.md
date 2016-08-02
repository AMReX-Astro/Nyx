# Nyx
*an adaptive mesh, cosmological simulation code*

`Nyx` is an N-body hydro cosmological simulation code with 
adaptive mesh refinement (AMR).  `Nyx` uses an unsplit PPM approach 
for the hydrodynamics, a particle representation of dark matter, 
and multigrid to solve the Poisson equation for self-gravity.

`Nyx` is parallelized with MPI + OpenMP.

More information on Nyx can be found here:

http://boxlib-codes.github.io/Nyx/


## Getting Started

To build `Nyx`, you need a copy of the `BoxLib` library:

https://github.com/BoxLib-Codes/BoxLib.git

There is a User's Guide in `Nyx/Docs/` (type `make` to build
from LaTeX source) that will guide you through running your first
problem.  

## Development Model:

New features are committed to the `development` branch.  Nightly
regression testing is used to ensure that no answers change (or if
they do, that the changes were expected).  No changes should ever
be pushed directly into `master`.

On the first workday of each month, we perform a merge of
`development` into `master`, in coordination with `BoxLib`. 
For this merge to take place, we need to be passing the regression tests.  
To accommodate this need, we close the
merge window into `development` a few days before the merge day.
While the merge window is closed, only bug fixes should be pushed into
`development`.  Once the merge from `development` -> `master` is done,
the merge window reopens.

## Mailing list

For questions on how to use `Nyx`, email Ann Almgren at ASAlmgren@lbl.gov

For more information about the science being done with `Nyx`, 
email Zarija Lukic at Zarija@lbl.gov
