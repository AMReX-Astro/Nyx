[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5059767.svg)](https://doi.org/10.5281/zenodo.5059767)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03068/status.svg)](https://doi.org/10.21105/joss.03068)
[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)


![Nyx](https://github.com/AMReX-Astro/Nyx/blob/development/Util/banner.jpeg)

*An adaptive mesh, massively-parallel, cosmological simulation code*

******

## About

Nyx code solves equations of compressible hydrodynamics on an adaptive grid
hierarchy coupled with an N-body treatment of dark matter. The gas dynamics in
Nyx uses a finite volume methodology on a set of 3-D Eulerian grids;
dark matter is represented as discrete particles moving under the influence of
gravity. Particles are evolved via a particle-mesh method, using Cloud-in-Cell
deposition/interpolation scheme. Both baryonic and dark matter contribute to
the gravitational field. In addition, Nyx includes physics needed to
accurately model the intergalactic medium: in optically thin limit and assuming
ionization equilibrium, the code calculates heating and cooling processes of the
primordial-composition gas in an ionizing ultraviolet background radiation field.
Additional physics capabilities are under development.

While Nyx can run on any Linux system in general, we particularly focus on supercomputer systems.
Nyx is parallelized with MPI + X, where X can be OpenMP on multicore architectures and
CUDA/HIP/DPC++ on hybrid CPU/GPU architectures.
In the OpenMP regime, Nyx has been successfully run at parallel concurrency
of up to 2,097,152 on NERSC's Cori-KNL. With Cuda implementation, it was run on up to
13,824 GPUs on OLCF's Summit.

More information on Nyx can be found at the [main web page](http://amrex-astro.github.io/Nyx/) and
the [online documentation](https://amrex-astro.github.io/Nyx/docs_html/).

## Standards and dependencies

To compile the code we require C++11 compliant compilers that support MPI-2 or
higher implementation.  If threads or accelerators are used, we require 
OpenMP 4.5 or higher, Cuda 9 or higher, or HIP-Clang.
To use Nyx, you also need [AMReX](https://github.com/AMReX-codes/amrex).

For example, to compile the Lyman alpha (LyA) executable on Summit:
```sh
$ module load gcc/6.4.0 cuda/11.0.3

$ git clone https://github.com/AMReX-Codes/amrex.git
$ git clone https://github.com/AMReX-astro/Nyx.git

$ cd Nyx/Exec/LyA
$ make -j 12 USE_CUDA=TRUE
```

See the [Getting Started section](https://amrex-astro.github.io/Nyx/docs_html/NyxGettingStarted.html) for more information.

## Development model

Please see CONTRIBUTING.md for details on how to contribute to AMReX.

## Outputs

Nyx outputs certain global diagnostics at each timestep and plot files at regular
intervals, or at user-specified redshifts. Visualization packages
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit),
[Paraview](https://www.paraview.org/),
[yt](http://yt-project.org/),
and [Amrvis](https://github.com/AMReX-Codes/amrvis)
have built-in support for the AMReX file format used by Nyx.

In addition, Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber
uses topological methods to construct merge trees of scalar fields, which is in
turn used to find halos.  Gimlet computes a variety of quantities
related to the Lyman-alpha forest science.  These suites are fully MPI-parallel and can
be run either *in situ* or in-transit, or with a combination of both
(see [Friesen et al. 2016](https://comp-astrophys-cosmol.springeropen.com/articles/10.1186/s40668-016-0017-2)).

## License

Nyx Copyright (c) 2017, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy).  All rights
reserved.

Details of the license can be found in [license.txt](license.txt) file.

If you have questions about your rights to use or distribute this software, 
please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

## Contact

For questions, comments, suggestions, contact Jean Sexton (JMSexton@lbl.gov)
or Zarija Lukic (zarija@lbl.gov).
