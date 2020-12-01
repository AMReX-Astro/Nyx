# Nyx

[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)

*An adaptive mesh, massively-parallel, cosmological simulation code*

******

### About

Nyx code solves equations of compressible hydrodynamics on an adaptive grid
hierarchy coupled with an N-body treatment of dark matter. The gas dynamics in
Nyx uses a finite volume methodology on an adaptive set of 3-D Eulerian grids;
dark matter is represented as discrete particles moving under the influence of
gravity. Particles are evolved via a particle-mesh method, using Cloud-in-Cell
deposition/interpolation scheme. Both baryonic and dark matter contribute to
the gravitational field. In addition, Nyx currently includes physics needed to
accurately model the intergalactic medium: in optically thin limit and assuming
ionization equilibrium, the code calculates heating and cooling processes of the
primordial-composition gas in an ionizing ultraviolet background radiation field.
Additional physics capabilities are under development.

While Nyx can run on any Linux system in general, we particularly focus on supercomputer systems.
Nyx is parallelized with MPI + X, where "X" can be OpenMP, CUDA, or HIP (DPC++ implementation
is ongoing). In the OpenMP regime, Nyx and has been successfully run at parallel concurrency
of up to 2,097,152 (on NERSC's Cori-KNL). With Cuda implementation, it was ran on up to
13,824 GPUs (on OLCF's Summit).

More information on Nyx can be found at the [main web page](http://amrex-astro.github.io/Nyx/) and
the [documentation section](https://amrex-astro.github.io/Nyx/docs_html/).


### Standards and dependencies

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

See the online documenation for further information.


### Development model

The `development` branch in also the main branch.  We use nightly
regression testing to ensure that no answers change (or if they do, that
the changes were expected). Contributions are welcomed and should be done via pull requests.
A pull request should be generated from your fork of Nyx and should target
the `development` branch.


### Outputs

Nyx outputs certain global diagnostics at each timestep and plot files at regular
intervals, or at user-specified redshifts. Visualization packages
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit),
[Paraview](https://www.paraview.org/)
and [yt](http://yt-project.org/)
have built-in support for the AMReX file format used by Nyx.

In addition, Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber
uses topological methods to construct merge trees of scalar fields, which is in
turn used to find halos.  Gimlet computes a variety of quantities
related to the Lyman-alpha forest science.  These suites are fully MPI-parallel and can
be run either "in situ" or "in-transit", or with a combination of both
(see [Friesen et al. 2016](https://comp-astrophys-cosmol.springeropen.com/articles/10.1186/s40668-016-0017-2)).


### License
Nyx is released under the LBL's modified BSD license, see the [license.txt](license.txt) file for details.


### Contact

For questions, comments, suggestions, contact Jean Sexton (JMSexton@lbl.gov)
or Zarija Lukic (zarija@lbl.gov).
