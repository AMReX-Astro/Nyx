*******
Preface
*******

Nyx is a adaptive mesh, N-body/hydro code for cosmological simulation on massively parallel
computers.  It couples the compressible hydrodynamic equations on a grid with a particle represenation
of dark matter.

The major capabilities:

  * 3-dimensional unsplit, 2nd-order hydrodynamics

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x between levels

  * full Poisson gravity (with triply periodic boundary conditions)

  * hybrid parallelization strategy with MPI + X, where X = OpenMP on multicore architectures
and CUDA/HIP/DPC++ on hybrid CPU/GPU architectures.

Nyx uses an Eulerian grid for the hydrodynamics solver and incorporates adaptive mesh refinement (AMR).
Our approach to AMR uses a nested hierarchy of logically-rectangular grids with simultaneous 
refinement in both space and time, utilizing the AMReX library.

