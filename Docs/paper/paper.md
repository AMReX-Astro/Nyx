---
title: 'Nyx: A Massively Parallel AMR Code for Computational Cosmology'
tags:
  - C++
  - cosmology
  - hydrodynamics
  - dark matter
  - N-body
authors:
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    affiliation: 1
  - name: Zarija Lukic
    orcid: 
    affiliation: 2
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Andrew Myers
    orcid: 0000-0001-8427-8330
    affiliation: 1
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Computational Cosmology Center, Lawrence Berkeley National Laboratory
    index: 2
date: November 2020
bibliography: paper.bib
---

# Summary
``Nyx`` is a highly parallel, adaptive mesh, finite-volume 
N-body compressible hydrodynamics solver for cosmological simulation.  
It has been used to simulate different cosmological scenarios with
a recent focus on the Lyman Alpha forest.
Together, Nyx, the compressible astrophysical simulation code, Castro [@castro], 
and the low Mach number code MAESTROeX [@maestroex], make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance-portable 
astrophysical simulation codes.

The core hydrodynamics solver in Nyx [@nyx] is based on the
directionally unsplit corner transport upwind method of @ctu with
piecewise parabolic reconstruction [@ppm].  In Nyx, we have
several modes of coupling the stiff heating-cooling terms to the hydro.  
The simplest method is the traditional operator splitting approach, 
using Strang splitting to achieve second-order in time.  However, 
this coupling can break down, and we have an alternative to Strang splitting
based on spectral deferred corrections (SDC), a method
that aims to prevent the hydro and stiff source terms from becoming decoupled.  
The simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the source terms and hydro, still to
second order [@simple_sdc].

The evolution of baryonic gas coupled with an N-body treatment of the dark
matter in an expanding universe. The mesh-based hydrodynamical baryonic gas
evolution is coupled through gravity to the particle-based representation of
dark matter. The dark matter particles are moved with a move-kick-drift algorithm
[@movekickdrift]. The Poisson equation for self-gravity of the baryonic gas and dark
matter is solved using geometric multigrid method. Nyx simulations can optionally
model neutrino particle effects and active galactic nuclei feedback.

Nyx is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is written in C++.
AMR levels are advanced at their own timestep (subcycling)
and jumps by factors of 2 and 4 are supported between levels.  We use
MPI to distribute AMR grids across nodes and use logical tiling with
OpenMP to divide a grid across threads for multi-core CPU machines
(exposing coarse-grained parallelism) or CUDA/HIP/DPC++ to spread the work across
GPU threads on GPU-based machines (fine-grained parallelism).  All of
the core physics can run on GPUs and have been shown to scale well.
For performance portability, we use the same source code
for both CPUs and GPUs, and implement our parallel loops in an abstraction
layer provided by AMReX. An abstract parallel for loop accepts as arguments
a range of indices and the body of the loop to execute for a given index,
and the AMReX backend dispatches the work appropriately (e.g., one cell per
GPU thread). This strategy is similar to the way the Kokkos [@Kokkos] and
RAJA [@RAJA] abstraction models provide performance portability in C++.

# Statement of Need

While there are a number of cosmological simulation codes, Nyx
offers a few unique features.  The original motivation for developing
Nyx was to build a simulation code based on a modern,
well-supported AMR library (BoxLib which evolved to AMReX), using
unsplit integration techniques.
The large developer community contributing to AMReX
(representing a large number of application codes across various domains)
results in Nyx continually gaining optimizations for new architectures.  

# Acknowledgements

The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.   
Nyx development was further supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.  

# References

