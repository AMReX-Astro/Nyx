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
    affiliation: 2
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Chris Daley
    affiliation: 3
  - name: Brian Friesen
    orcid: 0000-0002-1572-1631
    affiliation: 3
  - name: Andrew Myers
    orcid: 0000-0001-8427-8330
    affiliation: 1
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Computational Cosmology Center, Lawrence Berkeley National Laboratory
    index: 2
  - name: National Energy Research Scientific Computing Center (NERSC), Berkeley, CA, USA
    index: 3
date: February 2021
bibliography: paper.bib
---

# Summary
Nyx is a highly parallel, adaptive mesh, finite-volume
N-body compressible hydrodynamics solver for cosmological simulations.
It has been used to simulate different cosmological scenarios with
a recent focus on the intergalactic medium and Lyman alpha forest.
Together, Nyx, the compressible astrophysical simulation code, Castro [@castro],
and the low Mach number code MAESTROeX [@maestroex], make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance-portable
astrophysical simulation codes.
Other examples of cosmological simulation research codes include Enzo [@Bryan2014], Enzo-P/Cello [@Bordner2018; @Norman2018], RAMSES [@Teyssier2002], ART [@Kravtsov1997],
FLASH [@Fryxell2000], Cholla [@Villasenor2021], as well as Gadget [@Springel2020], Gasoline [@Wadsley2017], Arepo [@Weinberger2020], Gizmo [@Hopkins2014], and SWIFT [@Schaller2016].

The core hydrodynamics solver in Nyx [@nyx-paper1] is based on the
directionally unsplit corner transport upwind method of @ctu with
piecewise parabolic reconstruction [@ppm].  In Nyx, we have
several modes of coupling the stiff heating-cooling source terms to the hydro.
The simplest method is the traditional operator splitting approach,
using Strang splitting [@strang1968] to achieve second-order accuracy in time. However,
this coupling can break down, and we have an alternative to Strang splitting
based on spectral deferred corrections (SDC), a method
that aims to prevent the hydro and stiff source terms from becoming decoupled.
The simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the source terms and hydro, still to
second-order accuracy in time [@simple_sdc].

Nyx has a set of additional physics necessary to model the intergalactic medium
using heating-cooling source terms.
The code follows the abundance of six species: neutral and ionized hydrogen,
neutral, once and twice ionized helium, and free electrons. For these species,
all relevant atomic processes - ionization, recombination, and free-free transitions are
modeled in the code. Heating and cooling source terms are calculated using
a sub-cycled approach in order to avoid running the whole code on a short, cooling timescale.
Cosmological reionization is accounted for via a spatially uniform,
but time-varying ultraviolet background (UVB) radiation field,
inputted to the code as a list of photoionization and photoheating
rates that vary with redshift. Nyx also has the capability to model flash reionization
as well as inhomogeneous reionization [@onorbe2019].

The evolution of baryonic gas is coupled with an N-body treatment of the dark
matter in an expanding universe. The mesh-based hydrodynamical baryonic gas
evolution is coupled through gravity to the particle-based representation of
dark matter. The dark matter particles are moved with a move-kick-drift algorithm
[@movekickdrift]. The Poisson equation for self-gravity of the baryonic gas and dark
matter is solved using the geometric multigrid method. Nyx simulations can optionally
model neutrino particle effects and active galactic nuclei feedback.

Nyx is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is written in C++.
AMR levels are advanced at their own timestep (sub-cycling)
and jumps by factors of 2 and 4 are supported between levels.  We use
MPI to distribute AMR grids across nodes and use logical tiling with
OpenMP to divide a grid across threads for multi-core CPU machines
(exposing coarse-grained parallelism) or CUDA/HIP/DPC++ to spread the work across
GPU threads on GPU-based machines (fine-grained parallelism).  All of
the core physics can run on GPUs and have been shown to scale well.
For performance portability, we use the same source code
for both CPUs and GPUs; additionally, we implement our parallel loops
in an abstraction
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
helps in Nyx continually gaining optimizations for new architectures
and various operating systems.

At present, Nyx is mostly used to model the cosmological evolution of the
intergalactic medium (IGM) - a reservoir of low density gas that fills the space
between galaxies. Different physical effects, ranging from the nature of
dark matter to different reionization scenarios related to the radiation
from star-forming galaxies, set the observable
properties of IGM, making it a powerful probe of cosmology and astrophysics.
But in order to extract scientific insights, confronting observations of
the IGM (usually through the Lyman alpha forest) against simulated models is a necessity, and that is where Nyx steps in.
Incoming observations, for example, Dark Energy Spectroscopic Instrument (DESI) or
high-resolution spectrographs like the one mounted on the Keck telescope, are noticeably
reducing the statistical errors of measurements; the community needs tools
capable of producing model universes that are of the same or better accuracy as observations.
Nyx's scalability allows modeling of
representative cosmological volumes, while maintaining the resolution needed
to resolve small fluctuations in the intergalactic gas. Nyx includes physics to allow simulations
of different cosmological and reionization scenarios, enabling users to produce
mock universes for a variety of physically relevant models.

Our main targets are high-performance computer architectures and massively parallel simulations
needed for cosmological and astrophysical research.
Given these targets, Nyx's computational optimizations focus on large simulations on HPC systems,
although smaller simulations can be run on Linux distributions and macOS using AMReX's
build system support.


# Acknowledgements

The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.
Nyx development was further supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.

# References

