.. _developers-amrex-basics:

AMReX basics
================================

Nyx is built on the Adaptive Mesh Refinement (AMR) library `AMReX <https://github.com/AMReX-Codes/amrex>`__. This section provides a very sporadic description of the main AMReX classes and concepts relevant for Nyx.
For more details, please visit
the AMReX basics `doc page <https://amrex-codes.github.io/amrex/docs_html/Basics.html>`__, and the rest of the AMReX documentation.

* ``amrex::Box``: Dimension-dependent lower and upper indices defining a rectangular volume in 3D (or surface in 2D) in the index space. ``Box`` is a lightweight meta-data class, with useful member functions.

* ``amrex::BoxArray``: Collection of ``Box`` on a single AMR level. The information of which MPI rank owns which ``Box`` in a ``BoxArray`` is in ``DistributionMapping``.

* ``amrex::FArrayBox``: Fortran-ordered array of floating-point ``amrex::Real`` elements defined on a ``Box``. A ``FArrayBox`` can represent scalar data or vector data, with ``ncomp`` components.

* ``amrex::MultiFab``: Collection of `FAB` (= ``FArrayBox``) on a single AMR level, distributed over MPI ranks. The concept of `ghost cells` is defined at the ``MultiFab`` level.

* ``amrex::ParticleContainer``: A collection of particles; in Nyx the main ParticleContainer contains dark matter particles.
Particles in a ``ParticleContainer`` are organized per ``Box``. Particles in a ``Box`` are organized per tile
(this feature is off when running on GPU). 
Particles within a tile are stored in several structures, each being contiguous in memory: (i) an Array-Of-Struct (AoS) (often called `data`, they are the 3D position, the particle ID and the index of the CPU owning the particle), where the Struct is an ``amrex::Particle`` and (ii) Struct-Of-Arrays (SoA) for extra variables (often called ``attribs``).

The simulation domain is typically decomposed into multiple boxes, and each MPI rank owns, and performs operations on,
the fields and particles defined on the boxes assigned to that rank. However, every rank does have the full metadata, i.e
the list of all boxes and which ranks the data on those boxes is assigned to (the ``DistributionMapping``).
For convenience, AMReX provides iterators, to easily iterate over all the ``FArrayBox`` (with optional logical tiling)
in a ``MultiFab`` that are owned by that MPI rank (``MFIter``), or over all particles in a ``ParticleContainer`` on a per-box basis (``ParIter``).
These are respectively done in loops such as:

.. code-block:: cpp

   // mf is a MultiFab
   for ( amrex::MFIter mfi(mf, TilingIfNotGpu()); mfi.isValid(); ++mfi ) { ... }

and

.. code-block:: cpp

   // *this is a pointer to a ParticleContainer
   for (ParIter pti(*this, lev); pti.isValid(); ++pti) { ... }
