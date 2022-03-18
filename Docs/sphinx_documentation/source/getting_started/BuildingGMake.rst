Building Nyx with GNU Make
============================

Nyx is built on top of the AMReX framework so you must
download AMReX in order to build Nyx with GNU Make.

.. raw:: latex

   \vspace{.1in}

#. Clone/fork the AMReX repository:

   ::

       git clone https://github.com/AMReX-Codes/amrex

   You will want to periodically update AMReX by typing

   ::

       git pull

   in the ``amrex/`` directory.

.. note::
   When you check out AMReX (and Nyx), you will get the development
   branch.  Active development is done on this branch; monthly
   tags are made of this version that are compatible with the same
   monthly tag of AMReX itself.

#. Set the environment variable ``AMREX_HOME`` to point to
   your local copy of the AMReX repository.
   You can add this to your ``.bashrc`` as:

   ::

       export AMREX_HOME=/path/to/amrex/

   or to your ``.cshrc`` as:

   ::

       setenv AMREX_HOME /path/to/amrex/

   where ``/path/to/amrex/`` is the full path to the
   amrex directory.

#. Clone/fork the Nyx repository:

   ::

       git clone https://github.com/AMReX-Astro/Nyx

   As with AMReX development on Nyx is done in the
   ``development`` branch, so you should work there if you want
   the latest source.


#. Choose which executable to compile, for example MiniSB or LyA.

   MiniSB is a small version of the Santa Barbara problem, and LyA is a Lyman-:math:`\alpha` 
   forest simulation for investigating the large-scale structure formation of the universe.

   For some executables, namely LyA, AMR-density, AMR-zoom, and LyA_Neutrinos, Nyx uses a more complicated model for the heating-cooling.
   This requires you to also install a matching sundials installation to support the ODE solve. To install a stand-alone copy of Sundials, see :ref:`Compiling Nyx with SUNDIALS 5<sundials>`

#. From the directory in which you checked out Nyx, change directory to your build directory by typing for the small Santa-Barbara problem:

   ::

       cd Nyx/Exec/MiniSB

   or for the Lyman-:math:`\alpha` problem:

   ::

       cd Nyx/Exec/LyA

#. In your build directory, edit the GNUmakefile, and set

   ``COMP = your favorite compiler (e.g, gnu, Intel)``

   ``DEBUG = FALSE``

   We like ``COMP = gnu``.

   More information on the AMReX GNU Make setup can be found
   `here <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.

   All executables should work with MPI+CUDA by setting ``USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=TRUE``.

   HIP and DPC++ builds are under development and can be tested by compiling with ``USE_HIP=TRUE``  and ``USE_DPCPP=TRUE``  respectively.

   .. note::
      For executables with ``USE_HEATCOOL=TRUE`` in their GNUmakefile, a matching Sundials implementation is required. If Sundials is built with ``-DSUNDIALS_BUILD_PACKAGE_FUSED_KERNELS=ON``, Nyx should be built with ``USE_FUSED=TRUE``.
      The flag ``USE_FUSED`` tells the Nyx compile whether you compiled Sundials with fused cuda kernels. The default assumption is that non-cuda Nyx compiles set ``USE_FUSED=FALSE`` to match Sundials being built without fused cuda kernels.
      Starting with Sundials version 5.7.0, set ``USE_SUNDIALS_SUNMEMORY=TRUE`` to compile the optional Sundials SunMemory to AMReX Arena interface for GPU memory reuse.

#. Now type “make”. The resulting executable will look something like
   “Nyx3d.Linux.gnu.ex”, which means this is a 3-d version of the code,
   made on a Linux machine, with ``COMP = gnu``.

   Note that if you build with ``USE_MPI = TRUE`` in the GNUMakefile, then the
   name of the code will be something like “Nyx3d.Linux.gnu.MPI.ex”
