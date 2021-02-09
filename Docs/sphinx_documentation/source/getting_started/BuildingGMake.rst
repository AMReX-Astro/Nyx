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
   For some executables, namely LyA, AMR-density, AMR-zoom, and LyA_Neutrinos Nyx requires you to also install a matching sundials installation. To install a stand-alone copy of Sundials, see :ref:`sundials`

#. From the directory in which you checked out Nyx, type

   ::

       cd Nyx/Exec/MiniSB

   This will put you into a directory in which you can run a small
   version of the Santa Barbara problem. This will then be your compile directory.

   or

   ::

       cd Nyx/Exec/LyA

   This will put you into a directory in which you can run the Lyman-:math:`\alpha` problem. This will then be your compile directory.

#. In your compile directory, edit the GNUmakefile, and set

   ``COMP = your favorite compiler (e.g, gnu, Intel)``

   ``DEBUG = FALSE``

   We like ``COMP = gnu``.

   More information on the AMReX GNU Make setup can be found
   `here <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.


#. Now type “make”. The resulting executable will look something like
   “Nyx3d.Linux.gnu.ex”, which means this is a 3-d version of the code,
   made on a Linux machine, with ``COMP = gnu``.

   Note that if you build with ``USE_MPI = TRUE`` in the GNUMakefile, then the
   name of the code will be something like “Nyx3d.Linux.gnu.MPI.ex”
