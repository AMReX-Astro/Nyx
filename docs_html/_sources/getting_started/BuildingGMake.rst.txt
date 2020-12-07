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


#. From the directory in which you checked out Nyx, type

   ::

       cd Nyx/Exec/MiniSB

   This will put you into a directory in which you can run a small
   version of the Santa Barbara problem.

#. In Nyx/Exec/MiniSB, edit the GNUmakefile, and set

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
