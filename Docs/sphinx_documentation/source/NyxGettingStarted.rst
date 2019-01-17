***************
Getting Started
***************

Downloading the Code
====================

Nyx is built on top of the AMReX framework. In order to run
Nyx, you must download two separate git modules.

.. raw:: latex

   \vspace{.1in}

First, make sure that git is installed on your machine—we recommend version 1.7.x or higher.

.. raw:: latex

   \vspace{.1in}

#. Clone/fork the  repository:

   ::

       git clone https://github.com/AMReX-Codes/amrex

   You will want to periodically update AMReX by typing

   ::

       git pull

   in the ``amrex/`` directory.

   Note: when you check out AMReX (and Nyx), you will get the master
   branch. The ``master`` branch is guaranteed to be compatible
   with the ``master`` branch.
   Active development is done on the ``development`` branch
   in each repo, and merged into the master branch monthly.
   If you wish to use the Nyx ``development`` branch, then you
   should also switch to the ``development`` branch for AMReX.

#. Set the environment variable, ``AMREX_HOME``, on your
   machine to point to the path name where you have put .
   You can add this to your ``.bashrc`` as:

   ::

       export AMREX_HOME={\em /path/to/amrex/}

   or to your ``.cshrc`` as:

   ::

       setenv AMREX_HOME {\em /path/to/amrex/}

   where you replace ``/path/to/amrex/`` will the full path to the
   amrex/ directory.

#. Clone/fork the Nyx repository:

   ::

       git clone https://github.com/AMReX-Astro/Nyx

   As with AMReX development on Nyx is done in the
   ``development`` branch, so you should work there if you want
   the latest source.

Building the Code
=================

#. From the directory in which you checked out Nyx, type

   ::

       cd Nyx/Exec/LyA

   This will put you into a directory in which you can run a small
   version of the Santa Barbara test problem.

#. In Nyx/Exec/LyA, edit the GNUmakefile, and set

   COMP = your favorite compiler (e.g, gnu, Intel)

   DEBUG = FALSE

   We like COMP = gnu.

#. Now type “make”. The resulting executable will look something like
   “Nyx3d.Linux.gnu.ex”, which means this is a 3-d version of the code,
   made on a Linux machine, with COMP = gnu.

Running the Code
================

#. Type:

   ::

       Nyx3d.Linux.gnu.ex inputs

#. You will notice that running the code generates directories that look like
   plt00000, plt00020, etc.,
   and chk00000, chk00020, etc. These are “plotfiles” and
   “checkpoint” files. The plotfiles are used for visualization,
   and the checkpoint files for restarting the code.

Visualization
=============

There are several visualization tools that can be used for plotfiles. The standard tool used within the
AMReX community is Amrvis, a package developed and supported
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data.
Plotfiles can also be viewed using the VisIt, ParaView, and yt packages. Particle data can be viewed using ParaView.

Please see the Visualization chapter in the AMReX User’s Guide (available in amrex/Docs)
for more detail about using all of these visualization packages. The webpage details are available at :: https://amrex-codes.github.io/amrex/docs_html/Visualization.html
