.. _Chap:GettingStartedNew:

Getting Started
===============

The Nyx source code currently lives in a
`github repository <https://github.com/AMReX-Astro/Nyx.git>`_ that
can be cloned by using git:

.. code-block:: shell

   > git clone https://github.com/AMReX-Astro/Nyx.git --recursive

.. note::

   We recommend Git version 1.7.x or higher.

Once you have obtained the source code, the following sections describe the
source code contents, compiling, running a simple simulation, and visualizing
the simulations results.

.. toctree::
   :maxdepth: 1

   Source directory overview <getting_started/Structure>
   Building Nyx with GNU Make <getting_started/BuildingGMake>
   Building Nyx with CMake <getting_started/BuildingCMake>
   Running the code <getting_started/RunningTheCode>
   Compiling Nyx with SUNDIALS 5 <getting_started/NyxSundials>

Additional details about compiling Nyx in different environments can be found
in the `AMReX build documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html>`_ .
Containerized Ubuntu builds with CMake are run as 
`GitHub Actions <https://github.com/AMReX-Astro/Nyx/actions/workflows/linux.yml>`_
whose workflow files are stored in Nyx/.github/workflows .
