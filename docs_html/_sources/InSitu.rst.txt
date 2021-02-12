
 .. role:: cpp(code)
    :language: c++

.. _InSitu:

In situ Analysis
================

Nyx supports in situ visualization using Ascent, and can leverage AMReX's Sensei visualization pipeline. Additionally, Nyx supports in situ halo finding using Reeber2. This is useful both as
in situ analysis tool, and for subgrid models, like AGN (under development).

These insitu calls are controlled by input flags, in this example analysis happens starting after step 100, for step 199, 299, and so on and so forth::

  insitu.int = 100
  insitu.start = 100

Additionally, we can request analysis at specific redshifts, and Nyx will adjust the time-stepping to reach those red-shifts exactly::

  nyx.analysis_z_values = 7.0 6.0 5.0 4.0 3.0 2.0

Ascent visualization
--------------------

The primary example of this functionality is in `Nyx/Exec/LyA`. To compile with Ascent, add in GNUmakefile.summit::

  USE_ASCENT_INSITU = TRUE

The ``ascent_actions.yaml`` file will determine while the code is running what the Ascent publish action does. The ``ascent_actions_slicefile.yaml`` file gives an example of saving slice data, rather than an image file, while the simulation is running.

To build Ascent with Conduit for other machines and configurations, see `Building Ascent <https://ascent.readthedocs.io/en/latest/BuildingAscent.html>`_. For further information about using Ascent with AMReX-based applications, see `AMReX Blueprint Tutorial <https://amrex-codes.github.io/amrex/tutorials_html/Blueprint_Tutorial.html>`_ and `WarpX Ascent InSitu Documentation <https://warpx.readthedocs.io/en/latest/visualization/ascent.html>`_. 

Sensei interface
----------------

See AMReX documentation: `SENSEI <https://amrex-codes.github.io/amrex/docs_html/Visualization.html#sensei>`_

Halo finding
------------

To find halos in situ while Nyx is running we use the Reeber package.
To compile with Reeber, add in GNUmakefile::

  REEBER = TRUE

and set the location of Boost library::

  BOOST_INLUDE_DIR := ${BOOST_ROOT}/include/boost

GNUmake will default to::

  DIY_INCLUDE_DIR ?= ../../../diy/include
  REEBER_EXT_HOME ?= ../../../Reeber2

Note that these codes are in separate repositories and are not included in Nyx distribution.
If you intend to use in situ halo finding, you should clone Reeber from its
`GitHub page <https://github.com/mrzv/reeber>`_ and follow the installation instructions provided there.

In the Nyx inputs file, one should specify the time step interval for halo finder, fields which will be
used (to use the total density, one should specify both (gas) ``density`` and ``particle_mass_density``),
and thresholds of the boundary and extrema values::

  # Halo Finder
  reeber.halo_int = 1
  reeber.negate = 1
  reeber.halo_density_vars = density particle_mass_density
  reeber.halo_extrema_threshold = 20
  reeber.halo_component_threshold = 10
  
  # Call reeber insitu analysis
  insitu.reeber_int = 100


.. _note:
  These instructions are based on Reeber hash 8a274d35a415f7b15d8308a30763f52c4eeb7c7b, diy hash 88eca5107935b2d50eb352d99a6b0ed109b9c31c, and Nyx hash 33006ce18b1f945053c05a7cade0f4aba63378b5
