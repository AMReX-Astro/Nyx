
 .. role:: cpp(code)
    :language: c++

.. _InSitu:

In situ analysis
================

Currently, Nyx supports in situ halo finding. This is useful both as
in situ analysis tool, and for subgrid models, like AGN (under development).


Halo finding
------

To find halos in situ while Nyx is running we use the Reeber package.
To compile with Reeber, add in GNUmakefile::

REEBER = TRUE

and set the location of Boost library::

BOOST_INLUDE_DIR := ${BOOST_ROOT}/include/boost

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
