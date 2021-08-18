
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _PostProcessing:

Analyzing Outputs
=================

For analysis and visualizations purposes, Nyx outputs plotfiles with user-defined grid quantities and/or
particle data.  This file format is native to the AMReX framework, but there is also an option to use HDF5
outputs (this is implemented and is in a testing/optimization phase now).


Visualization
-------------

There are several visualization tools that can be used for AMReX plotfiles,
specifically VisIt, ParaView, and yt.

See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html>`_ for
guidance on how to use each of these tools.

The following inputs control when plotfiles are written

+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                         | Description                                                           |   Type      | Default   |
+=========================+=======================================================================+=============+===========+
| amr.plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                         | if -1 then no plotfiles will be written                               |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+

and the following control which variables appear in the plotfile

+----------------------------+---------------------------------------------------+------------+-----------+
|                            | Description                                       |   Type     | Default   |
+============================+===================================================+============+===========+
| amr.plot_vars              | Name of state variables to be in plotfiles        |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+
| amr.plot_dervive_plot_vars | Name of derived variables to be in plotfiles      |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+


Nyx also easily interfaces with two post-processing suites, Reeber used for halo finding
and Gimlet, used for calculating different summary statistics.


Reeber
------

Reeber uses topological methods to construct merge trees of scalar fields.
These trees are effectively parameter-independent and contain a complete
description of the field topology. In the context of Nyx, the field of interest
is usually the total matter density. Nyx then queries the merge tree with user-defined
runtime parameters in order to locate the volume-averaged center of dark matter
halos. The same tree can be queried with any number of such parameters to find
halos with different mass/density thresholds.  Reeber is publicly available on
`GitHub <https://github.com/mrzv/reeber>`_.


Gimlet
------

Gimlet computes a variety of quantities about the simulation, including optical
depths, Lyman-alpha fluxes, power spectra (both 1-D "line-of-sight" as well as
fully 3-D), and probability distribution functions. These suites are fully
MPI-parallel and can be run either "in situ" or "in-transit", or with a
combination of both. Gimlet code is not yet publicly available, but we are working
on releasing it.


Do It Yourself
--------------

Nyx and AMReX provide the capability for the user to execute an arbitrary
post-processing workflow.  In ``Util/Diagnostics/`` we provide a simple example
of opening an AMReX plotfile, reading and manipulating data in it.  That can be a
starting place for building analysis tools for any specific need.
