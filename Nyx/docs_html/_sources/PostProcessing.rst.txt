
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _PostProcessing:

PostProcessing
==============

Nyx interfaces with two post-processing suites, Reeber and Gimlet.

Reeber
------

Reeber uses topological methods to construct merge trees of scalar fields.
These trees are effectively parameter-independent and contain a complete
description of the field topology. In the context of Nyx, the field of interest
is the dark matter density. Nyx then queries the merge tree with user-defined
runtime parameters in order to locate the volume-averaged center of dark matter
halos. The same tree can be queried with any number of such parameters to find
halos with different mass/density thresholds.

Gimlet
------

Gimlet computes a variety of quantities about the simulation, including optical
depths, Lyman-alpha fluxes, power spectra (both 1-D ``line-of-sight'' as well as
fully 3-D), and probability distribution functions. These suites are fully
MPI-parallel and can be run either ``in situ'' or ``in-transit,'' or with a
combination of both. A detailed description of their usage is provided in the
Nyx User Guide.

Usage
------

Nyx can post-process with Gimlet alone, with Reeber alone, or with both
simultaneously. 

To compile with Gimlet, 

GIMLET = TRUE

to the GNUMakefile

To compile with Reeber, add 

REEBER = TRUE

Note that these codes are in separate repositories and are not included with Nyx.

Nyx and AMReX provide the capability for the user to execute an arbitrary
post-processing workflow in situ. An in situ workflow is one in which all 
MPI processes evolving the simulation stop at specified time steps and perform 
the post-processing before continuing with
the simulation.

