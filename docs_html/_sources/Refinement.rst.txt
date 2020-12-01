===================
Refinement Criteria
===================

The dynamic creation and destruction of grid levels is a fundamental part of `Nyx`'s capabilities. The
process for this is described in some detail in the `AMReX` documentation, but we summarize the key points
here.

At regular intervals (set by the user), each Amr level that is not the finest allowed for the run
will invoke a "regrid" operation.  When invoked, a list of error tagging functions is traversed. For each,
a field specific to that function is derived from the state over the level, and passed through a kernel
that "set"'s or "clear"'s a flag on each cell.  The field and function for each error tagging quantity is
identified in the setup phase of the code where the state descriptors are defined (i.e., in `Nyx_setup.cpp`).
Each function in the list adds or removes to the list of cells tagged for refinement. This final list of tagged
cells is sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm to create rectangular grids
which will define a new finer level (or set of levels).  State data is filled over these new grids, copying where
possible, and interpolating from coarser level when no fine data is available.  Once this process is complete,
the existing Amr level(s) is removed, the new one is inserted into the hierarchy, and the time integration
continues.

The traditional `AMReX` approach to setting up and controlling the regrid process involves explicitly
creating ("hard coding") a number of functions directly into `Nyx`'s setup code. (Consult the source code
and `AMReX` documentation for precisely how this is done).  `Nyx` provides a limited capability to augment
the standard set of error functions that is based entirely on runtime data specified in the inputs (ParmParse)
data.  The following example portion of a ParmParse'd input file demonstrates the usage of this feature:

::

      amr.refinement_indicators = flame_tracer lo_temp gradT

      amr.flame_tracer.max_level = 3
      amr.flame_tracer.value_greater = 1.e-6
      amr.flame_tracer.field_name = Y(H)

      amr.lo_temp.max_level = 1
      amr.lo_temp.value_less = 450
      amr.lo_temp.field_name = temp

      amr.gradT.max_level = 2
      amr.gradT.adjacent_difference_greater = 20
      amr.gradT.field_name = temp
      amr.gradT.start_time = 0.001
      amr.gradT.end_name = 0.002

Here, we have added three new custom-named criteria -- ``flame_tracer``: cells with the mass fraction of H greater than 1 ppm;
``lo_temp``: cells with T less than 450K, and ``gradT``: cells having a temperature difference of 20K from that of their
immediate neighbor.  The first will trigger up to Amr level 3, the second only to level 1, and the third to level 2.
The third will be active only when the problem time is between 0.001 and 0.002 seconds.

Note that these additional user-created criteria operate in addition to those defined as defaults.  Also note that
these can be modified between restarts of the code.  By default, the new criteria will take effect at the next
scheduled regrid operation.  Alternatively, the user may restart with ``amr.regrid_on_restart = 1`` in order to
do a full (all-levels) regrid after reading the checkpoint data and before advancing any cells.


