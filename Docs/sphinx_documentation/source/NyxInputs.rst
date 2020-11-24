.. role:: cpp(code)
  :language: c++
******
Inputs
******
.. toctree::
   :maxdepth: 1

  The  executable reads run-time information from an “inputs” file
   which you put on the command line. 

Problem Geometry
================

List of Parameters
------------------

+--------------------------+-----------------+-----------------+-------------+
| Parameter                | Definition      | Acceptable      | Default     |
|                          |                 | Values          |             |
+==========================+=================+=================+=============+
| **geometry.prob_lo**     | physical        | Real            | must be set |
|                          | location of low |                 |             |
|                          | corner of the   |                 |             |
|                          | domain          |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.prob_hi**     | physical        | Real            | must be set |
|                          | location of     |                 |             |
|                          | high corner of  |                 |             |
|                          | the domain      |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.coord_sys**   | coordinate      | 0 = Cartesian,  | must be set |
|                          | system          | 1 = r-z, 2 =    |             |
|                          |                 | spherical       |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.is_periodic** | is the domain   | 0 if false, 1   | 0 0 0       |
|                          | periodic in     | if true         |             |
|                          | this direction  |                 |             |
+--------------------------+-----------------+-----------------+-------------+

Examples of Usage
-----------------

-  **geometry.prob_lo** = 0 0 0
   defines the low corner of the domain at (0,0,0) in physical space.

-  **geometry.prob_hi** = 1.e8 2.e8 2.e8
   defines the high corner of the domain at (1.e8,2.e8,2.e8) in
     physical space.

-  **geometry.coord_sys** = 0
   defines the coordinate system as Cartesian

-  **geometry.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================

.. _list-of-parameters-1:

List of Parameters
------------------

+---------------+---------------------------------+-------------------+-------------+
| Parameter     | Definition                      | Acceptable Values | Default     |
+===============+=================================+===================+=============+
| **nyx.lo_bc** | boundary type of each low face  | 0,1,2,3,4,5       | must be set |
+---------------+---------------------------------+-------------------+-------------+
| **nyx.hi_bc** | boundary type of each high face | 0,1,2,3,4,5       | must be set |
+---------------+---------------------------------+-------------------+-------------+

[Table:BC]

Notes
-----

Boundary types are:

======================= ================
0 – Interior / Periodic 3 – Symmetry      
1 – Inflow              4 – Slip Wall     
2 – Outflow             5 – No Slip Wall  
======================= ================

Note – **nyx.lo_bc** and **nyx.hi_bc** must be consistent with
**geometry.is_periodic** – if the domain is periodic in a particular
direction then the low and high bc’s must be set to 0 for that
direction.

.. _examples-of-usage-1:

Examples of Usage
-----------------

-  **nyx.lo_bc** = 1 4 0

-  **nyx.hi_bc** = 2 4 0

-  **geometry.is_periodic** = 0 0 1

would define a problem with inflow (1) in the low-x direction,
outflow(2) in the high-x direction, slip wall (4) on the low and high
y-faces, and periodic in the z-direction.

Resolution
==========

.. _list-of-parameters-2:

List of Parameters
------------------

+---------------------------+-----------------+-----------------+-------------+
| Parameter                 | Definition      | Acceptable      | Default     |
|                           |                 | Values          |             |
+===========================+=================+=================+=============+
| **amr.n_cell**            | number of cells | Integer > 0     | must be set |
|                           | in each         |                 |             |
|                           | direction at    |                 |             |
|                           | the coarsest    |                 |             |
|                           | level           |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.max_level**         | number of       | Integer >= 0    | must be set |
|                           | levels of       |                 |             |
|                           | refinement      |                 |             |
|                           | above the       |                 |             |
|                           | coarsest level  |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio**         | ratio of coarse | 2 or 4          | must be set |
|                           | to fine grid    |                 |             |
|                           | spacing between |                 |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_int**        | how often to    | Integer > 0     | must be set |
|                           | regrid          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_on_restart** | should we       | 0 or 1          | 0           |
|                           | regrid          |                 |             |
|                           | immediately     |                 |             |
|                           | after           |                 |             |
|                           | restarting      |                 |             |
+---------------------------+-----------------+-----------------+-------------+

[Table:ResInputs]

Note: if **amr.max_level** = 0 then you do not need to set
**amr.ref_ratio** or **amr.regrid_int**.

.. _examples-of-usage-2:

Examples of Usage
-----------------

-  **amr.n_cell** = 32 64 64

   would define the domain to have 32 cells in the x-direction, 64 cells
   in the y-direction, and 64 cells in the z-direction *at the coarsest
   level*. (If this line appears in a 2D inputs file then the final
   number will be ignored.)

-  | **amr.max_level** = 2
   | would allow a maximum of 2 refined levels in addition to the coarse
     level. Note that these additional levels will only be created only
     if the tagging criteria are such that cells are flagged as needing
     refinement. The number of refined levels in a calculation must be
     :math:`\leq` **amr.max_level**, but can change in time and need not
     always be equal to **amr.max_level**.

-  | **amr.ref_ratio** = 2 4
   | would set factor 2 refinement between levels 0 and 1, and factor 4
     refinement between levels 1 and 2. Note that you must have at least
     **amr.>ax_level** values of **amr.ref_ratio** (Additional values
     may appear in that line and they will be ignored).

-  | **amr.regrid_int** = 2 2
   | tells the code to regrid every 2 steps. Thus in this example, new
     level-1 grids will be created every 2 level-0 time steps, and new
     level-2 grids will be created every 2 level-1 time steps.

Tagging
=======

.. _list-of-parameters-3:

List of Parameters
------------------

+-------------------------+------------------+------------------+---------+
| Parameter               | Definition       | Acceptable       | Default |
|                         |                  | Values           |         |
+=========================+==================+==================+=========+
| **nyx.allow_untagging** | are cells        | 0 or 1           | 0       |
|                         | allowed to be    |                  |         |
|                         | “untagged”       |                  |         |
+-------------------------+------------------+------------------+---------+

[Table:Tagging]

.. _notes-1:

Notes
-----

-  Typically cells at a given level can be tagged as needing refinement
   by any of a number of criteria, but cannot be “untagged”. That is,
   once tagged, no other criteria can untag them. If we set
   **nyx.allow_untagging** = 1 then the user is allowed to “untag” cells
   in the Fortran tagging routines.

Regridding
==========

Overview
--------

The details of the regridding strategy are described in Section
`5.5 <#subsec:grid-generation>`__, but first we cover how the input
parameters can control the gridding.

As described later, the user defines Fortran subroutines which tag
individual cells at a given level if they need refinement. This list of
tagged cells is sent to a grid generation routine, which uses the
Berger–Rigoutsos algorithm to create rectangular grids that contain the
tagged cells.

.. _list-of-parameters-4:

List of Parameters
------------------

+----------------------------+----------------+----------------+----------------+
| Parameter                  | Definition     | Acceptable     | Default        |
|                            |                | Values         |                |
+============================+================+================+================+
| **amr.regrid_file**        | name of file   | text           | no file        |
|                            | from which to  |                |                |
|                            | read the grids |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.grid_eff**           | grid           | Real > 0, < 1  | 0.7            |
|                            | efficiency at  |                |                |
|                            | coarse level   |                |                |
|                            | at which grids |                |                |
|                            | are created    |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.n_error_buf**        | radius of      | Integer >= 0   | 1              |
|                            | additional     |                |                |
|                            | tagging around |                |                |
|                            | already tagged |                |                |
|                            | cells          |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer > 0    | 128 in 2D, 32  |
|                            | of a grid in   |                | in 3D          |
|                            | any direction  |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer        | 128 in 2D, 32  |
+----------------------------+----------------+----------------+----------------+
| **amr.blocking_factor**    | grid size must | Integer > 0    | 2              |
|                            | be a multiple  |                |                |
|                            | of this        |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.refine_grid_layout** | refine grids   | 0 if false, 1  | 1              |
|                            | more if # of   | if true        |                |
|                            | processors     |                |                |
|                            | :math:`>` # of |                |                |
|                            | grids          |                |                |
+----------------------------+----------------+----------------+----------------+

[Table:GriddingInputs]

.. _notes-2:

Notes
-----

-  **amr.n_error_buf**, **amr.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **amr.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

.. _examples-of-usage-3:

Examples of Usage
-----------------

-  | **amr.regrid_file** = *fixed_grids*
   | In this case the list of grids at each fine level are contained in
     the file *fixed_grids*, which will be read during the gridding
     procedure. These grids must not violate the **amr.max_grid_size**
     criterion. The rest of the gridding procedure described below will
     not occur if **amr.regrid_file** is set.

-  | **amr.grid_eff** = 0.9
   | During the grid creation process, at least 90% of the cells in each
     grid at the level at which the grid creation occurs must be tagged
     cells. Note that this is applied at the coarsened level at which
     the grids are actually made, and before **amr.max_grid_size** is
     imposed.

-  | **amr.max_grid_size** = 64
   | The final grids will be no longer than 64 cells on a side at every
     level.

-  | **amr.max_grid_size** = 64 32 16
   | The final grids will be no longer than 64 cells on a side at level
     0, 32 cells on a side at level 1, and 16 cells on a side at level
     2.

-  | **amr.blocking_factor** = 32
   | The dimensions of all the final grids will be multiples of 32 at
     all levels.

-  | **amr.blocking_factor** = 32 16 8
   | The dimensions of all the final grids will be multiples of 32 at
     level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

   Having grids that are large enough to coarsen multiple levels in a
   V-cycle is essential for good multigrid performance in simulations
   that use self-gravity.

.. _subsec:grid-generation:

How Grids are Created
---------------------

The gridding algorithm proceeds in this order:

#. If at level 0, the domain is initially defined by **n_cell** as
   specified in the inputs file. If at level greater than 0, grids are
   created using the Berger–Rigoutsis clustering algorithm applied to
   the tagged cells, modified to ensure that the lengths of all new fine
   grids are divisible by **blocking_factor**.

#. Next, the grid list is chopped up if any grids have length longer
   than **max_grid_size**. Note that because **max_grid_size** is a
   multiple of **blocking_factor** (as long as **max_grid_size** is
   greater than **blocking_factor**), the **blocking_factor** criterion
   is still satisfied.

#. Next, if **refine_grid_layout** = 1 and there are more processors
   than grids at this level, then the grids at this level are further
   divided in order to ensure that no processor has fewer than one grid
   (at each level).

   -  if **max_grid_size** / 2 in the **BL_SPACEDIM** direction is a
      multiple of **blocking_factor**, then chop the grids in the
      **BL_SPACEDIM** direction so that none of the grids are longer in
      that direction than **max_grid_size / 2**

   -  If there are still fewer grids than processes, repeat the
      procedure in the **BL_SPACEDIM-1** direction, and again in the
      **BL_SPACEDIM-2** direction if necessary

   -  If after completing a sweep in all coordinate directions with
      **max_grid_size / 2**, there are still fewer grids than processes,
      repeat the steps above with **max_grid_size / 4**.

Simulation Time
===============

.. _list-of-parameters-5:

List of Parameters
------------------

+-----------------+--------------------------+--------------+---------+
| Parameter       | Definition               | Acceptable   | Default |
|                 |                          | Values       |         |
+=================+==========================+==============+=========+
| **max_step**    | maximum number           | Integer >= 0 | -1      |
|                 | of level-0 time          |              |         |
|                 | steps                    |              |         |
+-----------------+--------------------------+--------------+---------+
| **stop_time**   | final simulation         | Real >= 0    | -1.0    |
|                 | time                     |              |         |
+-----------------+--------------------------+--------------+---------+
| **nyx.final_a** | if                       | Real > 0     | -1.0    |
|                 | **nyx.use_comoving = t** |              |         |
|                 | and positive value       |              |         |
|                 | then this is             |              |         |
|                 | final value of a         |              |         |
+-----------------+--------------------------+--------------+---------+
| **nyx.final_z** | if                       | Real > 0     | -1.0    |
|                 | **nyx.use_comoving = t** |              |         |
|                 | and positive value       |              |         |
|                 | then this is             |              |         |
|                 | final value of z         |              |         |
+-----------------+--------------------------+--------------+---------+

[Table:TimeInputs]

.. _notes-3:

Notes
-----

To control the number of time steps, you can limit by the maximum number
of level-0 time steps (**max_step**), or the final simulation time
(**stop_time**), or both. The code will stop at whichever criterion
comes first. Note that if the code reaches **stop_time** then the final
time step will be shortened so as to end exactly at **stop_time**, not
pass it.

If running in comoving coordinates you can also set a final value of
:math:`a` by setting **nyx.final_a**, or a final value of :math:`z` by
setting **nyx.final_z**. You may only specify one or the other of these.
Once this value of :math:`a` or :math:`z` is reached in a time step, the
code will stop at the end of this full coarse time step. (Note it does
not stop at :math:`a` (or :math:`z`) exactly equal to the final value,
rather it stops once :math:`a` is greater than (or :math:`z` is less
than) this value.)

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **max_step** = 1000

-  **stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

-  If **nyx.do_hydro**\ :math:`= 1`, then typically the code chooses a
   time step based on the CFL number (dt = cfl \* dx / max(u+c) ).

-  If **nyx.do_hydro**\ :math:`= 0` and we are running with dark matter
   particles, then we use a time step based on the velocity of the
   particles, i.e., we set :math:`\Delta t` so that the particle goes no
   further than :math:`f \Delta t` in a coordinate direction where
   :math:`0 \leq f \leq 1.` The value for :math:`f` is currently
   hard-wired in Particles.H, but it will become an inputs parameter.

.. _list-of-parameters-6:

List of Parameters
------------------

+---------------------+----------------+----------------+----------------+
| Parameter           | Definition     | Acceptable     | Default        |
|                     |                | Values         |                |
+=====================+================+================+================+
| **nyx.cfl**         | CFL number for | Real > 0 and   | 0.8            |
|                     | hydro          | <= 1           |                |
|                     |                |                |                |
|                     |                |                |                |
+---------------------+----------------+----------------+----------------+
| **particles.cfl**   | CFL number for | Real > 0 and   | 0.5            |
|                     | particles      | <= 1           |
|                     |                |                |                |
|                     |                |              ` |                |
+---------------------+----------------+----------------+----------------+
| **nyx.init_shrink** | factor by      | Real > 0 and   | 1.0            |
|                     | which to       | <= 1           |                |
|                     | shrink the     |                |                |
|                     | initial time   |                |                |
|                     | step           |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.change_max**  | factor by      | Real >= 1      | 1.1            |
|                     | which the time |                |                |
|                     | step can grow  |                |                |
|                     | in subsequent  |                |                |
|                     | steps          |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.fixed_dt**    | level-0 time   | Real > 0       | unused if not  |
|                     | step           |                | set            |
|                     | regardless of  |                |                |
|                     | cfl or other   |                |                |
|                     | settings       |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.initial_dt**  | initial        | Real > 0       | unused if not  |
|                     | level-0 time   |                | set            |
|                     | step           |                |                |
|                     | regardless of  |                |                |
|                     | other settings |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.dt_cutoff**   | time step      | Real > 0       | 0.0            |
|                     | below which    |                |                |
|                     | calculation    |                |                |
|                     | will abort     |                |                |
+---------------------+----------------+----------------+----------------+

[Table:TimeStepInputs]

.. _examples-of-usage-5:

Examples of Usage
-----------------

-  | **nyx.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / umax_hydro.

-  | **particles.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / umax_particles where
     umax_particles is the maximum velocity of any particle in the
     domain.

-  | **nyx.init_shrink** = 0.01
   | sets the initial time step to 1% of what it would be otherwise.

-  | **nyx.change_max** = 1.1
   | allows the time step to increase by no more than 10% in this case.
     Note that the time step can shrink by any factor; this only
     controls the extent to which it can grow.

-  | **nyx.fixed_dt** = 1.e-4
   | sets the level-0 time step to be 1.e-4 for the entire simulation,
     ignoring the other timestep controls. Note that if
     **nyx.init_shrink** :math:`\neq 1` then the first time step will in
     fact be **nyx.init_shrink** \* **nyx.fixed_dt**.

-  | **nyx.initial_dt** = 1.e-4
   | sets the *initial* level-0 time step to be 1.e-4 regardless of
     **nyx.cfl** or **nyx.fixed_dt**. The time step can grow in
     subsequent steps by a factor of **nyx.change_max** each step.

-  | **nyx.dt_cutoff** = 1.e-20
   | tells the code to abort if the time step ever gets below 1.e-20.
     This is a safety mechanism so that if things go nuts you don’t burn
     through your entire computer allocation because you don’t realize
     the code is misbehaving.

Subcycling
==========

 supports a number of different modes for subcycling in time.

-  If **amr.subcycling_mode**\ :math:`=`\ Auto (default), then the code
   will run with equal refinement in space and time. In other words, if
   level :math:`n+1` is a factor of 2 refinement above level :math:`n`,
   then :math:`n+1` will take 2 steps of half the duration for every
   level :math:`n` step.

-  If **amr.subcycling_mode**\ :math:`=`\ None, then the code will not
   refine in time. All levels will advance together with a timestep
   dictated by the level with the strictest :math:`dt`. Note that this
   is identical to the deprecated command **amr.nosub = 1**.

-  If **amr.subcycling_mode**\ :math:`=`\ Manual, then the code will
   subcycle according to the values supplied by
   **subcycling_iterations**.

.. _list-of-parameters-7:

List of Parameters
------------------

+----------------+----------------+----------------+----------------+
| Parameter      | Definition     | Acceptable     | Default        |
|                |                | Values         |                |
+================+================+================+================+
| **amr.sub      | How shall we   | Auto, None or  | Auto           |
| cycling_mode** | subcycle       | Manual         |                |
+----------------+----------------+----------------+----------------+
| *              | Number of      | 1 or           | must be set in |
| *amr.subcyclin | cycles at each | ``ref_ratio``  | Manual mode    |
| g_iterations** | level          |                |                |
+----------------+----------------+----------------+----------------+

.. _examples-of-usage-6:

Examples of Usage
-----------------

-  | **amr.subcycling_mode**\ :math:`=`\ Manual
   | Subcycle in manual mode with largest allowable timestep.

-  | **amr.subcycling_iterations** = 1 2 1 2
   | Take 1 level-0 timestep at a time (required). Take 2 level-1
     timesteps for each level-0 step, 1 timestep at level 2 for each
     level-1 step, and take 2 timesteps at level 3 for each level 2
     step.

-  | **amr.subcycling_iterations** = 2
   | Alternative form. Subcycle twice at every level (except level 0).

Restart Capability
==================

|  has a standard sort of checkpointing and restarting capability. In
  the inputs file, the following options control the generation of
  checkpoint files (which are really directories):

.. _list-of-parameters-8:

List of Parameters
------------------

+----------------+----------------+----------------+----------------+
| Parameter      | Definition     | Acceptable     | Default        |
|                |                | Values         |                |
+================+================+================+================+
| **am           | prefix for     | Text           | “*chk*”        |
| r.check_file** | restart files  |                |                |
+----------------+----------------+----------------+----------------+
| **a            | how often (by  | Integer        | -1             |
| mr.check_int** | level-0 time   | :math:`> 0`    |                |
|                | steps) to      |                |                |
|                | write restart  |                |                |
|                | files          |                |                |
+----------------+----------------+----------------+----------------+
| **a            | how often (by  | Real           | -1.0           |
| mr.check_per** | simulation     | :math:`> 0`    |                |
|                | time) to write |                |                |
|                | restart files  |                |                |
+----------------+----------------+----------------+----------------+
| *              | name of the    | Text           | not used if    |
| *amr.restart** | file           |                | not set        |
|                | (directory)    |                |                |
|                | from which to  |                |                |
|                | restart        |                |                |
+----------------+----------------+----------------+----------------+
| **a            | should we      | 0 or 1         | 1              |
| mr.checkpoint_ | write          |                |                |
| files_output** | checkpoint     |                |                |
|                | files          |                |                |
+----------------+----------------+----------------+----------------+
| **amr.         | how parallel   | Integer        | 64             |
| check_nfiles** | is the writing | :math:`\geq 1` |                |
|                | of the         |                |                |
|                | checkpoint     |                |                |
|                | files          |                |                |
+----------------+----------------+----------------+----------------+
| *              | should we      | 0 or 1         | 0              |
| *amr.checkpoin | write a        |                |                |
| t_on_restart** | checkpoint     |                |                |
|                | immediately    |                |                |
|                | after          |                |                |
|                | restarting     |                |                |
+----------------+----------------+----------------+----------------+

.. _notes-4:

Notes
-----

-  You should specify either **amr.check_int** or **amr.check_per**. Do
   not try to specify both.

-  Note that if **amr.check_per** is used then in order to hit that
   exact time the code may modify the time step slightly, which will
   change your results ever so slightly than if you didn’t set this
   flag.

-  Note that **amr.plotfile_on_restart** and
   **amr.checkpoint_on_restart** only take effect if
   **amr.regrid_on_restart** is in effect.

-  See the Software Section for more details on parallel I/O and the
   **amr.check_nfiles** parameter.

-  If you are doing a scaling study then set
   **amr.checkpoint_files_output** = 0 so you can test scaling of the
   algorithm without I/O.

.. _examples-of-usage-7:

Examples of Usage
-----------------

-  **amr.check_file** = *chk_run*

-  **amr.check_int** = 10

   means that restart files (really directories) starting with the
   prefix “*chk_run*” will be generated every 10 level-0 time steps. The
   directory names will be *chk_run00000*, *chk_run00010*,
   *chk_run00020*, etc.

If instead you specify

-  **amr.check_file** = *chk_run*

-  **amr.check_per** = 0.5

   then restart files (really directories) starting with the prefix
   “*chk_run*” will be generated every 0.1 units of simulation time. The
   directory names will be *chk_run00000*, *chk_run00043*,
   *chk_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from *chk_run00061*,for example, then set

-  **amr.restart** = *chk_run00061*

.. _sec:PlotFiles:

Controlling PlotFile Generation
===============================

The main output from  is in the form of plotfiles (which are really
directories). The following options in the inputs file control the
generation of plotfiles

.. _list-of-parameters-9:

List of Parameters
------------------

+------------------+------------------+------------------+---------+
| Parameter        | Definition       | Acceptable       | Default |
|                  |                  | Values           |         |
+==================+==================+==================+=========+
| *                | prefix for       | Text             | “*plt*” |
| *amr.plot_file** | plotfiles        |                  |         |
+------------------+------------------+------------------+---------+
| **amr.plot_int** | how often (by    | Integer          | -1      |
|                  | level-0 time     | :math:`> 0`      |         |
|                  | steps) to write  |                  |         |
|                  | plot files       |                  |         |
+------------------+------------------+------------------+---------+
| **amr.plot_per** | how often (by    | Real :math:`> 0` | -1.0    |
|                  | simulation time) |                  |         |
|                  | to write plot    |                  |         |
|                  | files            |                  |         |
+------------------+------------------+------------------+---------+
| *                | name of state    | ALL, NONE or     | ALL     |
| *amr.plot_vars** | variables to     | list             |         |
|                  | include in       |                  |         |
|                  | plotfiles        |                  |         |
+------------------+------------------+------------------+---------+
| **amr.de         | name of derived  | ALL, NONE or     | NONE    |
| rive_plot_vars** | variables to     | list             |         |
|                  | include in       |                  |         |
|                  | plotfiles        |                  |         |
+------------------+------------------+------------------+---------+
| **amr.plo        | should we write  | 0 or 1           | 1       |
| t_files_output** | plot files       |                  |         |
+------------------+------------------+------------------+---------+
| **amr.plotf      | should we write  | 0 or 1           | 0       |
| ile_on_restart** | a plotfile       |                  |         |
|                  | immediately      |                  |         |
|                  | after restarting |                  |         |
+------------------+------------------+------------------+---------+
| **a              | how parallel is  | Integer          | 64      |
| mr.plot_nfiles** | the writing of   | :math:`\geq 1`   |         |
|                  | the plotfiles    |                  |         |
+------------------+------------------+------------------+---------+
| **ny             | Should we plot   | 0 or 1           | 0       |
| x.plot_phiGrav** | the              |                  |         |
|                  | gravitational    |                  |         |
|                  | potential        |                  |         |
+------------------+------------------+------------------+---------+
|                  | plot the         | 0 or 1           | 0       |
|                  | gravitational    |                  |         |
|                  | potential        |                  |         |
+------------------+------------------+------------------+---------+
| **particles.wri  | Should we write  | 0 or 1           | 0       |
| te_in_plotfile** | the particles in |                  |         |
|                  | a file within    |                  |         |
|                  | the plotfile     |                  |         |
+------------------+------------------+------------------+---------+
| **fab.format**   | Should we write  | NATIVE or IEEE32 | NATIVE  |
|                  | the plotfile in  |                  |         |
|                  | double or single |                  |         |
|                  | precision        |                  |         |
+------------------+------------------+------------------+---------+

All the options for **amr.derive_plot_vars** are kept in ``derive_lst``
in ``Nyx_setup.cpp``. Feel free to look at it and see what’s there.

.. _notes-5:

Notes
-----

-  You should specify either **amr.plot_int** or **amr.plot_per**. Do
   not try to specify both.

-  Note that if **amr.plot_per** is used then in order to hit that exact
   time the code may modify the time step slightly, which will change
   your results ever so slightly than if you didn’t set this flag.

-  See the Software Section for more details on parallel I/O and the
   **amr.plot_nfiles** parameter.

-  If you are doing a scaling study then set **amr.plot_files_output** =
   0 so you can test scaling of the algorithm without I/O.

-  **nyx.plot_phiGrav** is only relevant if **nyx.do_grav** = 1 

-  By default, plotfiles are written in double precision (NATIVE
   format). If you want to save space by writing them in single
   precision, set the fab.format flag to IEEE32.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **amr.plot_file** = *plt_run*

-  **amr.plot_int** = 10

   means that plot files (really directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. The
   directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.

If instead you specify

-  **amr.plot_file** = *plt_run*

-  **amr.plot_per** = 0.5

   then restart files (really directories) starting with the prefix
   “plt_run” will be generated every 0.1 units of simulation time. The
   directory names will be *plt_run00000*, *plt_run00043*,
   *plt_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

Screen Output
=============

.. _list-of-parameters-10:

List of Parameters
------------------

+----------------+----------------+----------------+----------------+
| Parameter      | Definition     | Acceptable     | Default        |
|                |                | Values         |                |
+================+================+================+================+
| **amr.v**      | verbosity of   | 0 or 1         | 0              |
|                | Amr.cpp        |                |                |
+----------------+----------------+----------------+----------------+
| **nyx.v**      | verbosity of   | 0 or 1         | 0              |
|                | Nyx.cpp        |                |                |
+----------------+----------------+----------------+----------------+
| **gravity.v**  | verbosity of   | 0 or 1         | 0              |
|                | Gravity.cpp    |                |                |
+----------------+----------------+----------------+----------------+
| **mg.v**       | verbosity of   | 0,1,2,3,4      | 0              |
|                | multigrid      |                |                |
|                | solver (for    |                |                |
|                | gravity)       |                |                |
+----------------+----------------+----------------+----------------+
| *              | verbosity of   | 0,1,2,3,4      | 0              |
| *particles.v** | pa             |                |                |
|                | rticle-related |                |                |
|                | processes      |                |                |
+----------------+----------------+----------------+----------------+
| **             | name of the    | Text           | not used if    |
| amr.grid_log** | file to which  |                | not set        |
|                | the grids are  |                |                |
|                | written        |                |                |
+----------------+----------------+----------------+----------------+
| *              | name of the    | Text           | not used if    |
| *amr.run_log** | file to which  |                | not set        |
|                | certain output |                |                |
|                | is written     |                |                |
+----------------+----------------+----------------+----------------+
| **amr.r        | name of the    | Text           | not used if    |
| un_log_terse** | file to which  |                | not set        |
|                | certain        |                |                |
|                | (terser)       |                |                |
|                | output is      |                |                |
|                | written        |                |                |
+----------------+----------------+----------------+----------------+
| **amr.         | if             |                |                |
| sum_interval** | :math:`> 0,`   |                |                |
|                | how often (in  |                |                |
|                | level-0 time   |                |                |
|                | steps)         |                |                |
+----------------+----------------+----------------+----------------+
|                | to compute and | Integer        | -1             |
|                | print integral |                |                |
|                | quantities     |                |                |
+----------------+----------------+----------------+----------------+
| **nyx.do_spe   |                | 0 or 1         | 1              |
| cial_tagging** |                |                |                |
+----------------+----------------+----------------+----------------+

.. _notes-6:

Notes
-----

-  **nyx.do_special_tagging** = 1 allows the user to set a special flag
   based on user-specified criteria. This can be used, for example, to
   calculate the bounce time in a core collapse simulation; the bounce
   time is defined as the first time at which the maximum density in the
   domain exceeds a user-specified value. This time can then be printed
   into a special file as a useful diagnostic.

.. _examples-of-usage-9:

Examples of Usage
-----------------

-  | **amr.grid_log** = *grdlog*
   | Every time the code regrids it prints a list of grids at all
     relevant levels. Here the code will write these grids lists into
     the file *grdlog*.

-  | **amr.run_log** = *runlog*
   | Every time step the code prints certain statements to the screen
     (if **amr.v** = 1), such as
   | STEP = 1 TIME = 1.91717746 DT = 1.91717746
   | PLOTFILE: file = plt00001
   | Here these statements will be written into *runlog* as well.

-  | **amr.run_log_terse** = *runlogterse*
   | This file, *runlogterse*, differs from *runlog* in that it only
     contains lines of the form
   | 10 0.2 0.005
   | in which “10” is the number of steps taken, “0.2” is the simulation
     time, and “0.005” is the level-0 time step. This file can be
     plotted very easily to monitor the time step.

-  | **nyx.sum_interval** = 2
   | if **nyx.sum_interval** :math:`> 0` then the code computes and
     prints certain integral quantities, such as total mass, momentum
     and energy in the domain every **nyx.sum_interval** level-0 steps.
     In this example the code will print these quantities every two
     coarse time steps. The print statements have the form
   | TIME= 1.91717746 MASS= 1.792410279e+34
   | for example. If this line is commented out then it will not compute
     and print these quanitities.

Gravity
=======

.. _list-of-parameters-11:

List of Parameters
------------------

+--------------------------+------------------+----------------+------------------+
| Parameter                | Definition       | Acceptable     | Default          |
|                          |                  | Values         |                  |
+==========================+==================+================+==================+
| **nyx.do_grav**          | Include          | 0 if false     | must be set if   |
|                          | gravity as a     | 1 if true      | USE_GRAV = TRUE  |
|                          | forcing term     |                | TRUE             |
+--------------------------+------------------+----------------+------------------+
| **gravity.no_sync**      | whether to       | 0 if false     |  0               |
|                          | perform the      | 1 if true      |                  |
|                          | “sync solve”     |                |                  |
+--------------------------+------------------+----------------+------------------+
| **gravity.no_composite** | whether to       | 0 if false     |  0               |
|                          | perform a        | 1 if true      |                  |
|                          | composite        |                |                  |
|                          | solve            |                |                  |
+--------------------------+------------------+----------------+------------------+

[Table:Gravity]

.. _notes-7:

Notes
-----

-  To include gravity you must set

   -  USE_GRAV = TRUE in the GNUmakefile

   -  **nyx.do_grav** = 1 in the inputs file

-  **gravity.no_sync** and **gravity.no_composite** are only relevant if
   USE_GRAV = TRUE; they both default to 0.

Physics
=======

.. _list-of-parameters-12:

List of Parameters
------------------

+---------------------------+-----------------+-----------------+-------------+
| Parameter                 | Definition      | Acceptable      | Default     |
|                           |                 | Values          |             |
+===========================+=================+=================+=============+
| **nyx.do_hydro**          | Time-advance    | 0 if false, 1   | must be set |
|                           | the fluid       | if true         |             |
|                           | dynamical       |                 |             |
|                           | equations       |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **nyx.do_react**          | Include         | 0 if false, 1   | must be set |
|                           | reactions       | if true         |             |
+---------------------------+-----------------+-----------------+-------------+
| **nyx.add_ext_src**       | Include         | 0 if false, 1   | 0           |
|                           | additional      | if true         |             |
|                           | user-specified  |                 |             |
|                           | source term     |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **nyx.use_const_species** | If 1 then read  | 0 or 1          | 0           |
|                           | h_species and   |                 |             |
|                           | he_species      |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **nyx.h_species**         | Concentration   | 0 :math:`<` X   | 0           |
|                           | of H            | :math:`<` 1     |             |
+---------------------------+-----------------+-----------------+-------------+
| **nyx.he_species**        | Concentration   | 0 :math:`<` X   | 0           |
|                           | of He           | :math:`<` 1     |             |
+---------------------------+-----------------+-----------------+-------------+

[Table:Physics]
