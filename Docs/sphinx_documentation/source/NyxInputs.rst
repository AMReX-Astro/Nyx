.. role:: cpp(code)
  :language: c++

******
Inputs
******
.. toctree::
   :maxdepth: 1

The Nyx executable reads run-time information from an “inputs” file which you put on the command line. 
This section describes the inputs which can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

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

+-----------------+---------------------------+--------------+---------+
| Parameter       | Definition                | Acceptable   | Default |
|                 |                           | Values       |         |
+=================+===========================+==============+=========+
| **max_step**    | maximum number of level 0 | Integer >= 0 | -1      |
|                 | time steps                |              |         |
+-----------------+---------------------------+--------------+---------+
| **stop_time**   | final simulation          | Real >= 0    | -1.0    |
|                 | time                      |              |         |
+-----------------+---------------------------+--------------+---------+
| **nyx.final_a** | final value of a          | Real > 0     | -1.0    |
+-----------------+---------------------------+--------------+---------+
| **nyx.final_z** | final value of z          | Real > 0     | -1.0    |
+-----------------+---------------------------+--------------+---------+

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
|                     | particles      | <= 1           |                |
|                     |                |                |                |
|                     |                |                |                |
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
| **nyx.dt_binpow**   | time step      | Real >=  0     | -1.0           |
|                     | chosen to be   |                |                |
|                     | a power of a   |                |                |
|                     | half times the |                |                |
|                     | comoving time  |                |                |
|                     | step           |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.relative_max**| max da/dt      | Real > 0       | 0.01           |
| **_change_a**       |                |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.absolute_max**| a_new-a_old    | Real > 0       | -1.0           |
| **_change_a**       |                |                |                |
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

-  | **nyx.dt_binpow** = 1.0
   | sets :math:`\mathrm{dt}=\left(\frac{1}{2}\right)^{n}\mathrm{dt}_{\mathrm{a}}|n:\mathrm{dt}_{\mathrm{cfl}}>\left(\frac{1}{2}\right)^{n}\mathrm{dt_{a}}`
     where :math:`\mathrm{dt}_{\mathrm{cfl}}` is determined by the more
     restrictive timestep of **nyx.cfl** and **particles.cfl**, and
     where :math:`\mathrm{dt}_{\mathrm{a}}` is determined by the
     **relative_max_change_a**, **absolute_max_change_a**, and the
     evolution of :math:`\frac{da}{dt}`

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
| **amr.subcycli | Number of      | 1 or           | must be set in |
| g_iterations** | cycles at each | ``ref_ratio``  | Manual mode    |
|                | level          |                |                |
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

+---------------------------------+----------------+----------------+----------------+
| Parameter                       | Definition     | Acceptable     | Default        |
|                                 |                | Values         |                |
+=================================+================+================+================+
| **amr.check_file**              | prefix for     | String         | “*chk*”        |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_int**               | how often (by  | Integer        | -1             |
|                                 | level-0 time   | :math:`> 0`    |                |
|                                 | steps) to      |                |                |
|                                 | write restart  |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_per**               | how often (by  | Real           | -1.0           |
|                                 | simulation     | :math:`> 0`    |                |
|                                 | time) to write |                |                |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.restart**                 | name of the    | String         | not used if    |
|                                 | file           |                | not set        |
|                                 | (directory)    |                |                |
|                                 | from which to  |                |                |
|                                 | restart        |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_files_output** | should we      | 0 or 1         | 1              |
|                                 | write          |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_nfiles**            | how parallel   | Integer        | 64             |
|                                 | is the writing | :math:`\geq 1` |                |
|                                 | of the         |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_on_restart**   | should we      | 0 or 1         | 0              |
|                                 | write a        |                |                |
|                                 | checkpoint     |                |                |
|                                 | immediately    |                |                |
|                                 | after          |                |                |
|                                 | restarting     |                |                |
+---------------------------------+----------------+----------------+----------------+

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

+-----------------------------+------------------+------------------+---------+
| Parameter                   | Definition       | Acceptable       | Default |
|                             |                  | Values           |         |
+=============================+==================+==================+=========+
| **amr.plot_file**           | prefix for       | String           | “*plt*” |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_int**            | how often (by    | Integer          | -1      |
|                             | level-0 time     | :math:`> 0`      |         |
|                             | steps) to write  |                  |         |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_per**            | how often (by    | Real :math:`> 0` | -1.0    |
|                             | simulation time) |                  |         |
|                             | to write plot    |                  |         |
|                             | files            |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_vars**           | name of state    | ALL, NONE or     | ALL     |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.derive_plot_vars**    | name of derived  | ALL, NONE or     | NONE    |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_files_output**   | should we write  | 0 or 1           | 1       |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plotfile_on_restart** | should we write  | 0 or 1           | 0       |
|                             | a plotfile       |                  |         |
|                             | immediately      |                  |         |
|                             | after restarting |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_nfiles**         | how parallel is  | Integer          | 64      |
|                             | the writing of   | :math:`\geq 1`   |         |
|                             | the plotfiles    |                  |         |
+-----------------------------+------------------+------------------+---------+
| **nyx.plot_rank**           | should we plot   | True / False     | False   |
|                             | the processor ID |                  |         |
|                             | in the plotfiles |                  |         |
+-----------------------------+------------------+------------------+---------+
| **fab.format**              | Should we write  | NATIVE or IEEE32 | NATIVE  |
|                             | the plotfile in  |                  |         |
|                             | double or single |                  |         |
|                             | precision        |                  |         |
+-----------------------------+------------------+------------------+---------+

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

Plotfile Variables
------------------

Native variables
^^^^^^^^^^^^^^^^

These variables come directly from the ``StateData``, either the
``State_Type`` (for the hydrodynamic variables), ``DiagEOS_Type``
(for the nuclear energy generation quantities). ``PhiGrav_Type`` and
``Gravity_Type`` (for the gravity quantities)

+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``density``                       | Baryonic mass density, :math:`\rho`               | M\ :math:`_\odot` / Mpc\ :math:`^3`  |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``xmom``                          | x-momentum, :math:`(\rho u)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``ymom``                          | y-momentum, :math:`(\rho v)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``zmom``                          | z-momentum, :math:`(\rho w)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_E``                         | Total energy density                              | :math:`{\rm erg~km^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_e``                         | Internal energy density                           | :math:`{\rm erg~km^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Temp``                          | Temperature                                       | :math:`{\rm K}`                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Ne``                            | Number density of electrons                       | dimensionless                        |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_X``                         | Mass density of species X (only valid for non-    | dimensionless                        |
| (where X is H or He, the species  | constant species)                                 |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiGrav``                       | Gravitational potential                           | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``grav_x``, ``grav_y``,           | Gravitational acceleration                        | :math:`{\rm km~s^{-2}}`              |
| ``grav_z``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+

Derived variables
^^^^^^^^^^^^^^^^^

+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| variable name                     | description                                       | derive routine              | units                                   |
+===================================+===================================================+=============================+=========================================+
| ``divu``                          | :math:`\nabla \cdot \ub`                          | ``derdivu``                 | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_e``                        | Specific internal energy computed from the        | ``dereint2``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | conserved :math:`(\rho e)` state variable as      |                             |                                         |
|                                   | :math:`e = (\rho e)/\rho`                         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_E``                        | Specific internal energy computed from the        | ``dereint1``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | total energy and momentum conserved state as      |                             |                                         |
|                                   | :math:`e=[(\rho E)-\frac{1}{2}(\rho \ub^2)]/\rho` |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``kineng``                        | Kinetic energy density,                           | ``derkineng``               | :math:`{\rm erg~km^{-3}}`               |
|                                   | :math:`K = \frac{1}{2} |(\rho \ub)|^2`            |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``logden``                        | :math:`\log_{10} \rho`                            | ``derlogden``               | M\ :math:`_\odot` / Mpc\ :math:`^3`     |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``MachNumber``                    | Fluid Mach number, :math:`|\ub|/c_s`              | ``dermachnumber``           | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``maggrav``                       | Gravitational acceleration magnitude              | ``dermaggrav``              | :math:`{\rm km~s^{-2}}`                 |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magmom``                        | Momentum density magnitude,                       | ``dermagmom``               | :math:`{\rm g~km^{-2}~s^{-1}}`          |
|                                   | :math:`|\rho \ub|`                                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvel``                        | Velocity magnitude, :math:`|\ub|`                 | ``dermagvel``               | :math:`\mathrm{km~s^{-1}}`              |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvort``                       | Vorticity magnitude, :math:`|\nabla\times\ub|`    | ``dermagvort``              | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``pressure``                      | Total pressure, including ions and electrons      | ``derpres``                 | :math:`{\rm dyn~km^{-2}}`               |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``soundspeed``                    | Sound speed                                       | ``dersoundspeed``           | :math:`\mathrm{km~s^{-1}}`              |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``H`` or ``He``                   | Mass fraction of species H or He                  | ``derspec``                 | --                                      |
|                                   |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``x_velocity``,                   | Fluid velocity,                                   | ``dervel``                  | :math:`\mathrm{km~s^{-1}}`              |
| ``y_velocity``,                   | :math:`\ub = (\rho \ub)/\rho`                     |                             |                                         |
| ``z_velocity``                    |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
	 
Screen Output
=============

.. _list-of-parameters-10:

List of Parameters
------------------

+----------------------------+------------------+----------------+----------------+
| Parameter                  | Definition       | Acceptable     | Default        |
|                            |                  | Values         |                |
+============================+==================+================+================+
| **amr.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | Amr.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **nyx.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | Nyx.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **gravity.v**              | verbosity of     | 0 or 1         | 0              |
|                            | Gravity.cpp      |                |                |
+----------------------------+------------------+----------------+----------------+
| **mg.v**                   | verbosity of     | 0,1,2,3,4      | 0              |
|                            | multigrid        |                |                |
|                            | solver (for      |                |                |
|                            | gravity)         |                |                |
+----------------------------+------------------+----------------+----------------+
| **particles.v**            | verbosity of     | 0,1,2,3,4      | 0              |
|                            | particle-related |                |                |
|                            | processes        |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.grid_log**           | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | the grids are    |                |                |
|                            | written          |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.run_log**            | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | certain output   |                |                |
|                            | is written       |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.run_log_terse**      | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | certain          |                |                |
|                            | (terser)         |                |                |
|                            | output is        |                |                |
|                            | written          |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.sum_interval**       | if               |                |                |
|                            | :math:`> 0,`     |                |                |
|                            | how often (in    |                |                |
|                            | level-0 time     |                |                |
|                            | steps)           |                |                |
+----------------------------+------------------+----------------+----------------+
|                            | to compute and   | Integer        | -1             |
|                            | print integral   |                |                |
|                            | quantities       |                |                |
+----------------------------+------------------+----------------+----------------+

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
| **nyx.do_grav**          | Include          | 0 if false     | must be set      |
|                          | gravity as a     | 1 if true      |                  |
|                          | forcing term     |                |                  |
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

.. _notes-6:

Notes
-----

-  To include gravity you must set **nyx.do_grav** = 1 in the inputs file

Physics
=======

.. _list-of-parameters-12:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **nyx.do_hydro**                 | Time-advance     | 0 if false, 1   | must be set |
|                                  | the fluid        | if true         |             |
|                                  | dynamical        |                 |             |
|                                  | equations        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.ppm_type**                 | Use PPM or       | 0 for PLM       | 1 (PPM)     |
|                                  | PLM for hydro    | 1 for PPM       |             |
|                                  | advance          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.enforce_min_density_type** | how to enforce   | "floor"         | "floor"     |
|                                  | rho greater than | "cons"          |             |
|                                  | small_dens       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.strang_split**             | Use strang       | 0 if false, 1   | 1           |
|                                  | splitting        | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.sdc_split**                | Use sdc          | 0 if false, 1   | 0           |
|                                  | splitting        | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.strang_grown_box**         | Use growntilebox | 0 if false, 1   | 1           |
|                                  | to avoid comms   | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.add_ext_src**              | Include          | 0 if false, 1   | 0           |
|                                  | additional       | if true         |             |
|                                  | user-specified   |                 |             |
|                                  | source term      |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.nghost_state**             | Set number of    | {1,2,3,4}       | 1           |
|                                  | ghost cells for  |                 |             |
|                                  | state variables  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.use_const_species**        | If 1 then read   | 0 or 1          | 0           |
|                                  | h_species and    |                 |             |
|                                  | he_species       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.h_species**                | Concentration    | 0 :math:`<` X   | 0           |
|                                  | of H             | :math:`<` 1     |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.he_species**               | Concentration    | 0 :math:`<` X   | 0           |
|                                  | of He            | :math:`<` 1     |             |
+----------------------------------+------------------+-----------------+-------------+


Cosmology
=========

.. _list-of-parameters-13:

List of Parameters
------------------

+----------------------------------+--------------------+-----------------+-------------+
| Parameter                        | Definition         | Acceptable      | Default     |
|                                  |                    | Values          |             |
+==================================+====================+=================+=============+
| **nyx.comoving_OmM**             | Relative (total)   |  0 :math:`<` X  | must be set |
|                                  | mass density       |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_OmB**             | Relative baryon    |  0 :math:`<` X  | must be set |
|                                  | density            |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_OmR**             | Relative           |  0 :math:`<` X  | must be set |
|                                  | radiation density  |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_h**               | Dimensionless      |  0 :math:`<` X  | must be set |
|                                  | Hubble parameter   |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.gamma**                    | Dimensionless      |  0 :math:`<` X  | :math:`5/3` |
|                                  | factor relating    |  :math:`<` 2    |             |
|                                  | :math:`p, \rho, e` |                 |             |
+----------------------------------+--------------------+-----------------+-------------+

Examples of Usage
-----------------

-  | **nyx.gamma** This changes :math:`\gamma` in the :math:`\gamma` law gas: :math:`p = (\gamma - 1) \rho e.`

Reionization models
===================

.. _list-of-parameters-14:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **uvb_rates_file**               | Name of the UVB  |  string         | must be set |
|                                  | (TREECOOL) file  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **uvb_density_A**                | Density-dependent|  real           | 1.0         |
|                                  | heating          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **uvb_density_B**                | Density dependent|  real           | 0.0         |
|                                  | heating          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_zHI_flash**       | Redshift of      |  0 :math:`<` X  | -1.0        |
|                                  | "flash" H reion. |  or -1 if off   |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_zHeII_flash**     | Redshift of      |  0 :math:`<` X  | -1.0        |
|                                  | "flash" He reion.|  of -1 if off   |             |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_reion**                 | Inhomogeneous    |  0 or 1         | 0           |
|                                  | reionization     |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_zhi_file**              | File with        |  string         | must be set |
|                                  | reionization map |                 | (if used)   |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_grid**                  | Size of the      |  integer        | must be set |
|                                  | reionization grid|                 | (if used)   |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_T_zHI**           | H reionization   |  real           | 2.0e4       |
|                                  | heat input       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_T_zHeII**         | He reionization  |  real           | 1.5e4       |
|                                  | heat input       |                 |             |
+----------------------------------+------------------+-----------------+-------------+



Multigrid Inputs
================

The following inputs can be set directly in the AMReX solver classes but we set them via the Nyx gravity routines.

These must be preceded by "gravity" in the inputs file:

+----------------------+---------------------------------------------------+-----------+--------------+
|                      | Description                                       | Type      | Default      |
+======================+===================================================+===========+==============+
| v                    |  Verbosity of Gravity class                       |  Int      |   0          |
+----------------------+---------------------------------------------------+-----------+--------------+
| ml_tol               |  Relative tolerance for multilevel solves         |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| sl_tol               |  Relative tolerance for single-level solves       |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| delta_tol            |  Relative tolerance for synchronization solves    |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| mlmg_agglomeration   |  Should we agglomerate deep in the V-cycle        |  Int      |   1          |
+----------------------+---------------------------------------------------+-----------+--------------+
| mlmg_consolidation   |  Should we consolidate deep in the V-cycle        |  Int      |   1          |
+----------------------+---------------------------------------------------+-----------+--------------+
| dirichlet_bcs        |  Should we use homogeneous Dirichlet bcs in the   |  Int      |   0          |
|                      |  gravity solves (used for testing only)           |           |              |
+----------------------+---------------------------------------------------+-----------+--------------+

These must be preceded by "mg" in the inputs file:

+----------------------+-----------------------------------------------------+-------------+--------------+
|                      | Description                                         |  Type       | Default      |
+======================+=====================================================+=============+==============+
| v                    |  Verbosity of multigrid solver                      |  Int        |   0          |
+----------------------+-----------------------------------------------------+-------------+--------------+
| bottom_solver        |  What is the bottom solver?                         |  String     |   "bicg"     |
|                      |  Options include "bicg", "smoother", "hypre", etc   |             |              |
+----------------------+-----------------------------------------------------+-------------+--------------+
| max_fmg_iter         |  Maximum number of F-cycles to do before            |  Int        |   0          |
|                      |  continuing with V-cycles in a multigrid solve      |             |              |   
+----------------------+-----------------------------------------------------+-------------+--------------+

There are a number of additional inputs that can be used to control the multigrid solver.  

See the `AMReX Multigrid documentation`_ for more details.

.. _AMReX Multigrid documentation: https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html

Memory Optimization
===================

.. _list-of-parameters-15:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **nyx.shrink_to_fit**            | Shrink Particle  | 0 if false, 1   |             |
|                                  | vector to save   | if true         |             |
|                                  | memory           |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.minimize_memory**          | Use less         | 0 if false, 1   |             |
|                                  | temporary scratch| if true         |             |
|                                  | memory in hydro  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_int**         | How often to     | Int < 0 if never| -1          |
|                                  | load-balance     | Int > 0         |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_start_z**     | Redshift to start| Real > 0        | 7.0         |
|                                  | load-balancing   |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_wgt_stategy** | Weight strategy  | {0, 1, 2}       | 0           |
|                                  | to load-balance  |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_wgt_nmax**    | Max ranks to     | 0 < Int < Ranks | -1          |
|                                  | load-balance     |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_stategy**     | Dmap strategy    | {KNAPSACK,      | SFC         |
|                                  | type for particle|  SFC,           |             |
|                                  | load-balancing   |  ROUNDROBIN     |             |
+----------------------------------+------------------+-----------------+-------------+
