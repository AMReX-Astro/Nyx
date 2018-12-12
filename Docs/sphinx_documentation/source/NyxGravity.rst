*******
Gravity
*******
Introduction
============

Integration Strategy
--------------------
Nyx uses subcycling to integrate levels at different timesteps.
The gravity algorithm needs to respect this. Self-gravity is computed
via multigrid. At coarse-fine interfaces, the stencil used in the
Laplacian understands the coarse-fine interface and is different than
the stencil used in the interior.

There are two types of solves that we discuss with AMR:

-  *composite solve* : This is a multilevel solve, starting at
   a coarse level (usually level 0) and solving for the potential on
   all levels up to the finest level.

-  *level solve* : This solves for the potential only on
   a particular level. Finer levels are ignored. At coarse-fine
   interfaces, the data from the coarse levels acts as Dirichlet
   boundary conditions for the current-level-solve.

Briefly:

-  At the beginning of a simulation, we do a multilevel composite
   solve (if ``gravity.no_composite`` = 0).

   We also do a multilevel composite solve after each regrid.

-  The old-time gravity on the coarse level is defined based on
   this composite solve, but we also do a level solve on the coarse
   level, and use it to determine the difference between the composite
   solve and the level solve, and store that in a MultiFab.

-  After the hydro advance on the coarse level, we do another level
   solve, and use the (level solve - compositive solve) as a lagged
   predictor of how much we need to add back to that level solve to get
   an effective new-time value for phi on the coarse level, and that’s
   what defines the phi used for the new-time gravity

-  Then we do the fine grid timestep(s), each using the same
   strategy

-  At an AMR synchronization step across levels, if we’re
   choosing to synchronize the gravitational field across levels
   (``gravity.no_sync`` = 0) we then do a solve starting from the coarse
   grid that adjusts for the mismatch between the fine-grid phi and
   the coarse-grid phi, as well as the mismatch between the fine-grid
   density fluxes and the coarse-grid density fluxes, and add the
   resulting sync solve phi to both the coarse and the fine level

   Thus, to within the gravity error tolerance, you get the same final
   result as if you had done a full composite solve at the end of the
   timestep (assuming ``gravity.no_sync`` = 0).

If you do ``gravity.no_composite`` = 1, then you never do a full
multilevel solve, and the gravity on any level is defined only by the
solve on that level. The only time this would be appropriate is if
the fine level(s) cover essentially all of the mass on the grid for
all time.

Controls
--------
In order to use gravity, we must always set::

   USE_GRAV=TRUE
   
in the GNUmakefile and ::
  
  nyx.do_grav= 1
  
in the inputs file.

Poisson Approximation
=====================

In Nyx we always compute gravity by solving a Poisson equation on the mesh hierarchy.
To define the gravitational vector we set

  .. math:: \mathbf{g}(\mathbf{x},t) = -\nabla \phi

  where

  .. math:: \mathbf{\Delta} \phi = \frac{4 \pi G}{a} (\rho - \overline{\rho}) \label{eq:Self Gravity}

where :math:`\overline{\rho}` is the average of :math:`\rho` over the entire domain if we assume triply periodic boundary conditions,
and :math:`a(t)` is the scale of the universe as a function of time.

Hydrodynamics Source Terms
==========================

There are several options to incorporate the effects of gravity into
the hydrodynamics system. The main parameter here is
``Nyx.grav_source_type``.

- ``Nyx.grav_source_type`` = 1 : we use a standard
  predictor-corrector formalism for updating the momentum and
  energy. Specifically, our first update is equal to :math:`\Delta t
  \times \mathbf{S}^n` , where :math:`\mathbf{S}^n` is the value of
  the source terms at the old-time (which is usually called time-level
  :math:`n`). At the end of the timestep, we do a corrector step where
  we subtract off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on
  :math:`\Delta t / 2 \times \mathbf{S}^{n+1}`, so that at the end of
  the timestep the source term is properly time centered.

.. ``Nyx.grav_source_type`` = 2 : we do something very similar
   to 1. The major difference is that when evaluating the energy source
   term at the new time (which is equal to :math:`\mathbf{u} \cdot
   \mathbf{S}^{n+1}_{\rho \mathbf{u}}`, where the latter is the
   momentum source term evaluated at the new time), we first update the
   momentum, rather than using the value of :math:`\mathbf{u}` before
   entering the gravity source terms. This permits a tighter coupling
   between the momentum and energy update and we have seen that it
   usually results in a more accurate evolution.

- ``Nyx.grav_source_type`` = 3 : we do the same momentum update as
  the previous, but for the energy update, we put all of the work
  into updating the kinetic energy alone. In particular, we explicitly
  ensure that :math:`(rho e)` maintains the same, and update
  :math:`(rho K)` with the work due to gravity, adding the new kinetic
  energy to the old internal energy to determine the final total gas
  energy. The physical motivation is that work should be done on the
  velocity, and should not directly update the temperature—only
  indirectly through things like shocks.

.. - ``castro.grav_source_type`` = 4 : the energy update is done in a
   “conservative” fashion. The previous methods all evaluate the value
   of the source term at the cell center, but this method evaluates the
   change in energy at cell edges, using the hydrodynamical mass
   fluxes, permitting total energy to be conserved (excluding possible
   losses at open domain boundaries). See
   :cite:`katzthesis` for some more details.
