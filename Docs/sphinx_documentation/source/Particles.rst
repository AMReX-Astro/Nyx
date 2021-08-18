*********************
Dark matter particles
*********************

Dark matter particles are included in the simulation by setting ::

    nyx.do_dm_particles = true

in the inputs file.

When dark matter particles are present, one has the option to run with or without the baryons; 
to build with no hydro capability set::

    NO_HYDRO=TRUE

It is also possible to build with ::

    NO_HYDRO=FALSE

but turn off hydro at run-time by setting ::

  nyx.do_hydro = 0

in the inputs file.

Our default is to use double precision for the mesh data and single precision for the particles; this is set in the 
GNUMakefile by::

  PRECISION = DOUBLE
  USE_SINGLE_PRECISION_PARTICLES = TRUE

We do not recommend using single precision for the mesh data, as this is not tested and will potentially degrade the quality of the hydrodynamical solve. To use single precision for the mesh use::

  PRECISION = FLOAT

Equations
=========

If we define :math:`{\mathbf x}_i` and :math:`{\bf u}_i` as the location and velocity of particle :math:`i`, respectively, then we wish
to solve

.. math::

   \begin{aligned}
   \frac{d {\mathbf x}_i}{d t} &=& \frac{1}{a} {\mathbf u}_i \\
   \frac{d (a {\mathbf u}_i) }{d t} &=& {\mathbf g}_i\end{aligned}

where :math:`{\mathbf g}_i` is the gravitational force evaluated at the location of particle :math:`i`, i.e.,
:math:`{\mathbf g}_i = {\mathbf g}({\mathbf x}_i,t).`


Particle Time Stepping: Move-Kick-Drift Algorithm
=================================================

In each time step:

-  Solve for :math:`{\mathbf g}^n` (only if multilevel, otherwise use :math:`{\mathbf g}^{n+1}` from previous step)

-  :math:`{\mathbf u}_i^{{n+\frac{1}{2}}} = \frac{1}{a^{{n+\frac{1}{2}}}} ( (a^n {\mathbf u}^n_i) + \frac{{\Delta t}}{2} \; {\mathbf g}^n_i )`

-  :math:`{\mathbf x}_i^{n+1 } = {\mathbf x}^n_i +  \frac{{\Delta t}}{a^{{n+\frac{1}{2}}}}  {\mathbf u}_i^{{n+\frac{1}{2}}}`

-  Solve for :math:`{\mathbf g}^{n+1}` using :math:`{\mathbf x}_i^{n+1}`

-  :math:`{\mathbf u}_i^{n+1} = \frac{1}{a^{n+1}} ( (a^{{n+\frac{1}{2}}} {\mathbf u}^{{n+\frac{1}{2}}}_i) + \frac{{\Delta t}}{2} \; {\mathbf g}^{n+1}_i )`

Note that at the end of the timestep :math:`{\bf x}_i^{n+1}` is consistent with :math:`{\bf g}^{n+1}` becasue
we have not advanced the positions after computing the new-time gravity. This has the benefit that
we perform only one gravity solve per timestep (in a single-level calculation with no hydro) because
the particles are only moved once.

Computing **g**
~~~~~~~~~~~~~~~

We solve for the gravitational vector as follows:

-  Assign the mass of the particles onto the grid in the form of density, :math:`\rho_{DM}`.
   The mass of each particle is assumed to be uniformly distributed over a cube of side :math:`\Delta x`,
   centered at what we call the position of the particle. We distribute the mass of each
   particle to the cells on the grid in proportion to the volume of the intersection of each cell
   with the particle’s cube. We then divide these cell values by :math:`\Delta x^3` so that the
   right hand side of the Poisson solve will be in units of density rather than mass.
   Note that this is the *comoving* density.

-  Solve :math:`\nabla^2 \phi = \frac{4 \pi G}{a} \rho_{DM}`.
   We discretize with the standard 7-point Laplacian (5-point in 2D)
   and use multigrid with Gauss-Seidel red-black relaxation to solve the equation for :math:`\phi` at cell centers.

-  Compute the normal component of :math:`{\bf g} = -\nabla \phi` at cell faces by differencing the adjacent values of :math:`\phi,`
   e.g. if :math:`{\bf g} = (g_x, g_y, g_z),` then we define :math:`g_x` on cell faces with a normal in the x-direction by computing
   :math:`g_{x,i-{\frac{1}{2}},j,k} = -(\phi_{i,j,k} - \phi_{i-1,j,k}) / \Delta x.`

-  Interpolate each component of :math:`{\bf g}` from normal cell faces onto each particle position using
   linear interpolation in the normal direction.

Output Format
=============

Checkpoint Files
----------------

  The particle positions and velocities are stored in a binary file in each checkpoint directory.
  This format is designed for being read by the code at restart rather than for diagnostics.
  We note that the value of :math:`a` is also written in each checkpoint directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

Particle Data in Plot Files
----------------------------

The particle positions and velocities will be written in binary files in each plotfile directory.
Dark matter particles will be in DM, active galactic nuclei particles will be in AGN,
neutrino particles will be in NPC.

  In addition, we can also
  visualize the particle locations as represented on the grid. There are multiple “derived quantities”
  which represent the particles. Including particle variables in the derived variables will make them
  be written as plotfile fields on the grid, i.e.::
  
    amr.derive_plot_vars = particle_count particle_mass_density 

  in the inputs file will generate plotfiles with only two variables.
  **particle_count** represents the number of particles in a grid cell;
  **particle_mass_density** is the density on the grid resulting from the particles.

  The same naming convention follows for particle velocities on the grid: **particle_x_velocity**,
  **particle_y_velocity**, **particle_z_velocity**

  Derived variables with **particle_** represent quantities from the Dark Matter Particle Container.
  Similar variables from the AGN particle container, and the Neutrino Particle Container
  are named **agn_** and **neutrino_**. Note these are particle fields written to the grid,
  which are distinct from the **density** field in the plotfile, which is baryonic density on the grid.

  We note that the value of :math:`a` is also written in each plotfile directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

