*********************
Dark matter particles
*********************

| For the moment, assume that we are running in comoving coordinates,
  with dark matter particles only (no hydro) and that the particles all exist at level 0. These assumptions are
  encapsulated in the following lines in the inputs file:

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

Initializing the Particles
==========================

There are several different ways in which particles can currently be initialized:

One must include the dark matter particles in the ``GNUmakefile`` by setting::

   USE_PARTICLES = TRUE
   nyx.do_dm_particles=1

And the particles can be initialized via ascii file, binary file, or other means.

Read from an ASCII file
-----------------------

To enable this option, set::
  
  nyx.particle_init_type = AsciiFile
  nyx.ascii_particle_file = *particle_file*

Here *particle_file* is the user-specified name of the file. The first line in this file is
assumed to contain the number of particles. Each line after that contains
x y z mass xdot ydot zdot
Note that the variable that we call the particle velocity, :math:`{\mathbf u} = a {\bf \dot{x}}`,
so we must multiply :math:`{\bf \dot{x}}`, by :math:`a` when we initialize the particles.

Read from a binary file
-----------------------

To enable this option, set::

  nyx.particle_init_type = BinaryFile
  nyx.binary_particle_file = *particle_file*
  
| As with the ASCII read, the first line in *particle_file* is
  assumed to contain the number of particles. Each line after that contains
| x y z mass xdot ydot zdot
| Note that the variable that we call the particle velocity, :math:`{\mathbf u} = a {\bf \dot{x}}`,
  so we must multiply :math:`{\bf \dot{x}}`, by :math:`a` when we initialize the particles.

Read from a binary "meta" file
------------------------------

| This option allows you to read particles from a series of files rather than
  just a single file. To enable this option, set::

    nyx.particle_init_type = BinaryMetaFile
    nyx.binary_particle_file =*particle file*
    
| In this case the *particle_file* you specify is an ASCII file specifying a
  list of file names with full paths. Each of the files in this list is assumed
  to be binary and is read sequentially (individual files are read in parallel) in
  the order listed.

Reading SPH particles
---------------------

For some applications it is useful to initialize the grid data with SPH-type
particles. To enable this option, you must set::

    nyx.do_santa_barbara = 1
    nyx.init_with_sph_particles =1

The SPH-type particles can then be read in by setting
where *sph_particle_file* is the user-specified name of the file
containing the SPH particles. The type of *sph_particle_file*
must be the same (Ascii, Binary or BinaryMeta) as the dark matter particle
file as specified by
The SPH particles will be discarded by the code once the grid data has been initialized.

Random placement
----------------

To enable this option, set::

  nyx.nyx.particle_init_type = Random
  
There are then a number of parameters to set, for example::
  
  nyx.particle_initrandom_count = 100000
  nyx.particle_initrandom_mass = 1
  nyx.particle_initrandom_iseed = 15

Cosmological
------------

Using cosmological initial conditions is a three step process:

#. Generating a transfer function (e.g. with ``camb``)

#. Generating an initial displacement field (with ``nyx-ic``)

#. Starting nyx

In the following we will look at each step a bit closer.

Generating a transfer function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transfer function is used in ``nyx-ic`` to generate the power
spectrum. The usual way is to use ``camb``\  [1]_
to calculate it for the desired universe. A sample ``camb.ini`` is
provided with ``nyx-ic``. The important options are:

-  **transfer_redshift(1) = 50**

-  **transfer_matterpower(1) = tf**

which determine the initial time for the simulation. You should make sure
that you catch all necessary wave numbers for the considered box length and
resolution.

From the ``camb`` output you have to note values for ``sigma_8``
for a redshift of zero and the initial redshift. We need this to compute
the right normalization.

Setting up the initial displacements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| We calculate the initial displacements with a stand-alone program called
  ``nyx-ic``. This takes a transfer function and some cosmological parameters
  as an argument and outputs an "init" directory which basically contains initial
  displacements for every grid point in an AMReX MultiFAB. Furthermore the mf
  contains a fourth field containing the density contrast as initial condition
  for the baryonic matter.
| ``nyx-ic`` is started with an “inputs“
  file similar to the one from Nyx. A sample one is provided. The options are

::

    #Omega_{Matter}
    cosmo.omegam = 0.272
    #Omega_{Lambda}
    cosmo.omegax = 0.728

    #equation of state paramater omega_{effective}
    cosmo.weff = -0.980

    #Omega_{baryon}*Hubble^2 
    cosmo.ombh2 = 0.0226
    #Hubble/100km/s
    cosmo.hubble = 0.704
    #scalar spectral index
    cosmo.enn = 0.963
    # initial z
    cosmo.z_init = 50

    #sidelength of the box (in Mpc)
    cosmo.boxside = 90.14
    #seed of the rng
    cosmo.isd = 100
    #resolution of the box
    cosmo.gridpoints = 256
    #the output file name
    cosmo.initDirName = init

    #choose the source of the transferfunction
    cosmo.transferfunction = CAMB

    #some tabulated transferfunction generated with camb (compare camb-ini-file)
    cosmo.tabulatedTk = tf
    # sigma8 for the input tf at z=0 and initial z (to calc the growthfactor)
    cosmo.init_sigma8_0 = 0.7891368
    cosmo.init_sigma8_init = 2.0463364E-02

The code solves the equation

.. math::

   \begin{aligned}
       P(k,a) = 2\pi^2\delta^2_H \frac{k^n}{H_0^{n+3}}T^2(k)\left( \frac{D(a)}{D(a=1)} \right)^2
       \end{aligned}

to calculate :math:`P` and from that gaussian distributed density perturbations
:math:`\delta` following that spectrum. Particle displacements are then calculated
as Zel’dovich displacements.

Non-gaussian effects as well as neutrino contributions are planned for the
future.

Using Nyx with cosmological initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  | **nyx.nyx.particle_init_type = Cosmological**
   | set the *right* init type

-  | **cosmo.initDirName = init**
   | set the name of the displacements directory (amrex format)

-  | **cosmo.particle_mass = 0.19178304E+10**
   | sets the mass [:math:`M_\odot`] of each particle

-  | **cosmo.omegam = 0.272**
   | set :math:`\Omega_{Matter}`

-  | **cosmo.omegax = 0.728**
   | set :math:`\Omega_\Lambda`

-  | **cosmo.hubble = 0.704**
   | set the reduced hubble constant :math:`h`

We will generate a particle of mass **particle_mass** in every grid cell
displaced from the center by the value found in the **initDirName** for
that cell. Velocities are calculated in the Zel’dovich approximation by

.. math::

   \begin{aligned}
           \vec{v} = \Delta{\vec{x}} \times 100 \text{km/s} \times a \sqrt{\Omega_M/a^3+\Omega_\Lambda} \times L_{\text{box}}
       \end{aligned}

where :math:`\Delta{\vec{x}}` is the displacement of the particle.

Time Stepping
=============

There are currently two different ways in which particles can be moved:

Random
------

| To enable this option, set
| Update the particle positions at the end of each coarse time step using a
  random number between 0 and 1 multiplied by 0.25 dx.

Motion by Self-Gravity
----------------------

| To enable this option, set

Move-Kick-Drift Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~

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

| The particle positions and velocities are stored in a binary file in each checkpoint directory.
  This format is designed for being read by the code at restart rather than for diagnostics.
| We note that the value of :math:`a` is also written in each checkpoint directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

Plot Files
----------

If **particles.write_in_plotfile =** 1 in the inputs file
then the particle positions and velocities will be written in a binary file in each plotfile directory.

| In addition, we can also
  visualize the particle locations as represented on the grid. There are two “derived quantities”
  which represent the particles. Setting
| in the inputs file will generate plotfiles with only two variables.
  **particle_count** represents the number of particles in a grid cell;
  **particle_mass_density** is the density on the grid resulting from the particles.

| We note that the value of :math:`a` is also written in each plotfile directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

ASCII Particle Files
--------------------

| To generate an ASCII file containing the particle positions and velocities,
  one needs to restart from a checkpoint
  file but doesn’t need to run any steps. For example, if *chk00350* exists, then one can set:
| = *chk00350*
| = 350
| = *particle_output*
| which would tell the code to restart from *chk00350*, not to take any further time steps, and to write an ASCII-format
  file called *particle_output*.
| This file has the same format as the ASCII input file:
| number of particles
| x y z mass xdot ydot zdot

Run-time Data Logs
------------------

| If you set
| in the inputs file, then at run-time the code will write out file
  *log_file* with entries every coarse
  grid time step, containing
| nstep time dt redshift a

and if **nyx.do_hydro** then also

max temp, rho-wgted temp, V-wgted temp, T @ :math:`\langle` rho :math:`\rangle`

Run-time Screen Output
----------------------

There are a number of flags that control the verbosity written to the screen at run-time. These are::

  amr.v
  nyx.v
  gravity.v
  mg.v
  particles.v
 
 These control printing about the state of the calculation (time, value of :math:`a`, etc) as well as
  timing information.


.. [1]
   See http://camb.info/
