Initial conditions
===================

There are two main ways in which initial conditions can be set in Nyx:
read from an ASCII file and read from a binary file(s). In both cases file should have 3 numbers
in the header: int DIM=3 (number of dimensions), int NX=4 (number of "extra" fields), and long NP
which is the total number of particles in the file.
This header is then folloed by the NP lines where each is:
x, y, z, mass, vx, vy, vz.
As said in the Units section of the documentation, the units are: Mpc, M\ :math:`_\odot`, and km/s,
and vx, vy, and vz are the components of the peculiar proper velocity.


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

This option allows you to read particles from a series of files rather than
just a single file. This is very convenient for large simulations.
To enable this option, set::

  nyx.particle_init_type = BinaryMetaFile
  nyx.binary_particle_file =*particle file*

In this case the *particle_file* you specify is an ASCII file specifying a
list of file names with full paths. Each of the files in this list is assumed
to be binary and is read sequentially (individual files are read in parallel) in
the order listed.

Since individual files are read sequentially, more particles should be read before
redistributing across MPI ranks. This is set by optimizing the maximum number of
readers and increasing the number of particles per read::

  amr.nreaders
  amr.nparts_per_read

Reading SPH particles
---------------------

The above initialization from a single particle specie assumes that baryons trace dark matter.
Baryon density and velocity is set by CIC-mapping particles onto the Eulerian grid.
Alternatively, one can initialize the baryonic gas from the SPH
particles. To enable this option, you must set::

    nyx.do_santa_barbara = 1
    nyx.init_with_sph_particles =1

The SPH-type particles can then be read in by setting
where *sph_particle_file* is the user-specified name of the file
containing the SPH particles. The type of *sph_particle_file*
must be the same (Ascii, Binary or BinaryMeta) as the dark matter particle
file as specified by
The SPH particles will be discarded by the code once the grid data has been initialized.


Test initial conditions
=======================

The following are used for code testing purposes and will not result in a meaningful cosmological simulation.


Random placement
----------------

To enable this option, set::

  nyx.particle_init_type = Random
  
There are then a number of parameters to set, for example::
  
  nyx.particle_initrandom_count = 100000
  nyx.particle_initrandom_mass = 1
  nyx.particle_initrandom_iseed = 15

Random placement (1 particle per grid cell)
-------------------------------------------

To enable this option, set::

  nyx.particle_init_type = RandomPerCell
  
Then only set the mass per particle::

  nyx.particle_initrandom_mass = 1

Note to increase the number of cells and keep the problem domain size 
and total mass fixed, the mass per particle must decrease proportionally.

Uniform placement
-----------------

To enable this option, set::

  nyx.particle_init_type = OnePerCell
  
There are then a number of parameters to set, for example::
  
  nyx.particle_inituniform_mass = 1
  nyx.particle_inituniform_vx = -1
  nyx.particle_inituniform_vy = 1
  nyx.particle_inituniform_vz = 1
