Initial conditions
===================

There are a couple main ways in which initial conditions can be set in Nyx:
using an ASCII particle file, using binary particle file(s), using
a uniform particle setup, using a random particle setup, using a binary mesh file setup, or
using an analytic setup.
As said in the Units section of the documentation, the units are: Mpc, M\ :math:`_\odot`, and km/s,
and particle velocities should be peculiar proper velocities.


Start from an ASCII file
------------------------

To enable this option, set::
  
  nyx.particle_init_type = AsciiFile
  nyx.ascii_particle_file = *particle_file*

Here *particle_file* is the user-specified name of the file. The first line in this file is
(long) assumed to contain the number of particles. Each line after that contains

x y z mass vx vy vz



Start from a binary file
------------------------

To enable this option, set::

  nyx.particle_init_type = BinaryFile
  nyx.binary_particle_file = *particle_file*
  
With binary file, the header should have 3 numbers:
(long) NP, which is the total number of particles in the file
followed by the (int) DM=3 (number of dimensions), and (int) NX=4 (number of "extra" fields).
Following this small header, 7 float numbers should be listed for each particle as before:
x y z mass vx vy vz.

The main difference between the ASCII and binary format thus amounts to a different header.
Here is an example C++ code which writes a Nyx-readable binary file::

      fwrite(&npart, sizeof(long), 1, outFile);
      fwrite(&DM, sizeof(int), 1, outFile);
      fwrite(&NX, sizeof(int), 1, outFile);
      for(i=0; i<Npart; i++) {
         fwrite(&x[i], sizeof(float), 1, outFile);
         fwrite(&y[i], sizeof(float), 1, outFile);
         fwrite(&z[i], sizeof(float), 1, outFile);
         fwrite(&mass[i], sizeof(float), 1, outFile);
         fwrite(&vx[i], sizeof(float), 1, outFile);
         fwrite(&vy[i], sizeof(float), 1, outFile);
         fwrite(&vz[i], sizeof(float), 1, outFile);
      }


Start from a binary "meta" file
-------------------------------

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


Start from a plotfile or checkpoint
-----------------------------------

To enable this option, set::

  nyx.particle_init_type = Restart
  nyx.restart_particle_file = *plot_file*

In this case the *plot_file* should contain particles in directory *DM*. Testing of this
functionality is mainly for the current default *Version_Two_Dot_Zero_single*.


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


Initial conditions for testing purposes
---------------------------------------

The following are used for code testing purposes and will not result in a meaningful cosmological simulation.


Random placement
^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = Random
  
There are then a number of parameters to set, for example::
  
  nyx.particle_initrandom_count = 100000
  nyx.particle_initrandom_mass_total = 100000
  nyx.particle_initrandom_iseed = 15
  nyx.fix_random_seed = 0


Random placement (1 particle per grid cell)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = RandomPerCell
  
Then only set the mass per particle::

  nyx.particle_initrandom_mass = 1

Note to increase the number of cells and keep the problem domain size 
and total mass fixed, the mass per particle must decrease proportionally.
An alternative is to set the total mass of all particles in the simulation.
This will cause Nyx to scale the mass per particle to fit the number of cells.
The following will all have the same total mass::

  nyx.particle_initrandom_mass = 1
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = 1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = -1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = 0.125
  amr.n_cell = 128 128 128

  nyx.particle_initrandom_mass = 0.125
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 128 128 128

  nyx.particle_initrandom_mass = -1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 128 128 128


Uniform placement
^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = OnePerCell
  
There are then a number of parameters to set, for example::
  
  nyx.particle_inituniform_mass = 1
  nyx.particle_inituniform_vx = -1
  nyx.particle_inituniform_vy = 1
  nyx.particle_inituniform_vz = 1

Initial Multifab-based setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = Cosmological
  nyx.do_readinics = 1

Then set the directory name of the MultiFab to restart the state variables from::

  nyx.readin_ics_fname = "mf"
  
Initial Analytic Problem Setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.do_santa_barbara = 0
  nyx.init_sb_vels = 0

For different executable directories, the ``Prob.cpp`` setup can be further customised
with ``prob.`` input flags. For the HydroTests directory, ``prob.prob_type=0`` corresponds
to Sod, StrongShockTube and DoubleRarefaction type tests, and ``prob.prob_type!=0``
corresponds to the Sedov type tests.
