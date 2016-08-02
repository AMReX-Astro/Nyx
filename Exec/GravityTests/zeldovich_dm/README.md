Nyx: Zeldovich Test Case
========================

This Exec sets up and runs the Zeldovich problem.

Create Initial Conditions
-------------------------

Generate the initial conditions by compiling and running
``generate_zeldovich_ics.cpp``. You may change some of the parameters in
that file, such as:

 - redshifts: initial and caustic
 - number of particles
 - box length
 - number of sheets (determines perturbation wavelength)
 - perturbation offset
 - sheet normal vector (determines direction of the perturbations)
 - Hubble rate (make sure it matches your cosmology in ``probin``!)

Then compile the ICs generator with:

    $ cd ics
    $ g++ generate_zeldovich_ics.cpp -o generate_zeldovich_ics.ex
    $ ./generate_zeldovich_ics.ex

This creates files ``zeldovich_particles.ascii`` and ``zeldovich_params.csv``.
The csv file contains the parameters used to generate the particle ascii file
and is used during analysis so you don't have to worry about adjusting params in
the analysis scripts manually.

Checking Initial Conditions
---------------------------

You might want to check the file with a quick ``head zeldovich_particles.ascii``
to make sure the number of particles in the first line is fine, and the particle
format ``x y z mass xdot ydot zdot`` is good.

You can also run

    $ cd analysis
    $ python check_ics.py

which will test ``ics/zeldovich_particles.ascii`` to make sure the sheet normal
is indeed a unit vector, that the velocities perpendicular to that are 0, and
that the density in the box matches the critical (matter) density given by the
cosmology. It also plots the particle distribution in comoving space (3d
projection) and comoving phase space for you to check by eye.

We might make this automatic in the future, but this is sufficient for now.

Run Nyx
-------

Compile the Nyx executable by checking the make options in ``build/GNUMakefile``
and running ``make`` in the build directory.

Check the ``inputs`` and ``probin`` files to make sure you match the conditions
set in your ICs. Check that ``prob_hi`` matches the size of the box and ``h`` is
the same.

If you want to run on hopper, run the normal ``qsub pbs_hopper``, adjusting for
however many cores you think your run requires. If you want to run locally, use
something like ``mpirun -np 4 ./Nyx3d.Darwin.gcc.gfortran.MPI.ex inputs >&
out``.

Analyze
-------

    $ cd analysis
    $ python process_zeldovich.py

will make a ``figs`` directory, open the ``particles.ascii`` file in each
directory beginning with ``plt``, and analyze the results. For each plt ascii
particle file, it will plot the phase space with zeldovich predictions, the 3d
projection of particle positions, and a histogram of velocities perpendicular to
the normal (should be close to 0).
