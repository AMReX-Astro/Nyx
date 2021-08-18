Building Nyx with CMake
=======================

CMake build is a two-step process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (``builddir``).
Next, the actual build is performed by invoking ``make`` from within ``builddir``.
If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_
from the HEP Software foundation is the perfect place to get started with it.

The CMake build process for Nyx is summarized as follows:

#. Create and enter the build directory:

   .. highlight:: console

   ::

       mkdir /path/to/builddir
       cd    /path/to/builddir

#. Perform the configuration step:

   .. highlight:: console

   ::

      cmake [nyx options] [dependencies options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] /path/to/Nyx/repo


   In the above snippet, ``[nyx options]`` indicates one or more of the options
   for the customization of the build listed in the :ref:`table <tab:nyxcmakeoptions>` below,
   whereas ``[dependecies options]`` are one or more of the CMake configuration options for the Nyx dependencies,
   i.e. `AMReX CMake options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_,
   and SUNDIALS. If the option ``CMAKE_BUILD_TYPE`` is omitted, ``CMAKE_BUILD_TYPE=Release`` is assumed.
   For example, to enable AMReX profiling capabilities in Nyx, configure as follows:

   .. code:: shell

             > cmake [nyx options] -DAMReX_TINY_PROFILE=yes -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] ..

#. Build the executable you are interested in (a small
   version of the Santa Barbara problem in this example):

   .. highlight:: console

   ::

      cd /path/to/builddir/Nyx/Exec/MiniSB
      make


   The resulting executable is ``nyx_MiniSB``.

.. note::
   **Nyx requires CMake 3.14 or higher.**


.. raw:: latex

   \begin{center}

.. _tab:nyxcmakeoptions:

.. table:: Nyx configuration options

   +-----------------+------------------------------+------------------+-------------+
   | Option name     | Description                  | Possible values  | Default     |
   |                 |                              |                  | value       |
   +=================+==============================+==================+=============+
   | CMAKE\_CXX\     | User-defined C++ flags       | valid C++        | None        |
   | _FLAGS          |                              | compiler flags   |             |
   +-----------------+------------------------------+------------------+-------------+
   | CMAKE\_CUDA\    | User-defined CUDA flags      | valid CUDA       | None        |
   | _FLAGS          |                              | compiler flags   |             |
   +-----------------+------------------------------+------------------+-------------+
   | BUILD\_SHARED\  | Build dependencies as shared | no/yes           | no          |
   | _LIBS           | libraries                    |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_MPI        | Enable build with MPI        | no/yes           | yes         |
   |                 |                              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_MPI\_THREAD| Concurrent MPI calls from    | no/yes           | yes         |
   | \_MULTIPLE      | multiple threads             |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_OMP        | Enable build with OpenMP     | no/yes           | no          |
   |                 |                              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_GPU\_      | On-node, accelerated GPU \   | NONE             | NONE,SYCL,\ |
   | BACKEND         | backend                      |                  | CUDA,HIP    |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_HYDRO      | Run with baryon hydro solve  | no/yes           | yes         |
   |                 | and fields                   |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_HEATCOOL   | Run with H and He heating-   | no/yes           | no          |
   |                 | cooling effects              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_CONST\_    | Don't evolve H and He, treat | no/yes           | no          |
   | SPECIES         | them as constant             |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_CGS        | Evolve quantities in CGS     | no/yes           | no          |
   |                 | units instead of code units  |                  |             |
   +-----------------+------------------------------+------------------+-------------+
.. raw:: latex

   \end{center}


Working with Git submodules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default Nyx CMake searches the system for existing installations of the required dependencies
(AMReX is always required, SUNDIALS may be required depending on the configuration options).
If the required dependencies are not found on the system, Nyx CMake will automatically checkout
and build those dependencies as part of its build process. To this end, Nyx CMake relies on git
submodules to checkout the AMReX and/or SUNDIALS git repositories. In what follows, we will
focus on the AMReX git submodule only, but the same concepts apply unchanged to the
SUNDIAL submodule as well.


If the AMReX submodule is not initialized, Nyx CMake will initialize it and checkout
the commit the submodule is pointing to. If instead the AMReX  submodule has already
been manually initialized and a custom commit has been checked out, that commit will
be used. For Nyx development or testing, you may need to build with a different
branch or release of AMReX.

The ``subprojects/amrex`` directory is a Git repo, so use all standard Git
commands. Either ``cd subprojects/amrex`` to run Git commands in the ``amrex``
directory, or use ``git -C subprojects/amrex`` in the Nyx repo. For
instance, to build with the ``my-amrex-branch`` branch of the AMReX repo:

.. code:: shell

    > git -C subprojects/amrex checkout my-amrex-branch
    > git status
    ...
    modified:   subprojects/amrex (new commits)

The branch ``my-amrex-branch`` will then be used when building Nyx.

To revert to the default version of the AMReX submodule, run ``git submodule
update``:

.. code:: shell

    > git submodule update
    > git status
    ...
    nothing to commit, working tree clean

You can edit, commit, pull, and push AMReX changes from ``subprojects/amrex``.
AMReX development is outside the scope of this document. Run ``git status`` in
the top-level Nyx repo to see whether the AMReX submodule has new commits,
modified files, or untracked files.

To update the AMReX submodule referenced by Nyx:

.. code:: shell

    > git -C subprojects/amrex checkout UPDATED_AMREX_COMMIT_SHA1
    > git add subprojects/amrex
    > git commit -m 'Updating AMReX version'

This will only update the AMReX SHA-1 referenced by Nyx. Uncommitted AMReX
changes and untracked AMReX files under ``subprojects/amrex`` are not added by
``git add subprojects/amrex``. (To commit to the AMReX repo, change directories
to ``subprojects/amrex`` and run Git commands there, before ``git add
subprojects/amrex``.)

.. note::

    Only update the AMReX submodule reference in coordination with the other
    Nyx developers!


.. _sec:build:external:

Using existing installations of required dependencies
-----------------------------------------------------

You may prefer to build Nyx against an AMReX and/or SUNDIALS installation already
present on your system. Unless these installations are located in standard system
paths, you need to tell Nyx CMake where to look for them.

.. code:: shell

    > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] [nyx options] -DAMReX_ROOT=/absolute/path/to/amrex/installdir /path/to/Nyx/repo

In the example above, ``-DAMReX_ROOT=/absolute/path/to/amrex/installdir`` instructs CMake to search
``/absolute/path/to/amrex/installdir`` before searching system paths for an available AMReX installation.
``AMReX_ROOT`` can also be set as an environmental variable instead of passing it as a command line option.
Similarly, you can define a ``SUNDIALS_ROOT`` variable, either via command line or the environment, to
teach CMake where to look for SUNDIALS. Choose one of the ``CMAKE_BUILD_TYPE`` to control the level of
optimization, the option ``-DCMAKE_BUILD_TYPE=Release`` will give the same defaults as GMake.


Few more notes on building Nyx
-----------------------------------

The system defaults compilers can be overwritten as follows:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=<c++-compiler> -DCMAKE_Fortran_COMPILER=<f90-compiler> [options]  ..

When building on a platform that uses the ``module`` utility, use either
the above command (with full path to the compilers) or the following:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn [options] ..

Nyx uses the same compiler flags used to build AMReX, unless
``CMAKE_Fortran_FLAGS``/``CMAKE_CXX_FLAGS`` is explicitly provided, or
the environmental variables ``FFLAGS``/``CXXFLAGS`` are set.


For GPU builds, Nyx relies on the `AMReX GPU build infrastructure <https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-with-cmake>`_
. The target architecture to build for can be specified via the AMReX configuration option ``-DAMReX_CUDA_ARCH=<target-architecture>``,
or by defining the *environmental variable* ``AMREX_CUDA_ARCH`` (all caps). If no GPU architecture is specified,
CMake will try to determine which GPU is supported by the system.


Building Nyx for Cori (NERSC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard build
--------------

For the Cori cluster at NERSC, you first need to load/unload modules required to build Nyx.

.. code:: shell

    > module unload altd
    > module unload darshan
    > module load cmake/3.14.0

The default options for Cori are the **Haswell** architecture and **Intel** compiler, if you want to compile with the **Knight's Landing (KNL)** architecture:

.. code:: shell

    > module swap craype-haswell craype-mic-knl

Or use the **GNU** compiler:

.. code:: shell

    > module swap PrgEnv-intel PrgEnv-gnu

Now Nyx can be built.

.. note::

    The load/unload modules options could be saved in the `~/.bash_profile.ext`


GPU build
---------

To compile on the GPU nodes in Cori, you first need to purge your modules, most of which won't work on the GPU nodes

.. code:: shell

    > module purge

Then, you need to load the following modules:

.. code:: shell

    > module load cgpu gcc/7.3.0 cuda/11.1.1 openmpi/4.0.3 cmake/3.14.4

Currently, you need to use OpenMPI; mvapich2 seems not to work.

Then, you need to use slurm to request access to a GPU node:

.. code:: shell

    > salloc -N 1 -t 02:00:00 -c 10 -C gpu -A m1759 --gres=gpu:8 --exclusive

This reservers an entire GPU node for your job. Note that you can’t cross-compile for the GPU nodes - you have to log on to one and then build your software.

Finally, navigate to the base of the Nyx repository and compile in GPU mode:

.. code:: shell

    > cd Nyx
    > mkdir build
    > cd build
    > cmake -DNyx_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=Volta -DCMAKE_CXX_COMPILER=g++ ..
    > make -j

For more information about GPU nodes in Cori -- `<https://docs-dev.nersc.gov/cgpu/>`_

Building Nyx for Summit (OLCF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the Summit cluster at OLCF, you first need to load/unload modules required to build Nyx.

.. code:: shell

    > module load gcc
    > module load cmake/3.14.0

To build Nyx for GPUs, you need to load cuda module:

.. code:: shell

    > module load cuda/11.0.3

To compile:

.. code:: shell

    > cd Nyx
    > mdkir build
    > cd build
    > cmake -DNyx_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=Volta -DCMAKE_C_COMPILER=$(which gcc)  -DCMAKE_CXX_COMPILER=$(which g++)   -DCMAKE_CUDA_HOST_COMPILER=$(which g++)  -DCMAKE_CUDA_ARCHITECTURES=70 ..
    > make -j

For more information about the Summit cluster: `<https://www.olcf.ornl.gov/for-users/system-user-guides/summit/>`_

.. note::

    The load/unload modules options could be saved in the `~/.profile`
