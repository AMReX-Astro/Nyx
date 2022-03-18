.. role:: cpp(code)
   :language: c++

.. _SUNDIALS:

Compiling Nyx with SUNDIALS 6
===============================

The following steps describe how to compile Nyx with
SUNDIALS 6 support. This library is necessary for non-adiabatic
heating-cooling Nyx runs, where ``USE_HEATCOOL=TRUE``, such as
the ``Nyx/Exec/LyA`` directory.

In order to use SUNDIALS:

#. We suggest using the Github mirror:
   https://github.com/LLNL/sundials and picking the type of
   parallelism that is appropriate for your architecture.

   To install with cuda and openmp support:
   
   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which gcc)  \
      -DCMAKE_CXX_COMPILER=$(which g++)   \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_FLAGS="-DSUNDIALS_DEBUG_CUDA_LASTERROR" \
      -DSUNDIALS_BUILD_PACKAGE_FUSED_KERNELS=ON \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"  \
      -DCUDA_ENABLE=ON  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=ON   \
      -DF2003_INTERFACE_ENABLE=OFF   \
      -DSUNDIALS_INDEX_SIZE:INT=32   \
      -DCUDA_ARCH=sm_70 ../
      make -j8
      make install -j8

   To install with openmp and no cuda support:
         
   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which gcc)  \
      -DCMAKE_CXX_COMPILER=$(which g++)   \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"  \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=ON   \
      -DF2003_INTERFACE_ENABLE=OFF   \
      -DSUNDIALS_INDEX_SIZE:INT=32 ../
      make -j8
      make install -j8

   To install with HIP support (with ROCm 4.5):

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/../instdir  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which clang) \
      -DCMAKE_CXX_COMPILER=$(which hipcc) \
      -DEXAMPLES_INSTALL_PATH=$(pwd)/../instdir/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_RELEASE=-O3 \
      -DCMAKE_CXX_FLAGS_RELEASE=-O3 \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=OFF   \
      -DF2003_INTERFACE_ENABLE=OFF \
      -DENABLE_HIP=ON \
      -DEXAMPLES_INSTALL=OFF ../
      make -j8
      make install -j8


   To install with SYCL support:

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/../instdir  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_CXX_COMPILER=$(which dpcpp)  \
      -DCMAKE_CXX_STANDARD=17 \
      -DEXAMPLES_INSTALL_PATH=$(pwd)/../instdir/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS_RELEASE=-O3  \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=OFF   \
      -DF2003_INTERFACE_ENABLE=OFF \
      -DENABLE_SYCL=ON ../
      make -j8
      make install -j8

#. Note that we give these examples for the gnu compiler or the appropriate parallel compiler.
   The compiler chosen needs to be consistent with Nyx's GNUMakefile
   variable COMP to ensure matching OMP runtime libraries for use with the OpenMP NVector. 

#. ``CUDA_ARCH`` must be set to the appropriate value for the GPU being targeted

#. For more detailed instructions for installing SUNDIALS with different flags and versions see
   the `SUNDIALS documentation <https://computing.llnl.gov/projects/sundials/sundials-software>`_.

#. In the ``GNUmakefile`` for the application which uses the interface to SUNDIALS, add
   ``USE_SUNDIALS = TRUE`` and ``SUNDIALS_ROOT=${INSTALL_PREFIX}``. Note that one must define the
   ``SUNDIALS_LIB_DIR`` make variable to point to the location where the libraries are installed
   if they are not installed in ``${INSTALL_PREFIX}/lib``. Note the default location
   for 64 is ``${INSTALL_PREFIX}/lib64``, which we override with ``-DCMAKE_INSTALL_LIBDIR=lib``.

#. If the application uses the SUNDIALS CVODE time integrator package, then the variable
   ``USE_CVODE_LIBS = TRUE`` should also be added in the ``GNUmakefile`` for the application.
   If the application used the SUNDIALS ARKode time integrator package, then the variable
   ``USE_ARKODE_LIBS = TRUE`` should be added.

Note that SUNDIALS can also be installed via Spack:

   ::
      
      spack install sundials+cuda+openmp
  
