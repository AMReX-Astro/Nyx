#!/bin/bash

set -ex
##Use the same make options and input file as Nyx/Exec/LyA

##Example test:
##clone and setup sundials
#git clone --recursive git@github.com:amrex-astro/Nyx.git

##Do sundials install
#Use existing documentation
if [ ! -d "../../subprojects/sundials/instdir" ]; then
    cd ../../subprojects/sundials
    mkdir builddir instdir
    cd builddir
    cmake \
	-DCMAKE_INSTALL_PREFIX="$(pwd)/../instdir"  \
	-DCMAKE_INSTALL_LIBDIR=lib \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DCMAKE_C_COMPILER="$(which gcc)"  \
	-DCMAKE_CXX_COMPILER="$(which g++)"   \
	-DCMAKE_CUDA_HOST_COMPILER="$(which g++)"    \
	-DEXAMPLES_INSTALL_PATH="$(pwd)/../instdir/examples" \
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
    cd ../../../Exec/HeatCoolTests
fi

export OMP_NUM_THREADS=4
export EXE='mpiexec -n 2 ./Nyx3d.gnu.TPROF.MPI.OMP.ex '

cd ../../
cd Exec/LyA
make -j8 SUNDIALS_ROOT=$(pwd)/../../subprojects/sundials/instdir
# Write to ascii if you want (format is a short header, then (i,j,k) data data data for each multifab)
#index mapping can be found in Nyx/Source/Hydro/IndexDefines.H or in the documentation
# The order in each Chunk.# is: old state data, old derived data, new state data, source terms from hydro, a correction term from enforcing conservation, the lagged IR term from old reactions
#${EXE} inputs max_step=1 nyx.hctest_example_write=1 nyx.v=2 fab.format=ASCII
${EXE} inputs max_step=1 nyx.hctest_example_write=1 nyx.v=2

cd ../HeatCoolTests
ln -s $(pwd)/../LyA/hctest .
make -j8 SUNDIALS_ROOT=$(pwd)/../../subprojects/sundials/instdir
${EXE} hctest/inputs.0
