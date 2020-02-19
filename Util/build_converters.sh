#/bin/bash
module purge

#Darwin platform
module load openmpi/3.1.3-intel_19.0.1
module load intel/19.0.1
module load python/3.5.1
module load cmake/3.15.3

export HDF5_DIR=/projects/exasky/hdf5-1.10.2-intel/

# Build blosc
if [ ! -d c-blosc ]; then
    echo "Building BLOSC ... "

    git clone https://github.com/Blosc/c-blosc.git
    cd c-blosc/
    git checkout v1.10.2
    mkdir install
    mkdir build
    cd build
    cmake .. -DCMAKE_C_FLAGS=-dynamic -DCMAKE_CXX_FLAGS=-dynamic -DCMAKE_INSTALL_PREFIX=../install
    make -j
    make install
    cd ..
    cd ..

    echo "Building BLOSC done!"
fi
BLOSC_DIR=$(pwd)/c-blosc

#INCLUDE_LOCATIONS += $(BLOSC_DIR)/install/include
#LIBRARY_LOCATIONS += $(BLOSC_DIR)/lib
#INCLUDE_LIBRARIES += -lblosc

cd convert_all_fields
make -j
cd ..

cd convert_nbody
make -j
cd ..

cd convert_lyaf
make -j
cd ..

