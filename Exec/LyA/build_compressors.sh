#/bin/bash

#PGI platform
#module purge

#module load python/3.5.1
#module load cmake/3.15.3
#module load pgi/19.7
#module load cuda/10.1
#module load openmpi/3.1.4-pgi_19.7
#export CUDA_HOME=/projects/opt/centos7/cuda/10.1/

if [ -z "$OLCF_C_BLOSC_ROOT" ]; then

  # Build blosc
  if [ ! -d c-blosc ]; then
      echo "Building BLOSC ... "

      git clone https://github.com/Blosc/c-blosc.git
      cd c-blosc/
      git checkout v1.10.2
      mkdir install
      mkdir build
      cd build
      #Need to patch ../blosc/blosc-export.h to add -- defined(__PGI)
      sed -i.bak '20s/$/ || defined(__PGI) /' ../blosc/blosc-export.h
      cmake .. -DCMAKE_INSTALL_PREFIX=../install -DBUILD_STATIC=ON
      make -j
      make install
      cd ..
      cd ..

      echo "Building BLOSC done!"
  fi
  #export BLOSC_DIR=$(pwd)/c-blosc/install
  #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/c-blosc/install/lib

else

  echo "Using olcf blosc"
  export BLOSC_DIR=$OLCF_C_BLOSC_ROOT
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BLOS_DIR/lib

fi
