#!/bin/bash
#BSUB -P PROJID
#BSUB -W 10 
#BSUB -nnodes 1
#BSUB -J test 
#BSUB -o test.%J
#BSUB -e test.%J

module list
set -x

module load gcc/10.2.0 cuda/11.2.0 hdf5 python

git clone --recursive -b development https://github.com/amrex-astro/Nyx.git

cd Nyx/Exec/LyA

make -j USE_CUDA=TRUE USE_OMP=FALSE SUNDIALS_ROOT=../../subprojects/sundials/instdir/ AMREX_HOME=../../subprojects/amrex/ USE_HEATCOOL=FALSE
jsrun -n 1 -a 1 -c 1 -g 1 ./Nyx3d.gnu.TPROF.MPI.CUDA.ex inputs.rt max_step=10
jsrun -n 1 -a 1 -c 1 -g 1 ./Nyx3d.gnu.TPROF.MPI.CUDA.ex inputs.rt max_step=10 amr.derive_plot_vars = particle_mass_density particle_x_velocity particle_y_velocity particle_z_velocity amr.plot_file=plt_exact

cd ../../subprojects/amrex/Tools/Plotfile/  
make -j programs=fcompare

cd ../../Util/Converters/assign_particle_vels/
cp ../../../subprojects/amrex/Tools/Plotfile/fcompare* .
#Current setup was getting strange behavior when using paths, so either copy or symlink
ln -s ../../../Exec/LyA/plt00010 .

##For example changes from original inputs file:
#diff inputs.rt ../../../Exec/LyA/inputs.rt
make -j AMREX_HOME=../../../subprojects/amrex/
jsrun -n1 ./WriteVelsDMPC3d.gnu.ex inputs.rt max_step=10
./fcompare.gnu.ex ../../../Exec/LyA/plt_exact00010 plt00010_gridDM

cd ../Plotfile2HDF5_grids
make -j AMREX_HOME=../../subprojects/amrex/ COMP=gnu HDF5_DIR=$OLCF_HDF5_ROOT
#First run does normal conversion and includes metadata
jsrun -n1 ./convert3d.gnu.PROF.MPI.ex input_path=../assign_particle_vels/plt00010 output_path=plt00010.hdf5
#Second run appends particle grid data as mesh datasets
jsrun -n1 ./convert3d.gnu.PROF.MPI.ex input_path=../assign_particle_vels/plt00010_gridDM output_path=plt00010.hdf5
h5dump -H plt00010.h5
