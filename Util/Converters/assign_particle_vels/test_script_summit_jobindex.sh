#!/bin/bash
#BSUB -P hep114
#BSUB -W 20
#BSUB -nnodes 1
#BSUB -J test[10]
#BSUB -o test.%J
#BSUB -e test.%J

module list
set -x

module load gcc/10.2.0 cuda/11.2.0 hdf5 python

git clone --recursive -b development https://github.com/amrex-astro/Nyx.git

export string_short_number="$LSB_JOBINDEX"
export string_number="$(echo `printf "%05.0f" $LSB_JOBINDEX`)"
export submit_dir="$LS_SUBCWD"

echo ${string_number}
echo ${string_short_number}

cd Nyx/Exec/LyA

make -j USE_CUDA=FALSE USE_OMP=FALSE SUNDIALS_ROOT=../../subprojects/sundials/instdir/ AMREX_HOME=../../subprojects/amrex/ USE_HEATCOOL=FALSE
jsrun -n 1 -a 1 -c 1 -g 1 ./Nyx3d.gnu.TPROF.MPI.ex inputs.rt max_step="${string_short_number}"
jsrun -n 1 -a 1 -c 1 -g 1 ./Nyx3d.gnu.TPROF.MPI.ex inputs.rt max_step="${string_short_number}" amr.derive_plot_vars = particle_mass_density particle_x_velocity particle_y_velocity particle_z_velocity amr.plot_file=plt_exact

cd ../../subprojects/amrex/Tools/Plotfile/  
make -j programs=fcompare

cd ../../../../Util/Converters/assign_particle_vels/
cp ../../../subprojects/amrex/Tools/Plotfile/fcompare* .
#Current setup was getting strange behavior when using paths, so either copy or symlink
ln -s ../../../Exec/LyA/plt${string_number} .

##For example changes from original inputs file:
#diff inputs.rt ../../../Exec/LyA/inputs.rt
make -j AMREX_HOME=../../../subprojects/amrex/ USE_CUDA=FALSE USE_MPI=TRUE
jsrun -n1 -a 1 -c 1 -g 1 ./WriteVelsDMPC3d.gnu.MPI.ex inputs.rt max_step="${string_short_number}" nyx.restart_particle_file=plt${string_number}
./fcompare.gnu.ex ../../../Exec/LyA/plt_exact${string_number} plt${string_number}_gridDM

./fcompare.gnu.ex ../../../Exec/LyA/plt_exact${string_number} plt${string_number}_gridDM &> $submit_dir/fcompare_results.txt

cd ../Plotfile2HDF5_grids
make -j AMREX_HOME=../../../subprojects/amrex/ COMP=gnu HDF5_DIR=$OLCF_HDF5_ROOT
#First run does normal conversion and includes metadata
jsrun -n1 ./convert3d.gnu.PROF.MPI.ex input_path=../assign_particle_vels/plt${string_number} output_path=plt${string_number}.hdf5
#Second run appends particle grid data as mesh datasets
jsrun -n1 ./convert3d.gnu.PROF.MPI.ex input_path=../assign_particle_vels/plt${string_number}_gridDM output_path=plt${string_number}.hdf5
h5dump -H "plt${string_number}.hdf5"
h5dump -H "plt${string_number}.hdf5" &> $submit_dir/h5dump-H_results.txt
