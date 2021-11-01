Use the same make options and input file as Nyx/Exec/LyA
Example test:
##clone and setup sundials
#git clone --recursive git@github.com:amrex-astro/Nyx.git

##Do sundials install
#Use existing documentation

export EXE="export OMP_NUM_THREADS=4;mpiexec -n 2 ./Nyx3d.gnu.TPROF.MPI.OMP.ex "
cd Nyx/Exec/LyA
make -j8 SUNDIALS_ROOT=$(pwd)/../../subproject/sundials/instdir
${EXE} inputs max_step=1 nyx.hctest_example_write=1 nyx.v=2

cd ../../Nyx/Exec/HeatCoolTests
ln -s $(pwd)/../../LyA/hctest .
make -j8 SUNDIALS_ROOT=$(pwd)/../../subproject/sundials/instdir
${EXE} hctest/inputs.0