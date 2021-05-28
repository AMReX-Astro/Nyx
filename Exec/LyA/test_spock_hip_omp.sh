#!/bin/bash                                                                                                       
#SBATCH -A hep114                                                                                                       
#SBATCH -J hip_omp
#SBATCH -o %x-%j.out                                                                                                       
#SBATCH -t 00:02:30                                                                                                       
#SBATCH -p ecp                                                                                                       
#SBATCH -N 2
# An example of running a 4 processes job
# modules needed for compilation
module load PrgEnv-gnu craype-accel-amd-gfx908 rocm

#compile Sundials using https://amrex-astro.github.io/Nyx/docs_html/getting_started/NyxSundials.html
#compile with jmsexton03/amrex hip_omp_flags_namespace branch
#cd ../../subprojects/amrex
#git remote add jmsexton03 https://github.com/jmsexton03/amrex
#git fetch jmsexton03
#git checkout hip_omp_flags_namespace
#cd ../../Exec/LyA
#make -j USE_HIP=TRUE USE_OMP=TRUE SUNDIALS_ROOT=../../subprojects/sundials/instdir-rocm/ AMREX_HOME=../../subprojects/amrex/ 
export omp=4; export OMP_NUM_THREADS=${omp}; srun -n 4 -c${omp} --gpus-per-task=1 --gpu-bind=closest ./Nyx3d.hip.x86-rome.TPROF.MPI.OMP.HIP.ex inputs nyx.binary_particle_file= ${WORLDWORK}/ast160/ICs/256sss_20mpc.nyx amr.n_cell=256 256 256 amr.max_grid_size=64 particles.n_readers=1 particles.nreaders=1 max_step=5 amrex.max_gpu_streams=${omp}

#srun -n 64 --ntasks-per-node=4 ./Nyx3d.hip.x86-rome.TPROF.MPI.HIP.ex inputs max_step=5 nyx.binary_particle_file= 1024s_20mpc.nyx amr.n_cell=1024 1024 1024 amr.max_grid_size=128 | tee out_1024_128_8.txt


#srun -n 16 --ntasks-per-node=4 ./Nyx3d.hip.x86-rome.TPROF.MPI.HIP.ex inputs.scaling.768 max_step=5 max_step=50 nyx.minimize_memory=1 nyx.shrink_to_fit=1 amrex.max_gpu_streams=8  amrex.the_arena_init_size=1000 | tee out_768_128_16_scal8.txt

#srun -n 16 --ntasks-per-node=4 ./Nyx3d.hip.x86-rome.TPROF.MPI.HIP.ex inputs.scaling.1024 max_step=5 max_step=50 nyx.minimize_memory=1 nyx.shrink_to_fit=1 amrex.max_gpu_streams=1 amr.regrid_on_restart=0 nyx.v=2 particles.v=2 amr.v=3 gravity.v=3 | tee out_1024_128_16_scaling.txt

#srun -n 16 --ntasks-per-node=4 ./Nyx3d.hip.x86-rome.TPROF.MPI.HIP.ex inputs.scaling.2048 max_step=5 max_step=50 nyx.minimize_memory=1 nyx.shrink_to_fit=1 amrex.max_gpu_streams=1 | tee out_2048_128_16.txt

#srun -n 16 --ntasks-per-node=4 ./Nyx3d.hip.x86-rome.TPROF.MPI.HIP.ex inputs max_step=5 nyx.binary_particle_file= 1024s_20mpc.nyx amr.n_cell=1024 1024 1024 amr.max_grid_size=128 nyx.minimize_memory=1 nyx.shrink_to_fit=1 amrex.max_gpu_streams=1 | tee out_1024_128_16.txt
