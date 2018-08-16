#!/bin/bash

#SBATCH -N 1
#SBATCH -C knl,quad,cache
#SBATCH -p debug
#SBATCH -A nyx
#SBATCH -J Nyx
#SBATCH -t 00:10:00
#SBATCH --export=ALL

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

cd $SLURM_SUBMIT_DIR

date
srun -n 8 -c 32 --cpu_bind=cores ./Nyx3d.intel.mic-knl.PROF.MPI.OMP.ex inputs
date
