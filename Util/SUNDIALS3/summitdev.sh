#!/bin/bash
#BSUB -P CSC308
#BSUB -W 8
#BSUB -nnodes 1
#BSUB -J Test
#BSUB -o Testo.%J
#BSUB -e Teste.%J

module load pgi/18.4
#module load cuda/9.1.85
module load cuda/9.0.69
module list
set -x

#export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
#export MPICH_RDMA_ENABLED_CUDA=1

omp=1
export OMP_NUM_THREADS=${omp}
EXE="main3d.pgi.CUDA.ex"
JSRUN="jsrun -n 1 -a 1 -g 1 -c 1 --bind=packed:${omp} "
#JSRUN="jsrun -n 4 -a 1 -c 1 -r 4 -d cyclic -b packed :${omp} js_task_info | sort"

rundir="${LSB_JOBNAME}-${LSB_JOBID}"
mkdir $rundir
cp $EXE $rundir
cp inputs $rundir
cp TREECOOL_middle $rundir
cp inputs_atomic $rundir
cd $rundir

# 0. Just run it
${JSRUN} --smpiargs="-gpu" ./${EXE} inputs

# 1. Run normally
#${JSRUN} --smpiargs="-gpu" cuda-memcheck ${EXE} inputs
#${JSRUN} --smpiargs="-gpu" cuda-memcheck ${EXE} inputs > memcheck${LSB_JOBID}.txt 

# 2. Run under nvprof and direct all stdout and stderr to nvprof.txt
#${JSRUN} --smpiargs="-gpu" nvprof --profile-child-processes ${EXE} inputs &> nvprof.txt
#${JSRUN} --smpiargs="-gpu" nvprof ${EXE} inputs &> nvprof${LSB_JOBID}.txt

# 3. Run under nvprof and store performance data in a nvvp file
# Can be converted to text using nvprof -i nvprof-timeline-%p.nvvp
#${JSRUN} --smpiargs="-gpu" nvprof --profile-child-processes -o nvprof-timeline-%p.nvvp ${EXE} inputs

# COLLECT PERFORMANCE METRICS - THIS IS MUCH SLOWER. Set nsteps=2 in the inputs files
# 4. Run under nvprof and collect metrics for a subset of kernels
#${JSRUN} --smpiargs="-gpu" nvprof --profile-child-processes --kernels '(deposit_current|gather_\w+_field|push_\w+_boris)' --analysis-metrics -o nvprof-metrics-kernel-%p.nvvp ${EXE} inputs

# 5. Run under nvprof and collect metrics for all kernels -- much slower!
#${JSRUN} --smpiargs="-gpu" nvprof --profile-child-processes --analysis-metrics -o nvprof-metrics-%p.nvvp ${EXE} inputs

#cp ../ElectromagneticPIC*.${LSB_JOBID} .
