#!/bin/bash

PROFILE_RANK=0
echo "hi from $PMIX_RANK"
if [ $PMIX_RANK == $PROFILE_RANK ]; then
    nvprof --openmp-profiling off --device-buffer-size 64 -f -o profile.${PMIX_RANK}.${LSB_JOBID}.nvvp "$@"
else
    "$@"
fi
