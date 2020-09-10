# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../../..

#Ascent location
USE_ASCENT_INSITU=TRUE 
ASCENT_HOME=/project/projectdirs/alpine/software/ascent/gnu/7.3.0/ascent-install/

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
#TINY_PROFILE = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = FALSE

USE_CONST_SPECIES = TRUE

USE_CXX_PROB = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/HydroTests/Make.Nyx
