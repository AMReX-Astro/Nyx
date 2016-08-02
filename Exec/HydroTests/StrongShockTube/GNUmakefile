# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= /project/projectdirs/nyx/src/BoxLib

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../../..

# compilation options
COMP    = gcc
FCOMP   = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
