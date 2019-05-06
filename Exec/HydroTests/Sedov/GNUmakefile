# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../../..

# compilation options
COMP    = pgi

USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = TRUE
TINY_PROFILE = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = FALSE

USE_OWN_BCS = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
