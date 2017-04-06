# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME = $(HOME)/BoxLib

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../../..

# compilation options
#COMP    = Cray
#FCOMP   = Cray
COMP    = gcc
FCOMP   = gfortran
USE_MPI = TRUE
USE_OMP = TRUE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = FALSE
USE_HEATCOOL = FALSE

#USE_CGS = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx

