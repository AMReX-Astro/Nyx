# AMREX_HOME defines the directory in which we will find all the BoxLib code
AMREX_HOME ?= ../../../../amrex
HDF5_DIR   ?=
 

# TOP defines the directory in which we will find Source, Exec, etc
TOP ?= ../../../../Nyx

EBASE     = WriteVelsDMPC

# compilation options
COMP    = gnu

USE_MPI = FALSE
USE_OMP = FALSE

PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE
USE_SINGLE_PRECISION_PARTICLES=TRUE
USE_CONST_SPECIES = TRUE

DIM      = 3

#include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
