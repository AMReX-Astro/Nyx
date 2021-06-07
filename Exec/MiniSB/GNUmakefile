# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../subprojects/amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    =  gnu

USE_MPI = TRUE
USE_OMP = FALSE

# Analysis
REEBER = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_CONST_SPECIES = TRUE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
