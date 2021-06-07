# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../subprojects/amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = gnu
USE_OMP = FALSE
USE_MPI = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = FALSE
DEBUG     = FALSE

# physics
DIM      = 3

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
